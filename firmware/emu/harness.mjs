#!/usr/bin/env node
/* array16 beam controller: the compiled RP2040 firmware in the rp2040js emulator.
 *
 * Boots the firmware image (UF2, Intel HEX or raw .bin) in rp2040js, plays scenario scripts
 * against it over emulated USB CDC, records the 74HC595 control pins with emulated timestamps,
 * decodes the register chain into latched 56-bit words and prints everything as JSON lines.
 * tests/test_firmware_emulated.py compares the result with array16/fwsim.py.
 *
 * rp2040js lacks or simplifies several things this firmware relies on. The harness fills them
 * in from outside (README.md has the details):
 *   - PIO timing: rp2040js runs the PIO in wall-clock batches and ignores delay cycles and the
 *     clock divider. Here the PIO is stepped once per clk_sys cycle in emulated time, honouring
 *     delays and the divider.
 *   - PIO FIFO join and clear: rp2040js ignores FJOIN_TX/FJOIN_RX, so the shifter's FIFO would be
 *     4 words instead of 8, and pio_sm_clear_fifos() (which toggles FJOIN_RX) would clear nothing.
 *   - Frequency counter FC0: absent (the SDK's frequency_count_khz() would spin forever). Modelled
 *     from the scenario's reference clocks, with the real 0.98 us * 2^interval measurement time.
 *   - GPOUT0..3 and GPIN0/1: the registers exist but no clock comes out. Modelled: the selected
 *     source divided by the DIV register toggles the GPOUT pin's input (the slot clock on GPIO21
 *     is read back by the firmware's ticker state machine).
 *   - USB CDC host: rp2040js's host model sends zero-length OUT packets whenever the device arms
 *     its endpoint. Here an armed OUT endpoint waits for data, and the script decides the packets.
 *   - ADC: CC1/CC2 voltages come from the scenario.
 *
 * Usage: node firmware/emu/harness.mjs [--image FILE] [--rp2040js DIR] (--script FILE|- | --demo)
 *                                      [--trace] [--main-c FILE]
 */
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));

// ---- command line ----------------------------------------------------------------------------
function parseArgs(argv) {
  const a = { trace: false, demo: false };
  for (let i = 2; i < argv.length; i++) {
    const k = argv[i];
    const val = () => {
      if (i + 1 >= argv.length) throw new Error(`${k} needs a value`);
      return argv[++i];
    };
    if (k === '--image') a.image = val();
    else if (k === '--rp2040js') a.rp2040js = val();
    else if (k === '--script') a.script = val();
    else if (k === '--main-c') a.mainC = val();
    else if (k === '--trace') a.trace = true;
    else if (k === '--demo') a.demo = true;
    else if (k === '--help' || k === '-h') a.help = true;
    else throw new Error(`unknown argument ${k}`);
  }
  return a;
}

const args = parseArgs(process.argv);
if (args.help) {
  const src = fs.readFileSync(fileURLToPath(import.meta.url), 'utf8');
  console.log(src.slice(src.indexOf('/*'), src.indexOf('*/') + 2));
  process.exit(0);
}
const RP2040JS = path.resolve(args.rp2040js || process.env.RP2040JS || path.join(os.homedir(), 'opt/src/rp2040js'));
const IMAGE = path.resolve(args.image || process.env.ARRAY16_FW_IMAGE || path.join(HERE, '../build/array16_controller.uf2'));
const MAIN_C = path.resolve(args.mainC ?? path.join(HERE, '../rp2040/main.c'));

function emit(obj) {
  process.stdout.write(JSON.stringify(obj) + '\n');
}

/** Exit once stdout has drained (writes to a pipe can still be queued on macOS). */
function finish(code) {
  process.exitCode = code;
  process.stdout.write('', () => process.exit(code));
}

function fail(message) {
  fs.writeSync(1, JSON.stringify({ ev: 'error', message }) + '\n'); // nothing else is queued yet
  process.exit(2);
}

if (!fs.existsSync(path.join(RP2040JS, 'src/rp2040.ts'))) fail(`rp2040js not found in ${RP2040JS}`);
if (!fs.existsSync(path.join(RP2040JS, 'node_modules/tsx'))) fail(`run npm install in ${RP2040JS}`);
if (!fs.existsSync(IMAGE) || !fs.statSync(IMAGE).isFile()) fail(`firmware image not found: ${IMAGE}`);

// rp2040js is TypeScript; tsx (one of its dev dependencies) compiles it on import. Its tsconfig
// matters: the classes rely on useDefineForClassFields = false.
const { register } = await import(pathToFileURL(path.join(RP2040JS, 'node_modules/tsx/dist/esm/api/index.mjs')).href);
register({ tsconfig: path.join(RP2040JS, 'tsconfig.json') });
const rp = (p) => import(pathToFileURL(path.join(RP2040JS, p)).href);
const { RP2040 } = await rp('src/rp2040.ts');
const { SimulationClock } = await rp('src/clock/simulation-clock.ts');
const { USBCDC } = await rp('src/usb/cdc.ts');
const { FIFO } = await rp('src/utils/fifo.ts');
const { GPIOPinState } = await rp('src/gpio-pin.ts');
const { bootromB1 } = await rp('demo/bootrom.ts');

// ---- pins: read from the firmware source so the harness follows it ---------------------------
const PINS = { SR_CLK: 2, SR_DATA: 3, SR_LATCH: 4, SR_OE_N: 5, TRIG: 6, REF40: 20, SLOTCLK: 21, REF10: 22, CC1: 26, CC2: 27 };
if (fs.existsSync(MAIN_C)) {
  for (const m of fs.readFileSync(MAIN_C, 'utf8').matchAll(/^#define\s+PIN_(\w+)\s+(\d+)/gm)) PINS[m[1]] = Number(m[2]);
}
const GPCK = 8; // GPIO function select: clock input/output
const GPOUT_PIN = [21, 23, 24, 25]; // silicon: GPOUT0..3
const GPIN_PIN = [20, 22]; // silicon: GPIN0, GPIN1

// ---- firmware image ------------------------------------------------------------------------------
function loadImage(mcu, file) {
  const FLASH = 0x10000000;
  const buf = fs.readFileSync(file);
  if (file.endsWith('.uf2')) {
    let blocks = 0;
    for (let off = 0; off + 512 <= buf.length; off += 512) {
      const b = buf.subarray(off, off + 512);
      if (b.readUInt32LE(0) !== 0x0a324655 || b.readUInt32LE(4) !== 0x9e5d5157 || b.readUInt32LE(508) !== 0x0ab16f30) continue;
      const addr = b.readUInt32LE(12), n = b.readUInt32LE(16);
      if (addr >= FLASH && addr + n <= FLASH + mcu.flash.length) mcu.flash.set(b.subarray(32, 32 + n), addr - FLASH);
      blocks++;
    }
    if (!blocks) throw new Error(`no UF2 blocks in ${file}`);
  } else if (file.endsWith('.hex')) {
    let hi = 0;
    for (const line of buf.toString('utf8').split(/\r?\n/)) {
      if (line[0] !== ':') continue;
      const n = parseInt(line.substr(1, 2), 16), off = parseInt(line.substr(3, 4), 16), type = parseInt(line.substr(7, 2), 16);
      if (type === 4) hi = parseInt(line.substr(9, 4), 16);
      if (type !== 0) continue;
      const addr = ((hi << 16) | off) >>> 0;
      for (let i = 0; i < n; i++) mcu.flash[addr - FLASH + i] = parseInt(line.substr(9 + 2 * i, 2), 16);
    }
  } else if (file.endsWith('.bin')) {
    mcu.flash.set(buf, 0);
  } else throw new Error(`unsupported image type: ${file}`);
}

// ---- PIO in emulated time --------------------------------------------------------------------------
// One tick() per clk_sys cycle. A state machine executes an instruction, then sits out its delay
// cycles; a stalled instruction (WAIT, blocking PULL) is retried every cycle. rp2040js resolves
// waits on GPIO, IRQ flags and FIFOs at the moment the condition changes, so the next instruction
// runs on the following tick. Timing is exact to within a couple of cycles (no input synchroniser).
class PioClock {
  constructor(pio) {
    this.pio = pio;
    this.busy = [0, 0, 0, 0];
    this.last = [0, 0, 0, 0];
    this.acc = [0, 0, 0, 0];
    for (const sm of pio.machines) {
      const i = sm.index;
      const restart = sm.restart.bind(sm);
      sm.restart = () => {
        restart();
        this.busy[i] = 0;
        this.last[i] = sm.cycles;
      };
      // SHIFTCTRL: FJOIN_TX (bit 30) / FJOIN_RX (bit 31). Changing either clears both FIFOs
      // (pio_sm_clear_fifos relies on it); a joined FIFO is 8 deep.
      const write = sm.writeUint32.bind(sm);
      sm.writeUint32 = (offset, value) => {
        const before = sm.shiftCtrl;
        write(offset, value);
        if (offset === 0x08 && (before ^ sm.shiftCtrl) & 0xc0000000) {
          const tx = (sm.shiftCtrl >>> 30) & 1, rx = sm.shiftCtrl >>> 31;
          sm.txFIFO = new FIFO(tx && !rx ? 8 : rx && !tx ? 0 : 4);
          sm.rxFIFO = new FIFO(rx && !tx ? 8 : tx && !rx ? 0 : 4);
          sm.updateDMATx();
          sm.updateDMARx();
          pio.checkInterrupts();
        }
      };
    }
  }

  /** No enabled state machine can do anything until something outside the PIO changes. */
  get idle() {
    for (const sm of this.pio.machines) {
      if (sm.enabled && (!sm.waiting || this.busy[sm.index] > 0 || sm.cycles !== this.last[sm.index])) return false;
    }
    return true;
  }

  tick() {
    const { machines } = this.pio;
    for (let i = 0; i < 4; i++) {
      const sm = machines[i];
      if (!sm.enabled) continue;
      if (sm.clockDivInt !== 1 || sm.clockDivFrac !== 0) {
        const div = (sm.clockDivInt || 0x10000) * 256 + sm.clockDivFrac;
        this.acc[i] += 256;
        if (this.acc[i] < div) continue;
        this.acc[i] -= div;
      }
      if (this.busy[i] > 0) {
        this.busy[i]--;
        continue;
      }
      // A wait released from outside (a GPIO edge, an IRQ flag) adds its delay cycles there.
      const owed = sm.cycles - this.last[i];
      if (owed > 0) {
        this.last[i] = sm.cycles;
        this.busy[i] = owed - 1;
        continue;
      }
      const before = sm.cycles;
      sm.step();
      const used = sm.cycles - before;
      this.last[i] = sm.cycles;
      if (used > 1) this.busy[i] = used - 1;
    }
    this.pio.checkChangedPins();
  }
}

// ---- 74HC595 chain, as seen on SR_CLK / SR_DATA / SR_LATCH / SR_OE_N -------------------------------
// Bit i of the control word is output Q(i % 8) of register i / 8; the first bit shifted in ends
// in the last register (see tests/test_firmware_core.py).
class Chain {
  constructor(bits = 56) {
    this.mask = (1n << BigInt(bits)) - 1n;
    this.bits = bits;
    this.sr = 0n;
    this.clocks = 0; // SR_CLK rising edges since the last latch
  }
  shift(bit) {
    this.sr = ((this.sr << 1n) | (bit ? 1n : 0n)) & this.mask;
    this.clocks++;
  }
  latch() {
    const c = this.clocks;
    this.clocks = 0;
    return { word: this.sr.toString(16).padStart(this.bits / 4, '0'), clocks: c };
  }
}

// ---- the emulated board --------------------------------------------------------------------------
const DEFAULT_ENV = {
  ref40_hz: 40e6, // on-board reference on GPIN0 (GPIO20); 0 = absent
  ref10_hz: 10e6, // HackRF 10 MHz on GPIN1 (GPIO22); 0 = absent
  cc1_mv: 1600, // USB-C CC1 across the 5.1k pull-down: 1.6 V = 3 A source
  cc2_mv: 0,
  // What FC0 does with a clock that never starts. The RP2040 documents FC0_STATUS.WAITING ("waiting
  // for test clock to start") but not whether the count ever finishes. 'done': DONE with 0 kHz and
  // DIED after the interval. 'wait': WAITING until the next measurement is started.
  fc_dead_clock: 'done',
};

class Board {
  constructor(scenario, emitEv, trace) {
    this.emit = emitEv;
    this.trace = trace;
    this.env = { ...DEFAULT_ENV, ...(scenario.env ?? {}) };
    this.warnings = new Map();
    const note = (kind) => (component, message) => {
      const key = `${kind} ${component}: ${message.replace(/0x-?[0-9a-f]+/g, '0x…')}`;
      this.warnings.set(key, (this.warnings.get(key) ?? 0) + 1);
    };
    this.clock = new SimulationClock();
    const mcu = (this.mcu = new RP2040(this.clock));
    mcu.logger = { debug() {}, info() {}, warn: note('warn'), error: note('error') };
    mcu.loadBootrom(bootromB1);
    loadImage(mcu, IMAGE);
    mcu.core.PC = 0x10000000; // boot2 at the start of flash, as rp2040js's demo does

    for (const pio of mcu.pio) pio.run = () => {}; // never on wall-clock timers
    this.pioClocks = mcu.pio.map((pio) => new PioClock(pio));
    this.installClockModels();
    this.installRecorders();
    this.host = new Host(this);
    this.applyEnv();
  }

  get t() {
    return Math.round(this.clock.nanos);
  }

  setEnv(change) {
    Object.assign(this.env, change);
    this.applyEnv();
    this.emit({ ev: 'env', t_ns: this.t, env: { ...this.env } });
  }

  applyEnv() {
    const adc = (mv) => Math.max(0, Math.min(4095, Math.round((mv * 4095) / 3300)));
    this.mcu.adc.channelValues[PINS.CC1 - 26] = adc(this.env.cc1_mv);
    this.mcu.adc.channelValues[PINS.CC2 - 26] = adc(this.env.cc2_mv);
    for (const g of this.gpout) this.retimeGpout(g);
  }

  // -- clock generators and the frequency counter (CLOCKS block) -------------------------------
  /** Frequency of a GPIN input: the scenario's reference, if its pin is set to the clock function. */
  gpinHz(n) {
    const hz = n === 0 ? this.env.ref40_hz : this.env.ref10_hz;
    return this.mcu.gpio[GPIN_PIN[n]].functionSelect === GPCK ? hz : 0;
  }

  installClockModels() {
    const { mcu } = this;
    const clocks = mcu.clocks;
    // GPOUTn generators: CTRL at 0x0c * n, DIV at 0x0c * n + 4. ENABLE bit 11, KILL bit 10,
    // AUXSRC bits 8:5. DIV: INT bits 31:8, FRAC bits 7:0.
    const AUX = [
      () => mcu.pllSys.frequency, () => this.gpinHz(0), () => this.gpinHz(1), () => mcu.pllUsb.frequency,
      () => mcu.roscFreq, () => mcu.xoscFreq, () => mcu.clkSys, () => 48e6, () => 48e6, () => 46875, () => clocks.refFreq,
    ];
    const AUX_NAME = ['pll_sys', 'gpin0', 'gpin1', 'pll_usb', 'rosc', 'xosc', 'clk_sys', 'clk_usb', 'clk_adc', 'clk_rtc', 'clk_ref'];
    this.gpout = [0, 1, 2, 3].map((n) => ({ n, pin: GPOUT_PIN[n], level: false, hz: 0, alarm: null, desc: '' }));
    for (const g of this.gpout) {
      g.alarm = this.clock.createAlarm(() => {
        g.level = !g.level;
        this.mcu.gpio[g.pin].setInputValue(g.level && this.mcu.gpio[g.pin].functionSelect === GPCK);
        if (g.hz > 0) g.alarm.schedule(1e9 / g.hz / 2);
      });
    }
    this.gpoutConfig = (g) => {
      const ctrl = clocks[`gpout${g.n}Ctrl`], div = clocks[`gpout${g.n}Div`];
      const aux = (ctrl >>> 5) & 0xf;
      const d = ((div >>> 8) || 0x1000000) + (div & 0xff) / 256;
      const on = (ctrl & (1 << 11)) && !(ctrl & (1 << 10));
      const src = AUX[aux] ? AUX[aux]() : 0;
      return { hz: on ? src / d : 0, src: AUX_NAME[aux] ?? `aux${aux}`, div: d, srcHz: src };
    };
    this.retimeGpout = (g) => {
      const c = this.gpoutConfig(g);
      const desc = JSON.stringify(c);
      if (desc === g.desc) return;
      g.desc = desc;
      g.hz = c.hz;
      g.alarm.cancel();
      if (c.hz > 0) g.alarm.schedule(1e9 / c.hz / 2); // keeps its level, new half period from now
      if (g.n === 0 || c.hz > 0) this.emit({ ev: 'gpout', t_ns: this.t, n: g.n, pin: g.pin, ...c });
    };
    // A GPIN only carries the reference once its pin is set to the clock function.
    for (const p of GPIN_PIN) {
      const pin = mcu.gpio[p];
      const check = pin.checkForUpdates.bind(pin);
      pin.checkForUpdates = () => {
        check();
        for (const g of this.gpout) this.retimeGpout(g);
      };
    }

    // FC0: REF_KHZ 0x80, MIN_KHZ 0x84, MAX_KHZ 0x88, DELAY 0x8c, INTERVAL 0x90, SRC 0x94,
    // STATUS 0x98, RESULT 0x9c (KHZ bits 29:5, FRAC 4:0). Writing a non-zero SRC starts a
    // measurement lasting 0.98 us * 2^INTERVAL. A clock that stops during it reads 0 kHz and DIED;
    // one that is not running at the start: see env.fc_dead_clock.
    const FC_SRC = [
      () => 0, () => mcu.pllSys.frequency, () => mcu.pllUsb.frequency, () => mcu.roscFreq, () => mcu.roscFreq,
      () => mcu.xoscFreq, () => this.gpinHz(0), () => this.gpinHz(1), () => clocks.refFreq, () => mcu.clkSys,
      () => mcu.clkPeri, () => 48e6, () => 48e6, () => 46875,
    ];
    const FC_NAME = ['null', 'pll_sys', 'pll_usb', 'rosc', 'rosc_ph', 'xosc', 'gpin0', 'gpin1', 'clk_ref', 'clk_sys', 'clk_peri', 'clk_usb', 'clk_adc', 'clk_rtc'];
    const fc = { ref: 0, min: 0, max: 0xffffffff, delay: 1, interval: 8, src: 0, done: 0, hz: 0, final: true };
    const finish = () => {
      if (fc.final || this.clock.nanos < fc.done) return;
      fc.final = true;
      if ((FC_SRC[fc.src] ? FC_SRC[fc.src]() : 0) === 0) fc.hz = 0; // stopped during the count
    };
    const read = clocks.readUint32.bind(clocks), write = clocks.writeUint32.bind(clocks);
    clocks.readUint32 = (offset) => {
      if (offset < 0x80 || offset > 0x9c) return read(offset);
      switch (offset) {
        case 0x80: return fc.ref;
        case 0x84: return fc.min;
        case 0x88: return fc.max;
        case 0x8c: return fc.delay;
        case 0x90: return fc.interval;
        case 0x94: return fc.src;
        case 0x98: {
          if (!fc.src) return 0;
          if (fc.hz === 0 && this.env.fc_dead_clock === 'wait') return 1 << 12; // WAITING
          if (this.clock.nanos < fc.done) return 1 << 8; // RUNNING
          finish();
          const khz = fc.hz / 1000;
          let s = 1 << 4; // DONE
          if (fc.hz === 0) s |= 1 << 28; // DIED
          if (khz < fc.min) s |= (1 << 20) | (1 << 16); // SLOW, FAIL
          else if (khz > fc.max) s |= (1 << 24) | (1 << 16); // FAST, FAIL
          else s |= 1; // PASS
          return s;
        }
        case 0x9c:
          finish();
          return fc.final ? Math.round((fc.hz / 1000) * 32) & 0x3fffffff : 0;
      }
      return 0;
    };
    clocks.writeUint32 = (offset, value) => {
      if (offset >= 0x80 && offset <= 0x9c) {
        switch (offset) {
          case 0x80: fc.ref = value & 0xfffff; break;
          case 0x84: fc.min = value & 0x1ffffff; break;
          case 0x88: fc.max = value & 0x1ffffff; break;
          case 0x8c: fc.delay = value & 7; break;
          case 0x90: fc.interval = value & 0xf; break;
          case 0x94:
            fc.src = value & 0xff;
            if (fc.src) {
              fc.hz = FC_SRC[fc.src] ? FC_SRC[fc.src]() : 0;
              fc.done = this.clock.nanos + 980 * 2 ** fc.interval;
              fc.final = false;
              this.emit({ ev: 'fc', t_ns: this.t, src: FC_NAME[fc.src] ?? `src${fc.src}`, khz: fc.hz / 1000 });
            }
            break;
        }
        return;
      }
      write(offset, value);
      if (offset <= 0x2c) for (const g of this.gpout) this.retimeGpout(g);
    };
  }

  // -- pins --------------------------------------------------------------------------------------------
  installRecorders() {
    const { mcu } = this;
    const level = (s) => (s === GPIOPinState.High ? 1 : s === GPIOPinState.Low ? 0 : null); // null = not driven
    this.chain = new Chain(56);
    const st = (this.timing = {
      words: 0, min_data_setup_ns: Infinity, min_clk_high_ns: Infinity, min_clk_low_ns: Infinity,
      min_latch_high_ns: Infinity, min_clk_to_latch_ns: Infinity, data_clk_same_instant: 0, clk_latch_same_instant: 0,
    });
    const last = { data: -Infinity, dataBefore: null, clkRise: -Infinity, clkFall: -Infinity, latchRise: -Infinity };
    const lv = (this.levels = {});
    // The board pulls SR_OE_N up: undriven means outputs off.
    const oeLevel = () => lv.SR_OE_N ?? 1;
    const watch = (name, fn) => {
      const pin = mcu.gpio[PINS[name]];
      lv[name] = level(pin.value);
      pin.addListener((state) => {
        const v = level(state), t = this.clock.nanos;
        const old = lv[name];
        lv[name] = v;
        if (this.trace || name === 'SR_OE_N' || name === 'TRIG') this.emit({ ev: 'pin', t_ns: Math.round(t), pin: name, level: v });
        fn?.(v, old, t);
      });
    };
    watch('SR_DATA', (v, old, t) => {
      if (last.data !== t) last.dataBefore = old;
      last.data = t;
    });
    watch('SR_CLK', (v, old, t) => {
      if (v === 1 && old !== 1) {
        // If SR_DATA changed in this same instant (a hold-time violation on a real chip), count it
        // and shift the level it had before.
        let data = lv.SR_DATA;
        if (last.data === t) {
          data = last.dataBefore;
          st.data_clk_same_instant++;
        } else if (level(mcu.gpio[PINS.SR_DATA].value) !== data) st.data_clk_same_instant++; // listener not run yet
        st.min_data_setup_ns = Math.min(st.min_data_setup_ns, t - last.data);
        if (last.clkFall > -Infinity) st.min_clk_low_ns = Math.min(st.min_clk_low_ns, t - last.clkFall);
        this.chain.shift(data === 1);
        last.clkRise = t;
      } else if (v !== 1 && old === 1) {
        st.min_clk_high_ns = Math.min(st.min_clk_high_ns, t - last.clkRise);
        last.clkFall = t;
      }
    });
    watch('SR_LATCH', (v, old, t) => {
      if (v === 1 && old !== 1) {
        if (level(mcu.gpio[PINS.SR_CLK].value) !== lv.SR_CLK || last.clkRise === t) st.clk_latch_same_instant++;
        st.min_clk_to_latch_ns = Math.min(st.min_clk_to_latch_ns, t - last.clkRise);
        const { word, clocks } = this.chain.latch();
        st.words++;
        last.latchRise = t;
        this.emit({ ev: 'latch', t_ns: Math.round(t), word, clocks, oe_n: oeLevel() });
      } else if (v !== 1 && old === 1) {
        st.min_latch_high_ns = Math.min(st.min_latch_high_ns, t - last.latchRise);
      }
    });
    watch('SR_OE_N');
    watch('TRIG');
  }

  // -- running ---------------------------------------------------------------------------------------------
  pioIdle() {
    return this.pioClocks[0].idle && this.pioClocks[1].idle;
  }

  /** Run until emulated time `endNs`, or until `until()` holds (checked between instructions). */
  runTo(endNs, until = null) {
    const { clock, mcu } = this;
    const core = mcu.core;
    const [p0, p1] = this.pioClocks;
    while (clock.nanos < endNs) {
      if (until && until()) return true;
      const cyc = 1e9 / mcu.clkSys;
      let n;
      if (core.waiting) {
        // WFI/WFE: nothing runs on the core until the next event; jump there unless the PIO is busy.
        const dt = clock.nanosToNextAlarm;
        n = this.pioIdle() ? Math.max(1, Math.min(Math.ceil((dt > 0 ? dt : 1000) / cyc), Math.ceil((endNs - clock.nanos) / cyc))) : 1;
      } else {
        n = core.executeInstruction();
      }
      // An alarm inside the window (a GPOUT edge, a DMA transfer) can wake a state machine: then
      // go cycle by cycle so that it reacts in the right cycle.
      const next = clock.nextAlarm;
      if (this.pioIdle() && !(next && next.nanos <= clock.nanos + n * cyc)) {
        clock.tick(n * cyc);
        p0.pio.checkChangedPins(); // pins written by the CPU through SMx_INSTR
        p1.pio.checkChangedPins();
      } else {
        for (let k = 0; k < n; k++) {
          clock.tick(cyc);
          p0.tick();
          p1.tick();
        }
      }
    }
    return until ? !!until() : true;
  }

  run(ns) {
    return this.runTo(this.clock.nanos + ns);
  }

  runUntil(pred, maxNs) {
    return this.runTo(this.clock.nanos + maxNs, pred);
  }
}

// ---- USB CDC host --------------------------------------------------------------------------------
// rp2040js's USBCDC enumerates the device and sets DTR. Its OUT side is replaced: an armed OUT
// endpoint waits until the script has a packet (as with a real host, which only sends when it has
// data), and the script decides where packets split.
class Host {
  constructor(board) {
    this.board = board;
    const usb = (this.usb = board.mcu.usbCtrl);
    this.cdc = new USBCDC(usb);
    this.connected = false;
    this.queue = []; // { step, data: Uint8Array }
    this.arms = []; // OUT buffer sizes the device has made available
    this.inflight = false; // rp2040js delivers one OUT buffer per read alarm: one packet at a time
    this.delivered = 0; // packets the device has received
    this.rx = [];
    this.rxTimes = [];
    this.replies = [];
    this.cdc.onDeviceConnected = () => {
      this.connected = true;
      board.emit({ ev: 'connected', t_ns: board.t });
    };
    this.cdc.onSerialData = (buf) => this.onRx(buf);
    const enumerate = usb.onEndpointRead;
    usb.onEndpointRead = (ep, size) => {
      if (ep !== this.cdc.outEndpoint) return enumerate(ep, size);
      this.arms.push(size);
      this.pump();
    };
  }

  send(step, data) {
    this.queue.push({ step, data });
    this.pump();
  }

  pump() {
    if (this.inflight || !this.arms.length || !this.queue.length) return;
    const size = this.arms.shift();
    const head = this.queue[0];
    let pkt = head.data;
    if (pkt.length > size) {
      head.data = pkt.subarray(size);
      pkt = pkt.subarray(0, size);
    } else this.queue.shift();
    const delay = this.usb.readDelayMicroseconds * 1000;
    this.inflight = true;
    this.usb.endpointReadDone(this.cdc.outEndpoint, pkt);
    // The device sees the packet when rp2040js's read alarm fires, `delay` from now; this alarm
    // fires in the same instant.
    const alarm = this.board.clock.createAlarm(() => {
      this.inflight = false;
      this.delivered++;
      this.board.emit({ ev: 'out', t_ns: this.board.t, step: head.step, hex: Buffer.from(pkt).toString('hex') });
      this.pump();
    });
    alarm.schedule(delay);
  }

  get idle() {
    return this.queue.length === 0 && !this.inflight;
  }

  onRx(buf) {
    const t = this.board.t;
    if (this.board.trace) this.board.emit({ ev: 'rx', t_ns: t, hex: Buffer.from(buf).toString('hex') });
    for (const b of buf) {
      this.rx.push(b);
      this.rxTimes.push(t);
    }
    // Reply frames: 0x5A cmd status len16 payload crc16. t_first_ns: when the first byte came in
    // (the firmware had handled the request by then).
    for (;;) {
      while (this.rx.length && this.rx[0] !== 0x5a) {
        this.rxTimes.shift();
        this.board.emit({ ev: 'junk', t_ns: t, byte: this.rx.shift() });
      }
      if (this.rx.length < 5) break;
      const n = 7 + (this.rx[3] | (this.rx[4] << 8));
      if (this.rx.length < n) break;
      const frame = this.rx.splice(0, n);
      const first = this.rxTimes.splice(0, n)[0];
      this.replies.push(frame);
      this.board.emit({ ev: 'reply', t_ns: t, t_first_ns: first, hex: Buffer.from(frame).toString('hex') });
    }
  }
}

// ---- scenarios -------------------------------------------------------------------------------------
function hexBytes(h) {
  if (!/^([0-9a-fA-F]{2})*$/.test(h)) throw new Error(`bad hex: ${h}`);
  return Uint8Array.from(Buffer.from(h, 'hex'));
}

function runScenario(sc) {
  const wall0 = Date.now();
  const name = sc.name ?? 'scenario';
  const ev = (o) => emit({ scenario: name, ...o });
  const board = new Board(sc, ev, args.trace || !!sc.trace);
  ev({ ev: 'start', t_ns: 0, image: IMAGE, pins: PINS, env: { ...board.env } });
  const host = board.host;
  if (!board.runUntil(() => host.connected, 300e6)) {
    ev({ ev: 'error', t_ns: board.t, message: 'USB CDC never connected', pc: board.mcu.core.PC.toString(16) });
    return false;
  }
  board.run((sc.settle_us ?? 2000) * 1000);
  let ok = true;
  (sc.steps ?? []).forEach((step, i) => {
    if (!ok) return;
    ev({ ev: 'step', t_ns: board.t, i, ...step });
    if (step.env) board.setEnv(step.env);
    if (step.dtr !== undefined) host.cdc.cdcSetControlLineState(step.dtr ? 3 : 2); // SET_CONTROL_LINE_STATE, RTS kept
    if (step.send !== undefined) {
      const data = hexBytes(step.send);
      const sizes = step.packets ?? [data.length];
      if (sizes.reduce((a, b) => a + b, 0) !== data.length) throw new Error(`step ${i}: packets do not add up`);
      let off = 0;
      sizes.forEach((n, j) => {
        const want = host.delivered + 1;
        host.send(i, data.subarray(off, off + n));
        off += n;
        if (step.gap_us !== undefined && j < sizes.length - 1) {
          // wait until the device has this packet (all of it, if it needed several USB packets)
          if (!board.runUntil(() => host.idle && host.delivered >= want, 50e6)) {
            ev({ ev: 'error', t_ns: board.t, message: `step ${i}: packet ${j} not taken by the device` });
            ok = false;
            return;
          }
          board.run(step.gap_us * 1000);
        }
      });
    }
    if (step.wait_us !== undefined) board.run(step.wait_us * 1000);
    if (step.wait_replies !== undefined) {
      const want = step.wait_replies;
      if (!board.runUntil(() => host.replies.length >= want, (step.timeout_us ?? 50000) * 1000)) {
        ev({ ev: 'timeout', t_ns: board.t, i, replies: host.replies.length, want });
      }
    }
  });
  const st = { ...board.timing };
  for (const k of Object.keys(st)) if (st[k] === Infinity) st[k] = null;
  ev({ ev: 'timing', t_ns: board.t, clk_sys_hz: board.mcu.clkSys, ...st });
  ev({ ev: 'warnings', t_ns: board.t, items: [...board.warnings].map(([k, n]) => ({ what: k, count: n })) });
  ev({ ev: 'end', t_ns: board.t, wall_ms: Date.now() - wall0, cpu_cycles: board.mcu.core.cycles });
  return ok;
}

// ---- built-in demo (frames as array16/protocol.py builds them) ----------------------------------
function crc16(bytes) {
  let crc = 0xffff;
  for (const b of bytes) {
    crc ^= b << 8;
    for (let i = 0; i < 8; i++) crc = crc & 0x8000 ? ((crc << 1) ^ 0x1021) & 0xffff : (crc << 1) & 0xffff;
  }
  return crc;
}

function request(cmd, payload = []) {
  const body = [0xa5, cmd, payload.length & 0xff, payload.length >> 8, ...payload];
  const crc = crc16(body);
  return Buffer.from([...body, crc & 0xff, crc >> 8]).toString('hex');
}

function word(v) {
  const out = [];
  for (let i = 0; i < 7; i++) out.push(Number((BigInt(v) >> BigInt(8 * i)) & 0xffn));
  return out;
}

function demoScript() {
  const PING = 0x01, SET_WORD = 0x10, LOAD_SEQ = 0x11, RUN = 0x12, STOP = 0x13, STATUS = 0x20;
  const seq = [0x0000000000ffffn, 0x00ffff0000ffffn, 0xff00000000ffffn].map(word);
  const load = [50, 0, 2, 0, seq.length & 0xff, seq.length >> 8, ...seq.flat()]; // slot 50 us, repeat 2, count u16
  return {
    scenarios: [
      {
        name: 'demo',
        steps: [
          { send: request(PING), wait_replies: 1 },
          { send: request(STATUS), wait_replies: 2 },
          { send: request(SET_WORD, word(0x80000000000001n)), wait_replies: 3 },
          { send: request(LOAD_SEQ, load), wait_replies: 4 },
          { send: request(RUN), wait_replies: 5 },
          { wait_us: 1000 },
          { send: request(STATUS), wait_replies: 6 },
          { send: request(STOP), wait_replies: 7 },
        ],
      },
    ],
  };
}

// ---- main ------------------------------------------------------------------------------------------
let script;
if (args.demo) script = demoScript();
else if (args.script) script = JSON.parse(args.script === '-' ? fs.readFileSync(0, 'utf8') : fs.readFileSync(args.script, 'utf8'));
else fail('give --script FILE (or - for stdin) or --demo');
let allOk = true;
for (const sc of script.scenarios ?? [script]) allOk = runScenario(sc) && allOk;
finish(allOk ? 0 : 1);
