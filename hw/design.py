"""The board's circuit, sheet by sheet. Run:  python3 hw/design.py

Writes the KiCad project to hw/kicad/ and a summary to hw/bom_summary.json.
First-pass schematics: every pinout marked "verify" below still needs a second
pair of eyes against its datasheet before layout.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import kicadgen as K                                   # noqa: E402
from array16 import board as B                         # noqa: E402

NC = K.NC
FP_R, FP_C, FP_L = "Resistor_SMD:R_0402_1005Metric", "Capacitor_SMD:C_0402_1005Metric", "Inductor_SMD:L_0402_1005Metric"
FP_UFL, FP_SMA = "Connector_Coaxial:U.FL_Hirose_U.FL-R-SMT-1_Vertical", "Connector_Coaxial:SMA_Samtec_SMA-J-P-H-ST-EM1_EdgeMount"

R, C, L = K.library("Device:R"), K.library("Device:C"), K.library("Device:L")
FB, XTAL, COAX = K.library("Device:FerriteBead"), K.library("Device:Crystal_GND24"), K.library("Connector:Conn_Coaxial")

# Pinouts from the datasheets (Infineon BGS12WN6 table 11, Skyworks SKY65404-31 table 1).
SW = K.custom("BGS12WN6", "U", [("5", "RFIN", "L"), ("6", "CTRL", "L"), ("3", "RF1", "R"), ("1", "RF2", "R"),
                                 ("4", "VDD", "T"), ("2", "GND", "B")], footprint="array16:TSNP-6-2")
LNA = K.custom("SKY65404-31", "U", [("3", "RF_IN", "L"), ("1", "V_ENABLE", "L"), ("4", "RF_OUT", "R"),
                                     ("6", "VCC", "T"), ("2", "GND", "B"), ("5", "GND", "B")], footprint="array16:QFN-6_1.5x1.5")
TL = K.custom("TLINE_70R7_QW", "TL", [("1", "A", "L"), ("2", "B", "R")], width=10.16, footprint="array16:TL_70R7_7.5mm")
PATCH = K.custom("PATCH", "AE", [("1", "FEED", "L")], width=10.16, footprint="array16:Patch_5G8_15.7x11.9")
PROBE = K.custom("PROBE_PAD", "AE", [("1", "FEED", "L")], width=10.16, footprint="array16:Probe_Pad_3mm")
TCXO = K.custom("TCXO_40MHz", "Y", [("1", "EN", "L"), ("4", "VDD", "T"), ("3", "OUT", "R"), ("2", "GND", "B")],
                footprint="Oscillator:Oscillator_SMD_Abracon_ASE-4Pin_3.2x2.5mm")
# verify: pinout taken from memory of the TI LMK1C1104 / CDCLVC1104 (TSSOP-8)
CLKBUF = K.custom("LMK1C1104", "U", [("1", "CLKIN", "L"), ("2", "1G", "L"), ("3", "Y0", "R"), ("5", "Y1", "R"), ("7", "Y2", "R"),
                                      ("8", "Y3", "R"), ("6", "VDD", "T"), ("4", "GND", "B")], footprint="Package_SO:TSSOP-8_4.4x3mm_P0.65mm")
CUSTOM = [SW, LNA, TL, PATCH, PROBE, TCXO, CLKBUF]


def r(sh, value, a, b):
    sh.add(R, "R", value, {"1": a, "2": b}, FP_R)


def c(sh, value, a, b="GND"):
    sh.add(C, "C", value, {"1": a, "2": b}, FP_C)


def ctl(sh, name):
    """RC filter where a control line enters an RF block; returns the filtered net."""
    r(sh, "1k", name, name + "_F")
    c(sh, "100p", name + "_F")
    return name + "_F"


def switch(sh, rfin, rf1, rf2, ctrl, decouple=True):
    """SPDT: CTRL low connects RFIN-RF1, CTRL high connects RFIN-RF2."""
    sh.add(SW, "U", "BGS12WN6", {"RFIN": rfin, "RF1": rf1, "RF2": rf2, "CTRL": ctrl, "VDD": "+3V3_RF", "GND": "GND"}, "array16:TSNP-6-2")
    if decouple:
        c(sh, "100p", "+3V3_RF")


def lna(sh, tag, rf_in, rf_out, enable):
    """SKY65404-31 with DC blocks and the datasheet's input match (figure 10); tune on the coupon."""
    c(sh, "5.6p", rf_in, tag + "_M")
    c(sh, "0.5p", tag + "_M")
    sh.add(L, "L", "0.6n", {"1": tag + "_M", "2": tag + "_IN"}, FP_L)
    sh.add(LNA, "U", "SKY65404-31", {"RF_IN": tag + "_IN", "RF_OUT": tag + "_OUT", "V_ENABLE": tag + "_EN", "VCC": "+3V3_RF", "GND": "GND"},
           "array16:QFN-6_1.5x1.5")
    r(sh, "100", enable, tag + "_EN")
    c(sh, "100n", "+3V3_RF")
    c(sh, "10p", "+3V3_RF")
    c(sh, "5.6p", tag + "_OUT", rf_out)


def wilkinson(sh, a, b, out):
    """Printed 2-way divider: two quarter-wave 70.7-ohm arms and a 100-ohm isolation resistor."""
    sh.add(TL, "TL", "70.7R quarter-wave", {"A": a, "B": out}, "array16:TL_70R7_7.5mm")
    sh.add(TL, "TL", "70.7R quarter-wave", {"A": b, "B": out}, "array16:TL_70R7_7.5mm")
    r(sh, "100", a, b)


def coax(sh, value, net, fp=FP_UFL):
    sh.add(COAX, "J", value, {"In": net, "Ext": "GND"}, fp)


# ---- element cell (x16) ----------------------------------------------------------------------
cell = K.Sheet("cell", "Element cell: T/R bypass, LNA, 2-bit phase shifter", ["RF", "RX", "P180", "P90"], "A3")
cell.add(PATCH, "AE", "Patch 15.7 x 11.9 mm", {"FEED": "ANT"}, "array16:Patch_5G8_15.7x11.9")
rx, p180, p90 = ctl(cell, "RX"), ctl(cell, "P180"), ctl(cell, "P90")
switch(cell, "ANT", "BYP", "LNA_A", rx)                       # RX high: through the LNA; low: bypass for transmit
lna(cell, "LNA", "LNA_A", "LNA_B", rx)
switch(cell, "PS_IN", "BYP", "LNA_B", rx, decouple=False)
switch(cell, "PS_IN", "S180_REF", "S180_LONG", p180)
switch(cell, "PS_MID", "S180_REF", "S180_LONG", p180, decouple=False)
switch(cell, "PS_MID", "S90_REF", "S90_LONG", p90)
switch(cell, "RF", "S90_REF", "S90_LONG", p90, decouple=False)
cell.notes = ["Net S180_LONG is 14.48 mm longer than S180_REF (180 deg at 5.8 GHz); S90_LONG is 7.24 mm longer than S90_REF (90 deg). 50-ohm traces are 0.37 mm wide on L4.",
              "BYP is the transmit bypass around the LNA; keep it the same electrical length in all 16 cells. No DC is allowed on the switch RF pins, hence the 5.6 pF blocks.",
              "RX high = receive through the LNA (and LNA enabled). All control lines low = safe power-up state."]

# ---- quadrant (x4) -----------------------------------------------------------------------------
qports = [f"{n}{i}" for n in ("RX", "P180_", "P90_") for i in range(4)] + ["MODE", "QRX", "LNA2_EN", "SHARED"]
quad = K.Sheet("quadrant", "Quadrant: four cells, 4-way combiner, mode switch, radio port", qports, "A2")
for i in range(4):
    quad.child(f"cell{i}", cell, {"RF": f"E{i}", "RX": f"RX{i}", "P180": f"P180_{i}", "P90": f"P90_{i}"})
wilkinson(quad, "E0", "E1", "P01")
wilkinson(quad, "E2", "E3", "P23")
wilkinson(quad, "P01", "P23", "QSUM")
mode, qrx = ctl(quad, "MODE"), ctl(quad, "QRX")
switch(quad, "QSUM", "QPORT", "SHARED", mode)                  # MODE low: own radio channel; high: shared combiner
switch(quad, "QPORT", "TX_SKT", "LNA_A", qrx, decouple=False)  # QRX high: receive through the port LNA
lna(quad, "LNA", "LNA_A", "RX_SKT", "LNA2_EN")
coax(quad, "RX U.FL", "RX_SKT")
coax(quad, "TX U.FL", "TX_SKT")
quad.notes = ["The four cell outputs E0..E3 must be equal in length to the combiner. All four quadrants are mirror images of each other."]

# ---- shared combiner and HackRF port -------------------------------------------------------------
hack = K.Sheet("hackrf_port", "Shared combiner and HackRF port", ["S0", "S1", "S2", "S3", "HRX", "LNA2_EN"])
wilkinson(hack, "S0", "S1", "A01")
wilkinson(hack, "S2", "S3", "A23")
wilkinson(hack, "A01", "A23", "HSUM")
hrx = ctl(hack, "HRX")
switch(hack, "HSUM", "BYP", "LNA_A", hrx)
lna(hack, "LNA", "LNA_A", "LNA_B", "LNA2_EN")
switch(hack, "HRF", "BYP", "LNA_B", hrx, decouple=False)
coax(hack, "HackRF SMA", "HRF", FP_SMA)
hack.notes = ["S0..S3 come from the four mode switches along equal-length lines. This block sits at the centre of the array; HRF runs 60 mm down to the edge connector."]

# ---- auxiliary patch and calibration probe ----------------------------------------------------------
aux = K.Sheet("aux_cal", "Auxiliary patch and calibration probe", ["CAL"])
aux.add(PATCH, "AE", "Auxiliary patch", {"FEED": "AUX"}, "array16:Patch_5G8_15.7x11.9")
coax(aux, "AUX U.FL", "AUX")
coax(aux, "CAL IN U.FL", "CAL_IN")
switch(aux, "CAL_IN", "CAL_TERM", "PROBE", ctl(aux, "CAL"))    # CAL low: terminated; high: to the centre probe
r(aux, "49.9", "CAL_TERM", "GND")
aux.add(PROBE, "AE", "Centre probe, 3 mm pad", {"FEED": "PROBE"}, "array16:Probe_Pad_3mm")
aux.notes = ["CAL IN takes a spare transmit channel. The host mutes that channel except while aligning the radio channels."]

# ---- 40 MHz reference ---------------------------------------------------------------------------------
clk = K.Sheet("clock", "40 MHz reference: one oscillator for both radios and the controller", ["CLK40"])
clk.add(FB, "FB", "600R@100MHz", {"1": "+3V3", "2": "VCLK"}, "Inductor_SMD:L_0603_1608Metric")
c(clk, "1u", "VCLK")
c(clk, "100n", "VCLK")
c(clk, "100n", "VCLK")
clk.add(TCXO, "Y", "40 MHz TCXO, 0.5 ppm", {"EN": "VCLK", "VDD": "VCLK", "OUT": "XO", "GND": "GND"}, "Oscillator:Oscillator_SMD_Abracon_ASE-4Pin_3.2x2.5mm")
clk.add(CLKBUF, "U", "LMK1C1104", {"CLKIN": "XO", "1G": "VCLK", "Y0": "Y0", "Y1": "Y1", "Y2": "Y2", "Y3": "Y3", "VDD": "VCLK", "GND": "GND"},
        "Package_SO:TSSOP-8_4.4x3mm_P0.65mm")
r(clk, "33", "Y0", "CLK40")
for i in (1, 2, 3):
    r(clk, "33", f"Y{i}", f"Y{i}_R")
    c(clk, "100n", f"Y{i}_R", f"CLK{i}")
    coax(clk, f"CLK{i} U.FL", f"CLK{i}")
clk.notes = ["CLK1 and CLK2 feed the two radio boards' external clock inputs; CLK3 is spare. Check the level each radio expects and pad if needed.",
             "verify: LMK1C1104 pinout and the TCXO's pin 1 function against their datasheets."]

# ---- beam controller --------------------------------------------------------------------------------------
bits = ([f"RX_{e}" for e in range(16)] + [f"P180_{e}" for e in range(16)] + [f"P90_{e}" for e in range(16)]
        + ["MODE"] + [f"QRX_{q}" for q in range(4)] + ["HRX", "CAL", "LNA2_EN"])
assert len(bits) == B.WORD_BITS and bits[B.BIT_MODE] == "MODE" and bits[B.BIT_QRX] == "QRX_0" and bits[B.BIT_LNA2] == "LNA2_EN"
cports = bits + ["CLK40", "USB_DP", "USB_DM", "CC1", "CC2"]
ctrl = K.Sheet("controller", "Beam controller: RP2040 and seven shift registers (56 control lines)", cports, "A1")
gpio = {"GPIO2": "SR_CLK", "GPIO3": "SR_DATA", "GPIO4": "SR_LATCH", "GPIO5": "SR_OE_N", "GPIO6": "TRIG", "GPIO7": "SYNC",
        "GPIO8": "EXP_IO0", "GPIO9": "EXP_IO1", "GPIO20": "CLK40", "GPIO22": "HCLK", "GPIO26_ADC0": "CC1", "GPIO27_ADC1": "CC2"}
rp = K.library("MCU_RaspberryPi:RP2040")
nets = {p["name"]: NC for p in rp.pins if p["name"].startswith("GPIO")}
nets.update(gpio)
nets.update({"IOVDD": "+3V3", "DVDD": "+1V1", "VREG_IN": "+3V3", "VREG_VOUT": "+1V1", "USB_VDD": "+3V3", "ADC_AVDD": "+3V3", "GND": "GND",
             "TESTEN": "GND", "XIN": "XIN", "XOUT": "XOUT", "RUN": "RUN", "SWCLK": "SWCLK", "SWD": "SWD", "USB_DP": "RP_DP", "USB_DM": "RP_DM",
             "QSPI_SS": "QSPI_SS", "QSPI_SCLK": "QSPI_SCLK", "QSPI_SD0": "QSPI_SD0", "QSPI_SD1": "QSPI_SD1", "QSPI_SD2": "QSPI_SD2", "QSPI_SD3": "QSPI_SD3"})
ctrl.add(rp, "U", "RP2040", nets, "Package_DFN_QFN:QFN-56-1EP_7x7mm_P0.4mm_EP3.2x3.2mm")
for _ in range(8):
    c(ctrl, "100n", "+3V3")
for v in ("100n", "100n", "1u"):
    c(ctrl, v, "+1V1")
c(ctrl, "1u", "+3V3")
r(ctrl, "27", "RP_DP", "USB_DP")
r(ctrl, "27", "RP_DM", "USB_DM")
r(ctrl, "10k", "RUN", "+3V3")
ctrl.add(XTAL, "Y", "12 MHz", {"1": "XIN", "3": "XTAL_B", "2": "GND", "4": "GND"}, "Crystal:Crystal_SMD_3225-4Pin_3.2x2.5mm")
r(ctrl, "1k", "XOUT", "XTAL_B")
c(ctrl, "15p", "XIN")
c(ctrl, "15p", "XTAL_B")
ctrl.add(K.library("Memory_Flash:W25Q32JVSS"), "U", "W25Q32JV", {"~{CS}": "QSPI_SS", "DO(IO1)": "QSPI_SD1", "IO2": "QSPI_SD2", "DI(IO0)": "QSPI_SD0",
                                                                  "CLK": "QSPI_SCLK", "IO3": "QSPI_SD3", "VCC": "+3V3", "GND": "GND"}, "Package_SO:SOIC-8_5.23x5.23mm_P1.27mm")
c(ctrl, "100n", "+3V3")
r(ctrl, "1k", "QSPI_SS", "BOOTSEL")
ctrl.add(K.library("Switch:SW_Push"), "SW", "BOOTSEL", {"1": "BOOTSEL", "2": "GND"}, "Button_Switch_SMD:SW_SPST_TL3342")
hc595, qn = K.library("74xx:74HC595"), ["QA", "QB", "QC", "QD", "QE", "QF", "QG", "QH"]
for k in range(7):
    n = {qn[i]: bits[8 * k + i] for i in range(8)}
    n.update({"SER": "SR_DATA" if k == 0 else f"SR_Q{k - 1}", "QH'": f"SR_Q{k}", "SRCLK": "SR_CLK", "RCLK": "SR_LATCH", "~{OE}": "SR_OE_N",
              "~{SRCLR}": "+3V3", "VCC": "+3V3", "GND": "GND"})
    ctrl.add(hc595, "U", "74HC595", n, "Package_SO:TSSOP-16_4.4x5mm_P0.65mm")
    c(ctrl, "100n", "+3V3")
r(ctrl, "10k", "SR_OE_N", "+3V3")
r(ctrl, "100", "TRIG", "TRIG_R")
coax(ctrl, "HackRF trigger U.FL", "TRIG_R")
r(ctrl, "100", "HCLK_SKT", "HCLK")
coax(ctrl, "HackRF 10 MHz in U.FL", "HCLK_SKT")
ctrl.add(K.library("Connector_Generic:Conn_02x05_Odd_Even"), "J", "Expansion",
         {"Pin_1": "+5V", "Pin_2": "+3V3", "Pin_3": "GND", "Pin_4": "GND", "Pin_5": "SR_Q6", "Pin_6": "SR_CLK", "Pin_7": "SR_LATCH",
          "Pin_8": "SYNC", "Pin_9": "EXP_IO0", "Pin_10": "EXP_IO1"}, "Connector_PinHeader_2.54mm:PinHeader_2x05_P2.54mm_Vertical")
ctrl.add(K.library("Connector_Generic:Conn_01x03"), "J", "SWD", {"Pin_1": "SWCLK", "Pin_2": "GND", "Pin_3": "SWD"},
         "Connector_PinHeader_2.54mm:PinHeader_1x03_P2.54mm_Vertical")
ctrl.notes = ["Control word: bit i is output Q(i mod 8) of register i div 8. Bits 0-15 RX, 16-31 P180, 32-47 P90 (by element), 48 MODE, 49-52 QRX, 53 HRX, 54 CAL, 55 LNA2_EN.",
              "SR_OE_N is pulled high so the registers stay off until the firmware has shifted in zeros. GPIO20 takes the 40 MHz reference, GPIO22 the HackRF's 10 MHz.",
              "The 12 MHz crystal keeps the RP2040's USB boot loader working; slot timing runs from the reference on GPIO20."]

# ---- USB-C input, hub and power ------------------------------------------------------------------------------
usb = K.Sheet("usb_power", "USB-C input, 4-port hub, 3.3 V supply", ["USB_DP", "USB_DM", "CC1", "CC2"], "A2")
usb.add(K.library("Connector:USB_C_Receptacle_USB2.0_16P"), "J", "USB-C",
        {"VBUS": "+5V", "GND": "GND", "CC1": "CC1", "CC2": "CC2", "D+": "UP_DP", "D-": "UP_DM", "SBU1": NC, "SBU2": NC, "SHIELD": "GND"},
        "Connector_USB:USB_C_Receptacle_GCT_USB4105-xx-A_16P_TopMnt_Horizontal")
r(usb, "5.1k", "CC1", "GND")
r(usb, "5.1k", "CC2", "GND")
# verify: FE1.1s strap pins (VBUSM, BUSJ, OVCJ, TESTJ) against its reference design
usb.add(K.library("Interface_USB:FE1.1s"), "U", "FE1.1s",
        {"VSS": "GND", "XIN": "HXIN", "XOUT": "HXOUT", "DMU": "UP_DM", "DPU": "UP_DP", "DM1": "A_DM", "DP1": "A_DP", "DM2": "B_DM", "DP2": "B_DP",
         "DM3": "USB_DM", "DP3": "USB_DP", "DM4": "S_DM", "DP4": "S_DP", "VDD5": "+5V", "VD33_O": "HUB_3V3", "VD33": "HUB_3V3", "VD18_O": "HUB_1V8",
         "VD18": "HUB_1V8", "REXT": "HUB_REXT", "~{XRSTJ}": "HUB_RST", "VBUSM": "HUB_VBUSM", "BUSJ": "HUB_BUSJ", "OVCJ": "HUB_OVCJ", "PWRJ": NC,
         "TESTJ/EESDA": NC, "DRV": NC, "LED1/EESCL": NC, "LED2": NC}, "Package_SO:SSOP-28_3.9x9.9mm_P0.635mm")
usb.add(XTAL, "Y", "12 MHz", {"1": "HXIN", "3": "HXOUT", "2": "GND", "4": "GND"}, "Crystal:Crystal_SMD_3225-4Pin_3.2x2.5mm")
c(usb, "15p", "HXIN")
c(usb, "15p", "HXOUT")
r(usb, "2.7k", "HUB_REXT", "GND")
r(usb, "10k", "HUB_RST", "HUB_3V3")
c(usb, "1u", "HUB_RST")
r(usb, "10k", "HUB_VBUSM", "+5V")
r(usb, "10k", "HUB_BUSJ", "GND")
r(usb, "10k", "HUB_OVCJ", "HUB_3V3")
for net, v in (("+5V", "10u"), ("HUB_3V3", "10u"), ("HUB_3V3", "100n"), ("HUB_1V8", "10u"), ("HUB_1V8", "100n")):
    c(usb, v, net)
usba = K.library("Connector:USB_A")
for name, tag in (("Radio A", "A"), ("Radio B", "B"), ("Spare", "S")):
    usb.add(usba, "J", f"USB-A: {name}", {"VBUS": "+5V", "D-": f"{tag}_DM", "D+": f"{tag}_DP", "GND": "GND", "Shield": "GND"},
            "Connector_USB:USB_A_Molex_67643_Horizontal")
    c(usb, "10u", "+5V")
usb.add(K.library("Regulator_Switching:AP63203WU"), "U", "AP63203 (3.3 V, 2 A)", {"IN": "+5V", "EN": "+5V", "GND": "GND", "SW": "BUCK_SW", "BST": "BUCK_BST", "FB": "+3V3"},
        "Package_TO_SOT_SMD:TSOT-23-6")
c(usb, "100n", "BUCK_BST", "BUCK_SW")
usb.add(L, "L", "4.7u 2A", {"1": "BUCK_SW", "2": "+3V3"}, "Inductor_SMD:L_Taiyo-Yuden_NR-40xx")
for v in ("10u", "22u", "22u"):
    c(usb, v, "+5V" if v == "10u" else "+3V3")
usb.add(FB, "FB", "600R@100MHz, 1 A", {"1": "+3V3", "2": "+3V3_RF"}, "Inductor_SMD:L_0603_1608Metric")
c(usb, "10u", "+3V3_RF")
c(usb, "100n", "+3V3_RF")
usb.notes = ["CC1/CC2 also go to two controller ADC pins, so the firmware can read how much current the port offers (1.5 A is needed for both radios).",
             "The hub's port 3 is the controller; ports 1 and 2 are the radio boards; port 4 is a spare socket.",
             "verify: FE1.1s strap pins against its reference design."]

# ---- top level ---------------------------------------------------------------------------------------------------
root = K.Sheet("root", "array16: 4x4 phased-array board, top level", [], "A1")
root.child("usb_power", usb, {p: p for p in usb.ports})
root.child("clock", clk, {"CLK40": "CLK40"})
root.child("controller", ctrl, {p: p for p in cports})
for q in range(4):
    el = [int(e) for e in B.QUAD_ELEMS[q]]
    n = {f"RX{i}": f"RX_{el[i]}" for i in range(4)}
    n.update({f"P180_{i}": f"P180_{el[i]}" for i in range(4)})
    n.update({f"P90_{i}": f"P90_{el[i]}" for i in range(4)})
    n.update({"MODE": "MODE", "QRX": f"QRX_{q}", "LNA2_EN": "LNA2_EN", "SHARED": f"S{q}"})
    root.child(f"quad{q}", quad, n)
root.child("hackrf_port", hack, {"S0": "S0", "S1": "S1", "S2": "S2", "S3": "S3", "HRX": "HRX", "LNA2_EN": "LNA2_EN"})
root.child("aux_cal", aux, {"CAL": "CAL"})
root.notes = ["Element e = 4*row + col, row 0 at the top and col 0 at the left seen from the front. Quadrants: 0 top-left, 1 top-right, 2 bottom-left, 3 bottom-right.",
              "Generated by hw/design.py. Edit the script, not these files, until the design settles."]


def build(out_dir=HERE / "kicad"):
    d = K.Design(root, out_dir)
    d.write()
    syms = "\n  ".join(K.dump([s.sx[0], K.q(s.name)] + s.sx[2:]) for s in CUSTOM)
    (out_dir / "array16.kicad_sym").write_text(f"(kicad_symbol_lib (version 20220914) (generator array16_kicadgen)\n  {syms}\n)\n")
    (out_dir / "sym-lib-table").write_text('(sym_lib_table\n  (lib (name "array16")(type "KiCad")(uri "${KIPRJMOD}/array16.kicad_sym")(options "")(descr "array16 custom symbols"))\n)\n')
    count = Counter()
    for (sheet, i, _), ref in d.refs.items():
        part = d.sheets[sheet].parts[i]
        count[part.value if part.prefix in ("U", "Y", "J", "AE", "TL", "SW", "FB") else {"R": "resistors", "C": "capacitors", "L": "inductors"}[part.prefix]] += 1
    summary = dict(parts=dict(sorted(count.items())), total=sum(count.values()), nets=len(d.netlist()),
                   pins=sum(len(d.sheets[s].parts[i].sym.pins) for (s, i, _) in d.refs))
    (HERE / "bom_summary.json").write_text(json.dumps(summary, indent=1))
    return d, summary


if __name__ == "__main__":
    _, s = build()
    print(json.dumps(s, indent=1))
