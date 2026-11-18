# Benches

A bench is a resource setup: the compute and the hardware an experiment may use. Every experiment script takes
`--bench NAME` and reads `benches/NAME.toml` (`array16/lab.py`). Without one it uses `local`: this machine, all
its cores but two, and no devices.

    name = "simulated"
    workers = 8                             # processes for compute; left out, all cores but two

    [device.array]                          # a device, under the role experiments ask for it by
    driver = "array16.sim:SimBackend"       # any callable, as module:attribute
    radio = "quad"                          # everything else is passed to the driver

An experiment says which roles it needs and what each must be able to do (`Experiment.needs`), and the harness
refuses a bench that cannot supply them, naming what is missing. A simulator and the hardware it stands in for
are two drivers for one role, so moving an experiment from simulation to the real bench is a change of bench
file and not of experiment. Each measurement is stored with the bench that made it, and with whatever a device
reports about itself through `identity()`: a serial number, a firmware version.

| Bench | What it is |
| --- | --- |
| `local` | compute only; needs no file |
| `simulated` | compute, and the simulated board, radios and room as the `array` |

A hardware driver is a class with the same methods as the simulator it replaces (`array16/backend.py` for the
`array`), named in a bench file here: a HackRF or a pair of bladeRFs as the `array`, a microcontroller as a
`controller`, a stepper stage as a `positioner`. Drivers build on open-source tools only.
