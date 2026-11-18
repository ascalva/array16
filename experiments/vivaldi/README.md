# experiments/vivaldi: openEMS run of the panel element

Everything here is written by `experiments/vivaldi_openems.py`. The study text is `experiments/vivaldi_study.md`.

## Check progress

```
cd /Users/ascalva/Documents/code/sdr_array/array16
.env/bin/python experiments/vivaldi_openems.py status
tail -f experiments/vivaldi/run.log
```

`status` prints the stage (validate, search, check, finished), the evaluations in the cache, the best design so
far with its parameters, the rate, and a rough ETA.

## Start, or resume after a crash, a kill or a reboot

The same command does both:

```
cd /Users/ascalva/Documents/code/sdr_array/array16
.env/bin/python experiments/vivaldi_openems.py launch
```

`launch` starts `all` in its own session under `caffeinate -i -s` (the Mac stays awake until it ends), appends to
`run.log`, and writes the Python process ID to `run.pid`. It refuses to start while that process is alive.
The plain form is equivalent:

```
nohup caffeinate -i -s .env/bin/python -u experiments/vivaldi_openems.py all >> experiments/vivaldi/run.log 2>&1 &
```

Resuming is safe at any point. Every finished simulation is one line in `evals.jsonl`; a design is looked up
there before it is simulated, and the optimiser is seeded, so a restart replays the finished evaluations from
the cache in seconds (`cached` lines in the log) and continues with the first one that is missing. Finished
stages (`validate.json`, `search.json`, `check.json` with `"complete": true`) are skipped.

To stop it: `kill $(cat experiments/vivaldi/run.pid)`.

To redo a stage, delete its `.json`; to redo simulations as well, delete their lines (or all of `evals.jsonl`).

## Results

| File | Content |
| --- | --- |
| `summary.md` | machine-written tables: validation, start against winner on both models, convergence, per-frequency active reflection, embedded pattern, compute used. Rewritten at the end of the search and of the check, and by `report` at any time |
| `*.png` | figures for `summary.md` |
| `evals.jsonl` | the cache: parameters, mesh, eps_r, the S-column of the excited element at 38 frequencies, far-field cuts for tile runs, solver time |
| `validate.json`, `search.json`, `check.json` | stage results and completion markers |
| `search_progress.json` | optimiser history so far |
| `run.log`, `run.pid` | log and process ID of the detached run |

`.env/bin/python experiments/vivaldi_openems.py report` rebuilds `summary.md` and the figures from these files
alone; it runs no simulation.

## Settings

Environment variables: `VIVALDI_ENGINE` (default `gpu`; `multithreaded` is the classic CPU engine, about 15 times
slower here; the engine is part of the cache key), `VIVALDI_THREADS`, `VIVALDI_JOBS`, `VIVALDI_SCRATCH` (solver
scratch folders; each is deleted after its result is stored, and kept only when a run fails).
