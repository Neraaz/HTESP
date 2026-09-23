# `tutorials/` — run every worked example, and say where it stopped

`examples/` ships 42 tutorials (`examples/QE/tutorial1..21` and
`examples/VASP/tutorial1..21`). Each one is a prose `README` telling you which
`mainprogram` commands to run, in which order, and which files to copy in from
an earlier tutorial first. This package turns that prose into a catalogue and
runs it: one command, launched from a batch script, walks the whole tree, waits
for the cluster jobs each step submits, checkpoints after every step — and when
it stops, tells you **exactly where**.

```
tutorials/
├── catalog.py            what the tutorials are (data only, no htesp import)
├── steps.py              the per-topic step lists the catalogue is built from
├── state.py              the checkpoint: state.json, save / load / resume
├── workdirs.py           seeding, input.in edits, squeue waiting, verification
├── runner.py             preflight + the execution engine
├── report.py             the stop report (report.md / report.json / console)
├── run_tutorials.py      main(argv) and the argument parser  ->  htesp-tutorials
├── selftest.py           unittest self-tests (stdlib only)
└── submit_tutorials.sh   the SLURM submission script
```

## The three modes

| mode | what it does | what it needs |
|------|--------------|---------------|
| *(default)* | the real campaign: submits jobs, polls `squeue`, waits, verifies | QE/VASP + SLURM |

Steps that write VASP inputs are recorded as *skipped* when pymatgen has no
POTCARs (`PMG_VASP_PSP_DIR` unset), naming
`htesp-check --config_vasp_pot /path/to/POT_GGA_PAW_PBE` as the fix: the inputs
they would write are unusable without one, and calling that done is a false
green. POTCARs are licensed, so they are never in `examples/`.

Steps that read the **output** of a real DFT run (`e0`, `2`, `4`, `phono2..4`,
`ev-collect`, `compute-elastic`, plotting…) are recorded as *skipped* with that
reason, instead of failing for a reason that is not their own.


## Where `examples/` is found

`examples/` is 185 MB of reference data and is **not** shipped inside the wheel,
so after an ordinary `pip install .` there is no `examples/` beside the
installed package. The runner looks for the tree in this order:

1. `$HTESP_EXAMPLES`
2. `./examples` under the current directory
3. beside the package (a source checkout or `pip install -e .`)

and `--examples DIR` overrides all three. `examples/` is read-only input —
`--workdir` is where runs are written, and a `--workdir` inside the example tree
is refused.

```bash
# on a laptop
python -m tutorials.run_tutorials --dry-run
python -m tutorials.run_tutorials --only QE/9,QE/12

# on a cluster
sbatch tutorials/submit_tutorials.sh --only QE
```

`htesp-tutorials` is the installed entry point (`pyproject.toml` registers
`htesp-tutorials = "tutorials.run_tutorials:main"`); `python -m
tutorials.run_tutorials` is the same thing from a checkout.

## Options

```
--workdir DIR        root for work directories, logs, checkpoint and report
--keep_output yes|no keep the generated tutorial_runs/ directories (default
                     yes).  'no' removes them at the end -- including failed
                     ones, so keep 'yes' while debugging.  Logs, the report
                     and the checkpoint are kept either way, and `--resume`
                     after a cleanup simply re-runs what was removed
--dry-run mode (default: the real thing)
--resume / --restart re-run only what is not done (default) / start over
--only  QE/9,VASP/14 run just these
--skip  QE/4,QE/5    leave these out
--from  <step-id>    start each selected tutorial at this step
--workers N          passed through to mainprogram --workers
--timeout HOURS      how long to wait for a step's cluster jobs (default 24)
--force              run even when preflight found errors
--list               print the catalogue (every tutorial, every step) and exit
-v/--verbose
```

Exit codes: `0` all finished · `1` something failed or was blocked · `2` a
preflight check failed before anything ran · `130` interrupted.

## Work directories: `examples/` is never written to

Each tutorial runs in `<workdir>/tutorial_runs/<CODE>/`, seeded in this order:

1. `examples/<QE|VASP>/` — the shared `batch.header`, `config.json`,
   `input.in` / `vasp.in`, and (QE) a symlink to the 81 `pp/*.upf`
   pseudopotentials, so 42 work directories do not cost 42 × 67 MB;
2. whatever the tutorials this one **depends on** produced — `R<mpid>-<name>/`,
   `scf_dir/`, `mpid.in`, `econv.csv`. Tutorial 9 (relaxation) is the hub most
   later tutorials start from;
3. the tutorial's own shipped files, which therefore win.

`README`, `log` and every `reference*.tar.gz` are skipped: those are the
*expected answer*, not input. The one exception is the `.cif` tutorial, whose
inputs genuinely live inside its reference archive and are unpacked from it.

## Resuming

Every step writes `<workdir>/state.json` the moment it finishes: its status
(`pending`/`running`/`done`/`failed`/`skipped`/`blocked`), exit code, duration,
log path, the artefacts expected and which were missing, and the job ids it
submitted. `--resume` (on by default) re-runs only what is not `done` or
`skipped`; `--restart` forgets the file. An interrupted run is saved too, so
Ctrl-C is a pause, not a loss.

## The stop report

On any stop — a non-zero exit, a step that exited 0 and produced nothing, a
timeout, a dependency that blocked a tutorial, or Ctrl-C — the driver prints a
summary and writes `<workdir>/report.md` and `<workdir>/report.json`:

* a one-screen table of all 42 tutorials: done / failed / blocked / skipped with
  step counts;
* **where it stopped**: tutorial code and title, step number and id, the exact
  command line, the working directory, the exit code, the last 40 lines of that
  step's log, the artefacts expected and which were missing, how long it ran;
* the exact retry line, e.g.
  `htesp-tutorials --resume --only QE/12 --from band-scf --workdir ...`;
* everything that was **blocked**, and by which dependency;
* any step that finished but could not be **verified** (no `squeue` on PATH, or
  no job ids recorded) — those are called out rather than counted as success.

Logs are at `<workdir>/logs/<CODE>/<NN>-<step-id>.log`, one per step, with the
command and working directory in the header.

## How a step is judged

A step is `done` when the command exits 0 **and** the artefacts the catalogue
declares for it exist. "Exited 0 and wrote nothing" is this package's most
common silent failure, so it is treated as a failure, and the report names the
glob that matched nothing.

After a step that submits, the runner reads the job ids the workflow layer
recorded in `<stage dir>/.htesp_job.json` (`{tag: [{"job": id, "time": epoch}]}`)
and polls `squeue -h -o %i -j <ids>` until none remain. The queue is never
grepped for a compound name. If `squeue` is missing, or no job ids were
recorded, the step is reported as *unverifiable* — never silently passed.

## Prerequisites (checked up front, all at once)

`preflight()` reports every problem before the first job is submitted:

* `python -m htesp` must start;
* `examples/<code>/batch.header` must exist, and QE needs `examples/QE/pp/*.upf`;
* **`MP_API_KEY` must be exported** for the database tutorials (QE/2, QE/3,
  QE/6, QE/16 and the VASP equivalents). The API key was removed from every
  shipped `config.json`, so these tutorials cannot run without it. Under
  `--dry-run` this is a warning and those steps are skipped; in the real mode it
  is an error;
* in the real mode: `sbatch` and `squeue` on PATH, plus `pw.x` (QE) or
  `vasp_std` (VASP);
* **enumlib** (`enum.x`/`multienum.x` and `makestr.x`/`makeStr.py`) for the
  magnetic-ordering tutorials, QE/21 and VASP/20. pymatgen's `EnumlibAdaptor`
  shells out to them and nothing pip-installs them — build from
  [enumlib](https://github.com/msg-byu/enumlib). This is a *warning*, not an
  error: the other forty tutorials run without it;
* **`PMG_VASP_PSP_DIR`** for anything that writes VASP inputs. POTCARs are
  licensed and not shipped; point pymatgen at yours with
  `htesp-check --config_vasp_pot /path/to/POT_GGA_PAW_PBE`.

Beyond that, the package's own optional dependencies must be installed for the
tutorials that use them — `ase` (OQMD, VASP structure handling), `pymatgen`
(AFLOW, convergence tests, elastic, magnetic enumeration), `bsym`
(substitutions), `phonopy` (phonons), `qmpy-rester` (OQMD), `ifermi` (Fermi
surfaces). A missing one shows up in the report as the traceback of the step
that needed it.

## Self-tests

```bash
python -m unittest tutorials.selftest -v      # or: pytest tutorials/selftest.py
```

30 tests, stdlib only, no `pymatgen` and no scheduler. They cover catalogue
integrity (dependencies resolve, no cycles, step ids unique, the QE/VASP
numbering offset, every example directory on disk), the checkpoint round trip,
`input.in` patching, artefact verification, job-id collection and the
missing-`squeue` path, the contents of a stop report for a synthesised failure,
and a full end-to-end run against a stub `mainprogram` — including a tutorial
that exits 0 and writes nothing, the tutorial blocked behind it, `--resume`,
`--restart` and `--from`.

## Known limitations

* VASP needs `POTCAR` files, which are licensed and are not shipped; configure
  pymatgen's `PMG_VASP_PSP_DIR` before running VASP tutorials for real.
* `examples/VASP/tutorial21` (3D Fermi surface) ships only `ifermi.tar.gz` — no
  `config.json`, no `input.in`, no `vasprun.xml`. It is flagged as a stub in the
  catalogue and cannot run as shipped; preflight warns about it, and
  `--skip` leaves it out if it is in the way.
* The relaxation loop (`2` → `3` → `e0` until `niteration < 3`) runs at most
  4 cycles by default; a system that has not converged by then is reported, not
  looped forever.
