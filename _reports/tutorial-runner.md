# The tutorial runner (`tutorials/`)

*What was asked for: "for tutorials, create a python files when run with batch
submission script, run all the calculations within examples. if stopped
somewhere, it will report where did it stop."*

## What was built

| file | lines | what it is |
|------|------:|------------|
| `tutorials/catalog.py` | 458 | the 42 tutorials as data: code, DFT code, title, example directory, dependencies, seeding rules, steps, convergence loop. Plus selection, topological ordering and `--list` formatting. |
| `tutorials/steps.py` | 407 | one builder per *topic* (not per tutorial number) returning the ordered `Step` list, the artefacts each step must produce, and which topics it must run after. |
| `tutorials/state.py` | 204 | the checkpoint. `RunState` / `TutorialState` / `StepState`, atomic save, versioned schema, `--resume` predicate, "where did it stop" query. |
| `tutorials/workdirs.py` | 242 | work-directory seeding, `input.in` patching, job-id collection from `.htesp_job.json`, `squeue` polling, artefact verification. |
| `tutorials/runner.py` | 498 | `preflight()` and `TutorialRunner`: subprocess execution, the relaxation loop, dependency blocking, per-step checkpointing. |
| `tutorials/report.py` | 320 | the stop report: summary table, `report.md`, `report.json`, console output, the retry command line. |
| `tutorials/run_tutorials.py` | 167 | the argument parser and `main(argv=None) -> int`, registered as `htesp-tutorials`. |
| `tutorials/selftest.py` | 490 | 30 `unittest` tests, stdlib only. |
| `tutorials/submit_tutorials.sh` | 118 | the SLURM submission script. |
| `tutorials/README.md` | 177 | how to run it. |
| `tutorials/__init__.py` | 15 | package docstring; nothing imports `htesp` at module scope. |

Nothing in `htesp/`, `tests/`, `docs/`, `utility/`, `legacy/` or `examples/` was
modified. `examples/` is read-only input: every tutorial runs in
`<workdir>/tutorial_runs/<CODE>/`.

## Catalogue coverage

All 42 tutorials have real, README-derived step lists; **one** is a stub.

* **QE/1 … QE/21** and **VASP/1 … VASP/20** — real command sequences taken from
  `examples/QE/README.txt`, `examples/VASP/README.txt` and the per-tutorial
  `README` files, cross-checked against `htesp/cli.py`'s `SPECIAL_COMMANDS`,
  `WORKFLOW_COMMANDS` and `NUMBERED` tables. Artefact globs were read off the
  shipped `reference*.tar.gz` archives and the writers in `htesp/workflow.py`.
* **VASP/21 (3D Fermi surface, IFermi) — STUB.** `examples/VASP/tutorial21`
  contains one file, `ifermi.tar.gz`: no `config.json`, no `input.in`, no
  `vasprun.xml`, no `README`. There is nothing to reconstruct a command sequence
  from beyond the index line in `examples/VASP/README.txt`, so its single
  `mainprogram fermisurface` step is a best-effort guess. It is flagged
  `stub=True`, called out in the summary table, and `--skip-stubs` drops it.

### The QE/VASP numbering offset

It is encoded, not hard-coded. `QE_TOPICS` is an ordered tuple of 21 topic keys;
`VASP_TOPICS` is derived from it by removing `elph` (QE/11, DFPT
electron-phonon — VASP has no counterpart) and appending `fermisurface`. So
VASP *n* ≡ QE *n* for *n* < 11 and VASP *n* ≡ QE *n+1* from 11 on, which
`vasp_number_to_qe_number()` states and a self-test checks for all 21 numbers.
A tutorial's dependencies are declared by *topic* (`("relax",)`) and resolved to
codes per tree, so `QE/12` depends on `QE/9` and `VASP/11` on `VASP/9` without
either list being written twice.

## How to run it

```bash
# laptop: exercises all input generation, needs no QE/VASP/SLURM
python -m tutorials.run_tutorials --dry-run
python -m tutorials.run_tutorials --dry-run --only QE/9,QE/12 --workdir /tmp/run

# cluster
sbatch tutorials/submit_tutorials.sh --code QE
sbatch tutorials/submit_tutorials.sh --no-dft        # prepare, submit nothing
sbatch tutorials/submit_tutorials.sh --resume        # carry on after a stop

python -m tutorials.run_tutorials --list             # the whole catalogue
python -m unittest tutorials.selftest -v             # the self-tests
```

Edit the four marked `#SBATCH` placeholders in `submit_tutorials.sh` first
(`<ACCOUNT>`, `<PARTITION>`, nodes, cpus, time). The site-specific
`--partition=dense -x dense001` that appears in 40 shipped example files is
deliberately **not** reproduced there. Set `HTESP_VENV` to activate an
environment; the script echoes its resolved configuration, tees the driver's
stdout to `<workdir>/driver-<jobid>.log`, and prints where the report is.

## The stop report

Printed and written to `<workdir>/report.md` + `report.json` on every stop —
a non-zero exit, a step that exited 0 and produced nothing, a timeout, a blocked
dependency, or Ctrl-C. It contains a one-screen table of all 42 tutorials, then
for the first failure and for every failure: tutorial code and title, step
number and id, the exact command line, the working directory, the exit code, the
duration, the last 40 lines of that step's log, the artefacts expected and which
were missing, and the exact retry line
(`htesp-tutorials --resume --only QE/12 --from band-scf --workdir …`). Blocked
tutorials name the dependency that blocked them; `--no-dft` adds everything that
would have been submitted; steps that finished but could not be verified get
their own section.

## Verification, waiting and checkpointing

* **Verification.** A step is `done` only when the command exits 0 *and* the
  artefact globs the catalogue declares for it match something. "Exited 0 and
  wrote nothing" is reported as a failure naming the glob that matched nothing.
* **Waiting.** After a submitting step the runner reads the job ids the workflow
  layer records in `<stage dir>/.htesp_job.json` and polls
  `squeue -h -o %i -j <ids>` (`--poll-interval`, `--job-timeout`). The queue is
  never grepped for a compound name. No `squeue`, or no recorded job ids ⇒ the
  step is reported *unverifiable*, never silently passed.
* **Checkpointing.** `state.json` is rewritten atomically after every step.
  `--resume` is the default; `--restart` forgets it.

## Limitations

1. **Nothing was executed for real.** This machine has no Quantum ESPRESSO,
   VASP, phonopy, SLURM, `pymatgen`, `ase`, `spglib`, `bsym` or `pytest`.
   Everything below the driver was proven with `--dry-run` against the real
   `mainprogram` and with a stub `mainprogram` in the self-tests. The real mode's
   `squeue` polling loop is unit-tested (including the missing-`squeue` path) but
   has never met a live queue.
2. **Artefact globs are conservative.** They were derived from reference
   archives and from the writers in `htesp/workflow.py`, not from a successful
   run. A real campaign may report a false failure if a step writes under a name
   this catalogue did not anticipate; the fix is one tuple in `steps.py`.
3. **`--no-dft` does not halt.** It runs each submitting step with
   `mainprogram … --dry-run`, records what would have been submitted, and
   carries on, so the whole pipeline's input generation is exercised in one go.
4. **The relaxation loop caps at 4 cycles** (`2` → `3` → `e0` until
   `econv.csv` reports `niteration < 3`). A system not converged by then is
   reported, not looped forever. The probe returns "cannot tell" when
   `econv.csv` is absent or unparsable, and the loop then does not repeat.
5. **`PYTHONPATH` is forced.** The `mainprogram` subprocess gets the repository
   root prepended to `PYTHONPATH`, so a checkout that was never `pip install`-ed
   works — and an installed copy elsewhere is shadowed by *this* tree.
6. **Seeding order is a judgement call.** Shared example configuration, then
   dependency output, then the tutorial's own shipped files, which win. A
   tutorial that ships a stale `Rmp-…/relax` therefore overrides the fresher one
   its dependency just produced.
7. **`data-combine` (QE/7, VASP/7)** is declared to depend on the Materials
   Project, OQMD and AFLOW tutorials and is seeded from all three. Whether that
   reproduces the intended three-database merge could not be confirmed, because
   none of the three can run here.

## Things in `examples/` that block a tutorial as shipped

| tutorial | problem |
|----------|---------|
| `VASP/tutorial21` | ships only `ifermi.tar.gz` — no `config.json`, `input.in`, `README` or `vasprun.xml`. Cannot run. Flagged as a stub. |
| `VASP/tutorial2` | has `config.json_search` and `config.json_download` but **no plain `config.json`**. The runner copies the first `config.json_*` it finds; without that, `mainprogram` falls back to the packaged default and the tutorial does not do what its index line says. |
| `QE/tutorial10` | ships `config.json` plus `config.json_ecut` / `config.json_kpoint`; only one of the two convergence modes is exercised per run. |
| `QE/tutorial11`, `QE/tutorial12` | ship no `config.json` and no `input.in` at all — their `README`s say "copy from tutorial9/tutorial1". The runner does that automatically via the dependency seeds. |
| `QE/tutorial15` | ships **no `README` at all** although it is the elastic-constants tutorial; its step list was reconstructed from `examples/VASP/tutorial14`'s README, which documents the same flow including the mid-tutorial `input.in` switch to `mpid-deformed.in`. |
| `VASP/tutorial1..6, 8, 9, 10, 12, 19, 21` | ship empty or absent `README` files; their step lists come from the QE twin plus `examples/VASP/README.txt`. |
| all database tutorials (2, 3, 6, and the hull tutorial) | the Materials Project API key was removed from every `config.json`, so they need `MP_API_KEY` in the environment. Preflight reports this once, up front, for all of them. |
| every VASP tutorial | `POTCAR` files are not shipped (licensed); `examples/VASP/README.txt` says so. pymatgen's `PMG_VASP_PSP_DIR` must be configured. |
| `QE/tutorial9` `run-scf.sh` and 40 other files | hard-code `--partition=dense -x dense001` and `mpirun -np 24`. Those are one site's and must be edited before any real run; the runner does not rewrite them. |

## What the dry run found in `htesp/` itself

The full `--dry-run` over all 42 tutorials (25 done, 15 failed, 2 blocked on
this machine) is dominated by uninstallable packages (`ase`, `pymatgen`,
`scipy`, `bsym`), which is expected here. One genuine defect surfaced:

* **`mainprogram 29` (`sitesub_scan`) exits 0 when every material failed.**
  `htesp/workflow.py`'s per-material guard logs `material mp-763-Mg1B2 failed`
  with a traceback and the process still returns 0, writing nothing. The
  tutorial runner catches it only because it verifies artefacts. The same
  pattern (`_guard` swallowing the per-material exception) affects
  `phono1`, `charge-input` and `download` for VASP. Worth a follow-up in the
  workflow layer: a scan whose materials all failed should not exit 0.
