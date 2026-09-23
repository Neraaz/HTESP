______________


                        ****    ****   ***********   *********   ********       ********                 
                        |  |    |  |       | |       | |____     | |            |  |  \ \                 
                        |  |____|  |       | |       |  ____|    | |*****       |  |__| |                       
                        |   ____   |       | |       | |_____          | |      |  _____/                       
                        |  |    |  |       | |       |_______|   ******| |      |_ |                       
                        |__|    |__|       |_| *********************************************                                          
                        **********************                                                               
                                      High Throughput Electron-Structure Package                              

                                              Program written by

                                Niraj K Nepal, PhD       &       Lin-Lin Wang, PhD                        
                          Email: nnepal@ameslab.gov                   Email: llw@ameslab.gov 
                                                               


_____________________________________________________________________________________________________________________________________
######################################################################################################################################

   # HTESP: DOCUMENTATION

##                                                                  Utilities:
                                

## Key functionalities:

a. Retrieving and Formatting Input Files from Materials Project, AFLOW, and OQMD Databases for Quantum Espresso (QE) and VASP Calculations.

b. Conducting Ground-State Calculations, including Structure Relaxation, Band Structure, and Density of States (DOS) Calculations, with Comprehensive Convergence Tests.

c. Performing Electron-Phonon Calculations and Investigating Superconductivity Utilizing Isotropic Eliashberg Approximation, with Spectral Function (α^2F) Plotting, Phonon Dispersion Analysis (with or without Atomic Projections).

d. Generating Input Files for Wannier90, EPW (Anisotropic Superconductivity), and WannierTools Calculations, with energies windows provided by users for wannierization.

e. Conducting Phonon and Thermodynamic Calculations Using the Phonopy Package.

f. Executing Ground-State Calculations to Construct Thermodynamic Phase Diagrams (Convex Hulls) with the Pymatgen Library.

g. Performing Fermi Surface Calculations Utilizing the IFERMI Package.

h. Computing Elastic Properties, Investigating Magnetic Ordering, and Other Related Analyses.

## Requirements

#### Current package is tested only for Linux Distribution, with Python and Bash languages.
Basic requirements

numpy, scipy, pandas, matplotlib

Pymatgen: https://pymatgen.org/

ASE: https://wiki.fysik.dtu.dk/ase/

mp-api (imported as `mp_api`): https://next-gen.materialsproject.org/api

`pyproject.toml` is the single source of truth for what is required. The full
required set is installed by `pip install .`:

numpy, scipy, pandas, matplotlib, pymatgen, mp-api, ase, spglib, PyYAML,
[bsym](https://bsym.readthedocs.io/) (site substitutions),
[lmfit](https://lmfit.github.io/lmfit-py/) (the SCDM fit for wannierisation) and
qmpy-rester (OQMD searches).

Optional features live in extras — see [Installation](#installation) below.
[IFermi](https://fermisurfaces.github.io/IFermi/) is one of them
(`htesp[fermisurface]`), not a requirement.

## Package structure

    HTESP/
      htesp/            the package: mainprogram, the workflow layer, the science modules
        data/           the packaged default config.json (every key, no API key)
      bin/              52 compatibility shims, one per former bash scan script
      tutorials/        the tutorial runner and its SLURM submission script
      tests/            the test suite (pytest or python -m unittest)
      tools/            check_names.py, gen_command_rst.py
      docs/             Sphinx documentation
      examples/         42 worked tutorials (21 QE + 21 VASP)
      utility/          templates and standalone helper scripts
      legacy/bash/      the original bash scan scripts, unmodified, for reference
      pyproject.toml    packaging and dependency metadata

## Installation

```bash
git clone https://github.com/Neraaz/HTESP.git
cd HTESP

conda create --name htesp python=3.11        # 3.10 or newer
conda activate htesp

pip install .                                 # or: pip install -e ".[test,docs]"
```

Optional extras, each pulling in only what that feature needs:

| extra | enables |
|---|---|
| `htesp[ml]` | `matminer` + `scikit-learn` for the machine-learning helpers |
| `htesp[fermisurface]` | `ifermi` + `plotly` for `mainprogram fermisurface` |
| `htesp[test]` | `pytest` |
| `htesp[docs]` | Sphinx and the theme |

Three things pip **cannot** install, because they are not Python packages:

| what | needed for | how |
|---|---|---|
| **phonopy** | `mainprogram phono1`…`phono5`, `phono-qha`, `eos-bm` | `htesp-check --install-phonopy` (conda-forge, or `--installer pip`) |
| **enumlib** | `mainprogram magenum` with `magmom.type: "ordering"` | `htesp-check --install-enumlib` (compiles from source; conda-forge has linux-64/osx-64 only) |
| **VASP POTCARs** | anything that writes VASP inputs | licensed — obtain them yourself, then `htesp-check --config_vasp_pot` |

## After installing

Run these once, in this order. Each prints what it did, and each is safe to
re-run.

```bash
# 1. does this machine have what it needs?
htesp-check                       # dependencies, architecture, page size
htesp-check --executables         # also pw.x, vasp_std, sbatch, enum.x, ...

# 2. Materials Project key -- stored, not exported
htesp-check --set_mp_api <your key>

# 3. VASP POTCARs (skip if you only run Quantum ESPRESSO)
htesp-check --config_vasp_pot /path/to/POT_GGA_PAW_PBE

# 4. phonopy, for the phonon commands
htesp-check --install-phonopy

# 5. enumlib, only for magnetic-ordering enumeration
htesp-check --install-enumlib

# 6. the compatibility shims, if you use the 1.x script names
export PATH="/path/to/HTESP/bin:$PATH"

# 7. a project configuration to edit
mainprogram config-init           # writes config.json here
mainprogram config-validate       # says which file is in use and what is wrong

# 8. a batch header for this cluster (see "Submitting to a scheduler")
mainprogram jobscript --init-header qe    # or: vasp
```

**Why `--set_mp_api` rather than `export MP_API_KEY=...`** — an exported
variable lives only in the shell that exported it, so a batch job, a `nohup`-ed
run or a new terminal loses it and every database command silently starts
skipping. `--set_mp_api` writes `~/.config/htesp/credentials` (mode `0600`) and
verifies the key against the API *before* storing it, so a typo never replaces a
working key. `$MP_API_KEY` still wins when set.

**Never put the key in `config.json`.** In 1.x a real key was committed to 218
tracked files; every shipped configuration now carries the placeholder
`use_your_API_KEY`. If you used 1.x, rotate that key — it remains in the
original repository's git history.

No `PYTHONPATH` juggling is needed: `htesp` is a real package and `mainprogram`
is a real entry point. Make sure the VASP and Quantum ESPRESSO executables are
on `PATH` as well.

`config.json` is optional. A fully-populated one ships inside the package, and
whatever you put in the working directory (or anywhere up to five directories
above it, or at `$HTESP_CONFIG`) is merged over it — so a `config.json` written
for an older version of HTESP still works, and you only need to write the keys
you want to change. `utility/input_files/config.json` is a copy to start from.

### Environment variables

| variable | effect |
|---|---|
| `MP_API_KEY` | Materials Project key; wins over the credentials file |
| `HTESP_CONFIG` | a `config.json` to use instead of searching |
| `HTESP_WORKERS` | default process-pool size (`1` disables multiprocessing) |
| `HTESP_EXAMPLES` | where `htesp-tutorials` looks for `examples/` |
| `HTESP_PYTHON` | the interpreter the `bin/` shims exec |
| `PMG_VASP_PSP_DIR` | pymatgen's POTCAR location (set by `--config_vasp_pot`) |

### Housekeeping

```bash
htesp-check --clean --dry-run     # list the build artifacts that would go
htesp-check --clean               # remove build/, *.egg-info/, __pycache__/, *.pyc
```

`--clean` touches only the checkout it is run from: it refuses an installed copy,
leaves `_removed/` alone (that is a deliberate archive), and never removes
enumlib, your API key or `PMG_VASP_PSP_DIR`.

If you copied this tree from macOS you may also have AppleDouble sidecars
(`._module.py`) beside every file — binary metadata, not code, which breaks any
tool that globs `*.py`. `--clean` reports them; remove them with
`find . -name '._*' -delete`, or avoid them with `rsync -a --exclude='._*'`.

## Running calculations

```bash
mainprogram 1                 # relax every material in the input.in range
mainprogram 1 --workers 16    # ... 16 materials at a time
mainprogram 4 --dry-run       # build all the inputs, submit nothing
mainprogram --list            # every process and command, one line each
```

Every per-material loop is a process pool. `--workers N` sets its size
(`$HTESP_WORKERS` sets the default, `min(cpu_count, 8)`; `1` disables it), each
material works in its own scratch directory, and results are ordered by material
index before anything is written, so `result.csv`, `econv.csv` and the
`mpid-*.in` lists are the same whatever the scheduling.

Exit status is `0` on success, `1` when one or more materials failed, `2` for a
bad command or a malformed `input.in`, and `130` on interrupt.

The 52 scan scripts that used to live in `src/bash` are now methods of
`HTESPWorkflow` in `htesp/workflow.py`. `bin/` holds a shim per script with the
same name and the same arguments (`start end trackfile [extra]`, `end`
exclusive), so existing habits and personal driver scripts keep working; the
originals are kept unmodified in `legacy/bash/`.

## Submitting to a scheduler

Every job HTESP submits is `batch.header` with a run command appended. The
header holds the `#SBATCH` directives and `module load` lines and nothing else.
The one in `examples/` names `--partition=dense` and a module path that exist on
one cluster and nowhere else, so copying it usually produces a job the scheduler
rejects before any calculation starts.

```bash
mainprogram jobscript --init-header qe     # or: vasp
mainprogram jobscript --init-header QuantumEspresso   # the same thing
mainprogram jobscript --init-header qe --force        # replace an existing header
```

This writes a header from what *this* machine reports: the partitions and their
cores per node (`sinfo`), the accounts you may charge (`sacctmgr`), whether the
cluster defines any generic resource at all (`scontrol show config`), the
`qe`/`quantum-espresso` or `vasp` modules Lmod offers — preferring the one Lmod
marks `(D)` — and whichever of `ibrun`, `srun` or `mpirun` is on `PATH`. On a
hierarchical Lmod site it also asks `module spider` which compiler/MPI the code
was built against and loads that chain first, because `module load qe` on its
own fails there with "these module(s) exist but cannot be loaded as
requested". The
other partitions and accounts it found are written as comments, so switching is
a matter of uncommenting a line. No run command is written into the header:
`mainprogram jobscript` appends that itself from `job_script.parallel_command`
and `job_script.nproc`, so the launcher it found is reported as a comment
naming the `config.json` values to set. How the process count is spelled
follows the launcher — `-np N` for the `mpirun` family, `-n N` for `srun` and
`aprun`, and nothing at all for `ibrun`, which runs the whole allocation.

`--init-header` prints a warning every time it writes one, and it is worth
heeding: read the file and make sure every module the build needs is loaded,
**dependencies included**. A chain one module short is accepted by the
scheduler and then fails inside the job with an error naming a shared library
rather than a module. Check it in a login shell first — `source batch.header
&& which pw.x` — and if that prints nothing the chain is incomplete.

It is a starting point, not a submit-ready job. Whatever the probes cannot
answer is left as a `# TODO` comment rather than guessed, and the node count and
wall time are placeholders — the generator knows nothing about the size of your
study. `--gres=gpu:N` appears only where SLURM actually defines GRES types; a
GPU cluster that schedules whole nodes reports `GresTypes = (null)` and gets no
`--gres` line, because one there would be rejected. Read the file before you
submit with it.

## The five commands

| command | purpose |
|---|---|
| `mainprogram` (also `htesp`) | the campaign driver — 58 named commands, 30 numbered processes |
| `htesp-check` | environment report and the one-off configuration steps |
| `htesp-tutorials` | run the 42 worked examples end to end |
| `htesp-workflow` | invoke one former scan script directly (51 names) |
| `bin/<name>` | 52 compatibility shims, `start end trackfile [extra]`, `end` exclusive |

```bash
mainprogram --list          # every process and command, one line each
mainprogram basicinfo       # the introduction
mainprogram process-info    # the numbered processes
htesp-check --help
htesp-tutorials --help
```

Full option-by-option reference, environment variables and exit codes:
[docs/usage.rst](docs/usage.rst#command-line-reference).

## Running the example tutorials

```bash
htesp-tutorials --list                    # the 42 tutorials and their steps
htesp-tutorials                           # every step, no QE/VASP/SLURM needed
htesp-tutorials --only QE/9,QE/11         # just these, with their dependencies
htesp-tutorials --only QE                 # one whole example tree
htesp-tutorials --jobs 4                  # four tutorials at once
htesp-tutorials --output                  # ... and list what it produced
htesp-tutorials --keep_output no          # drop the work dirs, keep logs+report
```

**It never runs a DFT calculation.** Every step is invoked as `mainprogram
<cmd> --dry-run`: all the input generation and file plumbing, nothing
submitted, nothing deleted. Driving a real campaign is `mainprogram`'s job.

Steps that cannot run here are *skipped with a reason* rather than failed —
reading the output of a real DFT run, a missing Materials Project key, VASP
inputs with no POTCARs configured, enumlib not on `PATH`. The first of those
also names the tutorial's own README, so you can run it for real. Steps that
need only a *relaxation* do run: its output is seeded from the tutorial
reference.

The runner checkpoints after every step, resumes with `--resume`, and when it
stops writes a report naming the tutorial, the step, the command, the working
directory, the exit code, the artifacts that were missing and the last lines of
the failing log. `--output` adds what each step wrote and what each file is
for. See `tutorials/README.md`.

`examples/` is found, not given: `$HTESP_EXAMPLES`, then `./examples`, then
beside the installed package.

## Documentation

The Sphinx sources are in `docs/`. **Two pages are generated — do not edit them
by hand:**

| page | generated from | regenerate / check |
|---|---|---|
| `docs/command.rst` | `htesp/help_text.py` | `python tools/gen_command_rst.py [--check]` |
| the `config.json` block in `docs/param.rst` | `htesp/data/config.json` | `python docs/gen_param_block.py [--check]` |

CI runs both with `--check` and fails if either is stale, so a hand edit is
silently reverted. `python docs/check_docs.py` additionally verifies that every
`:ref:` resolves, every JSON block parses, and every command named in the prose
actually exists.

## Tests

```bash
pytest tests/                                     # or:
python -m unittest discover -s tests -t .         # no pytest needed
python tools/check_names.py htesp tutorials       # undefined-name scan
```

See `tests/README.md` for what each file pins.

## Troubleshooting

| symptom | cause and fix |
|---|---|
| `PmgVaspPspDirError: PMG_VASP_PSP_DIR is not set` | POTCARs not configured — `htesp-check --config_vasp_pot /path/to/POT_GGA_PAW_PBE` |
| `EnumlibAdaptor requires the executables 'enum.x' ...` | `htesp-check --install-enumlib`, or set `magmom.order` to `["ferromagnetic"]` only |
| a search returns hundreds but writes 2 rows | the `ordering` filter — MP reports `Unknown` for materials with no magnetism calculation. The default is now `["NM", "Unknown"]`; a bare `"NM"` discards them |
| database steps skip themselves in a batch job | `$MP_API_KEY` was lost with the shell — store it with `htesp-check --set_mp_api` |
| `UnicodeDecodeError` / "null bytes" from a tool that globs `*.py` | macOS AppleDouble sidecars — `find . -name '._*' -delete` |
| `deltalake ... killed by SIGABRT` in `htesp-check` | a wheel built for 4 KiB pages on a 64 KiB-page kernel. Nothing required imports it; `htesp-check` prints the rebuild command if you want it gone |
| `Failed to resolve api.materialsproject.org` | transient DNS; re-run. `htesp-tutorials --resume` picks up where it stopped |

## Known limitations

* **POTCARs are never shipped** — they are licensed. Every VASP input-writing
  step is skipped, with that reason, until `--config_vasp_pot` is run.
* **`examples/VASP/tutorial21`** (3D Fermi surface) ships only `ifermi.tar.gz`
  and cannot run as shipped; preflight warns about it, and
  `htesp-tutorials --skip` leaves it out if it is in the way.
* **`htesp-tutorials` still uses the network** for the database tutorials: it
  runs no DFT, but it does query the Materials Project, OQMD and AFLOW. OQMD
  is the unreliable one — it gets a longer budget and a second attempt, and a
  timeout there is reported as skipped rather than failing the run.
* **enumlib and phonopy are not pip-installable**; see the table above.

### Contributors

Written and maintained by

#### Niraj K. Nepal (nnepal@ameslab.gov)

Senior Computational Scientist, Pittsburgh Supercomputing Center

#### Lin-Lin Wang

Staff Scientist, Ames National Laboratory 

## Citing HTESP

To support development activities, please cite the following paper and the papers referenced therein for calculations conducted.

N. K. Nepal, P. C. Canfield, and L.-L. Wang, HTESP (high-throughput electronic structure package): a package for the high-throughput ab initio calculations, Computational Materials Science, 244, 113247 (2024)

### Online Documentation

https://neraaz.github.io/HTESP/
