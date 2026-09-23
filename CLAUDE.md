# HTESP — repository context

High Throughput Electronic Structure Package: a driver that runs automated
Quantum ESPRESSO and VASP campaigns (database search → relaxation →
electron-phonon/Tc, bands, DOS, phonons, elastic constants, EPW/Wannier90,
WannierTools, phonopy, convex hulls, Fermi surfaces).

* Version in tree: **2.0.0** (`htesp/banner.py:__version__`, `pyproject.toml`).
* Authors: Niraj K. Nepal (nnepal@ameslab.gov), Lin-Lin Wang (llw@ameslab.gov).
* Paper to cite: Nepal, Canfield & Wang, *Comput. Mater. Sci.* **244**, 113247 (2024).
* Upstream: https://github.com/Neraaz/HTESP · docs https://neraaz.github.io/HTESP/

## Provenance — this is a working fork, not upstream

`git log` HEAD is upstream `aec0369` ("Update codebase for Python 3.11
compatibility"), i.e. HTESP 1.x. **The entire 2.0 rewrite is uncommitted in the
working tree**: 240 modified, 122 deleted, 19 untracked paths. `htesp/`,
`bin/`, `tests/`, `tools/`, `tutorials/`, `legacy/`, `pyproject.toml`,
`CHANGELOG.md`, `.gitignore`, `MANIFEST.in`, `_reports/`, `_removed/` are all
untracked; the old `src/` tree (109 files, including `src/bash/*`) is staged as
deleted. Do not assume `git diff` against HEAD is a small change, and do not
`git checkout`/`git clean` anything without checking — that would destroy the
rewrite.

`CHANGELOG.md` (2.0.0) is the authoritative narrative of what changed and why:
packaging, the removal of a committed Materials Project API key from 218 JSON
files, the bash→Python port, and a long list of fixed crashes and wrong numbers
(elastic constants 10× too large, wrong Fermi levels, SOC detection, k-paths,
INCAR `sed` collisions, etc.). `_reports/*.md` are the per-area work reports
behind it (`fixes-core`, `fixes-database`, `fixes-plot-epw`, `docs-update`,
`tutorial-runner`).

## Layout

```
htesp/            the package
  cli.py          `mainprogram` dispatcher (dict of commands, not if/elif)
  workflow.py     4.4k lines: all 52 former bash scan scripts as methods
  config.py       config.json search + deep-merge + API-key resolution
  inputin.py      the six-line input.in control file
  help_text.py    the four long help blocks (docs/command.rst generated from it)
  banner.py       logo, authors, citation, __version__
  check.py        `htesp-check` dependency/architecture report, plus
                  --set_mp_api, --config_vasp_pot, --install-enumlib,
                  --install-phonopy (conda-forge, `-p sys.prefix`;
                  phonopy has NO --version flag, probe with --help)
                  and --clean
  check_json.py   back-compat shim re-exporting htesp.config
  generate_submission.py  run-*.sh writer.  `launch()` spells the
                  process count per launcher: -np (mpirun family),
                  -n (srun/aprun), nothing (ibrun, which takes the
                  whole allocation).  Before this every branch wrote
                  `-np N`, so run-*.sh died on line 1 at TACC/Cray.
  batch_header.py `jobscript --init-header qe|vasp`: probes SLURM
                  (sinfo/sacctmgr/scontrol) and Lmod ($LMOD_CMD, and
                  note listings go to stderr and `(D)` may be padded
                  far from the name) to write a starting batch.header.
                  Lmod here is HIERARCHICAL: qe/7.3 lives under
                  /opt/apps/nvidia24/openmpi5/modulefiles and is
                  invisible until nvidia+openmpi are loaded, so
                  `module load qe` alone dies in a job script.
                  prerequisites() asks `module help` FIRST (the module
                  author's own "module load intel-oneapi QE/7.5-intel",
                  as on Bridges-2) and falls back to `module spider`
                  only when help names no prerequisite -- which is the
                  usual case, Vista included, where spider is the only
                  source of nvidia/cuda/openmpi.  Prose is rejected by
                  a module-spec regex: the Help block that follows the
                  hierarchy block would otherwise parse as modules.
                  Emitted unversioned, like the code module; the probe
                  uses the NEWEST build.  write() prints a warning to
                  check the file and its module dependencies.
  data/config.json  packaged default config (every key, placeholder API key)
  <the rest>      science modules, unchanged in purpose from 1.x
bin/              52 POSIX-sh shims, one per former bash script
tutorials/        htesp-tutorials.  It NEVER runs a DFT calculation:
                  every step is `mainprogram <cmd> --dry-run`, so nothing
                  is submitted and nothing deleted.  Real mode was deleted
                  along with --dry-run/--no-dft/--iter/--examples, the
                  sacct + convergence checks, squeue polling and
                  Step.job_dirs/check_converged (state.json SCHEMA_VERSION
                  2 -- older checkpoints are ignored, not repaired).
                  Steps needing only a RELAXATION do run: its output is
                  seeded from the QE/9 and VASP/9 references
                  (REFERENCE_OUTPUT -- a whitelist, deliberately excluding
                  econv.csv and scf-relax-*.in, which are the answers those
                  steps must produce).  Steps needing a real DFT run are
                  skipped naming the tutorial's README, falling back to the
                  QE<->VASP counterpart (12 ship none; only VASP/21 has
                  neither).  --jobs N runs a dependency level in a thread
                  pool; --timeout is MINUTES per tutorial (OQMD gets 100s
                  and 2 attempts, and a timeout there is SKIPPED, not
                  failed -- its client has no request timeout, so a socket
                  timeout is applied around the two qmpy_rester calls).
                  --output lists what each step wrote, from a before/after
                  file diff: that is how two wrong artefact declarations
                  and a pre-shipped "artefact" were found.
legacy/bash/      the 52 original bash scripts, unmodified, for reference
tests/            240 unittest tests (run under pytest or unittest)
tools/            check_names.py (undefined-name scan), gen_command_rst.py
  diagnostics/    MP search diagnostics + their README (not run by CI)
tutorials/        htesp-tutorials runner for the 42 worked examples
examples/         188 MB: QE/tutorial1..21 and VASP/tutorial1..21
utility/          input-file templates (input_files/) + standalone scripts
docs/             Sphinx sources + check_docs.py + gen_param_block.py
INSTALL/          README + requirements1.txt (pinned env incl. fermisurface)
_removed/         things pulled out of the tree (pyc, .DS_Store, standard.py …)
_reports/         narrative reports for the 2.0 work
build/, htesp.egg-info/   stale build artifacts (aarch64), not inputs
```

## Architecture

Three layers, top to bottom:

1. **`htesp/cli.py`** — `mainprogram` / `htesp` / `python -m htesp`.
   `Context` lazily builds config, `input.in` and the workflow object.
   `NUMBERED` maps processes 0–29 to workflow methods; `WORKFLOW_COMMANDS` maps
   named commands (with fixed extra args); `SPECIAL_COMMANDS` are the ones that
   call a science module directly. Exit status: `0` ok, `1` one or more
   materials failed, `2` bad command or malformed `input.in`, `130` interrupt.
   Global options: `--workers N`, `--dry-run`, `--root DIR`, `--config FILE`,
   `--force`, `-v`, `--list`, `--version`.
2. **`htesp/workflow.py`** — `HTESPWorkflow`, one method per former bash script,
   keeping each script's name and argument contract (`start end trackfile
   [extra]`, **`end` exclusive**). Also runnable directly:
   `python -m htesp.workflow relax-scan 1 5 mpid.in --workers 8`.
3. **science modules** — imported and called through `run_helper(module, entry,
   *argv)`, which temporarily installs the `sys.argv` those modules expect. No
   `os.system`, no interpreter restart per material.

Key workflow machinery worth knowing before editing:

* `Material` — one `v<N> <mpid> <compound>` track-file row; `material.dir` is
  `R<mpid>-<compound>/`, `material.sub(stage)` a stage below it. Stage
  directories in use: `relax`, `calc`, `bands`, `dos`, `phonon`, `phonopy`,
  `pressure`, `epw`.
* `Result` — what one material reports back; `.say()` buffers log lines that the
  parent replays **in material order**.
* `HTESPWorkflow.map(body, materials)` — the top-level loop. Process pool of
  `min(workers, len(materials))`; `--workers`/`$HTESP_WORKERS` control it
  (default `min(cpu_count, 8)`, `1` = serial). Results are sorted by index
  before anything is written, so `result.csv`, `econv.csv` and the `mpid-*.in`
  lists are deterministic. Bodies must be bound methods (they are pickled).
* `HTESPWorkflow.scratch(material, stage)` — a private
  `.htesp_scratch/<stage>-<mpid>-<compound>/` with a local `scf_dir`. Every
  fixed-name helper output (`mass.dat`, `qpoint.dat`, `kpoint.dat`, `BZ.pdf`,
  `kpathlines.dat`, `temp*.in`) goes there, then `collect_scratch()` moves
  artifacts out. **This is what makes the loop parallel-safe — do not write
  fixed-name files into the project root from a per-material body.**
* `pushd()` — every directory change restores the previous cwd and raises on a
  missing directory (the bash layer's 82 unguarded `cd`s were the most
  destructive bug class).
* `Scheduler` — `sbatch --parsable`, job ids recorded in
  `<stage>/.htesp_job.json` as `{tag: [{"job": id, "time": epoch}]}`;
  `is_running()` asks `squeue -j <ids>` instead of grepping the queue for a
  compound name. Job-name style comes from `CALC_VISIBLE_WITH_{ID,NAME,ID-NAME}`
  marker files in the project root.
* `self.remove()` / `self.remove_glob()` — deletion paths that become no-ops
  under `--dry-run` (used by `clean-scan` and `pressure-reset`).
* `QEText` — one place for parsing/editing QE input and output text (cards,
  namelists, `final coordinates`, Fermi level, nelec, pressure, volume).

## The two control files

**`input.in`** (`htesp/inputin.py`), six lines, created automatically if absent:

```
1              start index (1-based, inclusive)
30             end index   (EXCLUSIVE — loops run ii < end)
200 0          nkpt [kcut]
mpid-list.in   tracking file
phband dos     plot types, whitespace separated
DFT = QE       QE or VASP
```

`InputIn.parse/read` is strict (raises `InputInError` → exit 2);
`InputIn.load` is the tolerant variant the workflow layer uses and returns
`dft` lower-cased.

**Tracking file** — lines `v<N> <mpid> <compound>`; `R<mpid>-<compound>/` is the
per-material directory. Missing indices are warned about and skipped.

**`config.json`** (`htesp/config.py`) — searched at `$HTESP_CONFIG`, then the
working directory and up to `SEARCH_DEPTH = 5` parents; deep-merged over the
packaged `htesp/data/config.json`, which carries **every** key, so an old
config still works. Cached per resolved path (`clear_cache()` in tests).
`mainprogram` logs `configuration: <path>` and writes it into `log`.
Top-level sections: `job_script`, `mpi_key`, `download` (`element`/`chemsys`/
`oqmd`/`aflow`/`inp`), `conv_test`, `magmom`, `pseudo` (`pot` = VASP POTCAR
labels, `PSEUDO` = QE cutoffs), `substitute`, `pwscf_in`, `strain`,
`wanniertools_input`, `kptden`, `chull_cutoff`, `kpt_opt`, `elph_mode`, `plot`.

**API key** — never from `config.json` in practice: `$MP_API_KEY`, then
`~/.config/htesp/credentials`, then the file. `htesp/config.py::api_key()` is
the ONE resolver; anything asking `os.environ["MP_API_KEY"]` directly is a bug
(the tutorial runner did, so eight tutorials skipped on a machine configured
with `htesp-check --set_mp_api`). the shipped value is the
placeholder `use_your_API_KEY`, and `api_key()` returns `None` for it.
`require_api_key()` raises with instructions.

## Command surface

* Numbered processes `0`–`29` (no 22 in `NUMBERED`; 0, 19, 22 are special-cased
  in `run_numbered`). `mainprogram --list` prints them all; the canonical
  descriptions live in `htesp/help_text.py` (`PROCESS_SUMMARY`, `SUMMARY`).
* Named commands: `search`/`download`, `oqmd-search`/`oqmd-download`,
  `aflow-search`/`aflow-download`, `data-combine`, `jobscript` (with the
  global `--init-header qe|vasp`), `convtest`,
  `compound`, `checkph`, `checkfreq`, `change_k`, `magmom_extract`, `magenum`,
  `fermisurface`, `pd`, `primtoconv`, `e0`, `pressure-input`, `charge-input`,
  `elastic-input`, `compute-elastic`, `singlemode`, `history`,
  `config-init`, `config-validate`, `basicinfo`, `process-info`, `epw-info`,
  `wt-info`; the phonopy family (`phono1..5`, `phono-qha`, `eos-bm`,
  `eos-vinet`, `ev-collect`, `phono1..4-pressure`); the EPW/Wannier90 family
  (`epw1..5`, `qe-ph`, `wann-scdm|file|random`, `epw-scdm|file|random`); and
  WannierTools (`wt1`, `wt2`).
* Eight scan scripts have **no `mainprogram` verb** and are reachable only via
  `bin/<name>` or `python -m htesp.workflow <name>`: `atom-scan`,
  `band-distort-scan`, `energy-distort-scan`, `cancel_job`, `check_calc`,
  `sumpdos.sh`, `vasp-phonopy.sh`,
  `elph_finished_but_not_copied_to_completed_folder`.
* `bin/` shims are `exec python3 -m htesp.workflow <script-name> "$@"`; 51 are
  `/bin/sh`, only `jobscript.sh` is bash because it is meant to be sourced.

## Conventions this codebase holds itself to

* **Every fix carries an inline `# FIX(n): <what was wrong>` comment** — 199 of
  them across `htesp/*.py`, numbered 1–27 by work batch. Keep the convention
  when fixing something; several tests strip comments/strings via
  `code_only()` specifically so these comments do not defeat "old pattern is
  gone" assertions.
* Tests are `unittest` classes so they run under both `pytest tests/` and
  `python -m unittest discover -s tests -t .`, and the core ones need **no**
  scientific stack. Anything needing an optional package skips with a message.
  Every test pinning a fix names the original symptom.
* `tests/helpers.py` is the shared fixture layer: `TempProject` (throwaway cwd,
  clears `$HTESP_CONFIG`/`$MP_API_KEY`/`$HTESP_WORKERS`, clears the config
  cache), `python_sources()` (the *only* sanctioned way to enumerate modules —
  it excludes AppleDouble `._*.py` sidecars), and realistic QE fixtures
  (`VC_RELAX_OUT`, `RELAX_OUT`, `UNCONVERGED_OUT`, `SCF_IN`).
* Generated docs must not drift: `docs/command.rst` comes from
  `htesp/help_text.py` via `tools/gen_command_rst.py --check`, and a param block
  from `docs/gen_param_block.py --check`; `docs/check_docs.py` verifies every
  `:ref:`, every JSON block and every command named in prose. CI runs all three.
* No `os.system`, no bare `except:`, no `eval`, no shell string interpolation —
  `tests/test_packaging.py` enforces this.
* Architecture neutrality is enforced by `tests/test_portability.py`: no
  `platform.machine()` branch, no hardcoded `x86_64`/`aarch64`/`/opt/homebrew`,
  no compiled file inside `htesp/`, no `ctypes`, and every `bin/` shim must
  parse under POSIX `sh`.

## State on this machine (verified 2026-09-18)

**Two interpreters are reachable here, and they differ — check which one is
active before drawing conclusions.**

* Path: `/work/08910/nnepal/vista/HTESP_Fork/HTESP_claude` (TACC Vista,
  **aarch64**, kernel `...aarch64+64k` — 64 KiB pages, 144 CPUs).
* The **`htesp` conda env** — `/work/08910/nnepal/vista/anaconda3/envs/htesp/`,
  Python 3.11.16 — is the real one. Every required dependency imports
  (pymatgen, mp_api, ase, spglib, bsym, lmfit, qmpy_rester, …), `htesp` 2.0.0
  is pip-installed, and `mainprogram`, `htesp-workflow` and `htesp-tutorials`
  are on `PATH`. `sbatch`/`squeue` are present; no `pw.x`, `vasp_std` or
  `phonopy` on this login node.
* The **base anaconda** `python3` (3.12.7) has no scientific stack and no
  `htesp`; under it 29 science modules skip their import tests. If the test run
  shows a wall of "needs pymatgen" skips, that is the interpreter talking.
* The site-packages copy is a **non-editable** install from an earlier
  `pip install .`, so it lags the working tree. It still carries `doctor.py`
  and the old `htesp-doctor` script; `htesp-check` appears only after
  `pip install .` is run again. Use `python3 -m htesp.check` meanwhile.
* Two HTESP distributions are installed side by side — `htesp` 2.0.0 and the
  1.x `HTESP` — which `htesp-check` warns about:
  `pip uninstall -y HTESP` removes the old one.
* `python3 -m unittest discover -s tests -t .` → **Ran 240 tests, 1 failure**.
  The single failure is
  `tests.test_packaging.Hygiene.test_no_appledouble_sidecars`: 33 macOS
  AppleDouble `._*` files came along with the `.tar.gz` copy (24 in `htesp/`,
  plus `._pyproject.toml`, `._setup.py`, some in `tests/`, `_reports/`,
  `examples/`). They also break `tools/check_names.py`, which globs `*.py`
  directly instead of going through `tests.helpers.python_sources()`. Fix is
  `find . -name '._*' -delete`, or copy with `rsync --exclude='._*'`; they are
  not source.
* `deltalake` still aborts (SIGABRT) under the 64 KiB page size — its bundled
  jemalloc is built for 4 KiB pages. Nothing required imports it here, so it can
  be left alone; `htesp-check` otherwise recommends
  `JEMALLOC_SYS_WITH_LG_PAGE=16 pip install --no-binary deltalake deltalake`.

## Extras — there is no `htesp[oqmd]`

The declared extras are exactly **`ml`, `fermisurface`, `docs`, `test`, `all`**
(`pyproject.toml`, `[project.optional-dependencies]`). `qmpy-rester` is a
**required** dependency, not an extra: `oqmd-search` and `oqmd-download` are
ordinary commands, `oqmd_extract.py` imports it at module scope, and
`aflow_extract.py` imports `oqmd_extract`.

An `oqmd` extra was referenced in four places although it was never declared, so
`pip install ".[oqmd]"` — the one instruction a reader follows when
`import qmpy_rester` fails — could not work. Corrected in `README.md`,
`docs/usage.rst`, the `PROVIDED_BY` comment in `htesp/check.py` and two
`CHANGELOG.md` entries. The *runtime* behaviour was always right:
`PROVIDED_BY` is built from `EXTRAS`, which never had `oqmd`, so
`install_hint("qmpy_rester")` returns `pip install qmpy-rester`
(pinned by `tests/test_portability.py`).

`tests/test_packaging.py::Hygiene::test_no_document_offers_an_extra_that_does_not_exist`
now scans `README.md`, `CHANGELOG.md`, `INSTALL/README`, `docs/*.rst` and the
Python sources in `htesp/`, `tutorials/` and `tools/` for both spellings
(`htesp[...]` and `pip install ".[...]"`) and fails on any extra pyproject does
not declare. **If an extra is ever added, declare it in `pyproject.toml` first.**

## Materials Project API gotchas (both cost a day already)

Two `emmet-core` types do not stringify to what the code and `config.json`
expect. Both are normalised in `htesp/element_extract.py`; use the helpers
rather than `str()` anywhere an MP field reaches a file or a comparison.

* **`Ordering` is a plain `Enum`**, so `str(Ordering.NM)` is `'Ordering.NM'`
  and `Ordering.NM == 'NM'` is `False`. The `ordering` column of
  `download/data-<elm>.csv` therefore held `Ordering.NM`, and the `extract()`
  filter `data['ordering'] == 'NM'` matched **nothing** — a search returning
  465 compounds wrote 0 rows. In chemsys mode the same defect made
  `mag_logic` false for every entry whenever `chemsys.magnetic` was false.
  Use `plain_value()`. (Latent in 1.x too; the line is identical there.)
* **Material ids are migrating from `MPID` to `AlphaID`** (`mp-763` →
  `mp-bdj`). Both address the same material and the API accepts either, but
  only the legacy spelling matches `R<mpid>-<compound>/`, `scf-<mpid>.in` and
  the tracking files of an existing campaign. Two traps, one inside the other:
  `str(AlphaID)` gives the alphabetic form while `.string` gives the legacy
  one — **and `doc.dict()['material_id']` is already a plain `str` in the
  alphabetic form, with no `.string` left to recover it.** Attribute access
  (`doc.material_id`) keeps the `MPID` object. Read the attribute and pass it
  through `legacy_mpid()`, which also decodes a bare `'mp-bdj'` string via
  `AlphaID(...).string` and leaves non-MP ids (OQMD, AFLOW) alone.
  **This one was a regression the rewrite introduced**, not a pre-existing bug.

* **`ordering` is now `Unknown` for much of the database.** MP reports it for
  any material with no magnetism calculation — 202 of 465 in a boron-binary
  search, 158 of the 222 surviving the metal and formation-energy filters. The
  `extract()` filter was a bare `==` with no on/off switch, so `ordering: "NM"`
  discarded all of them: 465 hits, 2 rows written. This is MP's data changing,
  not a code defect — 1.x gives the same counts today. `filter_ordering()` now
  accepts `null` (no filter) and a list (`["NM", "Unknown"]`, usually what is
  wanted); a bare string still behaves as before.

* **Never request MP's nested `bandstructure` / `dos` documents.**
  `MpConnect.setting()` asked for `available_fields[:-29]`, a slice that still
  contains both. emmet-core validates every requested sub-document, and the
  API's payload for these no longer carries the fields the model declares
  (`dos.elemental.Ru.total.1.task_id ... Field required`,
  `bandstructure.setyawan_curtarolo.equivalent_labels ... Field required`).
  One unusable sub-document rejects the **whole** SummaryDoc, so
  `mainprogram download` raised `ValidationError` for every material that has
  band-structure or DOS data — 1 of 64 succeeded. `htepc.UNREQUESTABLE_FIELDS`
  now filters them out. Safe because HTESP *computes* bands and DOS; it never
  fetches them. The same slice was wrong the other way too: it dropped
  `ordering`, `total_magnetization` and `theoretical`, which the shipped
  config asks for and `setting()` then re-appends.

`tests/test_regressions.py::MagneticOrderingFilter` pins these, including a
check against the real `emmet.core` types that will speak up if upstream ever
makes `Ordering` a `StrEnum`.

## Pseudopotential tables

* **QE cutoffs** live in `pseudo.PSEUDO` of the packaged `config.json`, with
  `htepc.SSSP_EFFICIENCY` as the in-code fallback. The two are **identical by
  construction** and a test enforces it — edit the config table, then mirror it.
  103 elements, filled out from **SSSP 1.3.0 PBE efficiency**
  (`https://archive.materialscloud.org/api/records/rcyfm-68h65/files/SSSP_1.3.0_PBE_efficiency.json/content`
  — the `sssp.materialscloud.org` page is a JS app, this is the data it loads).
  `maxecut_sssp()` derives `ecutrho = 8 * ecutwfc`; SSSP's own ratio is 8 for
  PAW/ultrasoft but 4 for the ONCV norm-conserving ones (He, Ne, Ar, Kr, Xe),
  so those get a denser charge grid than needed — safe, just slower.
* **15 cutoffs predate SSSP 1.3.0 and disagree with it.** Thirteen of them —
  Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm — sit at exactly `30`,
  which looks like a placeholder rather than a measurement; SSSP gives 40–150
  (Fr 150, Rn 120, Ra 90, At 80). Also Cu (ours 55, SSSP 90) and O (ours 60,
  SSSP 50, so ours is the safer side). These were **left alone deliberately**;
  changing them changes the cutoffs of existing studies.
* **A VASP stage directory gets its POTCAR at submission time.** Only the
  download path (`htesp/vasp_input.py`) used to build one, so a directory that
  was *seeded* rather than downloaded -- every tutorial starting from a
  prepared `R<mpid>-<compound>/relax/` -- went to the scheduler with INCAR,
  KPOINTS and POSCAR only and VASP stopped on the first step.
  `_submit_vasp()` and the VASP branch of `stage_and_submit()` now call
  `write_potcar.stage_potcar()`, which never raises: with no POTCARs
  configured it explains how and the other inputs are still written.
  `htesp-check --config_vasp_pot` covers a tree already in
  `<functional>/<symbol>/POTCAR` layout; a raw VASP distribution needs
  `pmg config -p <src> <dst>` first, then
  `pmg config --add PMG_VASP_PSP_DIR <dst>`.
* **VASP POTCAR names** are `pseudo.pot`, used only by
  `write_potcar.poscar2potcar()`. All 96 elements through Cm are present and
  every name validates against the current `potpaw.64` PBE catalogue. The VASP
  wiki gives no machine-readable per-element recommendation — its advice is
  prose by element class — so this table is curated, not generated. It is
  harder than pymatgen's `MPRelaxSet` choices for 15 elements (`Ti_sv` vs
  `Ti_pv`, `Mo_sv` vs `Mo_pv`, `Fe` vs `Fe_pv` …), which is a deliberate
  difference, not drift.

## AFLOW / AFLUX query operators — do not "fix" these again

`htesp/aflow_extract.py` writes ranges as `nspecies(1*,*2)`, `natoms(1*,*N)`,
`enthalpy_formation_atom(-100*,*0)` — with a **comma**. The AFLUX docs describe
`,` as OR and `:` as AND inside a property's parentheses, which reads as though
a bounded range needs `1*:*2`. The rewrite changed it on that reasoning
(FIX(12)) and it was **wrong**; measured against the live API:

```
species('Mg','B'),nspecies(1*,*2)  -> 50 records, all binary   (B3Mg1, B4Mg2 …)
species('Mg','B'),nspecies(1*:*2)  -> 50 records, all TERNARY  (Ag2B1Mg2 …)
species('Mg','B'),nspecies(2)      -> 50 records, all binary
```

`examples/QE/tutorial5`'s reference is 64 binary Mg-B entries; the `:` spelling
returned 64 ternaries containing Ag. Reverted, and
`tests/test_regressions.py::AfluxRangeOperators` rejects reintroducing `:`.
**Verify any AFLUX operator change against the API and count species** — the
documentation alone is misleading here.

## Reading a QE input with pymatgen

`Structure.from_file()` must not be used for `scf-*.in`. pymatgen's handler
registry ends with a deliberately generic `("fleur-inpgen", "*.in*")` catch-all,
so every QE input matches Fleur and raises `ModuleNotFoundError: No module named
'pymatgen.io.fleur'` (an `ImportError`, so an `except (OSError, ValueError)`
fallback never fires — that is what killed `mainprogram data-combine`). Read QE
inputs with `PWInput.from_file(path).structure` directly.

## Other known inconsistencies (cosmetic)

* `pyproject.toml`'s `script-files` comment says "The 51 former bash scripts";
  there are 52, and `bin/`, `legacy/bash/` and `script-files` all agree on 52.
* `README.md` and `docs/` say `examples/` is 185 MB; `du -sh` reports 188 MB.
* README's 1.x-era "Basic requirements"/"Extra packages" prose lists `lmfit` and
  `bsym` under extras although both are required in `pyproject.toml`.

## Common commands

```bash
# run from the repo root (nothing is installed here)
python3 -m htesp --list                 # every process and command
python3 -m htesp basicinfo              # the introduction
python3 -m htesp config-validate        # which config.json, and what is wrong
python3 -m htesp 1 --workers 16         # relax the input.in range, 16 at a time
python3 -m htesp 4 --dry-run            # build all inputs, submit nothing
python3 -m htesp.workflow relax-scan 1 5 mpid.in --dry-run

python3 -m unittest discover -s tests -t .      # 239 tests, no pytest needed
python3 tools/check_names.py htesp tutorials tests tools
python3 tools/gen_command_rst.py --check
python3 docs/check_docs.py
python3 -m htesp.check                           # dependencies + architecture
python3 -m tutorials.run_tutorials --list        # the 42 tutorials and steps
python3 -m tutorials.run_tutorials --dry-run     # whole tree, no QE/VASP/SLURM
python3 -m tutorials.run_tutorials --list-iters  # the dependency iteration plan
python3 -m unittest tutorials.selftest           # 30 stdlib-only self-tests
```

Environment variables: `MP_API_KEY`, `HTESP_CONFIG`, `HTESP_WORKERS`,
`HTESP_EXAMPLES`, `HTESP_PYTHON` (used by the `bin/` shims).

## Tutorial runner

`tutorials/` turns the prose `README` of each of the 42 examples into a
catalogue (`catalog.py`, `steps.py`) and executes it (`runner.py`), checkpointing
to `state.json` after every step and writing `report.md`/`report.json` naming
the tutorial, step, command, working directory, exit code, missing artifacts and
the log tail wherever it stopped. Two modes only: `--dry-run` (laptop) and the
default real one (polls `squeue`).  `--no-dft` was removed as one mode more
than the runner needed explaining.
A step is `done` only when the command exits 0 **and** its declared artifact
globs match — "exited 0 and wrote nothing" is treated as failure. `examples/` is
read-only input; a `--workdir` inside it is refused. QE/VASP tutorial numbering
diverges from 11 onward (VASP has no DFPT el-ph tutorial, and adds an IFermi one
at 21); `vasp_number_to_qe_number()` is the mapping.
