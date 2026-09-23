# Changelog

## 2.0.0

A structural release. Every command keeps its name and its arguments —
`mainprogram 4`, `mainprogram phono1`, `relax-scan 1 30 mpid.in` all still work,
`end` is still exclusive, the tracking-file and `input.in` formats are unchanged
— but almost everything under them was rebuilt.

### Installation

* The package is now `htesp/`, built from `pyproject.toml`. `pip install .`
  produces a working `mainprogram`. In 1.x the entry point was
  `src.mainprogram:main` while every module imported flatly, so the installed
  command died with `ModuleNotFoundError` and the README told users to set
  `PYTHONPATH` by hand.
* `matplotlib`, `spglib` and `PyYAML` are declared; they were imported but
  missing from `install_requires`. `ifermi`, `matminer` and `scikit-learn` moved
  to the extras `htesp[fermisurface]` and `htesp[ml]` — they were hard
  requirements for features most users never run. `qmpy-rester` stays required:
  `oqmd-search` and `oqmd-download` are ordinary commands, `oqmd_extract.py`
  imports it at module scope, and `aflow_extract.py` imports `oqmd_extract`.
* The classifiers are valid (`Programming Language :: Python :: 3 + Bash
  Scripting` is not a trove classifier and PyPI would have rejected it), and the
  licence metadata points at `LICENSE`, which is an Iowa State University / DOE
  modified MIT rather than plain MIT.
* `.gitignore` and `MANIFEST.in` were added; the tracked `.pyc` files,
  `HTESP.egg-info/` and `.DS_Store` moved to `_removed/`.

### Security

* **The Materials Project API key is no longer stored in the repository.** A
  real key was committed in 218 tracked JSON files; all of them now carry the
  placeholder `use_your_API_KEY`. The key is read from `$MP_API_KEY`, then
  `~/.config/htesp/credentials`, then `config.json`.
  **Rotate the old key** — it remains in the git history of the original
  repository, and rewriting that history is a separate step.

### The bash layer

* All 52 scan scripts are methods of `HTESPWorkflow` (`htesp/workflow.py`). The
  originals are preserved unmodified in `legacy/bash/`; `bin/` holds a shim per
  script with the same name and arguments, so `export PATH=.../bin:$PATH`
  replaces `export PATH=.../src/bash:$PATH` and nothing else changes.
* **Every per-material loop is parallel.** `--workers N` on any command;
  `$HTESP_WORKERS` sets the default (`min(cpu_count, 8)`, `1` disables it).
* **No shared scratch.** `mass.dat`, `qpoint.dat`, `kpoint.dat`, `BZ.pdf`,
  `scf_dir/kpathlines.dat`, `temp*.in` and the rest are written in a private
  per-material directory instead of the project root. `qpoint.dat` was being
  *appended* to, so a second material gave `elph.py` a two-row array.
* **Results are ordered before anything is written**, so `result.csv`,
  `econv.csv`, `mpid-finished.in`, `mpid-list-elph-*.in`, `mpid-pressure-*.in`
  and `mpid-charge.in` no longer depend on submission order.
* **No unguarded `cd`.** 82 `cd` statements had no `|| exit`; when one failed the
  body ran in the project root and the paired `cd ../../` climbed above it.
  `clean-scan` then ran `rm -r ... slurm*` there, and `pressure-relax-scan`
  deleted every pressure point because one had failed.
* **Job ids are captured.** `sbatch --parsable`, recorded in
  `<stage>/.htesp_job.json`. `mainprogram checkph` asks `squeue -j` about those
  ids instead of grepping the queue for the compound name, which matched every
  job whose name contained `B`, `C` or `Si`.
* **`--dry-run`** builds every input file and submits nothing. It is also
  non-destructive: `mainprogram 20` (clean-scan) and `mainprogram 28`
  (pressure-reset) delete nothing under it.
* **Exit codes mean something.** `1` when materials failed, `2` for a bad
  command or `input.in`, `130` on interrupt. Previously every exit code was
  discarded, so a stage in which every material crashed exited `0` and the next
  stage ran on nothing.
* `phono1-pressure` … `phono4-pressure` work. They dispatched to `mainprogram
  vp-ph2` … `vp-ph5`, which are not commands, so all four printed "Bad input".
* `generate_submission_file` no longer requires a `submission here` line in
  `batch.header` — no shipped header has one, and the result was job scripts
  that allocated the job and ran nothing. Both conventions now work.

### Configuration

* Which `config.json` a run used is now recorded by the run. `mainprogram`
  logs `configuration: <path>` at startup and writes the same line into `log`,
  next to the resolved range. The file is searched for up the directory tree,
  so a stage running in `R<id>-<comp>/relax/` can pick up a different one from
  the file beside `input.in`; reconstructing that afterwards used to mean
  re-deriving the search by hand.
* **New command `mainprogram config-init`** writes the packaged default
  `config.json` into the project directory, the counterpart to
  `config-validate`: without one, a project runs on the packaged default and
  nothing says so. It refuses to overwrite an existing file unless `--force`
  (now a global option, since `rest` is `nargs="*"` and a per-command flag
  would be swallowed by the main parser), takes an optional target path and
  creates its parents, and points out that the API key is not read from the
  file.
* `mainprogram config-validate` also prints where it looked (`$HTESP_CONFIG`,
  then `config.json` in the working directory and up to 5 parents) and which
  packaged default the result was merged over, so `(packaged default)` is no
  longer a bare statement. With `-v` it dumps the merged configuration.

* `htesp/config.py` searches the working directory and up to five parents plus
  `$HTESP_CONFIG`, deep-merges what it finds over a packaged default carrying
  every key, and caches it. `check_json.config()` used to look in exactly two
  places, return `None` on failure (every caller then raised `TypeError` several
  frames later) and re-parse the file at each of its ~40 call sites.
* A `config.json` written for an older schema therefore keeps working: the 54
  keys the shipped example configs were missing now come from the default.
* New: `mainprogram config-validate`.
* New optional keys: `plot.a2f_smearing` (which electron-phonon smearing the
  a2F plot uses — it was a hard-coded line number) and `plot.dos_ylim` (the DOS
  height, which was wrongly taken from the energy `plot.ylim`), and
  `magmom.force_theorem`.
* `pseudo.PSEUDO` gained O, Se and Au, which were in the hard-coded fallback but
  not in the shipped table, so every oxide raised `KeyError`.

### Fixed — crashes

* `input.in` with four lines raised `IndexError` on `lines[4]`; with fewer it
  left `start`/`end` unbound and raised `NameError`. On the first run it wrote
  `plot_type` as the string `'phband'` and then iterated its letters, launching
  six plot jobs named `p`, `h`, `b`, `a`, `n`, `d`. Bare `mainprogram` raised
  `IndexError` on `sys.argv[1]`.
* `create_wt_inputs.py` read `config['wanniertool_input']` (no `s`) at import
  time, so `wt1` and `wt2` could never run.
* `create_epw_inputs.py` called `epw_write.write(...)` after its `with open(...)`
  block had closed — the EPW input writer raised on every run — and used `read()`
  without importing it.
* `plot_bandproj.py` used `procar_jband` in both `nspin == 2` branches with its
  assignment commented out; every spin-polarised band projection raised
  `NameError`.
* `element_extract.py` used `input_data`, which only `__main__` bound.
* `structure_group.py` called `PWInput(filepath)`; the constructor takes a
  `Structure`, so `data-combine` never worked for QE.
* `magnetic.py` passed a `monoclinic=` argument that does not exist.
* `wannier90.py` raised `KeyError: num_wann` with the shipped `wannier90.json`.
* `generate_submission.py` left `submission_files` unbound for
  `which_calc: "WANNIER"`.
* `oqmd_extract.py` looped forever for `limit < 10`; `aflow_extract.py` busy-
  looped on any HTTP error and parsed JSON with `ast.literal_eval`.
* `displace_phonopy.py` divided a nested list by a float.
* `mainprogram magenum` on a QE campaign printed "To be implemented" and exited
  `0`; it now fails with an explanation.

### Fixed — wrong numbers

* **Elastic constants.** VASP stress was converted with a factor of `1.0` where
  pymatgen expects GPa, so every modulus was 10× too large; the first,
  unrelaxed ionic step was used instead of the last. The QE path read the
  Ry/bohr³ columns and applied the Ry/Å³ factor (wrong by 6.75) without the sign
  flip, giving negative moduli and NaN velocities. Deformed cells were
  re-standardised, snapping the strains away. The CSV header had one more column
  than its rows.
* **Fermi levels.** QE took the *first* `Fermi` match — the first ionic step of a
  relaxation — and also matched `Fermi-Dirac`. VASP read one line per ionic step
  and then called `.item()` on it, raising for any relaxation. The VASP-5 and
  VASP-6 wordings were handled in two different code paths.
* **Spin-orbit.** `grep LSORBIT INCAR | wc -l` counted `LSORBIT = .FALSE.` as on,
  so the electron count was not halved and the wrong band was highlighted.
* **Magnetic QE inputs.** Spin-decorated structures keyed the pseudopotentials by
  `"Fe,spin=5"`; undecorated ones got `starting_magnetization = 0` for every
  species under `nspin = 2`, collapsing to non-magnetic. Site labels were read
  from the cell *before* standardisation and pasted onto the standardised one.
* **k-paths.** `generate_kpath` wrote fractional ASE points without `crystal`, so
  QE read them as `tpiba`. `create_matdyn` omitted `q_in_cryst_coord`, wrong for
  every non-cubic cell. `getkpt` and `setting_qeinput` symmetrised with different
  `symprec`, so the mesh was computed for a different cell than the one written.
  A `kcutoff` sliced the path *string* by characters, cutting `G1` and `K1` in
  half; `K|U` discontinuities became a spurious `K → U` segment.
* **INCAR edits.** `sed -i '/ENCUT/d'` also deleted `ENCUTGW`, and `NELM` deleted
  `NELMIN` and `NELMDL`. `vasp.in` keys and values were parsed into two
  independent lists, so a delete line shifted every later value onto the wrong
  key.
* **Relaxed-structure extraction** assumed the `vc-relax` block layout; for
  `calculation='relax'` it deleted the header *and* the first atoms, and an
  unconverged run silently produced an empty file.
* `ml_processing.py` had `recall` and `precision` labels swapped, and an
  `np.abs()` in the Allen–Dynes denominator that manufactured a finite Tc for
  λ below about 0.19.
* `site_subs.py` mutated pymatgen's `Element` enum singleton process-wide.
* `convergence_test.py` wrote a two-line INCAR for every iteration after the
  first.
* `freq_process.py` indexed the q-axis with the band index.

### Robustness

* 310 `os.system` calls became `subprocess.run([...])` with checked return codes
  and no shell. No element name or path is ever interpolated into a shell string
  again, and `eval()` is gone from `element_extract.py`.
* 20 bare `except:` clauses (which swallowed `KeyboardInterrupt`) are typed.
* Every `os.chdir` is in a `try/finally` or replaced by `cwd=`.
* Glob deletions are anchored: `rm scf_dir/kpoint-$A-$B*` also matched `MgB2O`
  while cleaning `MgB2`.
* `mp_api` is imported where a Materials Project query happens, not at module
  scope. `htepc.py`, `element_extract.py` and `ml_processing.py` each did
  `from mp_api.client import MPRester` at the top, so importing anything that
  reached `htepc` — `cif_to_gsinput`, `elastic`, `site_subs`, `poscar_to_vasp`,
  `qe_input` — loaded `mp_api`, and with it `deltalake` and its compiled Arrow
  extension. Where that extension cannot load, the interpreter does not raise:
  it aborts (`Fatal Python error: Aborted`, SIGABRT). Writing a QE input file
  now touches none of it; `htepc.mprester()` is the single lazy entry point and
  raises a plain `ImportError` when `mp_api` is absent.
* The import smoke test sweeps in a child interpreter. Run in-process, the
  abort above killed pytest itself part-way through `tests/test_imports.py`,
  losing the ~100 tests that had not run. The child now reports each module on
  its own flushed line, so a crash becomes one named test failure and the sweep
  restarts on the modules that were left.
* Every `*.py` sweep skips macOS AppleDouble sidecars. Copying the tree to a
  non-HFS filesystem (scp/rsync/tar to an HPC login node) writes one `._x.py`
  binary blob beside each `x.py`; they matched `PKG.glob("*.py")`, so
  `ast.parse` hit `UnicodeDecodeError`, `py_compile` hit "source code string
  cannot contain null bytes", and the import sweep tried `htesp.._module`.
  `tests/helpers.python_sources()` is now the one way the suite enumerates
  modules, and `Hygiene.test_no_appledouble_sidecars` reports the strays once,
  by name, with the `rsync --exclude='._*'` fix. `._*` is in `.gitignore`;
  `utility/sample_pp/._gbrv_1.5` moved to `_removed/appledouble/`.
* `drop_incar_keys` had a test that called it as `drop_incar_keys(lines, keys)`
  when the signature is `drop_incar_keys(keys, path="INCAR")`, so the
  ENCUT/ENCUTGW regression it was meant to pin was never exercised. It now
  writes a real INCAR to a temporary directory and checks the file afterwards.

### Architecture (x86_64 and arm64/aarch64)

* The package is pure Python and architecture-neutral, and `tests/test_portability.py`
  now keeps it that way: no `platform.machine()`/`uname` branch, no hardcoded
  `x86_64`/`aarch64`/`/usr/local`/`/opt/homebrew` path, no compiled file inside
  `htesp/`, no `ctypes`, and every `#!/bin/sh` shim in `bin/` parses under POSIX
  `sh` (dash on Linux, bash-as-sh on macOS) with no bashism. 51 of the 52 shims
  are `/bin/sh`; only `jobscript.sh`, which is meant to be sourced, is bash.
* **`htesp-check`** (`htesp/check.py`, also `python -m htesp.check`) reports
  what actually differs between machines: the interpreter, `platform.machine()`,
  byte order, CPU count, and the import status of every required dependency,
  every extra, and the compiled transitive ones (`pyarrow`, `deltalake`).
  Each is imported in a **child interpreter**, so a wheel built for the wrong
  architecture is reported as `ABORTED … killed by SIGABRT` instead of killing
  the check — which is exactly how `deltalake` behaved on an aarch64 login node.
  `--json` gives the report as data, `--executables` also looks for `pw.x`,
  `vasp_std` and `sbatch`. Exit status is `0` only when every required
  dependency imports.
* `htesp-check` separates `MISSING` (not installed) from `BROKEN` (installed, but
  its own import fails) and `ABORTED` (the interpreter died). The three need
  different answers and the first two were conflated: `from mp_api.client import
  MPRester` raising `ImportError: cannot import name 'BSPathType' from
  'emmet.core.electronic_structure'` was reported as MISSING, which advised
  reinstalling mp-api — useless, because mp-api was installed and the conflict
  was with emmet-core. A `BROKEN` row now says so and is left out of the
  install block. A `ModuleNotFoundError` naming a package *other* than the one
  probed counts as broken too.
* `INSTALL/requirements1.txt` explains why `mp-api` and `emmet-core` are pinned
  as a pair, verified against upstream: mp-api 0.46.0 declares
  `emmet-core>=0.86.3` with no ceiling and imports `BSPathType` from
  `emmet.core.electronic_structure`, which 0.86.3 has and 0.87.1 does not, so
  any emmet-core ≥0.87 breaks it. The earlier note in this file claimed 0.46.0
  introduced the `deltalake` dependency; that was wrong — 0.45.0 and 0.46.0 have
  neither `deltalake` nor `pyarrow`, and both arrive in a later release
  alongside `emmet-core>=0.87.1`.
* The third column is version information in every row: the installed
  version when a package imports, and `not installed, needs >=0.33` when it does
  not — read from the installed distribution's own metadata with
  `importlib.metadata.requires`, not from `pyproject.toml`, which is not shipped
  with the package. The pip commands moved out of the rows into one
  `To install what is missing:` block at the end, one line for the required
  packages and one for the extras.
* It also reports which HTESP is running and from where, and warns when both
  `htesp` and the 1.x `HTESP` distribution are installed or a stray top-level
  `src` package is importable — the two ways a reinstall appears to do nothing.
  `docs/usage.rst` gained a `Reinstalling from scratch` section.
* `htesp-check` probes `mp_api.client`, not `mp_api`. Importing the top-level
  package runs only its `__init__.py` and reported a false `ok` on the very
  machine where `mp_api.client` -- which is what `htepc.py` imports, and what
  reaches `deltalake` -- aborted.
* A missing package is reported with what to install, not just what failed to
  import. `No module named 'qmpy_rester'` is true and useless: the distribution
  is `qmpy-rester`. `htesp-check` now prints `pip install qmpy-rester`, and
  `scikit-learn` for `sklearn`, `PyYAML` for `yaml`, `mp-api` for `mp_api`.
  A package an extra provides is named as that extra instead
  (`pip install "htesp[ml]"`). A test checks every mapping against the names
  declared in `pyproject.toml`, so the map cannot drift from the dependency
  list.
* It also names the *kind* of failure rather than assuming a wrong-architecture
  wheel, and reports the kernel's memory page size alongside the machine. The
  first real case was an aarch64 wheel on an aarch64 node that still aborted,
  because its bundled jemalloc was built for 4 KiB pages and the kernel
  (`...aarch64+64k`) uses 64 KiB; reinstalling the same wheel cannot fix that,
  so the report gives `JEMALLOC_SYS_WITH_LG_PAGE=16 pip install --no-binary
  deltalake deltalake` instead. `Illegal instruction` and
  `incompatible architecture` get their own advice.
* `INSTALL/requirements1.txt` pinned `pymatgen-core==2026.4.7`; the distribution
  is `pymatgen` (`pymatgen-core` is a different, smaller package), so the pinned
  environment could not install as written. It was also **missing `spglib` and
  `PyYAML`**, both required, so the environment it produced could not import
  `crystal` or `displace_phonopy` — while pinning `ifermi` and `qmpy-rester`,
  which are extras. Both are fixed, `plotly` was added so the
  `htesp[fermisurface]` pin is complete, and the header now states the file's
  scope and that its pins are architecture-independent while the wheels behind
  them are not.
* Dependencies are declared by their **distribution** name, not their module
  name: `mp_api` became `mp-api` in `pyproject.toml`, both requirements files,
  `docs/usage.rst` and `README.md`. pip normalises the two to the same project
  so either installs, but writing the module name is how `pymatgen-core` got in
  — a plausible-looking name that is not the package meant. Two tests reject any
  declared or pinned name containing an underscore.
* `requirements.txt` is pinned, to the same versions as
  `INSTALL/requirements1.txt`, and both pin `emmet-core==0.86.3` alongside
  `mp-api==0.46.0` — without it an unpinned resolve installs emmet-core 0.87,
  which moved `BSPathType` and breaks `import mp_api.client`. A test asserts
  every line in both files is an exact `==` pin and that no version differs
  between them.
* The two requirements files do different jobs and now say so. `requirements.txt`
  mirrors `[project] dependencies` exactly — the required set, no extras — so
  `qmpy-rester` is in it (it is required, not an extra) while `matminer`,
  `scikit-learn`, `ifermi` and `plotly` are not. Four tests keep the pair
  honest: `requirements.txt`
  equals the declared set, `requirements1.txt` covers every required package,
  nothing is pinned that no extra declares, and every pin is exact. `INSTALL/README`
  and `docs/usage.rst` contradicted each other about which file was which; both
  were rewritten.

### Fixed — the database front ends

* **Material ids came back in the wrong form.** Materials Project is migrating
  `MPID` to `AlphaID` (`mp-763` becomes `mp-bdj`). Both address the same
  material, but only the legacy spelling matches `R<mpid>-<compound>/`,
  `scf-<mpid>.in` and the tracking files of an existing campaign, so a search
  wrote ids no later command could find. Two traps, one inside the other:
  `str(AlphaID)` gives the alphabetic form while `.string` gives the legacy
  one, *and* `doc.dict()['material_id']` is already a plain alphabetic `str`
  with no `.string` left to recover. `legacy_mpid()` reads the attribute and
  decodes a bare `'mp-bdj'` through `AlphaID(...).string`, leaving OQMD and
  AFLOW ids alone. This one was a regression the rewrite introduced.
* **`ordering` matched nothing.** `Ordering` is a plain `Enum` in emmet-core,
  so `str(Ordering.NM)` is `'Ordering.NM'` and `Ordering.NM == 'NM'` is
  `False`. The filter `data['ordering'] == 'NM'` therefore matched no row: a
  search returning 465 compounds wrote 0.
* **`ordering` is now `Unknown` for much of the database.** Materials Project
  reports it for any material with no magnetism calculation — 202 of 465 in a
  boron-binary search, 158 of the 222 that survived the metal and
  formation-energy filters. A bare `"NM"` discarded all of them: 465 hits, 2
  rows written, which reads as a broken search. `filter_ordering()` now takes
  `null` (no filter at all) and a list; the shipped default is
  `["NM", "Unknown"]`.
* **Requesting `bandstructure` or `dos` rejected the whole document.**
  `MpConnect.setting()` asked for `available_fields[:-29]`, a slice that still
  contained both. emmet-core validates every requested sub-document and the
  API's payload no longer carries the fields the model declares, so one
  unusable sub-document rejected the entire `SummaryDoc`: `mainprogram
  download` raised `ValidationError` for every material that has band-structure
  or DOS data, and 1 of 64 succeeded. Safe to drop, because HTESP *computes*
  bands and DOS and never fetches them.
* **The AFLUX range operator was changed on documentation alone, and was
  wrong.** The docs describe `,` as OR and `:` as AND inside a property's
  parentheses, which reads as though a bounded range needs `1*:*2`. Measured
  against the live API, `nspecies(1*,*2)` returns binaries and `nspecies(1*:*2)`
  returns *ternaries*; the reverted spelling is the comma. A test rejects
  reintroducing `:`.
* **OQMD prompted, and nobody was there.** `qmpy_rester.get_oqmd_phases()`
  defaults to `verbose=True`, which does not merely print — it calls
  `input('Proceed? [Y/n]:')`. HTESP is always run non-interactively, so that
  either blocked forever or raised "EOF when reading a line". Its client also
  passes no timeout to its `requests.Session`, so a stalled connection ran for
  nineteen minutes with two open sockets and no output; the two calls now run
  under a scoped socket timeout.

### Fixed — submission

* **Every run line spelled the process count `-np`, whatever the launcher.**
  That is mpirun syntax. `srun` spells it `-n` and rejects `-np`; `ibrun` takes
  no count at all and rejects both. On any site whose launcher is one of those
  — most TACC and many Cray machines — every generated `run-*.sh` died on its
  first line with a launcher usage message rather than anything from the DFT
  code. `launch()` now picks the flag from the launcher, and a value that
  already carries its own flags is used verbatim.
* **Thirteen tutorials submitted nothing and still passed.**
  `stage_and_submit()` copies `run-<stage>.sh` into the stage directory and
  submits *that*; when the script is absent it records the material as
  *skipped* and carries on. Thirteen tutorials declared submitting steps
  without ever running `mainprogram jobscript`, so every submission was
  skipped and the relaxation everything else depends on was never run.
* **A VASP stage directory reached the scheduler with no POTCAR.** Only the
  download path built one, so a directory that was *seeded* rather than
  downloaded went to `sbatch` with INCAR, KPOINTS and POSCAR only.
  `_submit_vasp()` now stages one, and never raises when POTCARs are not
  configured: it explains how and the other inputs are still written.
* **`phonopy -d` has not worked since phonopy 4.** The setup operations moved
  to `phonopy-init`, and phonopy prints that to stderr and *exits 0* — so
  `mainprogram phono1` reported "Number of supercells: 0" followed by "all
  done" while writing no `phonopy_disp.yaml` and no displaced supercells.
  Silent, and wrong on any current install. Displacement generation now goes
  to `phonopy-init` where it exists; phonopy 3 is unaffected.
* **`--dry-run` suppressed the displacement generation too.** That call is
  input generation — a symmetry analysis of the relaxed cell, with no DFT in
  it — and `--dry-run` promises to build every input file. Only the
  result-consuming phonopy calls are suppressed now.

### New — `htesp-check`

`htesp-doctor` was renamed `htesp-check` (module and command), and it took on
the one-off configuration steps, each of which had been prose in a README:

* `--set_mp_api KEY` verifies the key against the API *before* writing it to
  `~/.config/htesp/credentials` (mode 0600), so a typo never replaces a working
  key. An exported `$MP_API_KEY` lives only in the shell that exported it, so a
  batch job, a `nohup`-ed sweep or a new terminal silently loses it.
* `--config_vasp_pot DIR` points pymatgen at a POTCAR tree, accepting either
  `POT_GGA_PAW_PBE` or its parent — the setting wants the parent, and getting
  it wrong costs a "POTCAR not found" for every material. It proves the
  configuration by writing a POTCAR, and when that fails it spells out the
  `pmg config -p` reorganisation a raw VASP distribution needs.
* `--install-enumlib` builds enumlib from source: it is not on PyPI, has no
  Python packaging, and its conda-forge build is linux-64/osx-64 only, so on
  aarch64 a source build is the only route.
* `--install-phonopy` installs phonopy from conda-forge, pinned to the
  interpreter running `htesp-check` so it cannot land in `base`. It then runs
  the executable, because an installer exiting 0 is not a working command.
  phonopy has no `--version` flag — it exits 2 with a usage dump — so the
  probe is `--help`.
* `--clean` returns the checkout to its pre-build state. It refuses an
  installed copy, leaves `_removed/` alone (a deliberate archive), and never
  touches enumlib, an API key or `PMG_VASP_PSP_DIR`.
* `--executables` lists enumlib's binaries too; without them a successful
  `--install-enumlib` still looked like it had done nothing.

### New — `mainprogram jobscript --init-header`

The shipped `batch.header` says `--partition=dense` and loads no module. It was
written for one machine, so copying it produces a job the scheduler rejects
before any calculation starts. `--init-header qe|vasp` writes one from what
*this* machine reports: partitions and their cores per node (`sinfo`), the
accounts you may charge (`sacctmgr`), whether the cluster defines any generic
resource at all (`scontrol`), and the code's module.

* The module is loaded **unversioned** (`module load qe`), so Lmod resolves the
  site default and the header survives that build being retired; the versions
  found are listed in a comment for pinning by hand.
* On a hierarchical Lmod site the module is invisible until its compiler and
  MPI are loaded, so `module load qe` alone fails in a job script while working
  interactively. The toolchain is read from `module help` first — that is the
  module author speaking, and on Bridges-2 it names a runtime dependency Lmod
  does not model at all — and from `module spider` when help names none.
* `--gres` is emitted only where SLURM actually defines GRES types. A GPU
  cluster that schedules whole nodes reports `GresTypes = (null)`, and a
  `--gres` line there is rejected.
* No run command is written into the header: `mainprogram jobscript` appends
  that itself, and a command in the header would run first.
* Every generated header prints a warning to read it before submitting, with a
  one-line check (`source batch.header && which pw.x`).

`qe`, `QE`, `QuantumEspresso`, `quantum-espresso` and `espresso` all name the
same code.

### Changed — the tutorial runner

* **It never runs a DFT calculation.** Real mode and `--no-dft` were removed,
  along with `--iter`, `--examples`, the `sacct` and convergence checks and the
  `squeue` polling that existed to serve them. Every step is invoked as
  `mainprogram <cmd> --dry-run`: all the input generation, nothing submitted,
  nothing deleted. Driving a real campaign is `mainprogram`'s job.
  `state.json` is at schema 2; older checkpoints are ignored rather than
  repaired. `--timeout` went with them — see the step limits below.
* **Steps that need only a relaxation now run.** Their input — the finished
  relaxation — is seeded from the QE/9 and VASP/9 references, a narrow
  whitelist that deliberately excludes `econv.csv` and `scf_dir/scf-relax-*.in`
  because those are the answers the steps must produce. Eleven steps that had
  always been skipped now execute, and 42 of 42 tutorials complete where 18
  did.
* **`--output`** lists every file the run produced, grouped by the step that
  wrote it, with a line saying what each file is for. It compares the work
  directory before and after each step, so it reports what the step *did*
  rather than what its artefact globs matched — which is how two wrong artefact
  declarations were found (`pressure-input` claimed a directory that
  `mainprogram 26` creates later; `update-input` claimed a file written only
  when a structure is *not* relaxed) along with a tutorial whose declared
  artefact is shipped inside `examples/` itself.
* **`--jobs N`** runs a dependency level in a thread pool. Almost all of a
  sweep is spent waiting on other people's servers — in one measured run seven
  `search` steps accounted for six of eight minutes — and those tutorials are
  independent, so the waiting overlaps: a full sweep went from roughly fourteen
  minutes to three.
* **`--workers` defaults to 1, because it multiplies with `--jobs`.** A
  tutorial works on one or two materials, so `mainprogram`'s own default of
  `min(cpu_count, 8)` opened eight processes to do two things — and `--jobs 4`
  made that 36 processes before any library opened a pool of its own. `ulimit
  -u` caps processes *and threads* per user, and it is **100 on a TACC login
  node** against **16384 on a compute node**, so the same command succeeds on
  one and fails on the other. It failed badly: a measured sweep lost six
  tutorials to `BlockingIOError: [Errno 11]` from inside joblib and
  `ThreadPoolBuildError { ... WouldBlock }` from inside phonopy's rayon —
  three libraries deep, with nothing naming the real limit. Preflight now
  reads that limit at run time and warns, naming `ulimit -u`.
* **A child's libraries are told how many threads they may have** when the
  limit is tight (`RAYON_NUM_THREADS`, `OMP_NUM_THREADS`, `LOKY_MAX_CPU_COUNT`
  and the rest; a variable you set yourself is respected). Left alone they
  size themselves by `cpu_count` — 144 on this machine — so a single
  `phonopy-init` asks for 144 threads to displace a two-atom cell, which
  `--workers` cannot prevent because it reaches the pool `mainprogram` opens
  and nothing below it. Together with the default above, a full sweep went
  from 35 tutorials passing in 11m18s to **42 in 2m55s**: the fix is faster as
  well as more reliable, because the pools were never doing useful work.
* **Only the steps that name a limit are capped**, and only OQMD's do: six
  minutes for its `search`, three for its `download`. Everything else runs to
  completion. Picking one number for everything was tried three times and
  failed three times — six hours let the OQMD hang below go unnoticed for
  nineteen minutes, sixty seconds failed healthy searches, and a wall-clock
  budget per *tutorial* was worse still: it is contention-blind, so `--jobs`
  made every tutorial slower without moving its deadline and the same OQMD
  search measured 34.7s alone and 100.1s at `--jobs 4`. A cap on the one
  service that has ever hung is the shape that fits: a request which has
  produced nothing for six minutes is stuck however many tutorials are
  running, while one that is merely slow is left alone.
* **A skipped step says how to run it for real**, naming the tutorial's own
  README, falling back to the QE/VASP counterpart where a tutorial ships none
  (twelve do not; only VASP/21 has neither). Steps skipped for a missing
  POTCAR, API key or enumlib keep their own message, which is more actionable.
* **OQMD is treated as the unreliable service it is**: its own step limits, a
  second attempt when one runs out, and a timeout recorded as *skipped* rather
  than failing the run — the service is outside HTESP's control and a red
  sweep should mean something is actually wrong. `data-combine` then falls
  back to that tutorial's reference rather than being blocked by a database
  nobody can reach. Only a timeout is forgiven; a wrong answer still fails.
* **The API key is resolved the way HTESP resolves it everywhere** —
  `$MP_API_KEY`, then the credentials file, then `config.json`. The runner used
  to read only the environment variable, so eight tutorials skipped on a
  machine configured exactly as the README instructs.
* `--only` accepts a bare tree name (`--only QE`), which replaced `--code`.

### Packaging

* `pyproject.toml` pins exact versions instead of floors. `pip install -e .` is
  the documented way to install HTESP and it ignored both requirements files
  entirely, which is how an environment ended up several releases ahead of
  every pin with nothing to report the difference. `emmet-core` is now a
  declared dependency rather than only a pinned transitive one: mp-api declares
  it with no upper bound and imports a name that later versions removed, so an
  unpinned resolve breaks `import mp_api.client`.
* There is no OQMD extra, and there never was one — but four places named it,
  so the `pip install` line a reader follows when `import qmpy_rester` fails
  could not work. `qmpy-rester` is a required dependency: `oqmd-search` and
  `oqmd-download` are ordinary commands, and `aflow_extract.py` imports
  `oqmd_extract` too. (The literal spelling is left out here on purpose;
  `tests/test_packaging.py` scans every document for an extra pyproject does
  not declare, and would flag this line for quoting it.)
* The QE pseudopotential table covers all 103 elements, filled out from SSSP
  1.3.0 PBE efficiency. Fifteen cutoffs that predate it were left alone
  deliberately: changing them changes the cutoffs of existing studies.

### New

* **Tests.** 450 of them, as `unittest` classes so they run under `pytest tests/`
  and `python -m unittest discover -s tests -t .` alike, with no scientific stack
  required for the core. Each test pinning a fix names the original symptom.
  `tools/check_names.py` is an undefined-name scan for machines where `ruff`
  cannot be installed.
* **A tutorial runner.** `htesp-tutorials` turns the prose `README` of each of
  the 42 worked examples into a catalogue and executes it. It finds `examples/`
  through `$HTESP_EXAMPLES`, then `./examples`, then beside the package —
  `examples/` is 187 MB and is not in the wheel, so after `pip install .` the
  package-relative path names nothing and the only message was "the example
  tree .../site-packages/examples is missing". The failure now says where it
  looked. `--workdir` inside the example tree is refused (`examples/` is
  read-only input), and the work directory is not created until preflight
  passes. It checkpoints every step, resumes, and on any stop reports the
  tutorial, the step, the command, the working directory, the exit code, the
  missing artifacts and the tail of the failing log.
* **Documentation.** `docs/command.rst` is generated from `htesp/help_text.py`,
  which is what `mainprogram` prints, so the two can no longer drift — the old
  page told readers that partial DOS was `mainprogram 20`, which is `clean-scan`
  and deletes the run. New pages on the workflow layer, the tutorial runner, the
  examples and testing. `docs/check_docs.py` checks that every `:ref:` resolves,
  every JSON block parses and every command named in the prose exists.
