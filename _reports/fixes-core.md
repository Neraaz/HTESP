# Core science modules: confirmed-defect fixes

Scope: `htesp/htepc.py`, `elastic.py`, `magnetic.py`, `structure_group.py`,
`site_subs.py`, `convergence_test.py`, `vasp_process.py`, `poscar_to_vasp.py`.
Every change carries an inline `# FIX(<n>): <what was wrong>` comment.

Nothing could be executed end to end: `pymatgen`, `ase`, `spglib` and `pytest`
cannot be installed on this machine (proxy allowlist).  Verification was
`python3 -m py_compile htesp/*.py` plus unit tests of every shell-free helper
that does not need those packages (parsers, INCAR/QE-input rewriting, the
magnetic-input builder with a stubbed `Element`).

## Numbered fixes

| # | File | Function / place | What was wrong |
|---|------|------------------|----------------|
| 1a | htepc.py | `MpConnect.setting_qeinput` | `new_species`/`new_coord` were captured before `primitive=True` re-standardised the cell, so labels and coordinates were pasted onto the wrong sites.  The standardisation now runs first and everything below reads the cell that is written. |
| 1b | htepc.py | `MpConnect.setting_qeinput`, `_write_magnetic_input` | `pseudo` was keyed by the decorated species name (`"Fe,spin=5"`); pymatgen's `PWInput` validates `site.specie.symbol`.  Keyed by the bare symbol now.  As a consequence the ATOMIC_SPECIES block, `ntyp` and `nat` are rebuilt from the per-site type labels instead of being patched line-by-line with `sed -i`, which only lined up when pymatgen happened to emit one species line per decorated species. |
| 1c | htepc.py | `MpConnect._resolve_magnetization` (new) | Undecorated structures got `starting_magnetization(i) = 0` for every species under `nspin = 2`.  The per-element value from `config['magmom']['magmom']` is used instead, falling back to `DEFAULT_MAGMOM = 0.5` with a warning.  Values above 1 (Bohr magnetons in the config, fractions in QE) are clamped with a warning. |
| 2 | htepc.py | `MpConnect.getecut_sssp` | The `except KeyError` path re-raised.  Now: config `pseudo.PSEUDO` -> module-level `SSSP_EFFICIENCY` table (with a warning) -> a `KeyError` naming the element and telling the user to add it to `pseudo.PSEUDO`. |
| 3 | htepc.py | `INPUTscf.generate_kpath` | Wrote a bare `K_POINTS`, so QE read the fractional ASE `bandpath().kpts` as `tpiba`.  Now `K_POINTS crystal`, matching `htesp/kpath.py::printk`.  A `qualifier=` keyword allows `crystal_b`. |
| 4 | htepc.py | `INPUTscf.create_matdyn` | `q_in_cryst_coord=.true.` was missing; `htesp/matdyn.py:41` sets it.  Added. |
| 5 | htepc.py | `SYMPREC`, `INTERNATIONAL_MONOCLINIC` | `getkpt` used `symprec=0.1` and the default monoclinic setting, `setting_qeinput` used `0.01` with `international_monoclinic=False`, so the k-mesh belonged to a different cell than the one written.  Both go through the two module constants (`0.01`, `False`). |
| 6 | htepc.py | `MpConnect.setting_qeinput` | `control`/`system`/`electrons` were the module-global config sub-dicts, mutated in place; cutoffs and prefix leaked between materials.  Deep-copied per call. |
| 7 | htepc.py | `MpConnect.setting` | `kpt`, `evenkpt`, `kptype`, `kpshift`, `ecutwfc`, `ecutrho`, `comp_list`, `prefix`, `structure`, `data`, `mpid` are reset at the top. |
| 8 | htepc.py | `OUTPUTscf.extract_relax` | The `sed '$d' \| sed '1,4d' \| sed '5d'` pipeline assumed the vc-relax layout, so `calculation='relax'` lost its first two atoms and an unconverged run produced a silently empty file.  Parsed by keyword now, mirroring `workflow.py::QEText.final_coordinates`, with a last-`ATOMIC_POSITIONS` fallback and a warning when neither is present.  Returns the block. |
| 9 | htepc.py | `MpConnect.setting_qeinput` | `rm ... kpoint*` deleted other materials' files.  Only `temp.dat`, `kpoint-<mpid>.dat` and `kpath-<mpid>.dat` are removed. |
| 10 | elastic.py | `process_material` | VASP stress was `-1.0 * kBar` (pymatgen wants GPa) and read `ionic_steps[0]`, the unrelaxed step.  Now `-KBAR_TO_GPA (= -0.1) * ionic_steps[-1]['stress']`. |
| 11 | elastic.py | `qe_stress` (new) | The QE path took the Ry/bohr^3 columns `[:, :3]` and multiplied by 21798.7 (the Ry/Angstrom^3 -> kbar factor), wrong by 6.75x, and never flipped the sign.  Columns 3:6 (already kbar) are used with the same `-0.1` GPa conversion as VASP.  The `grep -A 3` shell-out is gone. |
| 12 | elastic.py | `deformation` (QE branch) | `getkpt()`/`setting_qeinput()` defaulted to `primitive=True`, which pushed the 0.5-1% strains back through `get_primitive_standard_structure()` and snapped them away.  Both are called with `primitive=False`. |
| 13 | elastic.py | `PROPNAME`, `process_material`, `main` | 18-column header vs 17-value rows.  `structure` is dropped from both, and rows are written in `PROPNAME` order (not dictionary order) so header and row cannot drift. |
| 14 | elastic.py | `main(argv=None)` | The `__main__` driver moved into `main(argv=None)` reading `input.in` through `htesp.inputin.InputIn`; `__main__` just calls it.  The old per-material `main(mpid, orig_prefix)` is now `process_material(mpid, orig_prefix, mode)`. |
| 15 | magnetic.py | `magnetic_structure` | `setting_qeinput(..., monoclinic=False, ...)` -- no such parameter, the call raised `TypeError`.  Replaced with `primitive=False`, matching the `obj.getkpt(primitive=False)` two lines above; the monoclinic setting is now `htesp.htepc.INTERNATIONAL_MONOCLINIC`. |
| 16 | magnetic.py | `magnetic_structure` | Hard-coded k-density 0.025.  Uses `config['kptden']` with 0.025 (`DEFAULT_KPTDEN`) as the fallback. |
| 17 | structure_group.py | `read_structure_file` | `PWInput(filepath)` -- the constructor takes a `Structure`.  Now `PWInput.from_file(filepath).structure`.  The surrounding bare `except:`s are typed. |
| 18 | site_subs.py | `substitution` (mode 2) | Assigning to `structure_sym[i].specie.symbol` mutates pymatgen's `Element` enum singleton process-wide.  Replaced by one `structure.replace_species(mapping)` call; names in `new_sub` that are absent from the structure are reported. |
| 19 | convergence_test.py | `main_vasp` (ecut) | `mv INCAR1` on the first iteration then `echo >>` afterwards produced a 2-line INCAR for every later cutoff.  A complete INCAR (`incar_base()` + ENCUT/ENAUG) is written every iteration. |
| 20 | convergence_test.py | `extract`, `main_qe`, `main_vasp` | `conv_test.ecut` is Ry for QE and eV for VASP and both landed in one `param-energy.txt`.  The numbers are still written verbatim; the file now starts with `# ecut unit: Ry` / `# ecut unit: eV` (`# kpoint unit: grid`) and the collected file is `convergence_result/<param>-<qe\|vasp>-<mpid>-<comp>.txt`. |
| 21 | convergence_test.py | `main_qe`, `main` | `rm scf_dir/temp*` removed other materials' scratch files.  `main_qe` returns the list of files it created and only those are removed. |
| 22 | convergence_test.py | `main(argv=None)`, `pushd` | The `__main__` driver moved into `main(argv=None)` using `InputIn`.  Every `os.chdir` chain is a `contextlib.contextmanager` `pushd()` with `try/finally`; the extraction path uses plain paths and no `chdir` at all. |
| 23 | vasp_process.py | `drop_incar_keys` | `sed -i '/{key}/d' INCAR` is a substring match: ENCUT also deleted ENCUTGW, NELM also deleted NELMIN/NELMDL.  Replaced by an in-Python rewrite matching the key at the start of a line followed by optional whitespace and `=`. |
| 24 | vasp_process.py | `parse_vasp_in` | Keys and values went into two independent lists, so a key-only (delete) line anywhere but at the end shifted every later value onto the wrong key.  One pass now yields ordered `(key, value_or_None)` pairs. |
| 25 | vasp_process.py | `vasp_process` | `m[element]` raised `KeyError` for any element missing from `config['magmom']['magmom']`.  Missing elements get `NONMAGNETIC_MAGMOM = 0.6` with a warning naming the element. |
| 26 | vasp_process.py | `band_phonopy`, `split_path_labels` | `,` was stripped from the ASE band path, joining two disconnected branches and interpolating a spurious segment.  Each comma-separated branch is now its own phonopy segment, in both `BAND =` and `BAND_LABELS =`. |
| 27 | poscar_to_vasp.py | `main` (QE branch) | `scf_dir` was never created on this branch, so the move failed and the input stayed behind as `scf-None.in`.  `os.makedirs("scf_dir", exist_ok=True)` added. |

## Cross-cutting changes (all eight files)

* **No `os.system` left.**  Replaced with `shutil.copy/move/copytree`,
  `os.makedirs(..., exist_ok=True)`, `Path.unlink`, pure-Python file rewriting,
  or `subprocess.run([...], check=True)`.  No element name or path is ever
  interpolated into a shell string.  `sbatch` uses `check=False` and logs a
  non-zero return code, because a queue refusing one job must not abort a scan.
* **No bare `except:`.**  Each is a specific exception (`KeyError`, `OSError`,
  `ValueError`, `TypeError`, `StopIteration`) or `except Exception` with a
  message; `KeyboardInterrupt` is no longer swallowed anywhere.
* **No `isdir`/`mkdir` TOCTOU.**  All `os.makedirs(..., exist_ok=True)`.
* **Every `os.chdir` is in `try/finally`** (`pushd()` in convergence_test.py,
  explicit `try/finally` in elastic/magnetic/site_subs where `poscar2potcar()`
  needs the working directory), or replaced by `subprocess.run(..., cwd=...)`.
* **API key**: `MpConnect.__init__` uses `htesp.config.api_key()` (non-raising,
  because most callers only write inputs) and `MpConnect.setting()` calls
  `require_api_key()` before contacting Materials Project.
* Public names and signatures preserved; `generate_kpath` gained a defaulted
  `qualifier=` keyword, `submission()` a defaulted `cwd=`, and `extract_relax`
  now returns the block it wrote (previously `None`).

## Not fixed / caveats

* **Nothing was executed.**  `pymatgen`/`ase`/`spglib`/`pytest` cannot be
  installed here, so `MpConnect.setting_qeinput`, `elastic.process_material`
  and `band_phonopy` were never run against real data.  Only `py_compile` and
  unit tests of the package-free helpers were possible.
* **`vasp.in` is no longer rewritten in place.**  `sed -i '/NSW/d' vasp.in`
  under `magmom.type == 'anisotropy'` destroyed the user's input file; NSW is
  now dropped from the parsed pairs instead.  The INCAR result is unchanged.
* **`starting_magnetization` clamping** (fix 1c) is a judgement call: config
  `magmom` values are Bohr magnetons and QE wants a fraction in [-1, 1], so
  anything larger is clamped to +/-1 with a warning rather than being written
  out and rejected by pw.x.
* **Fix 1b widened.**  Keying `pseudo` by the bare symbol changes what pymatgen
  writes into ATOMIC_SPECIES, so that block, `ntyp` and `nat` are rebuilt from
  the per-site labels.  This is the only way the two halves stay consistent,
  but it is more than a one-line change and deserves a real magnetic run.
* **`convergence_result` filenames changed** (fix 20) to include `qe`/`vasp`.
  Nothing in the repository reads them; only `docs/tutorial.rst` mentions the
  folder.
* `elastic.main` changed meaning (fix 14): it is now the CLI entry point, not
  the per-material routine.  Nothing imported the old one.
