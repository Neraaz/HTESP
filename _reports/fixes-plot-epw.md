# Fixes: plotting / EPW / WannierTools helper modules

Scope: `htesp/plot.py`, `htesp/plot_bandproj.py`, `htesp/kpath.py`,
`htesp/kpoint_path.py`, `htesp/wannier90.py`, `htesp/create_epw_inputs.py`,
`htesp/create_wt_inputs.py`, `htesp/projection_phband.py`,
`htesp/fitting_elph_smearing.py`, `htesp/displace_phonopy.py`,
`htesp/generate_submission.py`, `htesp/standard.py`, `htesp/freq_process.py`,
`htesp/elph.py`.

Every change carries a `# FIX(<n>): <what was wrong>` comment at the site.
`python3 -m py_compile htesp/*.py` passes. The modules could not be *executed*
here: `pymatgen`, `ase`, `spglib` and `pytest` are blocked by the proxy
allowlist, so everything below is verified by reading, by a module-level
undefined-name pass over the AST, and by running the pure-Python helpers
(path tokenising, label splitting, the eigenvector `np.array` fix) standalone.

---

## 1. `procar_jband` undefined in both `nspin == 2` branches
`htesp/plot_bandproj.py`, `PROCARProcessor.write_band_orbital` and
`PROCARProcessor.write_band_element`.

The assignment was commented out, so the very next line read an undefined name
(NameError) on every spin-polarised band projection. Restored to match what
the surrounding write statements index (`[0, j]` / `[1, j]`) and what the
`nspin == 1` branch does:

* `write_band_orbital`: `procar_jband = procar_iband[:, :, 0]` — after
  `data[spin, kpt, band, ion, col]` is sliced to `(2, nkpoint, nindices)`,
  taking the first matching orbital column gives `(2, nkpoint)`, mirroring the
  `[0]` the `nspin == 1` branch takes.
* `write_band_element`: `procar_jband = procar_iband` — the sum over the
  element's ions already yields `(2, nkpoint)`.

Both assignments were hoisted out of the k-point loop (they are loop-invariant).

## 2. `combine_orbital` defined twice
`htesp/plot_bandproj.py`, `DataProcessor`.

`grep -rn combine_orbital htesp/` found: def at ~434
(`(self, key, orb, if_gnu=False)`, unreachable), def at ~662 (`(self, orb)`,
the one that ran) and exactly one call site, `dp.combine_orbital(orb)` in
`main()`. Merged into a single
`combine_orbital(self, orb, key=None, if_gnu=False)`:

* `key=None` (the default, and what the only caller uses) keeps the behaviour
  that actually ran — sum `*-<orb>.dat` across every element into `<orb>.dat`;
* passing `key` restores the element-scoped variant — sum
  `<key>-*<orb>*.dat` into `<key>-<orb>.dat`, with `if_gnu` controlling the
  per-suborbital `.gnu` conversion.

Both branches now raise `FileNotFoundError` instead of `IndexError` when there
is nothing to combine, and the method returns the name it wrote.

## 3. `read` used but never imported in `ciftoxsf`
`htesp/create_epw_inputs.py`. Added `from ase.io import read`.

## 4. EPW writes after the `with` block closed
`htesp/create_epw_inputs.py`, `epw_sc`.

Three `epw_write.write(...)` calls sat *after* the `with open(out, "w")` block,
so they ran against a closed handle: `ValueError: I/O operation on closed
file` on every `epw` condition — the whole EPW input path was dead. The writes
moved inside the block. While there, `epw_sc_from_json("epw.json", out='epw.in')`
now passes the caller's `out` instead of the hard-coded `'epw.in'`, and the
local `config = json.load(f)` in `epw_sc_from_json` was renamed `epw_config`
because it shadowed the imported `htesp.check_json.config`.

## 5. `epsil` written for any non-`T` entry, at any q, magnetic or not
`htesp/create_epw_inputs.py`, `phonon_input`.

`docs/otherinput.rst:173` documents `epsil` as applying "for non magnetic and
at q = 0". The code now writes `epsil = .true.` only when all three hold:
not a metal, the q-point in `ph-q.in` is `0 0 0` (new `_is_gamma` helper), and
`config()['pwscf_in']['magnetic']` is false. Each skipped case prints why.
The `reduce_io`/namelist-terminator/q-line writes were de-duplicated out of the
two branches.

## 6. `wanniertool_input` key + config read at import time
`htesp/create_wt_inputs.py`.

The module read `json.load(...)['wanniertool_input']` at module scope; the real
key is `wanniertools_input`, and the surrounding `except FileNotFoundError`
does not catch the resulting `KeyError`, so **importing** the module raised and
`wt1`/`wt2` were unusable. Now:

* the defaults live in a module constant `DEFAULT_WT_INPUT` (no I/O at import);
* a new `read_wt_input()` calls `htesp.check_json.config()`, reads the correct
  `wanniertools_input` key and merges it over the defaults;
* `wt_body()` calls `read_wt_input()` itself.

## 7. Duplicate `kpoint_path` in `create_wt_inputs.py`
Deleted the local copy; `create_wt_inputs` now imports the shared
`htesp.kpoint_path.kpoint_path`. Signature difference reconciled by adding a
`slab=False` keyword to the shared function: the local copy also wrote
`wannier_kpath_slab.in` from `POSCAR-slab`, so `create_wt_inputs.main()` calls
`kpoint_path(file_name, slab=True)`. The shared version additionally brings the
VASP fallback and the `|` handling the copy lacked.

## 8. `config_settings["num_wann"]` KeyError
`htesp/wannier90.py`. The shipped `utility/input_files/wannier90.json` has no
`num_wann` key, so every VASP wannier90 input died. New `resolve_num_wann()`
resolves, in order: `config_settings['num_wann']` → the orbital count derived
from `projection.in` (`_count_projection_wannier`, s/p/d/f/sp…sp3d2 table) →
`config_settings['num_bands']` with a printed warning → otherwise a
`ValueError` naming exactly what to add to the JSON file.

## 9. `bands_plot == '.true.'` string comparison
`htesp/wannier90.py`. Added `_as_bool()` accepting `True`, `'.true.'`,
`'true'`, `'t'`, `'.t.'`, `'yes'`, `1`; the `Begin Kpoint_Path` branch uses it.

## 10. `K|U` emitted as a consecutive pair
`htesp/kpoint_path.py`. A combined discontinuity label was appended verbatim
(`band_dict['K|U']` → `KeyError`) and then paired with its neighbour, drawing a
bogus `K -> U` segment straight across the break. New `split_path_segments()`
splits on `|` and `,` and returns one list of labels per continuous stretch;
only pairs *within* a stretch are written. Unknown labels now raise a `KeyError`
that names them.

## 11. `kcutoff` slicing the path string by characters
`htesp/kpath.py`. `bandpath.path[:kcutoff]` cut `G1`/`K1` in half. New
`path_tokens()` / `cut_path()` tokenise the path into whole labels
(`[A-Za-z][0-9]*`, plus `|`/`,` separators) and slice the token list; a
trailing separator is dropped. Semantics preserved: `kcutoff = n` keeps the
first *n* labels, which is what the old character slice did for the
single-character labels it worked on. (Note: the `kpath` docstring claimed
"'1' means to remove one symmetry point"; the code never did that, and the
docstring was corrected to match the code rather than the other way round.)

## 12. `BZ.pdf` plotted on every call, default backend
`htesp/kpath.py`. `matplotlib.use("Agg")` is now set at import, before
`pylab`, and `kpath()` gained `plot_bz=False` / `bz_file="BZ.pdf"`. Only
`printk()` opts in (`plot_bz=True` by default) because `htesp/workflow.py:1441`
harvests `BZ.pdf` from that call. `kpoint_path`, `create_wt_inputs`, `plot` and
`htepc` no longer redraw and overwrite it on every call.
`matplotlib.use("Agg")` was also added to `create_epw_inputs.py` and
`fitting_elph_smearing.py`, which import `pyplot` without selecting a backend.

## 13. `scf_dir/` hard-coded in `printk`
`htesp/kpath.py`. `printk(out_dir="scf_dir", plot_bz=True)`; `kpathlines.dat`
and `kspecial-points.dat` are written under `out_dir`. The positional CLI
contract is kept and extended: `sys.argv[6]`, when present, overrides
`out_dir`. `htesp/workflow.py` calls
`run_helper("kpath", "main", "point", file, nkpt, kcut, 0)` — five positional
arguments, unchanged and still valid.

## 14. Fermi level: first match, `Fermi-Dirac`, `.item()`, split VASP wordings
`htesp/plot.py`, `plot()` and `band_wann_plot()`.

Old behaviour: QE — `grep Fermi scf.out | awk '{print $5}'` then take the
FIRST line (the first ionic step of a relaxation; it also matched
"Fermi-Dirac"); VASP — `grep E-fermi ../relax/OUTCAR` (one line per ionic step)
then `.item()`, which raises for any relaxation with more than one step; and
`band_wann_plot` grepped only the VASP-6 wording while `plot()` grepped only
VASP-5. Replaced by one set of helpers used by both:

* `fermi_from_qe_out()` — regex for `the Fermi energy is <x>`, LAST match;
  a spin-polarised run printing `the spin up/dw Fermi energies are <x> <y>` is
  recognised and the higher of the two is used, with a log line;
* `fermi_from_outcar()` — one regex matching **both** `E-fermi :` (VASP 5) and
  `Fermi energy:` (VASP 6), LAST match;
* `read_fermi(qe_out, outcar)` — dispatches on which file exists.

`NELECT` / `number of electrons` are read the same way
(`nelect_from_outcar`, `nelect_from_qe_out`), removing the
`fermi.dat` / `band_fermi.dat` scratch files and the `rm` that followed.
A missing output now raises `FileNotFoundError` instead of continuing with an
unbound `fermi`.

## 15. `grep LSORBIT INCAR | wc -l`
`htesp/plot.py`. Counted `LSORBIT = .FALSE.` as spin-orbit ON, so the band
index at the Fermi level was left un-halved. New `read_lsorbit()` parses the
value with `pymatgen.io.vasp.inputs.Incar`; the `lsorbit` scratch file is gone.

## 16. `Kpoints.from_file("KPOINTS").labels is None`
`htesp/plot.py`, `kptline()`. With the shipped `kpt_opt: true` the line-mode
path lives in `KPOINTS_OPT` and `KPOINTS` holds a plain mesh, so `.labels` is
None and `len(None)` raised `TypeError` — every VASP band plot failed with the
shipped config. New `kpt_labels()` tries, in order: `KPOINTS`, `KPOINTS_OPT`,
then the `label` file `band_wann_plot()` writes from `*_band.labelinfo.dat`;
if none of them yields labels it raises a `RuntimeError` that says which files
were tried and what to do. When the labels come from the `label` file the
segment length is inferred from the vasprun distances.
(`high_symm.in`, written by `vasp_process.py:302`, holds coordinates only — no
labels — so it is not usable as a label source.)

## 17. `a2f`: lambda from line 12, alpha2F from column 2
`htesp/plot.py`. One arbitrary smearing out of ten, with nothing tying the two
choices together. New `parse_lambda_out()` parses the whole
`lambda / omega_log / T_c` table after the `omega_log` header, and
`select_smearing()` picks a row. A new config key `plot.a2f_smearing`
(1-based smearing number) selects both the lambda row and the `alpha2F.dat`
column. **Leaving it unset reproduces exactly the old choice** (the row that
`lambda.out` line 12 resolves to, and column 2), so existing plots do not
change silently. The chosen smearing, lambda, omega_log and Tc are printed on
every run. An out-of-range column raises `IndexError` with the counts.

## 18. Four defects in the DOS plotting
`htesp/plot.py`.

* `.dos` header token `[8]` — right only for `nspin = 1`; the `nspin = 2`
  header has an extra column. Parsed by name now (`EFermi\s*=\s*<x>`), with a
  clear `ValueError` naming the offending header.
* `xlim1[0]` / `ylim1[0]` were indexed *before* the `is not None` tests, so a
  null `xlim` raised `TypeError`. Reordered: the bounds are only unpacked
  inside the `is not None` branch.
* the energy `ylim` was reused as a DOS-height limit (in `dos_plot_vasp`'s
  `nspin == 2` branch, twice, and in `dos_plot`). A separate `plot.dos_ylim`
  key now drives the DOS axis; when it is absent the previous data-driven /
  hard-coded defaults apply. **Behaviour change:** a config that set `ylim`
  and relied on it clamping the DOS height must now set `dos_ylim`.
* `ncol = int(ndos_data/4)` gave `legend(ncol=0)` (a `ValueError`) for fewer
  than four curves → `max(1, ...)`.

## 19. Second `grep > file` overwrote the first
`htesp/fitting_elph_smearing.py`, `fit_param`. Both the `*` (overflow) and
`NaN` counts were written to the same file `file`, so the `*` check never ran.
Both counts are taken in Python now and summed; the scratch file is gone.
`parse_lambda()` was rewritten to return the table as an array instead of
shelling out to `sed` into `lambda.txt` (and then `rm`-ing it), and the bare
`except:` became a typed, logged `except`.

## 20. Nested YAML list divided by a float
`htesp/displace_phonopy.py`, `applydisplace`. `self.eigen[i]` is
`[[re, im], [re, im], [re, im]]` straight out of YAML; `list / float` raised
`TypeError`. Wrapped in `np.array(..., dtype=float)`. Also added explicit
errors when `geteigenvec()` / `mass_extract()` were not called first.

## 21. `which_calc` accepted in one place, not the other
`htesp/generate_submission.py`. `main` dispatched on `'wannier'` exactly while
the generator accepted `'W'`/`'WANNIER'`, so `which_calc: "WANNIER"` fell
through every branch, `submission_files` stayed unbound and the next loop
raised `NameError`. Normalisation now happens in ONE place —
`normalise_calc()`, `str(x).strip().lower()` through an alias table
(`qe`, `epw`, `wannier`/`w`/`wan`, `vasp`) — and both `main` and
`generate_submission_files` dispatch on the normalised value. `main` picks the
command dict from a lookup table and raises a `ValueError` naming the bad value
instead of falling through.

## 22. `json.load('./config.json')` instead of `check_json.config()`
`htesp/generate_submission.py`. Replaced by `from htesp.check_json import
config` / `config()['job_script']`, so it gets the cwd+parents+`$HTESP_CONFIG`
search and the packaged-default merge. This also fixes the case where
`config.json` was not in the cwd and `parameters` stayed unbound.

## 23. `standard.py` removed
Verified with `grep -rn "standard" htesp/ legacy/ | grep -v standard.py` that
nothing imports it (all hits are `get_primitive_standard_structure` /
`get_conventional_standard_structure` / the unrelated `standardize` method on
`htepc.py`). It is also broken three ways: `Cell.get_bravais_lattice(Atoms)`
(called on an `Atoms`, not a `Cell`), `BravaisLattice[i][j]` (not subscriptable)
and `Fraction(coord).limit_denominator(10)` applied to Cartesian/scaled
coordinates. Moved with `mv htesp/standard.py _removed/standard.py`. **Not
deleted.**

## 24. `freq_process.py` and `elph.py`
### `freq_process.py`
* `row[i]` inside the band loop indexed the q-axis with the *band* index: every
  band got the same wrong abscissa, and the loop raised `IndexError` as soon as
  `nband > nqpoint` (the normal case for a dense q-mesh). Now `row[j]`.
* the output file was the fixed name `freq.plot`; an optional `sys.argv[2]`
  overrides it. `htesp/workflow.py` calls
  `run_helper("freq_process", "freq_process", comp)` — one argument, unchanged.

### `elph.py`
* `mass` / `nion` stayed unbound when `mass.dat` was absent, so every
  `write_file()` call raised `NameError` instead of simply omitting the
  `amass()` lines. Initialised up front.
* `qpt` stayed unbound in `write_file` and `irr_q` when no q-point file was
  present. New `read_qpoint_mesh()` raises a `FileNotFoundError` naming the
  files it looked for.
* `grep "irreducible representations" .../elph.out > irr.out` followed by
  `np.genfromtxt("irr.out")[:, 2]` — a fixed-name scratch file plus a fixed
  column index that silently yields NaN if QE's wording changes. Replaced by
  `count_irr(elph_out)`, a regex over the output; a missing `elph.out` or a
  short count now raises with a message.
* `parse_symmetry_analysis` called `.strip()` on `qpoints[i]`, which is a list
  of floats — `AttributeError`, so that function could never return. Removed.
* `os.system` calls (`rm`, `touch`, `echo >>`, `mkdir`) replaced with
  `os.remove`, `os.makedirs(exist_ok=True)` and plain file writes; the
  `PARALLEL_irr` file is now written from a single `with open(...)` block.
* `if "elph_mode" in input_data.keys()` → `input_data.get(...)`.

---

## Cross-cutting cleanups (all files in scope)

* **`os.system` removed everywhere in scope.** Replaced by `shutil.copy` /
  `shutil.move`, `os.remove`, `os.makedirs(exist_ok=True)`, plain file writes,
  or `subprocess.run([...])` with an argument list. No element name, compound
  name or path is interpolated into a shell string any more. Two calls survive
  as `subprocess.run`:
  - `plot.dos_plot`: `sumpdos.sh <element> <orbital>` — a shell script shipped
    with the package, run with `check=False` and a logged return code because a
    missing orbital file is tolerable there;
  - `create_epw_inputs.prepare_nscf`: `python -m htesp.kmesh_nscf n1 n2 n3 [wan]`
    with `check=True`, stdout redirected to the grid file by the `stdout=`
    argument rather than by a shell `>`.
* **Bare `except:` removed everywhere in scope** (`kpath`, `kpoint_path`,
  `plot.dos_plot_vasp`, `plot_bandproj.main`, `fitting_elph_smearing`) —
  replaced by the specific exceptions with a printed message.
* **`isdir`-then-`mkdir`** → `os.makedirs(..., exist_ok=True)` in `kpath`,
  `elph`, `create_epw_inputs`, `fitting_elph_smearing`.
* **Unbound-name guards**: `file_name` in `create_epw_inputs.main` and
  `create_wt_inputs.main`, `data` in `create_wt_inputs.wt_body` and
  `wannier90.epw_bandcheck`, `fermi`/`band_fermi` in `plot.plot`,
  `fileband`/label rows in `plot.band_wann_plot`, the cif glob in
  `create_epw_inputs.main`.
* **`os.chdir`**: none of the files in scope calls it, so no `try/finally`
  wrapping was needed.
* `projection_phband.reading_input` re-parsed the two header fields inside a
  loop over every line of `phonproj.in`, so comment stripping was decided by
  the LAST line rather than by lines 0 and 1. Fixed.
* `plot_bandproj._parse_proj_file` read `self.projdata[12]` for
  `natomwfc nkstot nbnd`. That line's position depends on the number of species
  and atoms in the projwfc.x header, so it was wrong for anything but the cell
  the code was written against — and wrong silently. New `parse_proj_sizes()`
  locates it by shape (three integers followed by the three `T`/`F` flags), and
  also records `lsda`, which is the `nspin = 2` marker in that file. The old
  fixed index is kept as a last resort with a warning.

## Public API changes other modules must know about

| symbol | before | after |
| --- | --- | --- |
| `htesp.kpath.kpath` | `kpath(filename, npoint, kcutoff)` | `kpath(filename, npoint, kcutoff, plot_bz=False, bz_file="BZ.pdf")` — **BZ.pdf is no longer written by default** |
| `htesp.kpath.printk` | `printk()` | `printk(out_dir="scf_dir", plot_bz=True)`, plus optional `sys.argv[6]` = out_dir |
| `htesp.kpoint_path.kpoint_path` | `kpoint_path(file_name)` | `kpoint_path(file_name, out="wannier_kpath.in", slab=False, slab_file="POSCAR-slab", slab_out="wannier_kpath_slab.in")` |
| `htesp.create_wt_inputs.kpoint_path` | local duplicate | **removed** — import from `htesp.kpoint_path` |
| `DataProcessor.combine_orbital` | two defs: `(key, orb, if_gnu)` and `(orb)` | one def: `(orb, key=None, if_gnu=False)` |
| `htesp.generate_submission` | — | new `normalise_calc(which_calc)`; `generate_submission_files` normalises its own `which_calc` |
| `htesp.freq_process.freq_process` | writes `freq.plot` | optional `sys.argv[2]` output name |
| `htesp.standard` | module | moved to `_removed/standard.py` |

New helpers other modules may reuse: `htesp.plot.read_fermi`,
`fermi_from_qe_out`, `fermi_from_outcar`, `read_lsorbit`, `kpt_labels`,
`parse_lambda_out`, `select_smearing`; `htesp.kpath.path_tokens` / `cut_path`;
`htesp.kpoint_path.split_path_segments`;
`htesp.wannier90._as_bool` / `resolve_num_wann`;
`htesp.create_wt_inputs.read_wt_input` / `write_wannier_centres`;
`htesp.elph.read_qpoint_mesh` / `count_irr`;
`htesp.plot_bandproj.parse_proj_sizes`.

New optional config keys (both default to the old behaviour when absent):
`plot.a2f_smearing` (1-based el-ph smearing for the a2f plot) and
`plot.dos_ylim` (DOS-axis limits, separate from the energy `ylim`).

## Not fixed / deliberately left alone

* **Nothing was executed.** `pymatgen`, `ase`, `spglib` and `pytest` cannot be
  installed here (proxy allowlist), so no module in scope could be imported or
  run. Only `py_compile`, an AST undefined-name pass and standalone runs of the
  pure-Python helpers were possible. The numerical claims about array shapes in
  fix (1) come from reading the reshape in `_parse_procar_data` together with
  the existing `[0, j]` / `[1, j]` indexing, not from a run.
* `plot.plot` still reads `../../input.in` by a fixed relative path and
  `inputline[2]` by index for `nkpoint`/`kcut`. Out of scope for the listed
  defects and shared with modules other agents own.
* `htesp/plot.py:'gammaband'` reads `lambda.dat` line 2 by index for `dosef`,
  and `plottype == 'phonband'`/`'gammaband'` read `sys.argv[4]`/`[5]` directly
  rather than through the function arguments. Same family as fix (17) but not
  listed; left unchanged to avoid changing published plots without a request.
* `create_epw_inputs.func_erfc` and `scdmfit` can leave `ydata` / `idx`
  unbound on a bad `entang`/`func_type` value. Pre-existing, not listed, and
  the callers only pass valid values.
* The `nspin == 2` branch of `PROCARProcessor._parse_procar_data` splits the
  PROCAR in half at `len(procar[1:]) / 2`, which assumes both spin blocks have
  exactly equal line counts. That looks fragile, but confirming it needs a real
  PROCAR, so it was left alone.
