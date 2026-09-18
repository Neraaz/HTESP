# Database / input-generation / QE-helper fixes

Scope: `htesp/element_extract.py`, `oqmd_extract.py`, `aflow_extract.py`,
`ml_processing.py`, `pymatgen_phase_diagram.py`, `qe_input.py`,
`cif_to_gsinput.py`, `vasp_input.py`, `write_potcar.py`, `crystal.py`,
`scftocif.py`, `checkfreq.py`, `lambda_in.py`, `q2r.py`, `matdyn.py`,
`matdyn_dos.py`, `band.py`, `dos.py`, `phonband.py`, `kmesh_nscf.py`,
`qe_axsf2cellpos.py`, `ifermi_plot.py`.

Every fix carries a `# FIX(<n>): ...` comment at the site.
`python3 -m py_compile htesp/*.py` passes.

---

## element_extract.py

**FIX(1) — `input_data` NameError on import.** `download()`, `create_input()`,
`extract()`, `download_by_entry()` and `main()` all read a module global that
was only bound inside `if __name__ == "__main__":`, so every one of them raised
`NameError` when the module was imported. Each now takes an optional
`input_data=None` argument and calls `config()` when it is omitted; `main()`
loads the configuration once at the top and threads it down. The `__main__`
block no longer assigns anything. The API-key lookup in `download()` now goes
through `require_api_key()` (`$MP_API_KEY` → `~/.config/htesp/credentials` →
`config.json`), replacing a branch that left `key` unbound when no
`config.json` was next to the working directory.

**FIX(2) — `eval(must_in)`.** `download_by_entry()` assembled a Python source
string out of `config.json` values (`"nelm < ntype_constraint and ('Mg' in
elm_list or 'B' in elm_list)"`) and `eval()`d it once per candidate entry, i.e.
config content was executed. Replaced by the explicit predicate
`nelm < ntype_constraint and any(e in elm_list for e in must_include)`.
Behaviour note: an **empty** `must_include` now selects nothing, where the old
code raised `SyntaxError` on `eval("... and ()")`. Every non-empty
`must_include` behaves identically.

**FIX(3) — `search['material_id'].string`.** An `MPDataDoc` is not
subscriptable. Made consistent with the generic branch two lines below:
`search.dict()[propty]`.

**FIX(4) — `summary.available_fields[:-29]`: NOT PRESENT in this file.**
`grep -n "available_fields" htesp/element_extract.py` returns nothing. The
positional slice exists only in `htesp/htepc.py`, which another agent owns —
flagged there, not fixed here.

Also in this file: `os.system("rm -r download_old")`, `os.system("mv download
download_old")` (×2 sites), `os.system("rm download/data*")` and
`os.system('echo "NA" > remove.list')` replaced with `shutil.rmtree`,
`shutil.move`, `glob`+`os.remove` and a plain file write; `isdir`-then-`mkdir`
replaced with `os.makedirs(..., exist_ok=True)`; the "up to 2 elements"
branch now raises `ValueError` instead of falling through with `data` unbound.

## oqmd_extract.py

**FIX(5) — infinite retry loop / unbound `data`.** `search()` decremented
`limit` by `int(0.1*limit)`, which is 0 for `limit < 10` (an infinite loop),
and on exhaustion fell out of the loop with `data` never assigned. The step is
now `max(1, int(0.1*orig_limit))`, the attempts are bounded by
`MAX_SEARCH_RETRIES = 8` with linear `RETRY_BACKOFF_SECONDS` backoff, the bare
`except:` is `except Exception as exc` with the reason printed, and exhaustion
raises `RuntimeError` naming the limits and the last error.

**FIX(6) — POSCAR species/position mismatch.** The counts line came from
`Counter(elements)` (grouped by species) while the positions were written in
the original OQMD site order. OQMD does not guarantee grouped sites, so for
sites `[Mg, B, Mg]` VASP read `Mg B / 2 1` against positions `(Mg, B, Mg)` and
gave the second Mg position to boron. Sites are now emitted grouped by species
in first-appearance order, so the counts and positions describe the same
structure. A regression comment marks it. `Counter` is no longer imported.

**FIX(7) — malformed `must_include` / `element_set`.** The loop's
`j < len(...) - 1 and j > 0` test never fired for `j == 0`, so no comma
followed the first element and two species were glued into one token:
`['Mg','B','C']` produced `(A-B),MgB,C`. Now `",".join(...)`.

**FIX(8) — KPOINTS computed for the wrong cell.** The k-mesh was computed from
the raw POSCAR three statements before `MPRelaxSet` overwrote POSCAR with the
standardised primitive cell, so the shipped KPOINTS described a different cell
from the shipped POSCAR. `pos_to_kpt()` is now called *after*
`relax_set.poscar.write_file("POSCAR")` in the VASP branch (and on the raw cell
in the QE branch, which keeps that cell).

**FIX(9) — driver moved out of `__main__`.** The whole driver (reading
`input.in` by hand, reading `config.json`, building `KWARGS`) is now
`build_kwargs(input_data=None)` + `main(argv=None)`; `__main__` is
`sys.exit(main())`. `input.in` is parsed with `htesp.inputin.InputIn.load` and
the configuration with `htesp.check_json.config()`. The `search` / `download`
argument contract is unchanged. The old "config.json not found" fallback
branch was dead and internally broken (it referenced an undefined `LIMIT`); it
is gone, because `config()` now always returns a fully populated dict.

Also: `os.system` `mkdir`/`mv`/`rm`/`vasp_process.py` replaced with
`os.makedirs(exist_ok=True)`, `shutil.move`, `shutil.copy`, `os.replace`,
`os.remove` and `subprocess.run([sys.executable, "-m", "htesp.vasp_process",
"POSCAR"], cwd=relax_dir, check=False)` (return code logged) — this also
removes the `os.chdir()` pair; bare `except:` in `download()` named and logged;
registry append replaced by `register_mpid()` is **not** applied here (the
per-material `mpid.in` write in `download()` still uses `register_mpid`, see
FIX(18)).

## aflow_extract.py

**FIX(10) — infinite busy loop with no timeout.** Both `urlopen()` call sites
(`search_data()` and `download()`) were `while not success: try: ... except:
continue`. Replaced by one `fetch_aflux(url)` helper with
`REQUEST_TIMEOUT = 60 s`, `MAX_REQUEST_RETRIES = 5`, linear backoff, named
exceptions (`HTTPError`, `URLError`, `TimeoutError`, `OSError`, `ValueError`)
logged per attempt, and a `RuntimeError` naming the URL on exhaustion.

**FIX(11) — `ast.literal_eval` on JSON.** AFLUX returns JSON, so any record
containing `null`/`true`/`false` raised `ValueError`. `fetch_aflux` uses
`json.loads`. The `ast` import is gone.

**FIX(12) — AFLUX range separator.** Inside a property's parentheses AFLUX
reads `,` as OR and `:` as AND (`!` = NOT, `*` = loose/wildcard). `natoms(1*,*4)`
therefore means "≥1 **or** ≤4", true for every entry, so every range filter
matched everything. Changed to `:` in all three range filters —
`nspecies(1*:*N)` (both the filtered and the unfiltered branch),
`natoms(1*:*N)`, `enthalpy_formation_atom(-100*:*0)` — with a comment citing
the operator semantics. `prop_string` is also initialised, so a config without
`prop` no longer raises `NameError` on `criteria + prop_string`.

Also: `os.system("mv ...")` ×2 → `os.replace`; `isdir`-then-`mkdir` →
`os.makedirs(exist_ok=True)`; bare `except:` in `download()` named and logged
with the structure index; `main(argv=None)` returns an exit code, parses
`input.in` with `InputIn`, and the registry append uses `register_mpid()`.

## ml_processing.py

**FIX(13) — swapped labels.** `class_accuracy()` formatted `(acc, precision,
recall)` against the labels `accuracy_score / recall_score / precision_score`,
so every reported recall was the precision and vice versa. Arguments reordered.

**FIX(14) — `np.abs()` around the Allen-Dynes denominator.** `ml_tc()` wrapped
`lambda - mu*(1 + 0.62 lambda)` in `np.abs()`, mirroring the pole at
`lambda ≈ 0.178` and manufacturing a finite, positive Tc for every lambda
below it, where Allen-Dynes has no superconducting solution. Extracted into a
new public `allen_dynes(wlog, lam, mustar=MUSTAR)` that uses the denominator
as-is and returns `NaN` where it is non-positive, with a `RuntimeWarning`
naming how many entries were affected. Verified numerically: values for
lambda > 0.178 are bit-identical to before; 0.05/0.15/0.1776 now give NaN
instead of 0.031/0.0/0.0.

**FIX(15) — optional imports at module scope.** `matminer` and `scikit-learn`
were imported at the top, so `import htesp.ml_processing` failed with
`ImportError` on a core-only install. They are now imported at the point of
use through a new `_require(module_name, extra="ml")` helper that raises
`ImportError("... pip install htesp[ml]")`. The module-scope
`json.load(config.json)` block (which left `input_data` unbound on
`FileNotFoundError`) is replaced by `config()` / `require_api_key()` inside
`jarvis_structure_feature()`.

Also: `os.system("sed -i '1d' <outfile>")` in `write_id_prop_csv()` replaced by
a pure-Python header strip.

## pymatgen_phase_diagram.py

**FIX(16) — CSV header/row column mismatch.** The header named six columns
(`...,e_above_hull_calc,decomposition`) while every row wrote five, because the
`get_decomposition()` write is commented out; pandas then read a shifted/NaN
last column. The header now names exactly the five fields written, with a
comment saying to re-enable both together.

**FIX(17) — `econv.csv` name: VERIFIED, no mismatch.** `grep -rn "econv"
htesp/ legacy/` shows `htesp/workflow.py:3330` writing `econv.csv` with header
`ID,comp,NIONS,energy,niteration`, and `legacy/bash/phonopy-scan` writing the
same name and header. `plot_phase()` reads `data.ID`, `data.comp` and
`data.energy`, all of which that header provides. Nothing to change.

Also: a material whose `relax/` has neither `POSCAR` nor `scf.in` now
`continue`s instead of falling through to `len(struc)` with `struc` unbound or
stale from the previous iteration.

## qe_input.py / cif_to_gsinput.py / vasp_input.py

**FIX(18) — non-atomic, non-idempotent `mpid.in` registry.** The
read-then-append pattern (`lines = open('mpid.in').readlines()`, then append
`v<len(lines)+1>`) appeared verbatim in `qe_input.py`, `vasp_input.py`,
`cif_to_gsinput.py`, `oqmd_extract.py` and `aflow_extract.py`. Two materials
written concurrently, or any hand-edited file, produced duplicate or skipped
`v<N>` numbers — and every consumer resolves a material with `grep "v$ii "`.
Three new functions in `cif_to_gsinput.py`:

* `read_mpid_entries(path="mpid.in") -> [(mpid, compound), ...]`
* `find_mpid(mpid, path="mpid.in") -> int | None` (1-based index)
* `register_mpid(mpid, compound, path="mpid.in") -> int`

`register_mpid` skips an mpid already present (returning its existing index),
renumbers `v<N>` densely from 1 on every write, and is atomic — the whole file
is written to a `tempfile.mkstemp` in the same directory and `os.replace`d.
Format is unchanged: `v<N> <mpid> <compound>` per line. All five call sites now
use it. In `qe_input.py` the "download only if new" behaviour is preserved via
`if find_mpid(mpid) is None:` (the old substring test `any(mpid in line ...)`
would also match `mp-12` inside `mp-123`; the new check is an exact field
match). Verified by test: three appends → v1 v2 v3; re-registering an existing
mpid is a no-op returning its index; a file hand-edited to `v1, v7` is
renumbered to `v1, v2, v3` on the next append.

**FIX(19) — duplicated `pos_to_kpt`.** NOT merged, as instructed. A module
docstring note was added at the top of `htesp/cif_to_gsinput.py` (and a
`.. note::` on the function) naming `htesp.htepc.pos_to_kpt`. **The reciprocal
note in `htesp/htepc.py` was NOT added — that file belongs to another agent.**
The differences, for whoever adds it:

| | `cif_to_gsinput.pos_to_kpt` | `htepc.pos_to_kpt` |
|---|---|---|
| signature | `(file, density, evenkpt=False)` | `(file, density)` |
| even mesh | rounds every odd division up | no |
| side effect | **writes `KPOINTS`** in cwd | none |
| returns | `kmesh` | `kmesh` |
| quirk | — | computes `kratio`/`klat`/`kmesh` twice (identical copy-paste; harmless) |

The numerical core is otherwise identical. Callers needing `KPOINTS` on disk
or an even mesh must use the `cif_to_gsinput` one.

**FIX(20) — `main()` guard and argument.** `cif_to_gsinput.py` did already have
`if __name__ == '__main__': main()`. What it lacked was an argument: the
calculation type was read from `sys.argv[1]` in the middle of `main()`. It is
now `main(calc_type=None, argv=None)` and falls back to the command line only
when `calc_type` is omitted, so `cif_to_gsinput.py <QE|VASP>` still works.

Also in these three files: `ciftoscf()`'s `try/except FileNotFoundError`
around `glob.glob` (which returns `[]` rather than raising, leaving
`file_path` unbound) now raises a real `FileNotFoundError`; `cif2cell` is run
with `subprocess.run([...], check=True)` instead of a shell string containing
the CIF path; the `sed '5 a ...' POSCAR > POSCAR_new; mv POSCAR_new POSCAR`
pair is an in-memory `lines.insert(5, ...)`; a cif2cell header with no
`order:` field raises instead of using a stale `index`; `mkdir`/`mv`/`cp`/`rm`
`os.system` calls became `os.makedirs(exist_ok=True)` / `shutil.move` /
`shutil.copy` / `os.remove`; `vasp_process.py POSCAR` became
`subprocess.run([sys.executable, "-m", "htesp.vasp_process", "POSCAR"],
cwd=relax_dir, check=False)` with the return code logged, removing the
`os.chdir` pair; the one remaining `os.chdir` (in `vasp_input.py`, needed
because `poscar2potcar()` reads `POSCAR` from the working directory) is now in
a `try/finally`; the `if os.path.isfile("config.json")` guards that left `d`
unbound are gone, since `config()` always returns a populated dict.

## Small QE helpers

**FIX(21) — callable top-level entry functions.** `htesp.workflow.run_helper`
imports the module and calls a named function with `sys.argv` installed.
Audited all fourteen; **every entry name the workflow layer uses already
existed at module level** — nothing had to be rescued out of a `__main__`
block. Verified against `grep -n "run_helper(" htesp/workflow.py`:
`band.main`, `phonband.phonband_in`, `dos.dos_in`, `q2r.q2r_in`,
`matdyn.matdyn_in`, `matdyn_dos.matdyn_dos`, `checkfreq.main`,
`lambda_in.main`, `crystal.main`, `scftocif.main`, `vasp_input.main`,
`qe_axsf2cellpos.main`.

What changed: each entry function now also accepts its values as real
arguments (falling back to `sys.argv` when they are omitted, so the CLI and
`run_helper` contracts are untouched), and modules whose entry name is not
`main` gained a `main = <entry>` alias so they can all be driven uniformly.
`kmesh_nscf.main` and `qe_axsf2cellpos.main` no longer call `sys.exit()` on the
success path (that would have killed an in-process caller); they return an exit
code and `__main__` does `sys.exit(main())`.

**FIX(22) — hard-coded scratch-file names.** Optional final parameters added,
with an optional extra `sys.argv` slot, defaulting to today's names, so a
caller can point them at a per-material file. Default behaviour unchanged.

* `matdyn.matdyn_in(mpid, compound, prefix2, mass_file=None)` — `mass.dat`
  (module constant `DEFAULT_MASS_FILE`)
* `matdyn_dos.matdyn_dos(mpid, compound, prefix2, mass_file=None)` — `mass.dat`
* `checkfreq.check_freq(filename, flag_file="freq.dat")` and
  `checkfreq.main(file_name, flag_file, argv)` — the soft-mode marker file

Not parameterised: `matdyn.py`'s `scf_dir/kpathlines.dat` (a path, not one of
the four named scratch files) and `lambda_in.py`'s `<compound>.dyn0` (already
per-material). `qpoint.dat` / `kpoint.dat` are read by `elph.py` and
`create_epw_inputs.py`, which are not in this file set.

**FIX(23) — optional third-party imports.**

* `ifermi_plot.py`: it turns out **nothing** third-party was imported here —
  the module only assembles the `ifermi ...` command line from `ifermi.json`.
  Added `check_ifermi()`, which looks the executable up with `shutil.which` at
  the point of use and raises
  `"pip install htesp[fermisurface] (or: pip install ifermi)"`, plus a
  docstring saying so. Also fixed: `command_params` was left unbound when
  `ifermi.json` was missing, so the loop below raised `NameError` instead of
  falling back to a bare command. Added `main(command, input_file, argv)`.
* `write_potcar.py`: its only third-party import is `pymatgen`, a core
  dependency imported by eight other modules — left at module scope, with a
  docstring note. Fixed instead: `pot1` was left unbound when no `config.json`
  sat beside the working directory (`NameError` one line later); now a missing
  `pseudo.pot` section, or a missing element in it, raises a `KeyError` that
  names what is missing. `poscar2potcar()` gained optional
  `poscar="POSCAR"`/`outfile="POTCAR"` arguments (all eight existing callers
  call it with no arguments).

Other helper fixes: `crystal.crystal_extract()` and `scftocif.scf_tocif()`
raise `ValueError` instead of printing and then using an unbound `cell_inp` /
a still-string `filename`; `scf_tocif` now compares the *basename* (and
accepts `CONTCAR` and any `*.in`), so `relax/scf.in` no longer falls through;
`crystal.py` dropped its unused `spglib`, `ase.io.espresso`, `ase.io.vasp` and
`ase.cell.Cell` imports (`spglib` cannot be installed in this environment, so
this also makes the module importable here); `lambda_in.py`'s three
`os.system` calls — a guarded `rm touch_list.txt` that never matched the file
it meant (the scratch file is `elph_list.txt`), a `touch`, and one
`grep 'Number of q in the star' <file> >> elph_list.txt` per `elph.out*` with
the filename interpolated into a shell string — became an in-memory scan; the
`elph_list.txt` scratch file is no longer written at all (verified unread by
the rest of the package and by `legacy/bash/`), so the trailing
`os.system("rm elph_list.txt")` is gone too; an empty `elph.out*` glob now
raises instead of `IndexError` on `.pop(0)`.

---

## Could not fix / deliberately out of scope

1. **FIX(4)** — the `summary.available_fields[:-29]` positional slice is not in
   `element_extract.py`; it lives in `htesp/htepc.py`, owned by another agent.
2. **FIX(19)** — the reciprocal docstring note in `htesp/htepc.py` was not
   added, for the same reason. The table above has what it needs to say.
3. **`setup.py` extras.** FIX(15) and FIX(23) point users at
   `pip install htesp[ml]` and `htesp[fermisurface]`, but `setup.py` declares
   no `extras_require` (it lists `ifermi` as a hard dependency and does not
   list `matminer` or `scikit-learn` at all). `setup.py` is outside this file
   set; someone should add
   `extras_require={"ml": ["matminer", "scikit-learn"], "fermisurface": ["ifermi"]}`
   and drop `ifermi` from `install_requires`.
4. **Nothing was executed.** `pymatgen`, `ase`, `spglib`, `mp_api`, `qmpy_rester`,
   `matminer` and `pytest` cannot be installed here (proxy allowlist), so
   verification was `python3 -m py_compile htesp/*.py` (passes), an AST check
   that every entry function is a module-level callable (passes), and direct
   execution of the dependency-free logic that was rewritten: `register_mpid`
   (dense/idempotent/atomic), the FIX(6) POSCAR grouping, and
   `allen_dynes` against the old `np.abs` expression.
5. **`ml_processing.ml_tc()` pre-existing type bug (NOT in the fix list, NOT
   fixed).** Its parameters are documented and defaulted as *file names*
   (`lam='lambda.csv'`), but the body immediately calls `lam.set_index('ID')`
   and `wlog['ID']`, i.e. it requires DataFrames. Calling it with its own
   defaults raises `AttributeError`. Flagging only.
