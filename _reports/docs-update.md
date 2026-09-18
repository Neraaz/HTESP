# `docs/` update for HTESP 2.0

`python3 -m unittest tests.test_docs` went from **4 failures** to **OK** (7 tests).
The full suite is `Ran 186 tests ... OK (skipped=1)`.

`sphinx-build` could not be verified on this machine: `pip3 install sphinx` is
refused by the proxy allowlist (`403 Forbidden` from the PyPI tunnel), and
`docutils` is not installed either. Verification was therefore done with
`tests/test_docs.py` plus a new `docs/check_docs.py` (below), and by reading the
RST for structural validity.

---

## The four test failures

1. **`test_no_document_uses_a_name_that_does_not_exist`**
   * `otherinput.rst` `projection.in`: `process = epw6-file or epw8-file` →
     ``mainprogram wann-file`` (wannierisation) / ``mainprogram epw-file`` (EPW).
   * `tutorial.rst` distortion section: `distort-extract.py` → `distort_extract.py`.
   * `tutorial.rst` submission section: `generate_submission.sh` →
     `generate_submission_file.sh`.
   * `tutorial.rst` convex-hull section: `econv_vasp.csv` → `econv.csv`
     (and the sentence now says it is written for both codes).
   * `elastic-compute`, `wanniertool_input` and `mainprogram inputinfo` were not
     present in `docs/` — they were in `docs/command.rst`, which someone else had
     already regenerated, and in `examples/VASP/tutorial14/README` (fixed below).

2. **`test_every_json_code_block_parses`** — the flagship block at
   `param.rst:92-473` had a stray `]` after `aflow.prop` and two missing closing
   braces. It is now **generated** from `htesp/data/config.json` by a new
   `docs/gen_param_block.py`, between `.. config-json-start` and
   `.. config-json-end` markers, so it cannot drift again. The `download`
   fragment further down was also unbalanced and is now closed.

3. **`test_no_duplicate_labels`** — `_pressure-label` was defined in
   `otherinput.rst` (the `pressure.in` file) and again in `tutorial.rst` (the
   pressure *campaign*). The tutorial one is now `_pressure-calc-label`, and the
   one `:ref:` that meant the campaign (the equation-of-state section) points at
   it. The `pressure.in` references keep `_pressure-label`.

4. **`test_every_ref_resolves`** — `tutorial.rst` referenced a label `command`
   that never existed; it now points at `command-label`, which `command.rst`
   defines.

---

## File by file

### `docs/usage.rst` — rewritten

* `python setup.py develop` and the two `PYTHONPATH` exports are gone.
  Installation is `pip install .`, development is `pip install -e ".[test,docs]"`.
* Every extra documented in a table: `ml`, `fermisurface`, `oqmd`, `test`,
  `docs`, `all`. Phonopy and cif2cell called out as *not* pip dependencies.
* The Python version conflict is resolved: `>=3.10`, with 3.11 named as the
  development version. The `conda create --name myenv python==3.11.0` pin is
  gone.
* The two contradictory requirements files are resolved: `pyproject.toml` is the
  source of truth; `requirements.txt` and `INSTALL/requirements1.txt` are
  described as pinned reproductions, not install instructions.
* `usage.rst:130`'s `~/src/bash` is gone; the section now documents
  `export PATH="path_to_HTESP/bin:$PATH"`, what the shims are, and that
  `legacy/bash/` holds the originals and is not on `$PATH`.
* New: `export MP_API_KEY=...` and the `~/.config/htesp/credentials` file.
* New: the global-option table (`--workers`, `--dry-run`, `--root`, `--config`,
  `-v`, `--list`, `--version`) and the exit-code table (0/1/2/130).

### `docs/param.rst`

* The ~316-key coverage is preserved — the generated block is the shipped file
  verbatim, so the key set matches by construction.
* `:68-90` English prose marked `code-block:: python` is now prose, with one
  small JSON block and one small Python block.
* New section "Where the file is read from": the working directory plus five
  parents, `$HTESP_CONFIG`, the deep merge over `htesp/data/config.json`, the
  cache, and `mainprogram config-validate`.
* `mpi_key` rewritten around `$MP_API_KEY`, the credentials file, and the
  three-step lookup order. It says plainly to leave the placeholder in the file.
  The top-level `mpi_key` bullet was updated to match.
* `:508-510` `only_init` is now a fourth peer value of `elph_mode`, with the
  two-pass instruction (`only_init` then `parallel_irr`).
* `:635` `which_calc` — states the value is normalised, so `"WANNIER"` and
  `"QE"` work as well as the lower-case spellings, and that anything else is
  rejected with a message.
* Nine `.. code-block:: json` directives holding bare `"key": {...}` fragments
  are now `text`.
* `:1176` `plot_type` → `proj_type`. The whole `bandproj` sub-list was
  re-indented (`colormap` was nested under `proj_type`; it is a peer) and `proj`,
  which was undocumented, is now documented.
* New `plot` keys documented: `a2f_smearing` (1-based smearing index, `null`
  keeps the old behaviour) and `dos_ylim` (DOS height; `ylim` no longer affects
  it). The section now carries a `_plot-label`.
* `:576-580` the `.. _ifermi:` label was inside a bullet list; it is now a
  sub-heading above the `ifermi.json` block, outside the list.
* `:832-838` chemsys and `:944-950` conv_test bullet damage (missing space after
  the hyphen, missing blank line before the list) fixed.
* `batch.header` section rewritten — see "code-side problems" below.

### `docs/tutorial.rst`

* `:613-614` PDOS was documented as `mainprogram 20`, which is `clean-scan`.
  Partial DOS is **18**; the whole band/DOS block is rewritten as one command
  per line with the scan-script name as a comment.
* `:622-623` the duplicated VASP DOS line is gone; VASP writes total and
  projected DOS from the same run, and that is now what it says.
* The range shorthand `8-12`, `13-15`, `16-17` is gone. The commands are written
  out and the text states there is no range syntax and that such an argument
  exits `2`.
* `mainprogram 20` is now described as destructive where it appears.
* `:921` **Magnetic force theorem** — written. What the force theorem is, the
  `"type": "anisotropy"` config, the five steps (`LCHARG`, `LSORBIT` + `ICHARG
  11` in `vasp.in`, `magenum`, submit, `e0`), what `magenum` actually copies and
  writes, and the fact that QE prints `To be implemented`.
* `:1001` **Wannier interpolated bandstructure** — written, as seven steps drawn
  from `examples/QE/tutorial19`: `4`, `epw1`, `wannier90.json` windows from
  `band_stat.csv` and the PDOS, `jobscript` with `which_calc: wannier`,
  `wann-file` / `wann-random` / `wann-scdm`, `pw2wan.in` / `ex.win`, plotting
  with `wann_band`, and `epw5` for the same-k-point comparison. Ends by pointing
  at the WannierTools and EPW pipelines that consume `ex_hr.dat`.
* `:23` the `submission here` claim — corrected, see below.
* RST damage fixed: numbered lists with no space after the number (`:48, 56, 58,
  60, 69, 571, 574, 583, 742, 762, 804`), the three stray leading `. ` bullets
  (`:229, 252, 267`), the Markdown link at `:788`. Three numbered "lists" that
  were really headings (`1. Element Replacement Mode`, `2. Dictionary
  Replacement Mode`, `A. Prepare input files`) became bold headings rather than
  enumerated lists that restart at 2 and at A.
* The "Preparing folder" section was restructured so the code blocks are inside
  their list items instead of terminating the list at each one.
* The `ph-q.in` bullet had a zero-indent code block inside it, which broke the
  list; re-indented.
* `mainprogram 26 : Perform SCF relaxation.` was written as if it were shell;
  now a command with a comment.
* `site_subs.py h` — that file is a module, not a script on `$PATH`; replaced by
  a pointer to the `substitute` config section.
* "look for ifermi-scan script inside bash" → the `ifermi-scan` shim in `bin/`.
* `useful_scripts` → `utility/usefull_scripts` (the spelling on disk).
* Fourteen "Worked example:" cross-links added, one per narrative section,
  naming the matching `examples/QE/tutorialN` and `examples/VASP/tutorialN`.

### `docs/otherinput.rst`

* `projection.in` command names (above).
* `:162-164` missing blank line before the bullet list after "Where:".
* `:259` Markdown-style link for `tot_charge` → RST.
* `mpid.in` pointed its `:ref:` at `pwd-label` (the working-folder listing)
  while claiming to reference `mpid-list.in`; rewritten without the misleading
  reference.

### `docs/utils.rst`

Regenerated from `htesp/`: 51 modules, each as `htesp.<module>` (they are no
longer importable as top-level names). `ifermi_plot`, `ml_processing`, `cli`,
`config`, `inputin`, `help_text`, `banner` and `workflow` added; `standard`
removed (it is in `_removed/`).

### `docs/conf.py`

* `sys.path.insert(0, '../src')` → the repository root, so `import htesp` works.
* `copyright` was "Niraj K Nepal"; now "2024, Iowa State University", matching
  `LICENSE` and `license.rst`, with a comment naming the DOE contract.
* `release` `v1.0` → `2.0.0`, plus `version = '2.0'`.
* `autodoc_mock_imports` added for pymatgen, ase, mp_api, ifermi &c., so the
  API build does not need the scientific stack installed.
* The theme import is now guarded, `exclude_patterns` is no longer empty, and
  `root_doc` is explicit.
* The dead `_ext` extension (commented out at `:16, 36-40, 60-61`) was retired:
  `docs/_ext/edit_on_github.py` → `_removed/docs/_ext/edit_on_github.py`.

### `docs/index.rst`

Rewritten intro (the old one was a single marketing paragraph), and four new
pages added to the toctree: `examples`, `workflow`, `tutorial_runner`,
`testing`.

### `docs/README.md`

Was 15 bytes ("# docu revised"). Now describes every page, how to build, which
two pages are generated and how to regenerate them, and the two check commands.

### `docs/cite.rst`, `contrib.rst`, `license.rst`

`cite.rst` already named the journal version; it now also carries a BibTeX entry
and says explicitly to cite the journal rather than the arXiv preprint.
`contrib.rst` points at the new testing page. `license.rst` gained one sentence
saying the text is `LICENSE` verbatim and is authoritative.

---

## New pages

| page | what it covers |
|---|---|
| `docs/workflow.rst` | `HTESPWorkflow`, the `bin/` shims and `legacy/bash/`, `--workers` and `$HTESP_WORKERS`, per-material scratch directories, deterministic output ordering, `--dry-run` and its two exceptions, job-id capture and what `checkph` now does, failure accounting and the exit codes |
| `docs/tutorial_runner.rst` | the three modes, selection options, how `examples/` is seeded without being written to, `state.json` and `--resume`, exactly what the stop report contains, how a step is judged, the preflight prerequisites, the known limitations, the self-tests. Drawn from `tutorials/README.md` and `_reports/tutorial-runner.md` |
| `docs/examples.rst` | the 42 worked tutorials, with a QE↔VASP mapping table: identical to 10, offset by one from 11 because QE/11 (DFPT el-ph) has no VASP twin, VASP/21 (IFermi) has no QE twin. Plus the four things to do before starting |
| `docs/testing.rst` | `pytest tests/` and `python -m unittest discover -s tests -t .`, a table of what each test file pins, `tools/check_names.py`, `tools/gen_command_rst.py --check`, `docs/gen_param_block.py --check`, `docs/check_docs.py`, the two generated pages, and the conventions — including that a test pinning a fix must name the original symptom |

## New scripts

* **`docs/check_docs.py`** — the documentation lint that needs no Sphinx. Checks
  JSON blocks parse, Python blocks compile, no duplicate labels, every `:ref:`
  and `:doc:` resolves, every toctree entry exists and every page is reachable,
  every `mainprogram <name>` is a real command or process according to
  `htesp/help_text.py`, every repository path in a literal exists, and no
  real-looking API key is present. Currently: `docs OK` over 14 pages, 7 JSON
  blocks, 2 Python blocks, 68 `:ref:`, 11 `:doc:`, 130 command mentions, 90
  paths. Documented in `testing.rst` and `docs/README.md`.
* **`docs/gen_param_block.py`** — regenerates the full `config.json` listing in
  `param.rst` from `htesp/data/config.json`; `--check` for CI.

## `examples/` corrections (the three permitted)

* `examples/VASP/tutorial14/README:53` — `mainprogram elastic-compute` →
  `mainprogram compute-elastic`.
* `examples/VASP/tutorial18/README:9` — "(tutorial12)" → "(tutorial11)": the
  VASP band+DOS tutorial is 11, not 12; 12 was copied from the QE twin.
* `examples/QE/input.in` — added the missing `DFT = QE` line.

---

## Things I could not fix

* **Sphinx could not be installed or run** here (`pip3 install sphinx` →
  `403 Forbidden` from the proxy). `docutils` is not installed either, so the RST
  was not parsed by a real parser. `docs/check_docs.py` and `tests/test_docs.py`
  cover the mechanical part; someone with a working Sphinx should still run
  `cd docs && make html` once and look at the warnings, particularly in the
  `bandproj` and `download` bullet trees of `param.rst`, which were deeply and
  inconsistently indented before this pass.
* **`docs/_ext/` is now an empty directory.** `rm`/`rmdir` are blocked in this
  folder, so the file was moved to `_removed/docs/_ext/edit_on_github.py` but
  the empty directory remains. Git will not track it; delete it when convenient.
* Two `.bak` files created by `sed -i` in `examples/VASP/tutorial14` and
  `tutorial18` were moved to `_removed/examples-README-bak/` rather than
  deleted, for the same reason.
* `docs/command.rst` is generated and I did not touch it. It is current
  (`tools/gen_command_rst.py --check` passes).

## Problems that are in the code, not the docs

1. **`generate_submission_file.sh` and `mainprogram jobscript` disagree about
   `batch.header`.** `HTESPWorkflow.generate_submission_file`
   (`htesp/workflow.py:4174`) still does
   `template.replace("submission here", cmd[name])`, while
   `htesp/generate_submission.py:204-217` copies the header and *appends* the
   command. None of the headers in `examples/` contains a `submission here`
   line, so `generate_submission_file.sh qe-elph mpirun 48` against a shipped
   header silently produces `run-*.sh` scripts that allocate the job and run
   nothing. The error message at `workflow.py:4131` even says the line is
   required — but only when `batch.header` is *missing*, which is the one case
   where it cannot help. Either make the substitution fall back to appending
   when the placeholder is absent, or make it refuse loudly. I documented both
   behaviours and recommended `mainprogram jobscript`, but this should be one
   behaviour, not two.

2. **`--dry-run` does not protect the two destructive processes.**
   `clean_scan` (`htesp/workflow.py:2397`, `_clean_one` at `:2411`) and
   `pressure_reset` (`:3131`, `_pressure_reset_one` at `:3139`) call `remove(...,
   recursive=True)` unconditionally; neither consults `self.dry_run`. Process 20
   deletes wavefunctions and `_ph0/` and moves the run to `completed/`, and
   process 28 removes the whole `pressure/` tree, so `mainprogram 20 --dry-run`
   destroys exactly what the user was checking. These are the two commands where
   a dry run matters most. I removed the `--dry-run` advice from the docs and
   said so explicitly, but the guard belongs in the code.

3. **`magnetic.py` QE path is a stub that exits 0.** `htesp/magnetic.py:217-232`
   prints `To be implemented` for `dft in ("QE", "qe")` and returns normally, so
   `mainprogram magenum` on a QE campaign reports success and writes nothing.
   It should be a non-zero exit or a raised error, so that the tutorial runner
   and any `&&` chain notice.

4. **The force-theorem step is commented out.** `htesp/magnetic.py:205` has
   `#rewrite_incar(..., add=(("ICHARG", 11),))` commented out, so the SAXIS
   directories produced by `magenum` are converged self-consistently rather than
   non-self-consistently on the copied `CHGCAR` — which is what the copied
   `CHGCAR` is there for. I documented the manual workaround (`ICHARG 11` in
   `vasp.in`), but a `magmom.force_theorem` boolean would be the right fix.

5. **`config.json` is still advertised as living in `utility/input_files/`** by
   `htesp/help_text.py` ("Look for config.json file in utility/input_files/"),
   although `htesp/data/config.json` is now the canonical, fully-populated
   default and the two files differ. `command.rst` is generated from that help
   text, so the stale sentence is reproduced in the documentation and I cannot
   fix it from `docs/`. Suggest pointing the help text at
   `htesp/data/config.json` (or at `mainprogram config-validate`).

6. Minor: `docs/usage.rst` used to tell users to add their conda
   `site-packages` to `PYTHONPATH`. That instruction is not just obsolete, it
   actively breaks an editable install. It is gone, but if it appears anywhere in
   `README.md` (which I do not own) it should go there too.
