# HTESP test suite

```bash
# with pytest (what CI uses)
pip install -e ".[test]"
pytest tests/ -q

# without pytest -- standard library only
python -m unittest discover -s tests -t . -v
```

The tests are `unittest` classes so they run under both runners. That matters
because most of them need nothing but the standard library: `mainprogram
--help`, the configuration loader, `input.in`, the tracking-file reader, the
Quantum ESPRESSO text handling and the whole command dispatcher are checked on
a machine with no `pymatgen`, no `ase` and no cluster. Anything that genuinely
needs a scientific package is skipped with a message naming it, rather than
erroring.

| file | what it pins |
|---|---|
| `test_config.py` | the search order, the deep merge over the packaged default, `MP_API_KEY` precedence, `config-validate` |
| `test_inputin.py` | `input.in` parsing, including the four ways the original parser failed on a short or `DFT`-less file |
| `test_qetext.py` | the QE reader: `vc-relax`, `relax` and timed-out outputs, cards, namelist edits |
| `test_workflow_core.py` | tracking files (`end` is exclusive), material paths, mesh files, `pushd`, job-id capture, failure accounting |
| `test_cli.py` | every command's dispatch, the plot-type loop, exit codes, the help blocks |
| `test_imports.py` | every module imports; nothing reads a file at import time |
| `test_packaging.py` | shims, entry points, no committed API key, no `os.system`, no bare `except:`, no `eval` |
| `test_regressions.py` | one test per defect from the 2026 review, each naming the original symptom |
| `test_docs.py` | the documentation cannot name a command that does not exist, cannot tell the reader to run the destructive process 20 for PDOS, and its JSON blocks must parse |
| `test_tutorials.py` | pulls in `tutorials/selftest.py` and ties the tutorial catalogue to the dispatcher |
| `test_portability.py` | the package stays architecture-neutral: no `platform.machine()` branch, no compiled file in `htesp/`, no `ctypes`, POSIX `sh` shims, and `htesp-check` reports a dependency that aborts instead of dying with it |

## Conventions

* Every test that pins a fix says, in its docstring or name, what the original
  symptom was. A test that only asserts current behaviour tells a later reader
  nothing about why it exists.
* `tests/helpers.py` has `TempProject`, which gives each test a throw-away
  project directory, makes it current, and clears `$HTESP_CONFIG`,
  `$MP_API_KEY` and the configuration cache — a stale value in a developer's
  shell would otherwise defeat the configuration tests.
* Negative source assertions go through `code_only()`, which strips comments
  and strings: every fix carries a `# FIX(n):` comment naming the old pattern,
  so a naive "the old pattern is gone" check would match the comment.

## Also worth running

```bash
python tools/check_names.py htesp tutorials tests    # undefined-name scan
htesp-check                                         # this machine's wheels (x86_64 / arm64)
python tools/gen_command_rst.py --check              # docs/command.rst is current
python -m unittest tutorials.selftest                # the tutorial runner alone
```
