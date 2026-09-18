# `docs/` — the HTESP documentation source

Sphinx sources for <https://neraaz.github.io/HTESP/>.

```bash
pip install -e ".[docs]"     # sphinx + the theme
cd docs && make html         # build/html/index.html
```

## What is here

| file | what it is |
|---|---|
| `index.rst` | the master `toctree` |
| `usage.rst` | installation, global options, exit codes |
| `param.rst` | the `config.json` reference — every key |
| `otherinput.rst` | the `.in` files: `input.in`, `vasp.in`, `kpoint.in`, … |
| `command.rst` | **generated** — every command and process |
| `tutorial.rst` | the narrative: one section per kind of campaign |
| `examples.rst` | the 42 worked tutorials in `examples/`, and their numbering |
| `workflow.rst` | `HTESPWorkflow`, `--workers`, `--dry-run`, job-id capture |
| `tutorial_runner.rst` | the `htesp-tutorials` driver |
| `testing.rst` | running the suite, the tooling, the conventions |
| `license.rst`, `contrib.rst`, `cite.rst`, `utils.rst` | short pages |
| `conf.py` | Sphinx configuration |
| `check_docs.py` | documentation lint that needs no Sphinx |
| `gen_param_block.py` | regenerates the `config.json` listing in `param.rst` |

## Two pages are generated — do not hand-edit them

* `command.rst` comes from `htesp/help_text.py`, which is also what
  `mainprogram` prints:

  ```bash
  python tools/gen_command_rst.py            # rewrite
  python tools/gen_command_rst.py --check    # CI: stale -> exit 1
  ```

* the full `config.json` listing in `param.rst`, between the
  `.. config-json-start` and `.. config-json-end` comments, comes from
  `htesp/data/config.json`:

  ```bash
  python docs/gen_param_block.py
  python docs/gen_param_block.py --check
  ```

## Checks

```bash
python docs/check_docs.py              # refs, labels, code blocks, names, paths
python -m unittest tests.test_docs     # the hard gate, run by the test suite
```

`tests/test_docs.py` is the specification: it fails if a document names a
command that does not exist, tells the reader to run the destructive process 20
for the partial DOS, contains a JSON block that does not parse, defines a label
twice, or points a `:ref:` at a label that is not there.
