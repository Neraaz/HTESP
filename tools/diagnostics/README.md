# `tools/diagnostics/` — the scripts behind the database-search investigation

Four small programs kept from the September 2026 investigation into
`mainprogram search` returning 465 compounds and writing 2 rows to
`mpid-list.in`. They are diagnostics, not tests: nothing here runs in CI, and
none of them modify a project.

**Interpreter.** These need the scientific stack, so use the `htesp`
environment, not the base one:

```bash
PY=/work/08910/nnepal/vista/anaconda3/envs/htesp/bin/python3
```

On this machine `python3` resolves to the base anaconda, which has no
`pymatgen`/`emmet-core` and will fail with `ModuleNotFoundError`.

| script | needs | what it answers |
|---|---|---|
| `search_funnel.py` | an MP API key | which filter drops the compounds |
| `mp_field_types.py` | emmet-core only | how MP ids and `ordering` stringify |
| `search_offline_stub.py` | no key | does the filter chain behave, on a laptop |
| `kmesh_symprec.py` | no key | does the symmetry tolerance move the k-mesh |

## `search_funnel.py`

`element_extract.extract()` applies five filters and then deletes
`download/data-<elm>.csv`, the pre-filter table — so a collapsed search leaves
no evidence of which filter fired. This reruns the query and prints the count
after each filter, plus the dtype and NaN count of every column a filter
touches. A `None` from the API becomes `NaN`, and `NaN < 0` is `False`, so rows
disappear silently.

```bash
cd <project directory with config.json>
export MP_API_KEY=<your key>
$PY <repo>/tools/diagnostics/search_funnel.py
```

The first line where the count falls off a cliff is the culprit. It writes
`download/data-<elm>.csv` into the working directory, so run it in a scratch
copy if you would rather not create that file.

**Answer from the September 2026 run** (`elm=['B']`, `ntype=[1,2]`): the cliff
is `ordering == 'NM'`. Materials Project now returns `ordering = 'Unknown'`
for 202 of the 465 hits (158 of the 222 surviving the metal and formation-energy
filters), where the 2024 reference data had none. This is a change in MP's data,
not a code defect — 1.x produces the same counts today. Note also that the
`ordering` filter had no on/off switch; it now accepts `null` (no filter) and
a list such as `["NM", "Unknown"]`. On that search, `"NM"` keeps 27
compounds, `["NM", "Unknown"]` keeps 123 and `null` keeps 129.

## `mp_field_types.py`

Two `emmet-core` types do not stringify to what the filters and `config.json`
expect. Verified output:

```
AlphaID(149, padlen=6, prefix='mp')   str()=mp-aaaaft   .string=mp-149
MPID('mp-763')                        str()=mp-763      .string=mp-763
Ordering.NM                           str()='Ordering.NM'  value='NM'  == 'NM' -> False
```

* **Material ids are migrating from `MPID` to `AlphaID`.** `str()` gives the new
  alphabetic spelling, `.string` the legacy `mp-<int>` one. 1.x wrote
  `…['material_id'].string`; the 2.0 rewrite dropped the `.string`, so
  `mpid-list.in` came back as `v1 mp-paf Yb1B2` instead of `v1 mp-10145 Yb1B2`,
  and every `R<mpid>-<compound>/` directory named from it changed with it.
  `htesp.element_extract.legacy_mpid()` restores the legacy spelling -- but
  note `doc.dict()['material_id']` is a plain `str` already in the alphabetic
  form, so the *attribute* `doc.material_id` must be read instead.
  **This one was a real regression, confirmed against a user's own run.**
* **`Ordering` is a plain `Enum`**, so `Ordering.NM != 'NM'`. This *looks* like
  it would break the `ordering` filter, and it does with a hand-built stub —
  but a live summary search returns `ordering` as a plain string, so the trap
  does not fire on the real path. `plain_value()` is a guard, not a fix.
  Trust a real run's `download.csv` over any reasoning about this.

## `search_offline_stub.py`

Drives `download()` and `extract()` against a synthetic 465-document result set
so the CSV writing and the filter chain can be exercised with no API key.
`--alpha-ids` makes the stub return `AlphaID`s the way current MP does, which is
how the `mp-paf` spelling was reproduced.

## `kmesh_symprec.py`

1.x symmetrised in `getkpt` with `symprec=0.1` + `international_monoclinic=True`
while `setting_qeinput` used `0.01` + `False` — the mesh was computed for a
different cell from the one written. 2.0 unifies both on `(0.01, False)`.

The useful finding is negative: across a near-cubic cell and a C2/m monoclinic
one, the **cell** differs visibly (space group, monoclinic setting) while the
**mesh is identical**, because `pos_to_kpt` depends only on reciprocal lattice
vector lengths. Run this before claiming a k-mesh changed.
