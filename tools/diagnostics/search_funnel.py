#!/usr/bin/env python
"""Where do the compounds go?  One MP query, then a count after each filter.

`mainprogram search` in element mode applies five filters inside
`element_extract.extract()` and then deletes `download/data-<elm>.csv`, the
pre-filter table -- so a search that returns 465 compounds and writes 2 rows
leaves no evidence of which filter fired.  This reruns the query and prints the
funnel, plus the dtype and NaN count of every column a filter touches (a `None`
from the API becomes NaN, and `NaN < 0` is False, which drops rows silently).

    cd <the project directory holding config.json>
    export MP_API_KEY=<your key>
    python <repo>/tools/diagnostics/search_funnel.py

NOTE: `download()` writes `download/data-<elm>.csv` into the working
directory (that is how the pre-filter table is produced). Run it in a scratch
copy of the project if you do not want that file created.
"""
import warnings

warnings.filterwarnings("ignore")
import _path  # noqa: F401,E402  (sys.path)

from htesp.check_json import config  # noqa: E402
import htesp.element_extract as ee  # noqa: E402


def main():
    cfg = config()
    d = cfg["download"]["element"]
    if cfg["download"].get("mode") != "element":
        print("note: download.mode is %r; this script follows the element-mode "
              "path.\n" % cfg["download"].get("mode"))
    print("config      :", d["elm"], "ntype", d["ntype"], "exclude", d["exclude"])
    print("filters     : metal=%s FE=%s thermo_stable=%s ordering=%r nsites<=%s\n"
          % (d["metal"], d["FE"], d["thermo_stable"], d["ordering"], d["nsites"]))

    data = ee.download(d["elm"][0], d["ntype"], d["exclude"], d["prop"], cfg)
    print("\n%-42s %5d" % ("returned by the MP query", len(data)))

    print("\ncolumn dtypes and how many are unusable (NaN):")
    for col in ("formation_energy_per_atom", "energy_above_hull", "ordering",
                "nsites", "band_gap"):
        if col in data:
            print("   %-42s %-10s NaN=%d"
                  % (col, data[col].dtype, data[col].isna().sum()))
        else:
            print("   %-42s MISSING FROM THE CSV" % col)
    print("\ndistinct ordering values stored:", list(data["ordering"].unique())[:8])

    step = data
    if d["metal"]:
        step = ee.metal_filter(step)
        print("%-42s %5d" % ("after metal_filter (band_gap<=1e-5)", len(step)))
    if d["FE"]:
        step = ee.stable(step)
        print("%-42s %5d" % ("after stable (formation_energy<0)", len(step)))
    if d["thermo_stable"]:
        step = ee.convexhull(step)
        print("%-42s %5d" % ("after convexhull (e_above_hull<0.001)", len(step)))
    step = step[step["ordering"] == d["ordering"]].reset_index(drop=True)
    print("%-42s %5d" % ("after ordering == %r" % d["ordering"], len(step)))
    step = step[step["nsites"] <= d["nsites"]].reset_index(drop=True)
    print("%-42s %5d" % ("after nsites <= %s" % d["nsites"], len(step)))
    if d["spacegroup"]:
        step = step[step["spacegroup"] == d["spacegroup"]].reset_index(drop=True)
        print("%-42s %5d" % ("after spacegroup", len(step)))

    print("\nsurvivors:", list(step["material_id"])[:20])
    print("\nThe first line where the count falls off a cliff is the culprit.")


if __name__ == "__main__":
    main()
