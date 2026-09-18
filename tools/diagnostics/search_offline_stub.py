#!/usr/bin/env python
"""Run the element-mode search end to end with a stub MPRester -- no API key.

Feeds `element_extract.download()`/`extract()` a synthetic result set so the
CSV writing and the whole filter chain can be exercised on a laptop.  Useful
for "did my change to the filters do what I think?" without spending a query.

    python <repo>/tools/diagnostics/search_offline_stub.py [--alpha-ids]

--alpha-ids makes the stub hand back emmet AlphaID material ids, the way
current Materials Project does, which is how the `mp-paf` / `mp-qd` spelling in
mpid-list.in was reproduced.

Writes into a temporary directory and removes it afterwards.
"""
import argparse
import json
import random
import shutil
import tempfile
import os
import warnings

warnings.filterwarnings("ignore")
import _path  # noqa: E402
from _path import REPO_ROOT  # noqa: E402

from pymatgen.core import Structure, Lattice  # noqa: E402
import htesp.element_extract as ee  # noqa: E402

N_DOCS = 465


def build_stub(alpha_ids):
    from emmet.core.electronic_structure import Ordering

    class Doc:
        """Mimics the MPDataDoc a summary search yields."""

        def __init__(self, i, nsites, ordering):
            a = 3.0 + 0.1 * i
            self.structure = Structure(Lattice.cubic(a), ["Mg", "B"],
                                       [[0, 0, 0], [.5, .5, .5]])
            if alpha_ids:
                from emmet.core.mpid import AlphaID
                mpid = AlphaID(700 + i, padlen=6, prefix="mp")
            else:
                mpid = "mp-%d" % (700 + i)
            self._d = {
                "material_id": mpid, "formula_pretty": "MgB2",
                "structure": self.structure,
                "formation_energy_per_atom": -0.3, "band_gap": 0.0,
                "energy_above_hull": 0.0, "total_magnetization": 0.0,
                "ordering": ordering,
                "total_magnetization_normalized_formula_units": 0.0,
                "num_magnetic_sites": 0, "theoretical": False, "nsites": nsites,
            }

        def dict(self):
            return dict(self._d)

    class Summary:
        def search(self, **kw):
            random.seed(0)
            return [Doc(i, random.choice([2, 4, 8, 12, 40]),
                        random.choice([Ordering.NM, Ordering.FM, Ordering.AFM]))
                    for i in range(N_DOCS)]

    class Materials:
        summary = Summary()

    class StubMPR:
        materials = Materials()

    return StubMPR()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha-ids", action="store_true",
                        help="hand back AlphaID material ids, as MP now does")
    args = parser.parse_args()

    stub = build_stub(args.alpha_ids)
    ee.mprester = lambda key: stub
    ee.require_api_key = lambda cfg=None: "stub"

    with open(REPO_ROOT / "htesp" / "data" / "config.json") as handle:
        cfg = json.load(handle)
    d = cfg["download"]["element"]

    previous = os.getcwd()
    workdir = tempfile.mkdtemp(prefix="htesp-stub-")
    try:
        os.chdir(workdir)
        data = ee.download(d["elm"][0], d["ntype"], d["exclude"], d["prop"], cfg)
        print("\npre-filter rows      :", len(data))
        print("material_id sample   :", list(data["material_id"][:4]))
        print("ordering values      :", list(data["ordering"].unique()))
        print("nsites <= %-2s keeps    : %d" % (d["nsites"],
                                                (data["nsites"] <= d["nsites"]).sum()))
        print("ordering == %-4r keeps: %d"
              % (d["ordering"], (data["ordering"] == d["ordering"]).sum()))
        kept = ee.extract(ntype=d["ntype"], properties=d["prop"], elm=d["elm"],
                          exclude_el=d["exclude"], nelm=1, metal=d["metal"],
                          neg_fe=d["FE"], thermo_stable=d["thermo_stable"],
                          ordering=d["ordering"], nsites=d["nsites"],
                          spacegroup=d["spacegroup"], input_data=cfg)
        print("surviving extract()  :", len(kept))
    finally:
        os.chdir(previous)
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
