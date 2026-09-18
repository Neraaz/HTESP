#!/usr/bin/env python
"""What the Materials Project API actually hands back, and how it stringifies.

Two emmet-core types do not stringify to what HTESP's filters and config.json
expect.  This prints both, so the question "is it still doing that?" is one
command rather than an afternoon.  Needs only emmet-core -- no API key, no
network.

    python <repo>/tools/diagnostics/mp_field_types.py
"""
import warnings

warnings.filterwarnings("ignore")
import _path  # noqa: F401,E402


def main():
    print("=== material ids: MPID vs AlphaID ===")
    from emmet.core.mpid import MPID, AlphaID

    alpha = AlphaID(149, padlen=6, prefix="mp")
    rows = [("AlphaID(149, padlen=6, prefix='mp')", alpha),
            ("AlphaID('mp-aaaaft')", AlphaID("mp-aaaaft")),
            ("MPID('mp-763')", MPID("mp-763"))]
    for label, value in rows:
        print("  %-38s str()=%-12s .string=%-12s"
              % (label, str(value), getattr(value, "string", "-")))
    print("  -> str() gives the NEW alphabetic spelling; .string the legacy one.")
    print("     htesp.element_extract.legacy_mpid() prefers .string, which is")
    print("     what 1.x wrote and what R<mpid>-<compound>/ directories use.\n")

    try:
        from htesp.element_extract import legacy_mpid, plain_value
    except ImportError as exc:                       # pragma: no cover
        print("  (could not import the helpers: %s)" % exc)
        return
    for _, value in rows:
        print("  legacy_mpid(%-12s) = %s" % (str(value), legacy_mpid(value)))

    print("\n=== magnetic ordering ===")
    try:
        from emmet.core.electronic_structure import Ordering
    except ImportError:
        print("  Ordering moved; skipping")
        return
    print("  base class:", Ordering.__mro__[1].__name__)
    for name in ("NM", "FM", "AFM"):
        member = getattr(Ordering, name)
        print("  Ordering.%-4s str()=%-14r value=%-6r == %r -> %s"
              % (name, str(member), member.value, name, member == name))
    print("  plain_value(Ordering.NM) =", plain_value(Ordering.NM))
    print("\n  NOTE: a *live* summary search returns `ordering` as a plain")
    print("  string, so this enum trap does not fire on the real path -- it was")
    print("  reproduced only with a hand-built stub.  plain_value() is a guard.")

    print("\n=== does SummaryDoc.dict() keep the enum? ===")
    from emmet.core.summary import SummaryDoc

    doc = SummaryDoc.model_construct(ordering=Ordering.NM)
    print("  SummaryDoc.model_construct(...).dict()['ordering'] = %r"
          % doc.dict().get("ordering"))
    print("  (model_construct skips validation; the real client path differs --")
    print("   trust a run's download.csv over this line.)")


if __name__ == "__main__":
    main()
