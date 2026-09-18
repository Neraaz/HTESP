#!/usr/bin/env python
"""Regenerate the full ``config.json`` listing in ``docs/param.rst``.

The listing used to be maintained by hand.  It drifted: a stray ``]`` after
the ``aflow.prop`` array and two missing closing braces meant the block anyone
would copy and paste was not valid JSON at all.  It is now produced from the
file the package actually ships, ``htesp/data/config.json``, which is also the
default every user configuration is merged over, so the two cannot disagree.

    python docs/gen_param_block.py            # rewrite the block
    python docs/gen_param_block.py --check    # exit 1 if it is out of date

The block lives between the ``.. config-json-start`` and ``.. config-json-end``
comments in ``docs/param.rst``.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
SOURCE = ROOT / "htesp" / "data" / "config.json"
TARGET = DOCS / "param.rst"

START = ".. config-json-start"
END = ".. config-json-end"
INDENT = "    "


def render() -> str:
    """Return the replacement text, including both marker lines."""
    with open(SOURCE) as handle:
        data = json.load(handle)
    body = json.dumps(data, indent=2)
    lines = [
        START,
        "   Generated from htesp/data/config.json by docs/gen_param_block.py.",
        "   Run that script after changing the shipped default; do not edit by hand.",
        "",
        ".. code-block:: json",
        "",
    ]
    lines += [(INDENT + line).rstrip() for line in body.splitlines()]
    lines += ["", END]
    return "\n".join(lines)


def splice(text: str, block: str) -> str:
    lines = text.splitlines()
    try:
        first = next(i for i, line in enumerate(lines) if line.strip() == START)
        last = next(i for i, line in enumerate(lines) if line.strip() == END)
    except StopIteration:
        raise SystemExit(
            f"{TARGET}: markers {START!r} and {END!r} not found")
    return "\n".join(lines[:first] + block.splitlines() + lines[last + 1:]) + "\n"


def main(argv: list[str]) -> int:
    check = "--check" in argv[1:]
    current = TARGET.read_text()
    updated = splice(current, render())
    if current == updated:
        if check:
            print(f"{TARGET.name} is up to date")
        return 0
    if check:
        print(f"{TARGET.name} is out of date; run python docs/gen_param_block.py")
        return 1
    TARGET.write_text(updated)
    print(f"rewrote the config.json block in {TARGET.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
