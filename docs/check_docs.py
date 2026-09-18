#!/usr/bin/env python
"""Documentation lint that does not need Sphinx installed.

``sphinx-build`` cannot be installed on every machine that touches this
repository, and the documentation had drifted a long way from the code while
nobody could build it: commands that do not exist, a flagship JSON block that
did not parse, a label defined twice, a ``:ref:`` to nothing.
``tests/test_docs.py`` pins the worst of that; this script is the wider sweep.

    python docs/check_docs.py          # report problems, exit 1 if any
    python docs/check_docs.py -v       # also say what was checked

Checks
------

1.  every ``.. code-block:: json`` body parses (when it is a whole document)
    and every ``.. code-block:: python`` body compiles;
2.  no label is defined twice, and every ``:ref:`` and ``:doc:`` target exists;
3.  every page named in a ``toctree`` exists, and every page is reachable from
    ``index.rst``;
4.  every ``mainprogram <name>`` in the prose is a command or process number
    that ``htesp/help_text.py`` knows about;
5.  every repository path written in a literal exists in the tree;
6.  no real-looking Materials Project API key is present.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
sys.path.insert(0, str(ROOT))

#: directories a documented path may start with
TREE_PREFIXES = (
    "htesp/", "bin/", "tests/", "tools/", "tutorials/", "examples/",
    "utility/", "legacy/", "docs/", "INSTALL/", "_reports/", "_removed/",
)
#: bare file names that must exist at the repository root
ROOT_FILES = {"pyproject.toml", "LICENSE", "README.md", "requirements.txt",
              "CHANGELOG.md", "MANIFEST.in", "requirement_docs.txt"}

#: an MP key is 32 alphanumerics; the placeholder is the only one allowed
API_KEY_RE = re.compile(r"\b[A-Za-z0-9]{32}\b")
API_KEY_PLACEHOLDER = "use_your_API_KEY"


def rst_files() -> list[Path]:
    return sorted(DOCS.glob("*.rst"))


def code_blocks(text: str):
    """Yield ``(language, line_number, body)`` for every code-block."""
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        match = re.match(r"\s*\.\.\s+code-block::\s*(\S+)\s*$", lines[index])
        if not match:
            index += 1
            continue
        language, start = match.group(1), index + 1
        index += 1
        while index < len(lines) and not lines[index].strip():
            index += 1
        body, indent = [], None
        while index < len(lines):
            line = lines[index]
            if not line.strip():
                body.append("")
                index += 1
                continue
            current = len(line) - len(line.lstrip())
            if indent is None:
                indent = current
            if current < indent:
                break
            body.append(line[indent:])
            index += 1
        yield language, start, "\n".join(body).strip("\n")


# --------------------------------------------------------------------------- #
#  the checks
# --------------------------------------------------------------------------- #
def check_code_blocks(problems: list[str], counts: dict) -> None:
    for path in rst_files():
        for language, line, body in code_blocks(path.read_text()):
            snippet = body.strip()
            if not snippet:
                continue
            if language == "json":
                counts["json"] += 1
                if not snippet.startswith("{"):
                    problems.append(
                        f"{path.name}:{line}: a json code-block that is not a "
                        f"whole document; declare it as 'text'")
                    continue
                try:
                    json.loads(snippet)
                except ValueError as exc:
                    problems.append(f"{path.name}:{line}: invalid JSON: {exc}")
            elif language == "python":
                counts["python"] += 1
                try:
                    compile(body, f"{path.name}:{line}", "exec")
                except SyntaxError as exc:
                    problems.append(
                        f"{path.name}:{line}: python block does not compile: {exc.msg}")


def collect_labels() -> dict[str, str]:
    labels: dict[str, str] = {}
    for path in rst_files():
        for number, line in enumerate(path.read_text().splitlines(), 1):
            match = re.match(r"\.\.\s+_([A-Za-z0-9_-]+):\s*$", line.strip())
            if match:
                labels.setdefault(match.group(1), f"{path.name}:{number}")
    return labels


def check_references(problems: list[str], counts: dict) -> None:
    seen: dict[str, str] = {}
    for path in rst_files():
        for number, line in enumerate(path.read_text().splitlines(), 1):
            match = re.match(r"\.\.\s+_([A-Za-z0-9_-]+):\s*$", line.strip())
            if not match:
                continue
            label, where = match.group(1), f"{path.name}:{number}"
            if label in seen:
                problems.append(
                    f"duplicate label {label!r} ({seen[label]} and {where})")
            seen[label] = where

    labels = set(seen) | {"search", "genindex", "modindex"}
    pages = {path.stem for path in rst_files()}
    for path in rst_files():
        text = path.read_text()
        counts["ref"] += len(re.findall(r":ref:`", text))
        for target in re.findall(r":ref:`[^<`]*<([^`>]+)>`", text):
            if target not in labels:
                problems.append(f"{path.name}: :ref: to unknown label {target!r}")
        for target in re.findall(r":ref:`([A-Za-z0-9_-]+)`", text):
            if target not in labels:
                problems.append(f"{path.name}: :ref: to unknown label {target!r}")
        counts["doc"] += len(re.findall(r":doc:`", text))
        for target in re.findall(r":doc:`[^<`]*<([^`>]+)>`", text):
            if target.lstrip("/") not in pages:
                problems.append(f"{path.name}: :doc: to unknown page {target!r}")
        for target in re.findall(r":doc:`([A-Za-z0-9_/-]+)`", text):
            if target.lstrip("/") not in pages:
                problems.append(f"{path.name}: :doc: to unknown page {target!r}")


def toctree_entries(text: str) -> list[str]:
    entries, lines, index = [], text.splitlines(), 0
    while index < len(lines):
        if lines[index].strip().startswith(".. toctree::"):
            index += 1
            while index < len(lines):
                stripped = lines[index].strip()
                if not stripped:
                    index += 1
                    continue
                if not lines[index].startswith((" ", "\t")):
                    break
                if not stripped.startswith(":"):
                    entries.append(stripped)
                index += 1
        else:
            index += 1
    return entries


def check_toctree(problems: list[str], counts: dict) -> None:
    pages = {path.stem for path in rst_files()}
    reachable = set()
    for path in rst_files():
        for entry in toctree_entries(path.read_text()):
            counts["toctree"] += 1
            name = entry.lstrip("/")
            if name not in pages:
                problems.append(f"{path.name}: toctree names missing page {entry!r}")
            reachable.add(name)
    for page in sorted(pages - reachable - {"index"}):
        problems.append(f"{page}.rst is not in any toctree")


def known_commands() -> tuple[set[str], set[str]]:
    from htesp.help_text import PROCESS_SUMMARY, SUMMARY
    return set(SUMMARY), set(PROCESS_SUMMARY)


def check_command_names(problems: list[str], counts: dict) -> None:
    try:
        commands, processes = known_commands()
    except Exception as exc:                      # pragma: no cover
        problems.append(f"could not import htesp.help_text: {exc}")
        return
    allowed = commands | processes | {"--list", "--version", "--help", "-h"}
    pattern = re.compile(r"mainprogram\s+([A-Za-z0-9_][A-Za-z0-9_.-]*)")
    for path in rst_files():
        if path.name == "command.rst":
            continue                              # generated from the same tables
        for number, line in enumerate(path.read_text().splitlines(), 1):
            for name in pattern.findall(line):
                counts["command"] += 1
                if name in allowed:
                    continue
                if re.fullmatch(r"\d+", name):
                    problems.append(
                        f"{path.name}:{number}: process {name!r} does not exist")
                else:
                    problems.append(
                        f"{path.name}:{number}: unknown command "
                        f"'mainprogram {name}'")


def check_paths(problems: list[str], counts: dict) -> None:
    literal = re.compile(r"``([^`\n]+)``")
    for path in rst_files():
        for number, line in enumerate(path.read_text().splitlines(), 1):
            for chunk in literal.findall(line):
                candidate = chunk.strip().rstrip(".,;:")
                if any(ch in candidate for ch in "{}<>*$ |"):
                    continue
                if not (candidate.startswith(TREE_PREFIXES)
                        or candidate in ROOT_FILES):
                    continue
                counts["path"] += 1
                if not (ROOT / candidate).exists():
                    problems.append(
                        f"{path.name}:{number}: path {candidate!r} is not in the tree")


def check_api_key(problems: list[str], counts: dict) -> None:
    for path in rst_files():
        text = path.read_text()
        for number, line in enumerate(text.splitlines(), 1):
            if API_KEY_PLACEHOLDER in line:
                continue
            for hit in API_KEY_RE.findall(line):
                if hit.isdigit() or hit.isalpha():
                    continue
                problems.append(
                    f"{path.name}:{number}: what looks like an API key: {hit}")
            counts["line"] += 1


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv[1:])

    problems: list[str] = []
    counts = {k: 0 for k in
              ("json", "python", "ref", "doc", "toctree", "command", "path", "line")}

    check_code_blocks(problems, counts)
    check_references(problems, counts)
    check_toctree(problems, counts)
    check_command_names(problems, counts)
    check_paths(problems, counts)
    check_api_key(problems, counts)

    if args.verbose:
        print(f"{len(rst_files())} pages: "
              f"{counts['json']} json and {counts['python']} python blocks, "
              f"{counts['ref']} :ref:, {counts['doc']} :doc:, "
              f"{counts['toctree']} toctree entries, "
              f"{counts['command']} command mentions, "
              f"{counts['path']} repository paths")

    if problems:
        print(f"{len(problems)} problem(s):")
        for problem in problems:
            print(f"  {problem}")
        return 1
    print("docs OK")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
