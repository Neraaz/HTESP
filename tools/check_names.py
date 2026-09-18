#!/usr/bin/env python
"""A small undefined-name checker, for machines where ruff cannot be installed.

It is not a replacement for ruff or pyflakes; it catches the one class of
defect that bit this package repeatedly -- a name used at module or function
scope that is never bound anywhere (``procar_jband``, ``qshift``, ``read``,
``input_data``, ``submission_files``) -- with nothing but the standard library.

    python tools/check_names.py htesp tutorials

Exit status is 1 when something is reported.
"""
from __future__ import annotations

import ast
import builtins
import sys
from pathlib import Path

# One definition of "a real Python source file".  Globbing ``*.py`` here meant
# this tool crashed with UnicodeDecodeError on macOS AppleDouble sidecars
# (``._module.py``, a binary blob beside every real module after a copy from a
# Mac) -- while the test suite, which already filtered them, ran clean.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tests.helpers import python_sources  # noqa: E402

BUILTIN_NAMES = set(dir(builtins)) | {"__file__", "__name__", "__doc__",
                                      "__spec__", "__package__", "__builtins__"}


class ScopeVisitor(ast.NodeVisitor):
    """Collects every name bound anywhere in a module, then every name used."""

    def __init__(self) -> None:
        self.bound: set[str] = set()
        self.used: list[tuple[str, int]] = []

    # -- bindings ----------------------------------------------------------
    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            self.bound.add((alias.asname or alias.name).split(".")[0])
        self.generic_visit(node)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        for alias in node.names:
            self.bound.add(alias.asname or alias.name)
        self.generic_visit(node)

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self.bound.add(node.name)
        for arg in [*node.args.args, *node.args.kwonlyargs, *node.args.posonlyargs]:
            self.bound.add(arg.arg)
        if node.args.vararg:
            self.bound.add(node.args.vararg.arg)
        if node.args.kwarg:
            self.bound.add(node.args.kwarg.arg)
        self.generic_visit(node)

    visit_AsyncFunctionDef = visit_FunctionDef        # type: ignore[assignment]

    def visit_Lambda(self, node: ast.Lambda) -> None:
        for arg in [*node.args.args, *node.args.kwonlyargs]:
            self.bound.add(arg.arg)
        self.generic_visit(node)

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        self.bound.add(node.name)
        self.generic_visit(node)

    def visit_ExceptHandler(self, node: ast.ExceptHandler) -> None:
        if node.name:
            self.bound.add(node.name)
        self.generic_visit(node)

    def visit_Global(self, node: ast.Global) -> None:
        self.bound.update(node.names)

    visit_Nonlocal = visit_Global                      # type: ignore[assignment]

    def visit_Name(self, node: ast.Name) -> None:
        if isinstance(node.ctx, (ast.Store, ast.Del)):
            self.bound.add(node.id)
        else:
            self.used.append((node.id, node.lineno))

    def visit_arg(self, node: ast.arg) -> None:
        self.bound.add(node.arg)
        self.generic_visit(node)

    def visit_comprehension(self, node: ast.comprehension) -> None:
        self.generic_visit(node)


def check(path: Path) -> list[str]:
    """Report names used in ``path`` that nothing in it ever binds."""
    try:
        tree = ast.parse(path.read_text())
    except SyntaxError as exc:
        return [f"{path}:{exc.lineno}: syntax error: {exc.msg}"]

    visitor = ScopeVisitor()
    visitor.visit(tree)
    problems = []
    seen: set[str] = set()
    for name, line in visitor.used:
        if name in visitor.bound or name in BUILTIN_NAMES or name in seen:
            continue
        seen.add(name)
        problems.append(f"{path}:{line}: undefined name {name!r}")
    return problems


def main(argv: list[str] | None = None) -> int:
    argv = list(argv if argv is not None else sys.argv[1:]) or ["htesp"]
    targets: list[Path] = []
    for entry in argv:
        path = Path(entry)
        targets.extend(python_sources(path, recursive=True)
                       if path.is_dir() else [path])

    problems: list[str] = []
    for target in targets:
        problems.extend(check(target))

    for problem in problems:
        print(problem)
    print(f"\n{len(targets)} file(s) checked, {len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
