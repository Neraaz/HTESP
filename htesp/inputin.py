#!/usr/bin/env python
"""Reading and writing ``input.in``.

``input.in`` is the six-line file every HTESP command starts from::

    1            <- start index (1-based, inclusive)
    30           <- end index   (EXCLUSIVE -- the loops run ii < end)
    200 0        <- nkpt [kcut]
    mpid-list.in <- tracking file
    phband dos   <- plot types, whitespace separated
    DFT = QE     <- QE or VASP

The original parser in ``mainprogram.py`` had three defects, all reproduced by
the tests in ``tests/test_inputin.py``:

* it guarded with ``len(lines) >= 4`` and then read ``lines[4]``, so a
  four-line file raised ``IndexError``;
* a shorter file left ``start``/``end``/``element``/``nkpt`` unbound, so the
  next statement raised ``NameError`` instead of reporting the real problem;
* on the first run it wrote ``plot_type`` as the string ``'phband'`` and then
  iterated over it, launching six plot jobs named ``p``, ``h``, ``b``, ``a``,
  ``n``, ``d``.

The sixth line is written unconditionally here.  The bash layer tested
``[ $dft == 'vasp' ]`` at 31 sites, which aborts with "unary operator expected"
when the line is missing.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_FILENAME = "input.in"

#: plot types understood by ``plot-scan``
PLOT_TYPES = ("phband", "band", "dos", "pdos", "a2f", "phdos", "bandproj", "phononproj")


class InputInError(ValueError):
    """Raised for an ``input.in`` that cannot be understood."""


@dataclass
class InputIn:
    """The parsed contents of ``input.in``."""

    start: int = 1
    end: int = 2
    nkpt: int = 200
    kcut: int = 0
    track: str = "mpid-list.in"
    plot_types: list[str] = field(default_factory=lambda: ["phband"])
    dft: str = "QE"
    path: Path | None = None

    # -- reading ---------------------------------------------------------- #
    @classmethod
    def parse(cls, text: str, path: os.PathLike | str | None = None) -> "InputIn":
        lines = [ln.rstrip("\n") for ln in text.splitlines()]
        if len(lines) < 4:
            raise InputInError(
                f"{path or DEFAULT_FILENAME} has {len(lines)} line(s); at least four are "
                "required (start, end, nkpt, tracking file).  "
                "Run 'mainprogram basicinfo' for the format."
            )

        def _int(value: str, what: str) -> int:
            try:
                return int(value.split()[0])
            except (ValueError, IndexError) as exc:
                raise InputInError(
                    f"{path or DEFAULT_FILENAME}: {what} must be an integer, got {value!r}"
                ) from exc

        start = _int(lines[0], "line 1 (start)")
        end = _int(lines[1], "line 2 (end)")
        kpt_tokens = lines[2].split()
        if not kpt_tokens:
            raise InputInError(f"{path or DEFAULT_FILENAME}: line 3 (nkpt) is empty")
        nkpt = _int(kpt_tokens[0], "line 3 (nkpt)")
        kcut = _int(kpt_tokens[1], "line 3 (kcut)") if len(kpt_tokens) > 1 else 0

        track_tokens = lines[3].split()
        if not track_tokens:
            raise InputInError(f"{path or DEFAULT_FILENAME}: line 4 (tracking file) is empty")
        track = track_tokens[0]

        plot_types = lines[4].split() if len(lines) > 4 and lines[4].strip() else ["phband"]

        dft = "QE"
        for line in lines[5:]:
            if "DFT" in line and "=" in line:
                dft = line.split("=", 1)[1].strip() or "QE"
                break

        if end <= start:
            raise InputInError(
                f"{path or DEFAULT_FILENAME}: end ({end}) must be greater than start "
                f"({start}); note that end is exclusive"
            )

        return cls(start=start, end=end, nkpt=nkpt, kcut=kcut, track=track,
                   plot_types=plot_types, dft=dft,
                   path=Path(path) if path else None)

    # -- tolerant loading (used by the workflow layer) -------------------- #
    @classmethod
    def load(cls, path: os.PathLike | str = DEFAULT_FILENAME,
             config: dict | None = None) -> "InputIn":
        """Best-effort parse that never raises.

        The workflow layer must keep going when a single field is malformed --
        the bash layer simply used an empty variable, which is how
        ``[ $dft == 'vasp' ]`` came to abort with "unary operator expected" at
        31 call sites.  Missing fields fall back to the documented defaults and
        ``dft`` falls back to ``download.inp.calc`` from ``config.json``.
        Unlike :meth:`parse`, ``dft`` is returned lower-cased.
        """
        def _int(value, fallback):
            try:
                return int(str(value).strip())
            except (TypeError, ValueError):
                return fallback

        try:
            lines = Path(path).read_text().splitlines()
        except OSError:
            lines = []

        obj = cls(start=1, end=2, nkpt=50, kcut=0, track="mpid.in",
                  plot_types=["phband"], dft="", path=Path(path))
        if len(lines) > 0:
            obj.start = _int(lines[0], obj.start)
        if len(lines) > 1:
            obj.end = _int(lines[1], obj.end)
        if len(lines) > 2:
            parts = lines[2].split()
            if parts:
                obj.nkpt = _int(parts[0], obj.nkpt)
            if len(parts) > 1:
                obj.kcut = _int(parts[1], 0)
        if len(lines) > 3 and lines[3].strip():
            obj.track = lines[3].split()[0]
        if len(lines) > 4 and lines[4].strip():
            obj.plot_types = lines[4].split()
        for line in lines:
            if line.strip().upper().startswith("DFT"):
                parts = line.replace("=", " = ").split()
                if len(parts) >= 3:
                    obj.dft = parts[2].strip().lower()
                break
        if not obj.dft:
            if config:
                obj.dft = str(config.get("download", {})
                              .get("inp", {}).get("calc", "qe")).lower()
            else:
                obj.dft = "qe"
        return obj

    @classmethod
    def read(cls, path: os.PathLike | str = DEFAULT_FILENAME) -> "InputIn":
        path = Path(path)
        return cls.parse(path.read_text(), path)

    @classmethod
    def load_or_create(cls, path: os.PathLike | str = DEFAULT_FILENAME,
                       dft: str | None = None) -> "InputIn":
        """Read ``input.in``, creating a sensible default when it is absent."""
        path = Path(path)
        if path.is_file():
            return cls.read(path)
        if dft is None:
            from htesp.config import config
            dft = str(config()["download"]["inp"]["calc"]).upper()
        created = cls(dft=dft, path=path)
        created.write(path)
        return created

    # -- writing ---------------------------------------------------------- #
    def render(self) -> str:
        return (
            f"{self.start}\n"
            f"{self.end}\n"
            f"{self.nkpt} {self.kcut}\n"
            f"{self.track}\n"
            f"{' '.join(self.plot_types)}\n"
            f"DFT = {self.dft}\n"
        )

    def write(self, path: os.PathLike | str | None = None) -> Path:
        target = Path(path or self.path or DEFAULT_FILENAME)
        target.write_text(self.render())
        self.path = target
        return target

    # -- convenience ------------------------------------------------------ #
    @property
    def is_vasp(self) -> bool:
        return self.dft.strip().lower() == "vasp"

    @property
    def is_qe(self) -> bool:
        return not self.is_vasp

    @property
    def count(self) -> int:
        """Number of materials the range covers (``end`` is exclusive)."""
        return max(0, self.end - self.start)

    def summary(self) -> str:
        return (
            f"start: {self.start}, end: {self.end}, tracking_file: {self.track}, "
            f"nkpoint: {self.nkpt}, DFT: {self.dft}"
        )

    def check_track_file(self, root: os.PathLike | str = ".") -> list[str]:
        """Warn when the range does not line up with the tracking file."""
        track = Path(root) / self.track
        if not track.is_file():
            return [f"tracking file {self.track!r} not found in {Path(root).resolve()}"]
        entries = [ln for ln in track.read_text().splitlines() if ln.strip()]
        problems = []
        if self.start < 1:
            problems.append(f"start ({self.start}) must be >= 1")
        if self.end - 1 > len(entries):
            problems.append(
                f"end ({self.end}) is past the {len(entries)} entries of {self.track}; "
                f"the last usable value is {len(entries) + 1} (end is exclusive)"
            )
        return problems
