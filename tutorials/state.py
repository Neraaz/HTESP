#!/usr/bin/env python
"""The tutorial runner's checkpoint file.

A campaign over ``examples/`` is hours of cluster time; it will be interrupted.
Everything the runner learns about a step is written to ``state.json`` in the
work directory *immediately after that step*, so ``--resume`` can pick the run
up exactly where it stopped and so the stop report can be rebuilt from the file
alone, without re-reading any logs.

The file is plain JSON with an explicit schema version, written atomically
(temporary file + :meth:`Path.replace`) so a checkpoint interrupted by SIGKILL
leaves the previous, complete checkpoint in place rather than a truncated one.
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

#: bumped whenever the on-disk layout changes incompatibly
SCHEMA_VERSION = 1

PENDING, RUNNING, DONE, FAILED, SKIPPED, BLOCKED = (
    "pending", "running", "done", "failed", "skipped", "blocked")

#: statuses that ``--resume`` does not re-run
TERMINAL = (DONE, SKIPPED)


@dataclass
class StepState:
    """What happened to one step of one tutorial."""

    step_id: str
    key: str
    index: int = 0
    status: str = PENDING
    command: list[str] = field(default_factory=list)
    workdir: str = ""
    exit_code: int | None = None
    duration: float = 0.0
    started: float | None = None
    finished: float | None = None
    log: str | None = None
    expected: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    reason: str = ""
    cycle: int = 1
    unverifiable: bool = False
    jobs: list[str] = field(default_factory=list)
    would_submit: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.status in TERMINAL


@dataclass
class TutorialState:
    """What happened to one tutorial."""

    code: str
    title: str = ""
    status: str = PENDING
    workdir: str = ""
    reason: str = ""
    blocked_by: str = ""
    started: float | None = None
    finished: float | None = None
    steps: dict[str, StepState] = field(default_factory=dict)

    def step(self, key: str) -> StepState | None:
        return self.steps.get(key)

    def counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for step in self.steps.values():
            out[step.status] = out.get(step.status, 0) + 1
        return out

    def first_failure(self) -> StepState | None:
        """The first step, in execution order, that failed."""
        for step in sorted(self.steps.values(), key=lambda s: (s.index, s.cycle)):
            if step.status == FAILED:
                return step
        return None


@dataclass
class RunState:
    """The whole checkpoint: every tutorial, every step, plus run metadata."""

    path: Path
    version: int = SCHEMA_VERSION
    mode: str = "dry-run"
    argv: list[str] = field(default_factory=list)
    started: float = field(default_factory=time.time)
    updated: float = field(default_factory=time.time)
    interrupted: bool = False
    tutorials: dict[str, TutorialState] = field(default_factory=dict)

    # -- construction ------------------------------------------------------- #
    @classmethod
    def load(cls, path: os.PathLike | str) -> "RunState":
        """Read a checkpoint; an absent or unreadable file yields a fresh one."""
        path = Path(path)
        state = cls(path=path)
        try:
            raw = json.loads(path.read_text())
        except (OSError, ValueError):
            return state
        if int(raw.get("version", 0)) != SCHEMA_VERSION:
            # An older checkpoint is kept on disk but not trusted: re-running is
            # cheap next to guessing what a changed field used to mean.
            return state
        state.version = int(raw["version"])
        state.mode = raw.get("mode", state.mode)
        state.argv = list(raw.get("argv", []))
        state.started = float(raw.get("started", state.started))
        state.updated = float(raw.get("updated", state.updated))
        state.interrupted = bool(raw.get("interrupted", False))
        for code, blob in (raw.get("tutorials") or {}).items():
            steps = {key: StepState(**step)
                     for key, step in (blob.get("steps") or {}).items()}
            payload = {k: v for k, v in blob.items() if k != "steps"}
            state.tutorials[code] = TutorialState(steps=steps, **payload)
        return state

    # -- persistence -------------------------------------------------------- #
    def to_dict(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "mode": self.mode,
            "argv": self.argv,
            "started": self.started,
            "updated": self.updated,
            "interrupted": self.interrupted,
            "tutorials": {
                code: {**{k: v for k, v in asdict(tut).items() if k != "steps"},
                       "steps": {key: asdict(step) for key, step in tut.steps.items()}}
                for code, tut in self.tutorials.items()},
        }

    def save(self) -> None:
        """Write the checkpoint atomically."""
        self.updated = time.time()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(self.path.suffix + ".tmp")
        tmp.write_text(json.dumps(self.to_dict(), indent=1, sort_keys=True))
        tmp.replace(self.path)

    # -- queries ------------------------------------------------------------ #
    def tutorial(self, code: str, title: str = "") -> TutorialState:
        """Get (creating if needed) the record for *code*."""
        tut = self.tutorials.get(code)
        if tut is None:
            tut = TutorialState(code=code, title=title)
            self.tutorials[code] = tut
        elif title and not tut.title:
            tut.title = title
        return tut

    def is_done(self, code: str, key: str) -> bool:
        """True when this step need not be re-run under ``--resume``."""
        tut = self.tutorials.get(code)
        step = tut.step(key) if tut else None
        return bool(step and step.ok)

    def reset(self, codes: Iterable[str] | None = None) -> None:
        """Forget the recorded progress (``--restart``)."""
        if codes is None:
            self.tutorials.clear()
            return
        for code in codes:
            self.tutorials.pop(code, None)

    def counts(self) -> dict[str, int]:
        """Tutorial-level status counts, for the summary line."""
        out: dict[str, int] = {}
        for tut in self.tutorials.values():
            out[tut.status] = out.get(tut.status, 0) + 1
        return out

    def stopped_at(self) -> tuple[TutorialState, StepState] | None:
        """The first failed step of the first failed tutorial, in run order.

        This is what "where did it stop" means: the earliest failure, since
        everything after it is either blocked by it or was never attempted.
        """
        ordered = sorted(self.tutorials.values(),
                         key=lambda t: (t.started is None, t.started or 0.0))
        for tut in ordered:
            step = tut.first_failure()
            if step is not None:
                return tut, step
        return None


def step_key(step_id: str, cycle: int = 1) -> str:
    """Checkpoint key for a step; loop cycles after the first get a suffix."""
    return step_id if cycle <= 1 else f"{step_id}#{cycle}"
