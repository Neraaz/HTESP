#!/usr/bin/env python
"""The stop report: *where* did the tutorial campaign stop, and what next.

The user's requirement for this driver was one sentence long -- "if stopped
somewhere, it will report where did it stop" -- so this module is the point of
the package, not an afterthought.  Every stop (a failed step, a step that
exited 0 and produced nothing, a timeout, a dependency that blocked a whole
tutorial, or Ctrl-C) produces the same three things:

* a one-screen summary table of all 42 tutorials, printed to the terminal;
* ``<workdir>/report.md`` -- the same table plus, for every stop, the tutorial,
  the step, the exact command line, the working directory, the exit code, the
  last 40 lines of that step's log, the artefacts that were expected and which
  were missing, how long it ran, and the exact command that retries from there;
* ``<workdir>/report.json`` -- the same content as data, for a wrapper script.

Nothing here runs anything; it reads a :class:`~tutorials.state.RunState`.
"""
from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path
from typing import Sequence

from tutorials.state import (BLOCKED, DONE, FAILED, PENDING, RUNNING, SKIPPED,
                             RunState, StepState, TutorialState)

#: how much of a failing step's log the report quotes
LOG_TAIL_LINES = 40

#: column order of the summary table
STATUSES = (DONE, FAILED, BLOCKED, SKIPPED, RUNNING, PENDING)

_MARK = {DONE: "ok", FAILED: "FAILED", BLOCKED: "blocked",
         SKIPPED: "skipped", RUNNING: "running", PENDING: "pending"}


def log_tail(path: str | None, lines: int = LOG_TAIL_LINES) -> list[str]:
    """Last *lines* lines of a step log, or a one-line explanation."""
    if not path:
        return ["(no log was written for this step)"]
    try:
        text = Path(path).read_text(errors="replace")
    except OSError as exc:
        return [f"(could not read {path}: {exc})"]
    tail = text.splitlines()[-lines:]
    return tail or ["(the log is empty)"]


def retry_command(code: str, step: StepState | None, mode: str,
                  workdir: Path | str) -> str:
    """The exact command line that retries the run from where it stopped."""
    parts = ["htesp-tutorials", "--resume", "--only", code]
    if step is not None and step.step_id:
        parts += ["--from", step.step_id]
    if mode == "dry-run":
        parts.append("--dry-run")
    parts += ["--workdir", str(workdir)]
    return " ".join(parts)


# --------------------------------------------------------------------------- #
#  the summary table
# --------------------------------------------------------------------------- #
def summary_rows(state: RunState, order: Sequence[str] | None = None
                 ) -> list[tuple[str, str, str, str, str]]:
    """``(code, status, steps, title, detail)`` per tutorial, in run order."""
    codes = list(order) if order else list(state.tutorials)
    rows = []
    for code in codes:
        tut = state.tutorials.get(code)
        if tut is None:
            rows.append((code, PENDING, "-", "", "never reached"))
            continue
        counts = tut.counts()
        done = counts.get(DONE, 0)
        total = len(tut.steps) or 0
        steps = f"{done}/{total}" if total else "-"
        detail = ""
        if tut.status == FAILED:
            step = tut.first_failure()
            detail = f"stopped at {step.step_id}" if step else tut.reason
        elif tut.status == BLOCKED:
            detail = f"blocked by {tut.blocked_by}"
        elif tut.reason:
            detail = tut.reason[:70]
        elif counts.get(SKIPPED):
            detail = f"{counts[SKIPPED]} step(s) skipped"
        rows.append((code, tut.status, steps, tut.title, detail))
    return rows


def summary_table(state: RunState, order: Sequence[str] | None = None) -> str:
    """The one-screen table that goes at the top of every report."""
    rows = summary_rows(state, order)
    width = max((len(r[3]) for r in rows), default=10)
    width = min(width, 58)
    head = f"{'tutorial':<10} {'status':<8} {'steps':<7} {'tutorial title':<{width}} notes"
    lines = [head, "-" * len(head)]
    for code, status, steps, title, detail in rows:
        lines.append(f"{code:<10} {_MARK.get(status, status):<8} {steps:<7} "
                     f"{title[:width]:<{width}} {detail}")
    counts: dict[str, int] = {}
    for _, status, *_rest in rows:
        counts[status] = counts.get(status, 0) + 1
    tally = ", ".join(f"{counts[s]} {s}" for s in STATUSES if counts.get(s))
    lines += ["-" * len(head), f"{len(rows)} tutorials: {tally or 'nothing run'}"]
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
#  the stop sections
# --------------------------------------------------------------------------- #
def stop_block(tut: TutorialState, step: StepState) -> list[str]:
    """The "where it stopped" block for one failed step, as markdown lines."""
    expected = ", ".join(step.expected) or "(none declared)"
    missing = ", ".join(step.missing) or "(none)"
    command = " ".join(step.command) or "(no command was started)"
    out = [
        f"### {tut.code} -- {tut.title}",
        "",
        f"* **stopped at step** {step.index} of {len(tut.steps)}: "
        f"`{step.step_id}`" + (f" (loop cycle {step.cycle})" if step.cycle > 1 else ""),
        f"* **reason**: {step.reason or 'unknown'}",
        f"* **command**: `{command}`",
        f"* **working directory**: `{step.workdir}`",
        f"* **exit code**: {step.exit_code if step.exit_code is not None else 'n/a'}",
        f"* **ran for**: {step.duration:.1f} s",
        f"* **log**: `{step.log or '(none)'}`",
        f"* **expected artefacts**: {expected}",
        f"* **missing artefacts**: {missing}",
    ]
    if step.jobs:
        out.append(f"* **cluster jobs**: {', '.join(step.jobs)}")
    if step.unverifiable:
        out.append("* **note**: this step could not be verified against the queue")
    out += ["", f"Last {LOG_TAIL_LINES} lines of the log:", "", "```"]
    out += log_tail(step.log)
    out += ["```", ""]
    return out


def blocked_block(tut: TutorialState) -> list[str]:
    return [f"### {tut.code} -- {tut.title}", "",
            f"* **blocked by**: {tut.blocked_by}",
            f"* **reason**: {tut.reason}", ""]


# --------------------------------------------------------------------------- #
#  assembling
# --------------------------------------------------------------------------- #
def build_markdown(state: RunState, order: Sequence[str] | None = None,
                   workdir: Path | str = ".") -> str:
    """The whole ``report.md``."""
    failed = [t for t in state.tutorials.values() if t.status == FAILED]
    blocked = [t for t in state.tutorials.values() if t.status == BLOCKED]
    lines = [
        "# HTESP tutorial run",
        "",
        f"* mode: `{state.mode}`",
        f"* started: {time.ctime(state.started)}",
        f"* updated: {time.ctime(state.updated)}",
        f"* work directory: `{workdir}`",
        f"* command line: `{' '.join(state.argv)}`",
    ]
    if state.interrupted:
        lines.append("* **the run was interrupted with Ctrl-C**")
    lines += ["", "## Summary", "", "```", summary_table(state, order), "```", ""]

    stopped = state.stopped_at()
    if stopped is None and not blocked:
        lines += ["## Where it stopped", "",
                  "It did not stop: every selected tutorial ran to the end.", ""]
    else:
        lines += ["## Where it stopped", ""]
        if stopped is not None:
            tut, step = stopped
            lines += [f"**The first thing that went wrong was {tut.code}, step "
                      f"{step.index} (`{step.step_id}`).**", "",
                      "Retry from exactly there with:", "",
                      "```",
                      retry_command(tut.code, step, state.mode, workdir),
                      "```", ""]
    if failed:
        lines += ["## Failures", ""]
        for tut in sorted(failed, key=lambda t: t.started or 0.0):
            step = tut.first_failure()
            if step is not None:
                lines += stop_block(tut, step)
            else:
                lines += [f"### {tut.code} -- {tut.title}", "",
                          f"* **reason**: {tut.reason}", ""]
            lines += ["Retry with:", "", "```",
                      retry_command(tut.code, step, state.mode, workdir), "```", ""]
    if blocked:
        lines += ["## Blocked", "",
                  "These never started because a tutorial they build on did not "
                  "finish.", ""]
        for tut in sorted(blocked, key=lambda t: t.code):
            lines += blocked_block(tut)
    unverifiable = [(t, s) for t in state.tutorials.values()
                    for s in t.steps.values() if s.unverifiable]
    if unverifiable:
        lines += ["## Unverified steps", "",
                  "These finished, but the driver could not confirm it:", ""]
        for tut, step in unverifiable:
            lines += [f"* `{tut.code}` step {step.index} (`{step.step_id}`): "
                      f"{step.reason}"]
        lines.append("")
    return "\n".join(lines) + "\n"


def build_json(state: RunState, order: Sequence[str] | None = None,
               workdir: Path | str = ".") -> dict:
    """The machine-readable twin of :func:`build_markdown`."""
    stopped = state.stopped_at()
    payload = {
        "mode": state.mode,
        "argv": state.argv,
        "started": state.started,
        "updated": state.updated,
        "interrupted": state.interrupted,
        "workdir": str(workdir),
        "summary": [
            {"code": c, "status": s, "steps": n, "title": t, "note": d}
            for c, s, n, t, d in summary_rows(state, order)],
        "counts": state.counts(),
        "stopped_at": None,
        "failures": [],
        "blocked": [],
        "would_submit": [
            {"tutorial": t.code, "step": s.index, "step_id": s.step_id,
             "lines": s.would_submit}
            for t in state.tutorials.values() for s in t.steps.values()
            if s.would_submit],
    }
    if stopped is not None:
        tut, step = stopped
        payload["stopped_at"] = {
            "tutorial": tut.code, "title": tut.title, "step": step.index,
            "step_id": step.step_id, "reason": step.reason,
            "command": step.command, "workdir": step.workdir,
            "exit_code": step.exit_code, "duration": step.duration,
            "log": step.log, "log_tail": log_tail(step.log),
            "expected": step.expected, "missing": step.missing,
            "retry": retry_command(tut.code, step, state.mode, workdir)}
    for tut in state.tutorials.values():
        if tut.status == FAILED:
            step = tut.first_failure()
            payload["failures"].append({
                "tutorial": tut.code, "title": tut.title, "reason": tut.reason,
                "step": asdict(step) if step else None,
                "log_tail": log_tail(step.log) if step else [],
                "retry": retry_command(tut.code, step, state.mode, workdir)})
        elif tut.status == BLOCKED:
            payload["blocked"].append({"tutorial": tut.code, "title": tut.title,
                                       "blocked_by": tut.blocked_by,
                                       "reason": tut.reason})
    return payload


def write_report(state: RunState, workdir: Path,
                 order: Sequence[str] | None = None) -> tuple[Path, Path]:
    """Write ``report.md`` and ``report.json``; returns both paths."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    md_path, json_path = workdir / "report.md", workdir / "report.json"
    md_path.write_text(build_markdown(state, order, workdir))
    json_path.write_text(json.dumps(build_json(state, order, workdir), indent=1))
    return md_path, json_path


def console_report(state: RunState, workdir: Path,
                   order: Sequence[str] | None = None) -> str:
    """What the driver prints when it finishes or stops."""
    lines = [summary_table(state, order), ""]
    stopped = state.stopped_at()
    blocked = [t for t in state.tutorials.values() if t.status == BLOCKED]
    if stopped is None and not blocked:
        lines.append("Every selected tutorial ran to the end.")
    if stopped is not None:
        tut, step = stopped
        lines += [
            "WHERE IT STOPPED",
            f"  tutorial : {tut.code} -- {tut.title}",
            f"  step     : {step.index}/{len(tut.steps)}  {step.step_id}"
            + (f"  (loop cycle {step.cycle})" if step.cycle > 1 else ""),
            f"  command  : {' '.join(step.command) or '(none)'}",
            f"  cwd      : {step.workdir}",
            f"  exit     : {step.exit_code if step.exit_code is not None else 'n/a'}"
            f"   after {step.duration:.1f}s",
            f"  reason   : {step.reason}",
            f"  missing  : {', '.join(step.missing) or '(nothing)'}",
            f"  log      : {step.log}",
            "",
            "  last lines of that log:",
        ]
        lines += [f"    {line}" for line in log_tail(step.log, 12)]
        lines += ["", "  retry from there with:",
                  f"    {retry_command(tut.code, step, state.mode, workdir)}"]
    if blocked:
        lines += ["", "BLOCKED (a tutorial they build on did not finish)"]
        lines += [f"  {t.code:<9} blocked by {t.blocked_by}"
                  for t in sorted(blocked, key=lambda t: t.code)]
    lines += ["", f"full report: {Path(workdir) / 'report.md'}",
              f"             {Path(workdir) / 'report.json'}"]
    return "\n".join(lines)
