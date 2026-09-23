#!/usr/bin/env python
"""``htesp-tutorials`` -- run every calculation in ``examples/`` end to end.

This is the driver a batch submission script launches (see
``tutorials/submit_tutorials.sh``).  It walks the catalogue in
:mod:`tutorials.catalog`, runs each tutorial's ``mainprogram`` steps in its own
work directory, waits for the cluster jobs each step submits, checkpoints after
every step -- and, when it stops, says exactly where.

    htesp-tutorials                           # the whole tree
    htesp-tutorials --only QE/9,QE/12         # the real thing
    htesp-tutorials --resume                  # carry on where it stopped
    htesp-tutorials --list                    # print the catalogue and exit

Exit codes: ``0`` everything finished, ``1`` something failed or was blocked,
``2`` a preflight check failed before anything ran, ``130`` interrupted.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from tutorials import report as report_mod
from tutorials import catalog as catalog_mod
from tutorials.catalog import (CATALOG, format_catalog, parse_codes,
                               select)
from tutorials.runner import (RunOptions, TutorialRunner,
                              preflight)
from tutorials.state import BLOCKED, FAILED, RunState

LOG = logging.getLogger("htesp.tutorials")


def build_parser() -> argparse.ArgumentParser:
    """The command line.  Kept in one place so the self-tests can exercise it."""
    parser = argparse.ArgumentParser(
        prog="htesp-tutorials",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Run every worked example under examples/ and report where "
                    "the run stopped.",
        epilog="Tutorial codes look like QE/9 or VASP/14; run --list to see them.")
    parser.add_argument("--workdir", default="tutorial_runs_root", type=Path,
                        help="root for the work directories, logs, checkpoint and "
                             "report (default: ./tutorial_runs_root)")

    resume = parser.add_mutually_exclusive_group()
    resume.add_argument("--resume", action="store_true", default=True,
                        help="re-run only what is not already done (the default)")
    resume.add_argument("--restart", action="store_true",
                        help="forget the checkpoint and run everything again")

    parser.add_argument("--keep_output", choices=("yes", "no"), default="yes",
                        help="keep the generated tutorial_runs/ directories "
                             "when the run ends (default: yes).  'no' removes "
                             "them, including those of failed tutorials, so "
                             "keep 'yes' while debugging.  Logs, the report "
                             "and the checkpoint are kept either way")
    parser.add_argument("--only", default=None,
                        help="comma-separated tutorial codes to run, e.g. "
                             "QE/9,VASP/14.  A bare tree name means all of it: "
                             "--only QE")
    parser.add_argument("--skip", default=None,
                        help="comma-separated tutorial codes to leave out; a "
                             "bare tree name works here too")
    parser.add_argument("--from", dest="from_step", default=None,
                        help="start each selected tutorial at this step id "
                             "(the retry line in the report uses this)")
    parser.add_argument("--workers", type=int, default=1,
                        help="passed through to mainprogram --workers "
                             "(default: 1).  A tutorial works on one or two "
                             "materials, so a bigger pool buys nothing and "
                             "multiplies with --jobs: four tutorials at "
                             "mainprogram's own default of 8 is 36 processes, "
                             "and a login node allows 100")
    parser.add_argument("--jobs", type=int, default=1, metavar="N",
                        help="run N tutorials at once (default: 1).  Most of a "
                             "sweep is spent waiting on the Materials Project, "
                             "OQMD and AFLOW servers, and those tutorials are "
                             "independent, so the waiting overlaps.  Keep it "
                             "modest: the same APIs rate-limit")
    parser.add_argument("--force", action="store_true",
                        help="run even when preflight reports errors")
    parser.add_argument("--output", action="store_true",
                        help="after the run, list every file it produced and "
                             "what that file is for, grouped by the step that "
                             "wrote it.  A step that wrote nothing is named as "
                             "such -- which is how two wrong artefact "
                             "declarations were found")
    parser.add_argument("--list", action="store_true",
                        help="print the catalogue and exit")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    return parser


def _configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")


def main(argv: list[str] | None = None) -> int:
    """Entry point for ``htesp-tutorials``; returns the process exit code."""
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    _configure_logging(args.verbose)

    # `examples/` is found, never given: $HTESP_EXAMPLES, then ./examples,
    # then beside the installed package.  The flag is gone because the answer
    # is the tree this package ships with -- see catalog.find_examples.
    catalog = CATALOG
    examples = catalog_mod.EXAMPLES

    if args.list:
        print(format_catalog(catalog))
        return 0

    try:
        only, skip = parse_codes(args.only), parse_codes(args.skip)
    except ValueError as exc:
        LOG.error("%s", exc)
        return 2

    codes = select(only=only, skip=skip, catalog=catalog)
    if not codes:
        LOG.error("no tutorials selected (--only %s --skip %s)",
                  args.only, args.skip)
        return 2

    # not created until preflight has passed: `--workdir examples/` should not
    # leave a directory behind in the reference tree before being rejected
    workdir = Path(args.workdir).resolve()
    options = RunOptions(
        workdir=workdir, resume=not args.restart,
        workers=args.workers,
        jobs=args.jobs,
        from_step=args.from_step, verbose=args.verbose, examples=examples,
        keep="all" if args.keep_output == "yes" else "none",
    )

    problems = preflight(codes, options, catalog=catalog)
    for problem in problems:
        (LOG.error if problem.fatal else LOG.warning)("preflight: %s", problem)
    if any(p.fatal for p in problems) and not args.force:
        LOG.error("%d preflight error(s); nothing was run.  Fix them, or pass "
                  "--force to run anyway.", sum(p.fatal for p in problems))
        return 2

    workdir.mkdir(parents=True, exist_ok=True)
    state = RunState.load(workdir / "state.json") if options.resume else RunState(
        path=workdir / "state.json")
    state.argv = ["htesp-tutorials", *argv]
    if not options.resume:
        state.reset(codes)

    LOG.info("running %d tutorial(s) under %s", len(codes), workdir)
    runner = TutorialRunner(codes, options, state=state, catalog=catalog)
    try:
        state = runner.run()
    except Exception as exc:                        # noqa: BLE001 - reported below
        LOG.exception("the driver itself failed: %s", exc)
        state.save()

    md_path, _json_path = report_mod.write_report(state, workdir, codes)
    print()
    print(report_mod.console_report(state, workdir, codes))
    if args.output:
        from tutorials import manifest

        print(manifest.render(state))
    LOG.info("report written to %s", md_path)

    if state.interrupted:
        return 130
    bad = [t for t in state.tutorials.values() if t.status in (FAILED, BLOCKED)]
    return 1 if bad else 0


if __name__ == "__main__":                          # pragma: no cover
    raise SystemExit(main())
