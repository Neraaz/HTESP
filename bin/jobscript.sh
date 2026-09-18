#!/bin/bash
# HTESP compatibility shim -- sourced, not executed.
#
# The original src/bash/jobscript.sh defined qerun/vasprun for the 21 scan
# scripts that sourced it.  Those scans are now htesp/workflow.py methods and
# submit through htesp.workflow.Scheduler, so this file exists only for
# personal scripts that still source it.  Three things differ from the original
# (kept verbatim at legacy/bash/jobscript.sh):
#
#   * sbatch is called with --parsable and the job id is recorded in
#     .htesp_job.json next to the script, in the same format the Python
#     Scheduler writes, so `mainprogram checkph` can ask squeue about *this*
#     job instead of grepping the queue for the compound name (which matched
#     every job containing "B", "C", "Si", ...);
#   * a failed sbatch is reported instead of being discarded;
#   * SLEEP defaults to 1 second rather than 0.01, which was no throttle at all.
#     Override it by exporting SLEEP before sourcing this file.

: "${SLEEP:=1}"

_htesp_record_job() {
    # $1 = job id, $2 = tag
    [ -n "$1" ] || return 0
    "${HTESP_PYTHON:-python3}" - "$1" "$2" <<'PY' 2>/dev/null || true
import json, sys, time, pathlib
job, tag = sys.argv[1], sys.argv[2]
store = pathlib.Path(".htesp_job.json")
try:
    data = json.loads(store.read_text())
except Exception:
    data = {}
data.setdefault(tag, []).append({"job": job, "time": time.time()})
store.write_text(json.dumps(data, indent=1))
PY
}

_htesp_submit() {
    # $1 = script file, $2 = tag
    local script="$1" tag="$2" jid
    if ! jid=$(sbatch --parsable "$script"); then
        echo "htesp: sbatch failed for $script in $PWD" >&2
        return 1
    fi
    echo "$jid"
    _htesp_record_job "$jid" "$tag"
    sleep "$SLEEP"
}

_htesp_job_name() {
    # $1 = mpid, $2 = compound, $3 = suffix (may be empty)
    local a="$1" b="$2" suffix="$3" base
    if   [ -f ../../CALC_VISIBLE_WITH_ID ];      then base="$a"
    elif [ -f ../../CALC_VISIBLE_WITH_NAME ];    then base="$b"
    elif [ -f ../../CALC_VISIBLE_WITH_ID-NAME ]; then base="$a-$b"
    else base=""
    fi
    if [ -z "$base" ]; then
        printf '%s' ""
    elif [ -n "$suffix" ]; then
        printf '%s-%s' "$base" "$suffix"
    else
        printf '%s' "$base"
    fi
}

vasprun() {
    local a="$1" b="$2" name
    name=$(_htesp_job_name "$a" "$b" "")
    if [ -n "$name" ]; then
        mv -f run.sh "$name.sh" && _htesp_submit "$name.sh" vasp
    else
        _htesp_submit run.sh vasp
    fi
}

qerun() {
    local a="$1" b="$2" stage="$3" name
    name=$(_htesp_job_name "$a" "$b" "$stage")
    if [ -n "$name" ]; then
        mv -f "run-$stage.sh" "$name.sh" && _htesp_submit "$name.sh" "$stage"
    else
        _htesp_submit "run-$stage.sh" "$stage"
    fi
}
