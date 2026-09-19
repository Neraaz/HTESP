#!/bin/bash
# -----------------------------------------------------------------------------
# submit_tutorials.sh -- run every HTESP tutorial from a batch job.
#
#   sbatch tutorials/submit_tutorials.sh --dry-run --only QE
#   sbatch tutorials/submit_tutorials.sh --only QE/9,QE/12
#   bash   tutorials/submit_tutorials.sh --dry-run          # no scheduler needed
#
# Every argument is passed straight through to the driver
# (htesp-tutorials / python -m tutorials.run_tutorials), so anything the driver
# accepts works here.
#
# EDIT BEFORE THE FIRST SUBMISSION -- the four #SBATCH lines below are
# placeholders.  They deliberately do NOT carry the partition and node names
# from the shipped examples (--partition=dense, -x dense001): those belong to
# one site and already appear in 40 example files.  Replace <ACCOUNT>,
# <PARTITION>, the wall time and the core count with your own, or delete a line
# entirely if your site does not use it.
#
# Optional environment:
#   HTESP_VENV   path to a virtualenv/conda prefix to activate (or its
#                bin/activate script).  Left unset, the current python is used.
#   HTESP_ROOT   the HTESP checkout.  Default: the parent of this script.
#   HTESP_OUT    where to write this run's output.  Default: $HTESP_ROOT/tutorial_runs_root
# -----------------------------------------------------------------------------
#SBATCH --job-name=htesp-tutorials
#SBATCH --account=<ACCOUNT>            # <-- EDIT (delete if your site has none)
#SBATCH --partition=<PARTITION>        # <-- EDIT
#SBATCH --nodes=1                      # <-- EDIT
#SBATCH --cpus-per-task=8              # <-- EDIT (matches --workers below)
#SBATCH --time=2-00:00:00              # <-- EDIT (the whole campaign is long)
#SBATCH --output=htesp-tutorials-%j.out
#SBATCH --error=htesp-tutorials-%j.err

set -o pipefail

# --- where things are --------------------------------------------------------
script_path="${BASH_SOURCE[0]}"
if [ -n "${SLURM_JOB_ID:-}" ] && [ ! -f "$script_path" ]; then
    # sbatch copies the script into its spool directory; recover the original
    found="$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
             | awk -F= '/Command=/{print $2; exit}')"
    [ -n "$found" ] && [ -f "$found" ] && script_path="$found"
fi
script_dir="$(cd "$(dirname "$script_path")" && pwd)"
HTESP_ROOT="${HTESP_ROOT:-$(dirname "$script_dir")}"
HTESP_OUT="${HTESP_OUT:-$HTESP_ROOT/tutorial_runs_root}"
mkdir -p "$HTESP_OUT"

# --- python environment ------------------------------------------------------
if [ -n "${HTESP_VENV:-}" ]; then
    if [ -f "$HTESP_VENV/bin/activate" ]; then
        # shellcheck disable=SC1091
        . "$HTESP_VENV/bin/activate"
    elif [ -f "$HTESP_VENV" ]; then
        # shellcheck disable=SC1091
        . "$HTESP_VENV"
    else
        echo "HTESP_VENV=$HTESP_VENV has no bin/activate -- ignoring it" >&2
    fi
fi
PYTHON="${PYTHON:-python3}"
export PYTHONPATH="$HTESP_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1

# --- the driver's own stdout goes somewhere predictable ----------------------
driver_log="$HTESP_OUT/driver-${SLURM_JOB_ID:-local}.log"

echo "=============================================================="
echo " HTESP tutorial runner"
echo "   host            : $(hostname)"
echo "   date            : $(date)"
echo "   job id          : ${SLURM_JOB_ID:-<not a slurm job>}"
echo "   HTESP_ROOT      : $HTESP_ROOT"
echo "   work directory  : $HTESP_OUT"
echo "   python          : $($PYTHON -V 2>&1) at $(command -v "$PYTHON")"
echo "   virtualenv      : ${VIRTUAL_ENV:-${HTESP_VENV:-<none>}}"
echo "   driver log      : $driver_log"
echo "   arguments       : $*"
echo "=============================================================="

if ! command -v squeue >/dev/null 2>&1; then
    echo "note: squeue is not on PATH -- only --dry-run can be trusted"
fi

# --- run ---------------------------------------------------------------------
# --workers defaults to the cores this job was given.
workers="${SLURM_CPUS_PER_TASK:-4}"
"$PYTHON" -m tutorials.run_tutorials \
    --workdir "$HTESP_OUT" \
    --workers "$workers" \
    "$@" 2>&1 | tee "$driver_log"
status=${PIPESTATUS[0]}

echo
echo "=============================================================="
echo " driver exited with status $status"
echo "   0 = everything finished, 1 = something failed or was blocked,"
echo "   2 = a preflight check failed, 130 = interrupted"
echo
echo " WHERE IT STOPPED is at the end of:"
echo "   $HTESP_OUT/report.md"
echo "   $HTESP_OUT/report.json   (the same thing as data)"
echo "   $HTESP_OUT/state.json    (the checkpoint --resume reads)"
echo "   $HTESP_OUT/logs/<tutorial>/<NN>-<step>.log"
echo
echo " resume with:"
echo "   sbatch $script_path --resume $*"
echo "=============================================================="
exit "$status"
