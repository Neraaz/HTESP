#!/bin/sh
# HTESP compatibility shim.
#
# `vasp-phonopy.sh` used to be a bash script in src/bash/.  Its logic now lives in
# htesp/workflow.py as HTESPWorkflow.vasp_phonopy_template(), which is parallel by default,
# never leaves the process in the wrong directory and records the SLURM job id.
# The argument contract is unchanged: start end trackfile [extra], `end`
# exclusive.  The original script is kept, unmodified, at legacy/bash/vasp-phonopy.sh.
#
# Extra options this shim accepts: --workers N, --dry-run, --root DIR, -v.
exec "${HTESP_PYTHON:-python3}" -m htesp.workflow "vasp-phonopy.sh" "$@"
