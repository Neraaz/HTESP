#!/bin/bash
cat > vasp-phonopy.in <<_EOF_
PREC Accurate
IBRION -1
ISMEAR 0
SIGMA 0.05
IALGO 38
LREAL Auto
LWAVE .FALSE.
LCHARG .FALSE.
NSW 0
ISIF
EDIFFG
_EOF_

exit 0
