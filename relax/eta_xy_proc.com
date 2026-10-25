#!/bin/tcsh
# Original nico_nnh_tr_noSE experiment only, never one of its TopSpin split copies.
# Tune these common phases on the reference plane; apply them to all four planes.
# TopSpin phase numbers are not automatically transferable to NMRPipe.
set H_P0 = 0
set H_P1 = 0
set N_P0 = 0
set N_P1 = 0
set H_ZF = 4096
set N_ZF = 1024

set SCRIPT_PATH = "$0"
if ("$SCRIPT_PATH" !~ /*) set SCRIPT_PATH = "$cwd/$SCRIPT_PATH"
python3 "$SCRIPT_PATH:h/eta_xy.py" process --experiment "$cwd" \
    --h-p0 "$H_P0" --h-p1 "$H_P1" --n-p0 "$N_P0" --n-p1 "$N_P1" \
    --h-zf "$H_ZF" --n-zf "$N_ZF"
exit $status
