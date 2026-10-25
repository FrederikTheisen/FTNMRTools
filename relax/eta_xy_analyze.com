#!/bin/tcsh
# Run from the original nico_nnh_tr_noSE experiment folder.
# One list, exported from a reference plane. Two- and three-coordinate lists work.
set SPARKY_PEAKS = "*.list"
set DATA_FILE = "eta_xy_processed/m.ft3"
set NOISE_FILE = "eta_xy_processed/vEstNoise.txt"
set OUTPUT_PREFIX = ""

set SCRIPT_PATH = "$0"
if ("$SCRIPT_PATH" !~ /*) set SCRIPT_PATH = "$cwd/$SCRIPT_PATH"
python3 "$SCRIPT_PATH:h/eta_xy.py" analyze --experiment "$cwd" \
    --peaks "$SPARKY_PEAKS" --data "$DATA_FILE" --noise "$NOISE_FILE" \
    --prefix "$OUTPUT_PREFIX"
exit $status
