#!/bin/tcsh
# Run from the original nico_nnh_tr_noSE experiment folder.
# One list, exported from a reference plane. Two- and three-coordinate lists work.
set SPARKY_PEAKS = "*.list"
set DATA_FILE = "eta_xy_processed/m.ft3"
set NOISE_FILE = "eta_xy_processed/vEstNoise.txt"
set OUTPUT_PREFIX = ""

# Use the shared FTNMRTools checkout; keep experiment settings in this script.
set SCRIPT_PATH = "$0"
if ("$SCRIPT_PATH" !~ */*) then
    if (! -f "$SCRIPT_PATH") set SCRIPT_PATH = "`which "$SCRIPT_PATH"`"
endif
if ("$SCRIPT_PATH" !~ /*) set SCRIPT_PATH = "$cwd/$SCRIPT_PATH"
set TOOLS_DIR = "$HOME/Documents/Projects/FTNMRTools/relax"
if ($?FTNMRTOOLS_ROOT) then
    set TOOLS_DIR = "$FTNMRTOOLS_ROOT/relax"
endif
set HELPER = "$TOOLS_DIR/eta_xy.py"
if (! -f "$HELPER") then
    echo "Shared FTNMRTools helper not found: $HELPER"
    echo "Set FTNMRTOOLS_ROOT to your FTNMRTools checkout."
    exit 2
endif
python3 "$HELPER" analyze --experiment "$cwd" \
    --peaks "$SPARKY_PEAKS" --data "$DATA_FILE" --noise "$NOISE_FILE" \
    --prefix "$OUTPUT_PREFIX"
exit $status
