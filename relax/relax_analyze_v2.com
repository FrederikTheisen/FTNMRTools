#!/bin/tcsh

### RELAXATION DATA ANALYSIS V2 — editable experiment settings ###
# Quote a filename or pattern. Exactly one peak list must match.
set SPARKY_PEAKS = "*.list"
set DATA_FILE = ./m.ft3
set DATA_FMT = nmrpipe
set PEAK_FMT = sparky3d # sparky3d or sparky2d
set RLX_DFIT_HEADER = relax_header.txt
set RLX_DFIT_ESD = vEstNoise.txt
# Empty preserves the existing prefix (filename before the first dot).
set OUTPUT_PREFIX = ""
set EARLY_POINT_POLICY = auto # auto, report (keep full-fit rates), off, exclude_shortest

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
set HELPER = "$TOOLS_DIR/relax_tools.py"
if (! -f "$HELPER") then
    echo "Shared FTNMRTools helper not found: $HELPER"
    echo "Set FTNMRTOOLS_ROOT to your FTNMRTools checkout."
    exit 2
endif
python3 "$HELPER" analyze --peaks "$SPARKY_PEAKS" --data "$DATA_FILE" \
    --format "$DATA_FMT" --peak-format "$PEAK_FMT" \
    --header "$RLX_DFIT_HEADER" --noise "$RLX_DFIT_ESD" --prefix "$OUTPUT_PREFIX" \
    --early-point-policy "$EARLY_POINT_POLICY"
exit $status
