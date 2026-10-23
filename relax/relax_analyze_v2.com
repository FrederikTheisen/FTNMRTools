#!/bin/tcsh

### RELAXATION DATA ANALYSIS V2 — editable experiment settings ###
# Quote a filename or pattern. Exactly one peak list must match.
set SPARKY_PEAKS = "*.list"
set DATA_FILE = ./m.ft3
set DATA_FMT = nmrpipe
set DATA_TYPE = r1 # r1, r2, noe
set PEAK_FMT = sparky3d # sparky3d or sparky2d
set RLX_DFIT_HEADER = relax_header.txt
set RLX_DFIT_ESD = vEstNoise.txt
# Empty preserves the existing prefix (filename before the first dot).
set OUTPUT_PREFIX = ""

# Keep relax_tools.py beside this script. It stages and validates the complete run.
set HELPER = "$0:h/relax_tools.py"
python3 "$HELPER" analyze --peaks "$SPARKY_PEAKS" --data "$DATA_FILE" \
    --format "$DATA_FMT" --mode "$DATA_TYPE" --peak-format "$PEAK_FMT" \
    --header "$RLX_DFIT_HEADER" --noise "$RLX_DFIT_ESD" --prefix "$OUTPUT_PREFIX"
exit $status
