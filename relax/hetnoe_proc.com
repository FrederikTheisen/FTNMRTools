#!/bin/tcsh
# Run from the ORIGINAL interleaved HetNOE experiment, not a TopSpin split copy.
# Tune phases on the reference plane, then use the SAME settings for both planes.
# Zero phase is a starting setting, not a validated phase for every dataset.
set H_P0 = 0
set H_P1 = 0
set N_P0 = 0
set N_P1 = 0
set H_ZF = 4096
set N_ZF = 1024
# N correction runs with 15N active, before transposing back to H.
set N_BASELINE_ORDER = 4
set H_BASELINE_ORDER = 4

set SCRIPT_PATH = "$0"
if ("$SCRIPT_PATH" !~ */*) then
    if (! -f "$SCRIPT_PATH") set SCRIPT_PATH = "`which "$SCRIPT_PATH"`"
endif
if ("$SCRIPT_PATH" !~ /*) set SCRIPT_PATH = "$cwd/$SCRIPT_PATH"
# Resolve entrypoint symlinks so helpers can stay in the source repository.
set HELPER = `python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve().with_name("hetnoe_proc.py"))' "$SCRIPT_PATH"`
if (! -f "$HELPER") then
    echo "Missing helper: $HELPER"
    exit 2
endif
python3 "$HELPER" process --experiment "$cwd" \
    --h-p0 "$H_P0" --h-p1 "$H_P1" --n-p0 "$N_P0" --n-p1 "$N_P1" \
    --h-zf "$H_ZF" --n-zf "$N_ZF" \
    --n-baseline-order "$N_BASELINE_ORDER" --h-baseline-order "$H_BASELINE_ORDER" $argv:q
exit $status
