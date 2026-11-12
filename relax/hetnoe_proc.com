#!/bin/tcsh
# Copy into the ORIGINAL interleaved Bruker experiment and run:
#     tcsh hetnoe_proc.com
# Applies ordinary Bruker split 2 semantics before converting each 2D dataset.
# Preserved Bruker split datasets: split/{reference,saturated}/
# Standalone 2D planes: ft/test001.ft2 reference, ft/test002.ft2 saturated.
# Tune phases on the reference; both planes use the same settings.

set H_P0 = 0
set H_P1 = 0
set N_P0 = 0
set N_P1 = 0
set BASELINE_ORDER = 4
set OUTDIR = .

set TOOLS_DIR = "$HOME/Documents/Projects/FTNMRTools/relax"
if ($?FTNMRTOOLS_ROOT) then
    set TOOLS_DIR = "$FTNMRTOOLS_ROOT/relax"
endif
set PROCESSOR = "$TOOLS_DIR/hetnoe_proc.py"
if (! -f "$PROCESSOR") then
    echo "Shared HetNOE processor not found: $PROCESSOR"
    echo "Set FTNMRTOOLS_ROOT to your FTNMRTools checkout."
    exit 2
endif

python3 "$PROCESSOR" process --experiment "$cwd" --output-dir "$OUTDIR" \
    --h-p0 "$H_P0" --h-p1 "$H_P1" --n-p0 "$N_P0" --n-p1 "$N_P1" \
    --h-baseline-order "$BASELINE_ORDER" --n-baseline-order "$BASELINE_ORDER" $argv:q
exit $status
