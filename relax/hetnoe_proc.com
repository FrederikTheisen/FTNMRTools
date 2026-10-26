#!/bin/tcsh
# Process the pseudo-3D NMRPipe FIDs made by the Bruker `bruker` conversion.
# Run from the original experiment directory after running fid.com.
# H and 15N phases are shared by both HetNOE planes; tune on the reference.

set H_P0 = 0
set H_P1 = 0
set N_P0 = 0
set N_P1 = 0
set BASELINE_ORDER = 4
set OUTDIR = hetnoe_processed

if (! -e fid/test001.fid) then
    echo "Missing fid/test001.fid. Run the Bruker-generated fid.com first."
    exit 1
endif

mkdir -p "$OUTDIR"
if ($status != 0) exit 1

xyz2pipe -in fid/test%03d.fid \
| nmrPipe -fn SOL \
| nmrPipe -fn SP -off 0.5 -end 0.95 -pow 2 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -zf 1 -auto \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 $H_P0 -p1 $H_P1 -di \
| nmrPipe -fn EXT -left -sw \
| nmrPipe -fn TP \
| nmrPipe -fn SP -off 0.5 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -zf 1 -auto \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 $N_P0 -p1 $N_P1 -di \
| nmrPipe -fn POLY -auto -ord $BASELINE_ORDER \
| nmrPipe -fn TP \
| nmrPipe -fn POLY -auto -ord $BASELINE_ORDER \
| pipe2xyz -out "$OUTDIR/test%03d.ft3" -ov

if ($status != 0) then
    echo "NMRPipe processing failed. Check the first command that reported an error."
    exit 1
endif

echo "Processed NMRPipe plane series: $OUTDIR/test%03d.ft3"
echo "15N baseline correction is applied before the final transpose."
