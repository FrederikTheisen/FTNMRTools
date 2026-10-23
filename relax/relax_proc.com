#!/bin/csh

# Editable sequence-specific timing and delay-list settings.
set delay_block = 0.01
set delay_constant = 0.01
set path2vc = "vclist"
set path2vd = "vdlist"
set path2ncyc = "ncyc_list"
set delay_source = auto # auto, vc, vd, ncyc; auto requires exactly one list
# ncyc_list has no implicit units. Configure seconds per cycle before selecting it.
set cycle_seconds = ""
set cycle_constant = 0

# Prefix bare/relative filenames before :h; tcsh leaves a bare name unchanged.
set SCRIPT_PATH = "$0"
if ("$SCRIPT_PATH" !~ /*) set SCRIPT_PATH = "$cwd/$SCRIPT_PATH"
set HELPER = "$SCRIPT_PATH:h/relax_tools.py"
if (! -f "$HELPER") then
    echo "Missing helper: $HELPER"
    echo "Keep relax_tools.py beside this script (and any entrypoint symlink)."
    exit 2
endif
if ($#argv == 0) then
    python3 "$HELPER" run-proc --script "$SCRIPT_PATH" --source "$cwd"
    exit $status
endif
if ($#argv != 2) then
    echo "Usage: csh relax_proc.com"
    exit 2
endif
if ("$argv[1]" != "--worker") then
    echo "Usage: csh relax_proc.com"
    exit 2
endif
if (! $?FTNMR_RELAX_STAGE) then
    echo "Internal worker must be started by the controller."
    exit 2
endif
# tcsh may use a /home alias while Python uses /data for the same directory.
# Check directory identity, preserving the guard against a different location.
python3 -c 'import os, sys; sys.exit(not os.path.isdir(sys.argv[1]) or not os.path.samefile(sys.argv[1], os.getcwd()))' "$FTNMR_RELAX_STAGE"
if ($status != 0) then
    echo "Worker directory differs from its temporary stage."
    exit 2
endif
set SOURCE = "$argv[2]"
# The controller starts this worker in a new temporary directory. A pipeline
# fails if ANY command fails, even when its final command returns success.
set anyerror
onintr failed

python3 "$HELPER" delays --base "$SOURCE" --vc "$path2vc" --vd "$path2vd" \
    --ncyc "$path2ncyc" --source "$delay_source" --block "$delay_block" \
    --constant "$delay_constant" --cycle-seconds "$cycle_seconds" --cycle-constant "$cycle_constant"
if ($status != 0) goto failed
set vList = (`cat delaylist`)
if ($status != 0) goto failed
set n = $#vList

# Bruker conversion settings: edit for the acquisition, as before.
bruk2pipe -verb -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 2760 -dspfvs 20 -grpdly 67.9857330322266  \
  -xN              1280  -yN              2200  \
  -xT               640  -yT              1100  \
  -xMODE            DQD  -yMODE  Complex  \
  -xSW         7246.377  -ySW         1515.152  \
  -xOBS         600.483  -yOBS          60.853  \
  -xCAR           4.771  -yCAR         118.577  \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -ov

if ($status != 0) goto failed
python3 "$HELPER" check-fid test.fid delaylist
if ($status != 0) goto failed
mkdir fid ft
if ($status != 0) goto failed
set cList = ()
foreach v ($vList)
    set cList = ($cList 0.0)
end
set i = 1
foreach v ($vList)
    set cList[$i] = 1.0
    set outName = (`printf fid/test%03d.fid $i`)
    nmrPipe -in test.fid \
    | nmrPipe -fn COADD -time -axis Y -cList $cList -verb \
    | nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr \
      -out $outName -ov
    if ($status != 0) goto failed
    set cList[$i] = 0.0
    @ i++
end
python3 "$HELPER" series fid/test*.fid
if ($status != 0) goto failed

# Processing choices are unchanged; the same scheme is used for every plane.
set out = proc.ucsf
set i = 1
foreach tau ($vList)
    set inName = (`printf fid/test%03d.fid $i`)
    set outName = (`printf ft/test%03d.ft2 $i`)
  nmrPipe -in $inName \
  | nmrPipe  -fn SOL                             \
  | nmrPipe  -fn SP -off 0.4 -end 0.98 -pow 2 -c .5    \
  | nmrPipe  -fn ZF -size 4096                          \
  | nmrPipe  -fn FT                                     \
  | nmrPipe  -fn EXT -left -sw                               \
  | nmrPipe  -fn PS -p0 0 -p1 0 -di                     \
  | nmrPipe  -fn TP                                     \
  | nmrPipe  -fn SP -off 0.4 -end 0.98 -pow 2 -c .5    \
  | nmrPipe  -fn ZF -size 1024                           \
  | nmrPipe  -fn FT                                      \
  | nmrPipe  -fn PS -p0 0 -p1 0 -di                 \
  | nmrPipe  -fn TP                                     \
  | nmrPipe  -fn POLY -auto -ord 4                      \
     -verb -ov -out $outName

    if ($status != 0) goto failed
    sethdr $outName -tau $tau
    if ($status != 0) goto failed
    @ i++
end
python3 "$HELPER" series ft/test*.ft2
if ($status != 0) goto failed
xyz2pipe -in ft/test%03d.ft2 > m.ft3
if ($status != 0) goto failed
pipe2ucsf m.ft3 $out
if ($status != 0) goto failed
ucsfdata -a1 15N -o1 10. -sw1 100. -f1 100. -a2 15N $out
if ($status != 0) goto failed

foreach spec (ft/test*.ft2)
    sethdr $spec -ndim 2 -zN 1 -zT 1 -zMODE Real
    if ($status != 0) goto failed
    setfdata $spec -fileCount 1
    if ($status != 0) goto failed
end

# Estimate each processed plane separately, in the same order as the delays.
foreach spec (ft/test*.ft2)
    specStat.com -in $spec -stat vEstNoise > noise_stat.txt
    if ($status != 0) goto failed
    python3 "$HELPER" noise-value noise_stat.txt >> noise_values.txt
    if ($status != 0) goto failed
end
set noiseList = (`cat noise_values.txt`)
if ($status != 0) goto failed
echo "$n $noiseList" > vEstNoise.txt
if ($status != 0) goto failed
python3 "$HELPER" validate-proc
if ($status != 0) goto failed
exit 0

failed:
echo "Processing failed; previous outputs have not been replaced."
exit 1
