#!/bin/tcsh
# Run from the experiment directory after processing and exporting peak positions.
# test001 = reference; test002 = saturated. Intensities are read from both.
set PEAKS = "*.list"
set PEAK_FORMAT = sparky3d # Your export has w1 (plane), w2 (N), w3 (H).
set REF_FILE = ft/test001.ft2
set NOE_FILE = ft/test002.ft2
set OUTPUT_PREFIX = ""
set TOOLS_DIR = "$HOME/Documents/Projects/FTNMRTools/relax"
if ($?FTNMRTOOLS_ROOT) then
    set TOOLS_DIR = "$FTNMRTOOLS_ROOT/relax"
endif
set HELPER = "$TOOLS_DIR/relax_tools.py"

if (! -f "$HELPER" || ! -f "$REF_FILE" || ! -f "$NOE_FILE") then
    echo "Need the shared FTNMRTools analyzer and $REF_FILE, $NOE_FILE."
    exit 1
endif

# Assemble a temporary series; keep the displayed spectra as independent 2Ds.
set work = `mktemp -d /tmp/hetnoe-analysis.XXXXXX`
if ($status != 0) exit 1
onintr failed
cp "$REF_FILE" "$work/test001.ft2"
if ($status != 0) goto failed
cp "$NOE_FILE" "$work/test002.ft2"
if ($status != 0) goto failed
set noiseList = ()
foreach plane (001 002)
    sethdr "$work/test${plane}.ft2" -ndim 2 -zN 1 -zT 1 -zMODE Real
    if ($status != 0) goto failed
    setfdata "$work/test${plane}.ft2" -fileCount 1
    if ($status != 0) goto failed
    specStat.com -in "$work/test${plane}.ft2" -stat vEstNoise >! "$work/noise.txt"
    if ($status != 0) goto failed
    set sigma = `python3 "$HELPER" noise-value "$work/noise.txt"`
    if ($status != 0) goto failed
    set noiseList = ($noiseList $sigma)
end
echo "2 0 1" >! "$work/relax_header.txt"
echo "2 $noiseList" >! "$work/vEstNoise.txt"
python3 "$HELPER" series "$work/test001.ft2" "$work/test002.ft2"
if ($status != 0) goto failed
foreach plane (001 002)
    sethdr "$work/test${plane}.ft2" -zFD 2
    if ($status != 0) goto failed
end
xyz2pipe -in "$work/test%03d.ft2" >! "$work/m.ft3"
if ($status != 0) goto failed

# The shared analyzer extracts at identical coordinates in both spectra,
# calculates saturated/reference with noise errors, and reports overlap flags.
python3 "$HELPER" analyze --peaks "$PEAKS" --peak-format "$PEAK_FORMAT" \
    --data "$work/m.ft3" --format nmrpipe --mode noe \
    --header "$work/relax_header.txt" --noise "$work/vEstNoise.txt" --prefix "$OUTPUT_PREFIX"
if ($status != 0) goto failed
/bin/rm -rf "$work"
exit 0

failed:
echo "Analysis failed; temporary inputs retained in $work"
exit 1
