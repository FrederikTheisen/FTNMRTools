#!/bin/tcsh
# Run from the experiment directory after processing and exporting peak positions.
# test001 = reference; test002 = saturated. Intensities are read from both.
set PEAKS = "*.list"
set PEAK_FORMAT = sparky3d # Your export has w1 (plane), w2 (N), w3 (H).
set TOOLS_DIR = "$HOME/Documents/Projects/FTNMRTools/relax"
if ($?FTNMRTOOLS_ROOT) then
    set TOOLS_DIR = "$FTNMRTOOLS_ROOT/relax"
endif
set HELPER = "$TOOLS_DIR/relax_tools.py"

if (! -f "$HELPER" || ! -f ft/test001.ft2 || ! -f ft/test002.ft2) then
    echo "Need the shared FTNMRTools analyzer and ft/test001.ft2, ft/test002.ft2."
    exit 1
endif

# Assemble a temporary series; keep the displayed spectra as independent 2Ds.
set work = `mktemp -d /tmp/hetnoe-analysis.XXXXXX`
if ($status != 0) exit 1
onintr failed
cp ft/test001.ft2 ft/test002.ft2 "$work/"
if ($status != 0) goto failed
set noiseList = ()
foreach plane (001 002)
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
xyz2pipe -in "$work/test%03d.ft2" >! "$work/m.ft3"
if ($status != 0) goto failed

# The shared analyzer extracts at identical coordinates in both spectra,
# calculates saturated/reference with noise errors, and reports overlap flags.
python3 "$HELPER" analyze --peaks "$PEAKS" --peak-format "$PEAK_FORMAT" \
    --data "$work/m.ft3" --format nmrpipe --mode noe \
    --header "$work/relax_header.txt" --noise "$work/vEstNoise.txt"
if ($status != 0) goto failed
/bin/rm -rf "$work"
exit 0

failed:
echo "Analysis failed; temporary inputs retained in $work"
exit 1
