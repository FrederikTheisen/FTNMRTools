Relaxation processing and analysis
=================================

Available scripts:

    relax_analyze.com          Original analyzer, restored unchanged.
    relax_analyze_v2.com       Revised analyzer with validation and diagnostics.
    relax_proc.com             Updated processing script; scientific settings retained.
    relax_proc_original.com    Original processing template, preserved unchanged.

The original scripts are byte-for-byte copies from commit cbbd3e9, before the
relaxation fixes. They remain independent of relax_tools.py. The descriptions
below apply to updated relax_proc.com and relax_analyze_v2.com.

Keep the updated scripts and relax_tools.py together. Run the .com files
from the experiment directory:

    csh /path/to/FTNMRTools/relax/relax_proc.com
    tcsh /path/to/FTNMRTools/relax/relax_analyze_v2.com

Bare filenames, relative paths and absolute paths are supported, including
paths with spaces. For a local copied script, `tcsh relax_proc.com` or
`./relax_proc.com` also works. Keep relax_tools.py beside each copied script
or entrypoint symlink. A missing helper stops before processing begins.

Requirements: Python 3.8+, tcsh, the existing NMRPipe commands, data_xtr,
rlx_dfit and gnuplot with pdfcairo support. Python uses only its standard
library. The controller inherits your configured command PATH and NMRTXT.

These scripts do not automatically install themselves in dataset folders.
Copying or using an updated script does not run it. Running relax_proc.com
processes spectra; running relax_analyze_v2.com extracts intensities and fits.
Review each script's editable settings before running it on an acquisition.
The two analyzers use the same output filenames; changing the entrypoint name
does not create a separate set of result filenames.

Processing settings and delays
------------------------------

The default Bruker conversion and interleaved splitting retain the original
template. Solvent suppression, apodization, Fourier transformation,
extraction, phasing, zero filling and POLY -auto settings are unchanged.
These settings are acquisition-specific. The former series.com calls are
expanded into the same sethdr settings, with each update checked individually.

The original conversion-output selection is preserved:

    If conversion writes test.fid, the original COADD selection and
    bruk_ranceY.M shuffling split the 2D interleaved FID into planes.

    If conversion writes fid/test%03d.fid directly, the script uses those
    planes without COADD, additional shuffling or FID header rewriting.
    The converter must arrange H/N in X/Y and a real delay axis in Z, and
    must handle any gradient shuffling.

No new layout setting is needed. Producing both output layouts, or neither,
is an error. Both routes create fid/ and ft/ before conversion. Direct-plane
inputs are checked for exact filenames/counts against the delay list,
compatible headers, untransformed X/Y axes, a real Z axis and complete payload
sizes before Fourier processing. Delays keep acquisition order and repeats.

For example, a 3D conversion with -yLAB Tau -zLAB 15N can transpose its stream
into HN/15N/Tau planes before pipe2xyz -x -out ./fid/test%03d.fid. Keep that
conversion's MULT, TP and ZTP operations. Echo-AntiEcho/Rance-Kay conversion
already invokes the appropriate gradient shuffling; do not apply it again.
See the NMRPipe conversion/shuffling documentation:
https://www.nmrscience.com/ref/nmrpipe/coadd.html

Configure path2vc, path2vd and path2ncyc at the top of relax_proc.com. Paths
are relative to the experiment directory or absolute; patterns are accepted
only when they identify one file. For example, path2vd may be
"lists/vd/15N_T1rho_9_1-200.vd". delay_source=auto requires exactly one
available source. Set delay_source to vc, vd or ncyc to select explicitly.
No source or multiple matches is an error; old delaylist files are not used.

vdlist values accept bare numbers or s for seconds, m/ms for milliseconds,
and u/us/µs for microseconds, including scientific notation. Blank lines
and whole-line # comments are ignored. Order and repeated delays are kept.
Values must be finite and nonnegative.

vclist retains the original formula:
    delay_constant + count * delay_block

ncyc_list requires an explicit cycle_seconds setting and uses:
    cycle_constant + count * cycle_seconds

Counts must be nonnegative integers. Never infer a cycle duration from an
instrument name or folder name. Match it to the pulse sequence.

Processing still writes delaylist, relax_header.txt, fid/, ft/, m.ft3,
proc.ucsf and vEstNoise.txt. Per-plane noise is estimated with specStat.com
and retains the plane/delay order. A single peak-list noise comment is not
substituted for the per-plane noise vector. NMRPipe series header metadata,
plane counts, delay counts and positive finite noise SDs are checked.

Analysis and peak identity
--------------------------

Set SPARKY_PEAKS to one filename or a pattern matching exactly one list.
DATA_TYPE accepts r1, r2 or noe; PEAK_FMT accepts sparky3d or sparky2d.
The default output prefix retains the original filename-before-first-dot
rule. Set OUTPUT_PREFIX explicitly for filenames containing other dots.

Peak lists require assignment, coordinate and overlap columns. Comments,
including compact noise/source_folder metadata, are ignored during table
extraction. Input peak lists are never rewritten. Each peak receives a
unique temporary tag accepted by data_xtr. Matching uses these tags, not
row order or assignment alone. Duplicate labels at different coordinates
remain separate; an identical assignment AND coordinate row is rejected.

R1 and R2 keep rlx_dfit -mode exp2 -rms_noise and their original starting
rates. The fitter still controls parameter estimation and uncertainty
rescaling. NOE keeps -mode noe -rms_noise; supply the type header required
by that fitter mode, rather than treating its type codes as decay times.
No exponential fit or decay residual is generated for NOE.

Existing result names remain:
    <prefix>_results.out
    <prefix>_results_peak_overlap.out
    <prefix>_results_filtered.out
    data_xtr.out, rlx_dfit.inp, rlx_dfit.out, intensity.out, plot.gnu
    <peak-list-filename>_intensity.pdf

The three-column result files retain residue, fitted value and error.
Filtered results exclude overlaps and invalid numerical fits. Duplicate
labels and other quality warnings are flagged without automatic exclusion.
The fitter output rlx_dfit.out is preserved verbatim with its temporary
unique tags; <prefix>_fit_diagnostics.tsv maps them to original assignments
and coordinates. data_xtr.out and intensity.out restore original labels.

The diagnostics TSV also contains amplitude/error, chi-square and its
reported threshold, overlap, validity and warning flags. Native warning
text is preserved, including any malformed number formatting. A rescaling
warning is not silently interpreted as an exact numeric scale factor.

Plots use the reported R and I0 without a second fit. Measured intensities,
input noise bars and predicted intensities are divided by fitted I0.
The residual panel is (measured - predicted) / input noise SD. These bars
condition on fitted I0; they do not include uncertainty in the fitted
amplitude or quantify overlap/processing systematics. Invalid or zero I0
receives a flagged raw-intensity panel without a normalized curve or
residuals. The rate/NOE overview shows valid fitted peaks and marks overlaps.

intensity.out retains time and plotted intensity as its first two columns.
For decay data its additional columns are plotted noise SD, plotted
prediction, plotted residual and residual/input-noise-SD. Invalid curves
use NaN for prediction and residuals. For NOE the columns are plane index,
raw intensity and per-plane input noise SD.

Failure handling
----------------

Each invocation uses a new hidden temporary directory beside its outputs.
The two scripts share a .relax.lock to prevent overlapping runs in the same
folder. Inputs stay at their original paths; the processing stage contains a symlink
to ser, which the conversion reads while writing its outputs in that stage.

Every external command and every stage of the script's pipelines is checked.
Missing, empty or inconsistent products fail the run. Old fit files and
leftover temporary directories cannot supply inputs to a new run. Outputs
are published only when the entire run, including plotting, succeeds.
Publication errors roll back replacements. Output symlinks are rejected.
The lock file may remain after a run; the OS releases its lock automatically.
The worker verifies that its current directory is the controller's temporary
stage using directory identity, so /home and /data aliases of the same folder
are accepted. A different directory, missing path or file remains an error.

Tests without processing data
-----------------------------

    PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s relax/tests -v

Run this from the FTNMRTools repository root. All processing, extraction,
fitting and plotting commands in the integration tests are stubs operating
only in temporary folders on synthetic files. The fixture native_exp2.out
contains a small excerpt of an existing fitter output for parser validation;
it does not trigger a fit. No dataset or CCPN project is opened or modified.
