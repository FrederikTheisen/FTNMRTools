Relaxation processing and analysis
=================================

Keep relax_proc.com, relax_analyze.com and relax_tools.py together. Run the
.com files from the experiment directory, as before:

    csh /path/to/FTNMRTools/relax/relax_proc.com
    tcsh /path/to/FTNMRTools/relax/relax_analyze.com

Requirements: Python 3.8+, tcsh, the existing NMRPipe commands, data_xtr,
rlx_dfit and gnuplot with pdfcairo support. Python uses only its standard
library. The controller inherits your configured command PATH and NMRTXT.

These scripts do not automatically install themselves in dataset folders.
Copying or using an updated script does not run it. Running relax_proc.com
processes spectra; running relax_analyze.com extracts intensities and fits.
Review each script's editable settings before running it on an acquisition.

Processing settings and delays
------------------------------

The Bruker conversion, splitting, solvent suppression, apodization, Fourier
transformation, extraction, phasing, zero filling and POLY -auto settings
remain the same as the original template. They are acquisition-specific.
The former series.com calls are expanded into the same sethdr settings,
with each header update checked individually.

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

Tests without processing data
-----------------------------

    PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s relax/tests -v

Run this from the FTNMRTools repository root. All processing, extraction,
fitting and plotting commands in the integration tests are stubs operating
only in temporary folders on synthetic files. The fixture native_exp2.out
contains a small excerpt of an existing fitter output for parser validation;
it does not trigger a fit. No dataset or CCPN project is opened or modified.


Automatic shortest-delay test (v2)
---------------------------------

relax_analyze_v2.com now sets EARLY_POINT_POLICY = auto. The helper accepts
--early-point-policy auto|report|off. auto can select a shortened fit; report
runs the same checks and candidate fit but keeps all-point rates; off keeps
all-point rates without attempting a candidate. NOE bypasses this feature.
The independent original relax_analyze.com remains unchanged.

The shortest numerical delay is tested, not the first acquisition plane.
All replicates at that delay are excluded together. This is a single,
experiment-wide decision, with no iterative trimming or per-peak fallback.
All fits continue to use rlx_dfit exp2. The helper adds no Python dependencies.

The assessment cohort is frozen from positive, finite all-point fits without
overlap flags and with exactly one peak per residue. At least 20 residues
and five distinct remaining delays are required. Repeats are averaged within
each residue/delay with inverse-noise-variance weights. The median across
residues is normalized by the fixed all-point I0. Input noise is never
replaced by the fitter's rescaled uncertainties for selection decisions.

The shortest-delay residual must exceed both 0.5 percent of I0 and three
times the median normalized noise, with 80 percent residual-sign agreement.
The candidate is accepted only if the RMS structure score on the SAME
retained delay bins improves at least 50 percent, every retained median
lies within its limit, and median absolute fractional rate change is at
most 20 percent. These named defaults live in EARLY_DEFAULTS and are recorded
in every audit. They are diagnostic heuristics, not evidence of a physical
cause or a validated temperature/sequence correction.

Existing three-column result filenames now carry comment headers naming
policy, fit basis, decision, excluded delay and original 1-based planes,
and original/retained counts. The diagnostics TSV records both rates/errors,
all fitted parameters/warnings, rate change, and selected fit basis. The
machine-readable <prefix>_early_point_summary.json records thresholds,
cohort IDs, residual bins, comparison scores and reasons. Native artifacts
are validated numerically without hashing.

Native artifacts for each run:
  rlx_dfit_all_points.inp/.out and rlx_dfit_all_points_noise.txt
  rlx_dfit_without_shortest.inp/.out and rlx_dfit_without_shortest_noise.txt
  rlx_dfit.inp/.out and rlx_dfit_noise.txt (the selected, matching trio)
If no candidate is attempted, its three files contain an explicit comment
instead of leaving a previous run's candidate looking current. Native fitter
outputs for attempted fits are retained verbatim. The processing noise file
vEstNoise.txt and relaxation headers keep ALL acquisition observations.
For manual reproduction use the selected trio, not the original full-length
noise file with a shortened fitter input.

data_xtr.out and intensity.out retain every observation in acquisition
order. For v2 decay analyses intensity.out appends a seventh used_in_fit
column (1/0); the first six columns retain their meaning. The PDF records the
selection, compares aggregate residuals on retained observations, and shows
excluded points with open orange symbols. Original-noise goodness-of-fit
statistics use only retained points and N_used-2 degrees of freedom; these
are conditional diagnostics after a data-driven selection.

CSV schema 2 is self-contained: it includes the original observations,
per-observation mask, both fits/warnings and the decision. Existing schema-1
CSVs remain supported. Both CSV and --saved-prefix report refresh validate
and retain the recorded selection without extraction or refitting. Inconsistent masks, decisions, fitted
values, observations or archived numerical inputs are rejected. A report refresh
never decides to exclude another delay.

Insufficient data or a rejected candidate keeps the full fit and records a
reason. A malformed candidate or failed external command aborts publication,
leaving previous outputs intact. Invalid non-cohort candidate fits remain
invalid if the experiment-wide candidate is selected; they never silently
borrow an all-point rate.
