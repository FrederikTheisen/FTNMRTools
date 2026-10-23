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
The v2 analyzer infers R1, R1rho (T1rho), R2 or NOE from the peak-list filename;
include exactly one of those type tokens in the name. PEAK_FMT accepts
sparky3d or sparky2d. A legacy --mode argument is accepted only as a check
against the filename inference. The original relax_analyze.com remains an
independent legacy script with its own DATA_TYPE setting.
The default output prefix retains the original filename-before-first-dot
rule. Set OUTPUT_PREFIX explicitly for filenames containing other dots.

Peak lists require assignment, coordinate and overlap columns. Comments,
including compact noise/source_folder metadata, are ignored during table
extraction. Input peak lists are never rewritten. Each peak receives a
unique temporary tag accepted by data_xtr. Matching uses these tags, not
row order or assignment alone. Duplicate labels at different coordinates
remain separate; an identical assignment AND coordinate row is rejected.

R1, R1rho and R2 use rlx_dfit -mode exp2 -rms_noise. Filename inference picks
the starting rate: 1.0 for R1/R1rho and 10.0 for R2. The fitter still controls
parameter estimation and uncertainty rescaling. NOE keeps -mode noe -rms_noise;
supply the type header required by that fitter mode, rather than treating its
type codes as decay times.
No exponential fit or decay residual is generated for NOE.

The analyzer publishes three user-facing files:
    <peak-list-filename>_intensity.pdf
    <prefix>_numbers.csv
    <prefix>_support.zip

The PDF contains the plots and dataset assessment. The numbers CSV has one row
per peak, including assignments, coordinates, overlap/validity/warning flags,
selected rates and errors, fit statistics, any all-points versus
shortest-delay-excluded comparison, and every acquisition plane's delay,
intensity, input noise and fit-use flag. Repeated delays stay in acquisition
order. The support ZIP groups the detailed report table, diagnostics, saved
fitter inputs and outputs, and plotting records. Native fitter files are
preserved byte-for-byte inside that archive instead of being scattered beside
the PDF and numbers table.

Plots use the reported R and I0 without a second fit. Measured intensities,
input noise bars and predicted intensities are divided by fitted I0. The
residual panel is (measured - predicted) / input noise SD. These bars condition
on fitted I0; they do not include uncertainty in the fitted amplitude or
quantify overlap/processing systematics. Invalid or zero I0 receives a flagged
raw-intensity panel without a normalized curve or residuals. The rate/NOE
overview shows valid fitted peaks and marks overlaps.

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

Detailed diagnostics, the full measurement table, the selection audit, and
compatibility information are grouped inside the support ZIP rather than
published as separate files. The archive retains the fit comparisons, warning
text, thresholds, cohort and residual-bin records used to review a shortened
fit. Native artifacts are validated numerically without hashing.

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

The support ZIP keeps the original observations, per-plane fit mask, both
candidate fits/warnings and the selection audit. Report refresh validates the
saved selection and numerical contents without extraction or refitting. A
refresh never decides to exclude another delay.

Insufficient data or a rejected candidate keeps the full fit and records a
reason. A malformed candidate or failed external command aborts publication,
leaving previous outputs intact. Invalid non-cohort candidate fits remain
invalid if the experiment-wide candidate is selected; they never silently
borrow an all-point rate.
