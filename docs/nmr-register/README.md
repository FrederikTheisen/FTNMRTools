# NMR register

The pilot keeps sample and experiment metadata in one private Grist document and
uses a Linux hourly scan to reconcile transferred Bruker data. Researchers edit
samples and assessments in Grist from either Mac or Linux. The collector records
actual acquisition facts in the canonical experiment fields after a unique,
exact match. Planned or unmatched registrations keep their entered values;
identity and sample-tag evidence stays separate. The collector never changes
assessments, notes, or sample links. After validating the exact observed
`AcquisitionKey` and all identity, fingerprint, and copy guards, it fills blank
`SourceID`, `Folder`, and `Expno` components on new or uniquely matched rows.
These identity fields are fill-only, not source-owned overwrite fields: nonblank
researcher values remain unchanged, and disagreements require reconciliation.
An explicitly observed empty folder represents a source-root acquisition. The
hourly cached-observation check also repairs blank identities without a rescan.

The authoritative pilot is [the team register](https://theisen.getgrist.com/3oWrbb4osxG8/registernmrtest):
[Electrostatics](https://theisen.getgrist.com/3oWrbb4osxG8/registernmrtest/p/12),
[Samples](https://theisen.getgrist.com/3oWrbb4osxG8/registernmrtest/p/9), and
[Experiments](https://theisen.getgrist.com/3oWrbb4osxG8/registernmrtest/p/8).

## Existing-data reconciliation

Instrument files define acquired experiments and acquisition metadata. The
workbook is mainly evidence for sample registration and explicit sample links;
its dates, EXPNOs, types and statuses are not an authoritative experiment catalog.
Workbook typos and duplicate rows alone do not block an identified acquisition.
Preserve current human assessments and frozen conditions; retain genuine sample,
source-identity and assessment conflicts for review.

The 2026-11-04 desktop scan found 1,081 acquisitions; all have linked Grist
records. Reviewed migrations reconciled 79 duplicate workbook entries with
those linked records, preserving their assessments, notes and sample links.
After the 2026-11-05 cleanup, the register contains 1,128 experiment records,
including 47 unlinked workbook entries. Their `MigrationIssue` begins with a specific action shown
in Needs attention. These are unresolved workbook identities or conflicts,
not a hidden set of desktop acquisitions waiting for import.

Keep using Grist for samples, registration and assessments. For new instrument
data, copy the complete experiment folders into the configured desktop data
root. The hourly collector discovers them while the desktop session is running;
`nmr-register sync` runs the same process immediately. Registering in Grist
first also works: enter the exact folder and EXPNO so synchronization can join
the planned entry to its data. An acquisition title can supply the sample ID.

Do not rerun the original workbook import to resolve these duplicate rows.
The originals remain read-only migration references. Reconciliation uses explicit
reviewed source/target pairs, a full backup, fresh field guards, and one atomic
Grist action bundle. Conflicting effective assessments and sample conditions
block a merge. The compact operational report is
`~/.local/share/nmr-register/migration-reconciliation.json`; Grist remains the
working register. Full backups and internally stored attachment bytes were
verified by a restore after this pass.

## Daily use

The Grist document has four working views: Experiments with a card for one
record; Samples with series and one tube-rack view (numbered rows 1–12 and
lettered columns A–H);
Needs attention; and Overview. Use the Overview for the current register scan,
and use Needs attention to resolve machine observations that need a person.

On Samples, the main table and selected-sample card remain open. The bar above
them opens Register samples, Series, Storage grid, Sample experiments, and
More sample details. The last card shows DNA, TFE and any saved other
conditions or notes, and lets you edit preparation date/month and storage;
DNA and TFE remain available as table filters. Preparation date and storage
are display fields; their underlying values stay editable in the register.

Tube status is tracked per physical sample: `Usable` is the neutral working
default (not a positive QC finding), `Suspected degradation` records concern,
`Degraded` records a degraded tube, and `Dead` means permanently unavailable,
including used-up or discarded tubes. Dead samples remain in their recorded
rack slots until the position is cleared or moved; status alone does not prove
physical removal. Older `Retired` or `Unavailable` values are read as `Dead`.
Tube status is independent of experiment assessment.

`Series (optional)` is collapsed by default. Use it only to reuse shared
preparation defaults; it can be ignored for one-off samples or bulk tube entry.

Create a sample for each physical tube using the embedded sample registration
form, then optionally choose its rack position visually (numbered rows 1–12,
lettered columns A–H). Bulk entry previews the exact consecutive positions and
stops if any are occupied; it never skips to a later slot. Existing rack
locations keep their current A–H/1–12 storage values. Use the experiment registration form to register an acquisition
and link it to its sample. Leave a numeric condition blank when it is unknown;
enter `0` only when zero was measured or deliberately set. Blank/unknown values
inherit a selected series default on samples, while an explicit zero overrides
it. Record numeric salt, DNA, and TFE values in their stated units; zero is a
valid value. Do not put a sample temperature on the sample record. Temperature
belongs to the experiment.

Sample IDs are readable labels built from the preparation date, isotope,
protein, distinguishing condition, and sample number. If only the month is
known, enter it as `YYYY-MM`; the register keeps month precision without adding
a day.

The sample form suggests protein names already used for samples, series, or
experiment snapshots, and also accepts new names; saving adds a new spelling to
the sample and series choices without changing existing records.

The collector reads only the acquisition's original `EXPNO/pdata/1/title`.
It never enumerates processing versions or reads their other title files.
Missing original titles and absent sample tags do not create warnings or Needs
attention items. Previously collected comments are retained when the original
title is missing or unreadable; a readable empty title can clear them.
Unreadable or unsafe original titles and explicit conflicting or malformed
sample tags still require attention.

### Link an acquisition title to a tube

Register the physical tube in Grist first. In TopSpin, copy that tube's exact
Sample ID into the acquisition's original `EXPNO/pdata/1/title`, on its
own line with no added spaces or quotes:

```text
SAMPLE_ID=20251017_HeTDR1_15N_NaCl50
```

Use the same ID for acquisitions made from that same tube. Do not reuse it for a
different tube. The marker is explicit text, not a cryptographic signature. A
per-acquisition `SAMPLE_ID` applies only to that acquisition. For sibling
EXPNOs of the same tube in one immediate folder, add this line to the title of
one acquisition in that folder:

```text
SAMPLE_ID_FOLDER=20251017_HeTDR1_15N_NaCl50
```

Copy the tag-holding acquisition with its siblings into a monitored source; the
collector reads only `EXPNO/pdata/1/title`, never other processing versions. The folder default is
limited to that source and its immediate folder and does not descend into
nested folders. An individual `SAMPLE_ID` takes precedence over a folder
default. Keep the title files with the transferred acquisitions.

After a title edit, the next scan observes the change and may fill an unassigned
Sample link only when the marker matches exactly one registered Sample ID and
the acquisition identity is unique and clear. A folder default can fill an
unassigned sibling only when there is one clear folder marker; conflicting
folder defaults are flagged for review. Both markers are optional; an absent marker produces no warning or absence evidence. A missing marker leaves the current
Sample link alone. It does not create a Sample record or replace an existing
human-owned link. Conflicting markers, unreadable titles, duplicate Sample IDs,
or ambiguous acquisition identity appear for review; resolve the source or the
human-owned link in Grist rather than editing observed fields. If the sample is
registered or its ID corrected after the tag was observed, choose the matching
sample manually in the experiment's Attention details. Manually clearing a
Sample link is respected; a later change to the title or acquisition identity
can cause the exact tag match to be checked again. A title change does not
create a new acquisition or copy its assessment.

The reviewed register contains 42 sample records and 757 experiment records.
Keep source relationships unresolved when the workbook does not identify which
physical tube belongs to an experiment; do not infer links from a series label.
When evidence establishes a relationship, set the sample on the experiment and
resolve its editable `MigrationIssue` in the Attention details card.

The collector has 757 exact targets; 753 are currently observed and four are
absent. The current sync has no identity conflicts, metadata warnings, pending
work, or pending results. Sixteen internal result attachments are present.

The Experiments page uses its standard pinned **Experiment type** filter.
**Calibration** and **1D** are unchecked in the saved default view. Select them
in that filter whenever you want to include those EXPNOs.

An EXPNO without raw data is a setup placeholder: changing its pulse program
or parameters, or copying setup templates between EXPNOs, is normal. These
changes update the same record without replacement or copy warnings when the
collector has verified clean raw absence and has never observed raw data there.
Existing reviews or results require reconciliation before a placeholder changes.

The collector distinguishes a clean setup gaining raw data from an acquired
dataset being replaced. Metadata finalization during an acquisition requires
the same start and raw-file generation, with matching audit progress and terminal
events when the pulse program changes. Proven transitions preserve the existing
record and human review; genuine replacements remain flagged for reconciliation.
A raw calibration record whose acquisition and setup pulse programs disagree is
flagged for review rather than silently hidden or reclassified.

Register experiments from the Experiments page. New registrations start as
`Planned`; enter any known acquisition date, instrument, type, pulse program, or
settings directly in their named fields. When the collector finds one exact
acquisition match, it fills actual values into those same fields; until then,
planned values remain. The Assess experiment tray shows collected acquisition
title comments as plain text. The experiment keeps a frozen sample/condition
snapshot from when its sample link was set; later sample edits do not rewrite it.
The raw experiment `Temperature_C` value is retained. Summary and coverage
views show a nominal `Temperature group (°C)`: within 0.2 °C of a whole degree
it groups at that integer; otherwise, within 0.02 °C of a half degree it groups
at that half degree; values outside those bounds remain exact. For example,
readings near 24.85 °C, including records near 298 K, group at
25 °C without changing the original reading; 5.5 °C and 16.0 °C remain distinct
groups.
Review the Needs attention page for unregistered or ambiguous acquisitions,
scan problems, stopped acquisitions, and finished experiments awaiting assessment
for more than a day when a nonblank Study is assigned. Assessments are optional
outside a Study; existing ratings and notes remain visible. Every assigned Study,
including Other or a custom label, enables this reminder. Calibration and Dummy
entries remain exempt; identity, sample, and acquisition problems still require
attention regardless of Study. Missing or ambiguous log evidence shows `Unknown`. A raw
directory being present does not mean acquisition or processing is complete.
The detector accepts reviewed TopSpin 3.6.5, 3.7, 4.4.1, and 5.0.0 audit
formats, including inter-record `$$` comments. Completion still requires a
matching acquisition UUID, start time, and explicit terminal event; copied-row
processing history and an “in progress” entry do not prove completion.
Collector reconciliation notices have a separate local state: clearing a notice
only in Grist does not acknowledge it in the collector, and it may reappear
after the acquisition changes. A supported local acknowledgement workflow is
still pending; retain these notices for manual review until it is available.

## Linux collector

The public commands are:

```text
nmr-register sync [--dry-run]
nmr-register status
nmr-register backup
```

Run `sync --dry-run` after configuring local transferred-data roots to preview
what the next sync would reconcile. A normal hourly sync updates observations;
repeat scans are idempotent. Exact registered acquisitions are matched to their
registration. Ambiguous copies and reruns are surfaced for human review. The
collector prunes protected subtrees before walking them and scans all configured
data types. It is a local user service and emits no operating-system
notifications.

The collector reuses the existing local API credential at
`~/.config/nmr-register/grist.json`. It verifies the file permissions and that
its Grist host and document ID match the configured document. The key stays
local; do not commit it or paste it into chat. `nmr-register credential --test`
checks the existing enrollment with a read-only request. The scheduled service
uses the same credential without a GUI prompt.

From the repository root, install into a private user environment and expose
the command where the supplied systemd unit expects it. Then initialize its
non-secret configuration and verify the existing local credential:

```sh
python3 -m venv ~/.local/share/nmr-register/venv
~/.local/share/nmr-register/venv/bin/pip install '.[migration]'
mkdir -p ~/.local/bin
ln -s ~/.local/share/nmr-register/venv/bin/nmr-register ~/.local/bin/nmr-register
nmr-register config --init --source nmrdata=/absolute/path/to/transferred/data
nmr-register credential --test
nmr-register sync --dry-run
```

The initial configuration uses the approved 11-acquisition pilot allowlist.
Repeat `--source ID=/absolute/path` for every transferred-data root that should
be monitored; keep each source ID stable. Review `nmr-register config` before
starting the timer. After pilot acceptance, `nmr-register config --init
--all-acquisitions` explicitly broadens the scanner. Install the supplied
`units/nmr-register.service` and `units/nmr-register.timer` under
`~/.config/systemd/user/`, then reload and enable the hourly user timer:

```sh
mkdir -p ~/.config/systemd/user
cp units/nmr-register.service units/nmr-register.timer ~/.config/systemd/user/
systemctl --user daemon-reload
systemctl --user enable --now nmr-register.timer
```

It runs while the user service manager is active. Check it with:

```sh
systemctl --user status nmr-register.timer
nmr-register status
```

On the pilot Linux machine, the command is installed at
`~/.local/bin/nmr-register`, the private environment is
`~/.local/share/nmr-register/venv`, and the hourly timer is enabled. The current
status is 753 observations and zero pending work; four configured targets are
currently absent. User lingering is disabled, so collection runs only while the
user session is active. The existing local credential is in use.

## Results and backups

Attach only the existing current result CSV selected by an explicit manifest
and verified checksum. Do not choose the newest CSV by timestamp or upload a
whole experiment folder. Download attachments from Grist to the laptop when
needed.

`nmr-register backup` creates a complete Grist document backup, including
attachments, and retains the seven most recent daily backups. Before a schema
migration, run `nmr-register backup --premigration` and retain that separate
internal Grist backup. Excel export is a read-only data export, not an editable
master or a full document backup. To verify a local backup before restoring,
run `nmr-register backup --validate FILE` and optionally `--restore-to NEWFILE`
to validate a restored local copy. Full document restore uses Grist's internal
restore workflow into a separate document; verify records, links, views, and
attachments there before reconnecting any collector.

The first-time schema/layout setup is an explicit Grist write. Preview it with
`nmr-register setup`; apply it only to the named empty private pilot document
with `nmr-register setup --apply`. Pilot import plans can be reviewed with
`nmr-register import --workbook FILE --reviews FILE --source-id ID` and optional
catalog/status-workbook evidence inputs. Selected pilot rows are applied only
with the explicit `--apply` option. Use
`nmr-register review-snapshot --source-reviews FILE --source-id ID --output
FILE` to export Grist assessment decisions for the legacy review workflow; this
writes a new snapshot file and leaves the source TSV unchanged. The CLI does
not apply a full migration. After the QC schema is installed, the importer
refuses plans that would write a legacy `Assessment` value into the computed
`Assessment` column. Do not import legacy decisions as-is. Any later experiment
import needs row-by-row reviewed QC choices; old status text is never copied
into the formula field. Sample-only imports are separate.

The two reviewed acquisition batches added 745 exact-key records (50 in the
first batch and 695 in the expansion). The guarded importer rechecked manifest
and original source-input hashes, every selected row against a fresh full
inventory, and current Grist identities before inserting. It did not update
existing experiments or infer Sample, Series, Study, QC, or result links.
Observed acquisition facts take precedence; reviewed workbook facts fill only
unknown fields, while observed `Unknown` status remains `Unknown`. Preparation
snapshots stay blank when no exact physical tube is linked; source conditions
remain in `SourceEvidence` for review. A fixed registration time was reused for
preview and retries. After an uncertain response, the importer rereads Grist and
reconciles exact identities instead of blindly repeating a write. Another 414
scan datasets remain held for review; these are not verified independent
acquisitions. The legacy report bridge has not yet been cut over.

Restore into a new or empty Grist document using the Grist internal restore
workflow and verify tables, attachments, and views before reconnecting the
collector. The collector's local state backup alone cannot restore Grist data.
Keep a backup copy away from the host running the collector.

## Rollout boundary

This is a live pilot for the named private document. The Mac browser walkthrough
passed on 2026-11-02. The two reviewed batches added 745 experiments, bringing
the working register to 42 samples and 757 experiments. The document has 16
internal attachments; 11 newly selected CSVs were attached, and the full restore
verified selected CSV bytes against their checksums. Existing experiment rows,
human fields, assessments and prior result pointers were preserved. The package
suite passed 228 tests (four optional skips).

The hourly collector watches 757 exact targets; 753 are currently observed and
four are absent. Another 414 scan datasets remain held for review and are not
verified independent acquisitions. Unknown tube links and legacy QC remain
unresolved. Replacement and copy cases remain excluded for review; do not clear
machine reconciliation fields by hand. The legacy report review-writing cutover
is still incomplete.

## Verification

Run the register acceptance suite from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

Automated tests use temporary synthetic directories and fake transports. See
[ACCEPTANCE.md](ACCEPTANCE.md) for the live checks and remaining rollout limits.

### Electrostatics coverage

Classify each experiment with the human-owned Study field: Electrostatics, DNA,
TFE, Other, or blank when unknown. A blank DNA or TFE concentration does not
establish membership, and the zero-dose DNA control remains a DNA experiment.
On the Electrostatics page, choose a protein, then set and save its field, NaCl,
and experiment-type axes. The matrix has experiment rows and nominal field
columns; each cell lists the expected NaCl concentrations separately. Counts
include only identity-trusted acquisitions with overall `Use`. Blue means one
`Use`; green means two or more in that protein, temperature, field, experiment
and salt cell. Sample or protocol differences remain available in the details
and do not change this colour threshold. Acquisition status remains separate.
`Do not use` acquisitions do not count and appear in the discarded subtitle.
Limited, unassessed, pending, and identity-review records remain visible
separately. All visible salt cells share the tallest required height and
realign when the available width changes. The default types are T1,
T1rho, hetNOE, and etaXY; water T1 is optional. Select one temperature and click a
NaCl line to show its individual acquisitions. A zero count does not establish
that an experiment was never run. Sample and protocol differences remain in
acquisition details.

Analysis status is a human-owned current state. A selected CSV and an Analyzed
status report an output, not a validation or a usability decision. The current
human assessment is shown separately. The acquisition details card contains the
full source identity, assessment reason, and all sample-condition fields.

To review the narrow page installation and its existing-pilot field updates:

```sh
nmr-register electrostatics --catalog /path/to/acquisition-catalog.csv
```

After reviewing the output, use the same command with `--apply`. It writes a
validated backup first, then installs only the Electrostatics page, experiment
form and the reviewed blank pilot fields. If Grist returns HTTP 429, wait until
the site quota or rate limit is restored, then rerun the same command; the
installation can resume. No acquisitions are imported by this command.
Experiment details show instrument, type, date, acquisition status, field,
temperature, pulse program, D1 and delay values, and spin-lock field and offset.
Unknown settings remain blank; entered zeroes remain zero. The collector updates
these same fields only after an exact acquisition match. The title comments are
shown as text in the Assess experiment tray; they are not interpreted as instructions.

Sample-link alerts follow the current Sample assignment immediately. An
unassigned acquisition shows “Sample not assigned”; assigning a tube clears
that alert. Redundant imported “Sample link unresolved” clauses and their
adjacent workbook conditions remain stored as migration evidence but do not
clutter attention messages. Genuine sample conflicts, title mismatches, and
other acquisition or review warnings remain visible.

### Assessment checklist

In Experiments, select a row in the main list and open Assess experiment in the
tray beside Register experiment. The tray edits that same experiment. Enter the
six ratings and optional note, save, then close the tray. Overall `Assessment`
remains visible in the main table; the tray shows the selected acquisition's
source and title comments. Electrostatics also has an Assess experiment tray:
select a map cell, choose a row in Individual acquisitions, then open the tray.
Both pages edit the same experiment and calculated assessment.
The six choices are human-entered; the collector,
CSV presence, and `AnalysisStatus` do not automatically mark a criterion as good.

For a discovered acquisition with a verified sample and a unique, consistent
folder/EXPNO match, saving a changed checklist rating also confirms registration
at the save time. Existing registration dates are preserved. Note-only edits
and unresolved acquisition or sample identities do not confirm registration.

| Criterion | Choices, from acceptable to failing |
| --- | --- |
| Acquisition coverage | Complete; Sufficient; Insufficient |
| Settings | Appropriate; Usable with restrictions; Invalid |
| FIDs and spectrum | Clean; Usable with exclusions; Unusable |
| SNR | Adequate; Marginal; Inadequate |
| Sample during acquisition | Acceptable; Concern; Unacceptable |
| Analysis quality | Reliable; Limited; Unreliable |

The calculated `Assessment` is read-only. Any failing choice yields `Do not
use`, even if another criterion is still unassessed. With no failure, any
`Not assessed` choice leaves the canonical assessment blank for counts and exports. Once all criteria
are assessed, zero or one caution-level choice yields `Use`; two or more
yield `Limited use`. The explanation identifies an accepted single caution. Acquisition coverage `Complete` and `Sufficient` are both
acceptable. This is a decision checklist, not an average or score.

The main table's `Assessment` display uses the computed `AssessmentDisplay`
column. An incomplete checklist shows progress such as `Use so far · analysis
missing` or `Limited use so far`. These are provisional;
they do not count as approved `Use` decisions. Only a sole missing analysis
check adds ` · analysis missing`; other missing checks stay in the detailed
explanation. A failing choice remains decisive as `Do not use`. With no rated checks or imported
decision, the display says `Not assessed`. Existing imported decisions remain
visible when all six checklist choices are blank. Completing the checklist
removes the provisional wording.

`AssessmentExplanation` is read-only and explains failures, limited choices, or
remaining unassessed criteria together, including known cautions and failures on partial checklists. That calculated explanation satisfies the reason
requirement for `Limited use` or `Do not use`; add optional context in the shared
manual `AssessmentReason` field when useful. The migration leaves checklist
items `Not assessed` unless existing evidence supports a reviewed mapping. It
preserves existing notes and reasons and does not infer all-pass ratings.

`AnalysisStatus` separately tracks processing progress (`Not started`, `In
progress`, or `Analyzed`). It is not a quality rating and does not change the
checklist or calculated assessment.


Optional instrument frequency profiles in collector `config.json` identify an
instrument from its raw observed proton frequency (MHz). Universal defaults are
empty. For this lab, the confirmed profiles can be configured as:

```json
"instrument_profiles": [
  {"id": "600_1300", "proton_mhz": 599.8, "nominal_mhz": 600},
  {"id": "600Cryo", "proton_mhz": 600.48, "nominal_mhz": 600},
  {"id": "700", "proton_mhz": 700.25, "nominal_mhz": 700},
  {"id": "850", "proton_mhz": 850.13, "nominal_mhz": 850},
  {"id": "950", "proton_mhz": 950.23, "nominal_mhz": 950}
]
```

Each profile has a unique nonempty `id`, positive finite `proton_mhz` and
`nominal_mhz`, and optional positive finite `tolerance_mhz` (default `0.05`).
`nominal_mhz` describes the nominal field; matching always uses `proton_mhz`
against the unrounded `ObservedField_MHz`. Exactly one profile must match within
its inclusive tolerance. Unknown or overlapping matches leave the instrument
unresolved; overlapping profiles are allowed so ambiguity never selects an
arbitrary instrument. A configured source `instrument` takes precedence.
Inferred matches append their frequency and tolerance to
`ObservedProtocolEvidence`, preserving existing evidence and the raw field.

Confirmed TopSpin split datasets and explicit single-FID extractions are excluded from the acquisition inventory before folder sample defaults are resolved. Exclusion uses the explicit `created by split N: ... from EXPNO N` audit statement or the collector's existing machine evidence, never experiment numbers, folder names, pulse programs, titles, or free-form notes. A single-FID extraction requires a complete structured `mvrow`/`rser ... from` audit origin; any recorded `goN` acquisition prevents that exclusion. Experiment numbers, file sizes, Dummy status, comments, and ordinary mentions do not establish extraction. Original acquisitions and acquisitions with unknown metadata remain eligible.

On observation or flush, local state removes confirmed split observations and pending experiment writes, including caches from older collectors. It removes only their exact targets from generated matching-metadata warnings; other copy targets and replacement or ambiguous warnings remain. This does not delete source files or existing Grist rows. Existing split rows require a separate reviewed cleanup. Batch migration also rejects a split retained in a stale observation inventory.
