# Acceptance record

The automated suite checks safety properties using synthetic files and fake
Grist transports. The independent acceptance tests use temporary data and an
injected transport; the full suite also exercises a local HTTP server.

On 2026-10-30, all 74 package tests passed in the installed Python 3.12
environment, including the localhost HTTP test, and passed again after the
central integration was installed. The 12 independent acceptance tests also
passed with `/opt/ccpnmr/miniconda/bin/python`. The existing
FTNMRTools relaxation suite was run in an isolated copy before and after
replacing only `relax/condition_overview.py` with the staged integration. Both
runs had the same 23 failure/error names and order (188 tests, 5 failures,
18 errors, 1 skipped); no new failures appeared.

Live pilot checks also passed on 2026-10-30. The private document has one
series, seven samples, and eleven experiments. Eight acquisitions are linked
to collector observations, seven experiment records have a physical sample
link, and four experiment records still need sample links. Four known-DNA
titration samples remain unlinked because the source workbook does not identify
their experiment relationships. No identity conflicts were reported. Three
metadata differences remain visible for review; these block legacy review
export until resolved. Two current result attachments downloaded and matched
their recorded SHA-256 values. Two consecutive syncs after the first produced
no creates, updates, or attachments, and left no pending work. Human fields
survived the initial sync unchanged. The read-only field audit found one
sample-level DNA value of zero and three blank/unknown values; those remain
distinct in the register. The hourly user timer is installed and enabled on the
pilot Linux machine. Its service completed successfully (exit status 0); local
status showed eight observations, zero pending work, and no last error. User
lingering is disabled, so the timer runs only during the active user session.

The read-only document audit found all four working pages, both full-access
registration forms, the expected SampleID/SeriesID reference labels, the
storage grid (rows A–H by columns 1–12), and no formula error cells. The XLSX
and CSV exports matched the live register counts. A live internal backup was
restored to a temporary document in the same workspace; all seven tables,
record fields, references, view/filter metadata, and both attachment hashes
matched. The temporary document was deleted after those checks. The Mac browser
walkthrough remains pending.

| Requirement | Local acceptance evidence | Still requires |
|---|---|---|
| Protected-path pruning | Synthetic protected directory and a sentinel child whose metadata access fails if traversed | Review configured roots on Linux |
| Repeated imports and scans | Synthetic retry/idempotence tests; two live follow-up syncs created no records or attachments | Review replacement/copy reconciliation before broad rollout |
| Human-field preservation | Synthetic sync tests and live sync preserved all seeded human fields | Mac browser review |
| Timeout after remote commit | Fake HTTP transport commits then times out; retry reconciles without duplication | No induced remote timeout; normal sync completed and retried idempotently |
| Zero-valued conditions and frozen snapshot | Synthetic trigger-contract test; live zero remains distinct from blank, and no formula errors were found | Mac browser review of editing and recalculation |
| Partial dates | Synthetic month-precision test; live table audit found no formula errors | Mac browser review |
| Pilot import conflicts | Synthetic conflict tests; live sync reported zero identity conflicts and four missing experiment Sample links | Review the four links and three metadata differences |
| Result identity | Synthetic manifest/hash tests; both live result attachments downloaded and matched their SHA-256 metadata | Confirm result meaning on the laptop |
| Backup restore | Synthetic SQLite/attachment test and live restore copy matched seven tables, records, layout metadata, refs, and attachment hashes | Mac browser review |
| Export parity | Live XLSX and CSV exports contained the same sample, experiment, and series counts as Grist | Mac browser review |

The independent acceptance suite was run with:

```sh
PYTHONDONTWRITEBYTECODE=1 /opt/ccpnmr/miniconda/bin/python -m unittest discover -s tests -p test_acceptance.py -v
```

All 12 tests passed. The full package suite was run with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

The existing FTNMRTools relaxation suite is independent and can be run from
that repository with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s relax/tests -v
```

As of the 2026-10-30 acceptance record, full migration remained disabled until
the user completed the Mac browser walkthrough. Acquisition replacement/copy cases remain flagged for review until
the collector acknowledgement workflow is implemented and reviewed; do not
clear machine-generated reconciliation fields by hand.

## Current update — 2026-11-02

The Mac browser walkthrough passed. Two reviewed batches added 745 exact-key
experiments, bringing the register to 42 Samples, 757 Experiments and 16 internal
attachments. All 62 experiments present before the expansion retained their
stored values, human fields, assessments and prior result pointers. Eleven newly
selected CSVs were attached; the full restore verified 14 selected CSV byte
hashes. No sample links, studies, or QC ratings were inferred.

The collector has 757 exact targets, with 753 currently observed and four
absent. The current sync created no records and reported zero conflicts,
metadata warnings, pending work, or pending results. A repeated import plan
proposed zero inserts. Another 414 scan datasets remain held for review; these
are not verified independent acquisitions. The complete restore and package
suite passed (228 tests, four optional skips), and the hourly timer was restored.

The historical notes above describe the 2026-10-30 state. Broad migration,
legacy-QC reconciliation and the report review-writing cutover remain incomplete;
replacement/copy candidates stay excluded pending review.

## Current update — 2026-11-04

The current full desktop scan found 1,081 acquisitions, all linked in Grist.
The 79 reviewed workbook/acquisition merges preserved every nonformula field
according to the conflict and ownership contract. All prior nonblank effective
assessments survived, sample links and notes transferred without replacing
conflicting human values, and selected result pointers stayed on their linked
acquisitions. The current document has 1,132 experiment records and 54 internal
attachments. A full backup restored with the experiment records, layout and all
attachment bytes unchanged; original workbook bytes and collector state were
unchanged. The hourly collector timer was restored.

There are 51 visible, unlinked workbook rows: 35 missing EXPNOs,
four effective-assessment conflicts, seven identities needing
correction, two setup-only matches, and three absent acquisitions. Each has a
specific migration action in Needs attention. These are reconciliation cases;
there is no undisclosed hold on the 1,081 desktop acquisitions. The earlier
historical counts and rollout holds above do not describe the current state.

The reconciliation regression tests cover real zero values, unset references,
source-fact ownership, assessment conflicts and masked imported decisions,
selected-result conflicts, reference rewrites, fresh-edit guards, atomic retry
recognition, and removal of obsolete generated migration warnings.

The user-confirmed follow-up reconciled eight 2026-04-22 HeTDR1 entries as
30 bp DNA and five 2026-03-02 SeV Ntail entries as 250 mM NaCl. The original
workbook explicitly names the existing 30 bp tube for the eight DNA acquisitions;
those now link to that sample with frozen preparation conditions. The salt
acquisitions retain the existing 250 mM tube and frozen snapshots. The separate
50 mM sample record was unchanged. All current assessments, complementary noise
notes, and the existing shared-fingerprint warning were preserved. A fresh full
backup retained all 54 attachment byte hashes. The original workbook was unchanged.
