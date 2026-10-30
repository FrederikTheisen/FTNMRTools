# NMR register pilot

The pilot keeps sample and experiment metadata in one private Grist document and
uses a Linux hourly scan to reconcile transferred Bruker data. Researchers edit
samples and assessments in Grist from either Mac or Linux. The collector records
what it observed in separate `Observed…` fields; it does not take ownership of
human-entered metadata, assessment, notes, or sample links.

## Daily use

The Grist document has four working views: Experiments with a card for one
record; Samples with series and storage views (separate storage boxes, each
with rows A–H and columns 1–12);
Needs attention; and Overview. Use the Overview for the current register scan,
and use Needs attention to resolve machine observations that need a person.

Create a sample for each physical tube using the embedded sample registration
form, then assign its storage box and position (each box has rows A–H and
columns 1–12). Use the experiment registration form to register an acquisition
and link it to its sample. Leave a numeric condition blank when it is unknown;
enter `0` only when zero was measured or deliberately set. Blank/unknown values
inherit a series default on samples, while an explicit zero overrides it. Put
shared preparation conditions on the series and
override them for an individual tube only when needed. Record numeric salt, DNA,
and TFE values in their stated units; zero is a valid value. Do not put a sample
temperature on the sample record. Temperature belongs to the experiment.

Sample IDs are readable labels built from the preparation date, isotope,
protein, distinguishing condition, and sample number. If only the month is
known, enter it as `YYYY-MM`; the register keeps month precision without adding
a day.

The reviewed pilot import contains seven sample records and eleven experiment
records. Four DNA titration samples remain unlinked because the source workbook
does not map those tubes to particular experiments. Their source evidence stays
with the records for review; do not infer an experiment link from a titration
series label. When source evidence identifies a relationship, set the physical
sample on the experiment and resolve its editable `MigrationIssue` in the
Attention details card.

The current pilot sync matched eight acquisitions to collector observations
and found no identity conflicts. Four experiment rows still need a sample link,
and three metadata differences need review; those differences block legacy
review export. Two current result attachments are present and match their
recorded checksums.

Register experiments from the Experiments page. Choose the sample and set the
current assessment to `Use`, `Limited use`, or `Do not use`. Give a reason for
the latter two. The experiment retains a frozen sample/condition snapshot from
when its sample link was set; later sample edits do not rewrite that snapshot.
Review the Needs attention page for unregistered or ambiguous acquisitions,
scan problems, and experiments awaiting assessment for more than a day. A raw
directory being present does not mean acquisition or processing is complete.
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
`~/.local/share/nmr-register/venv`, and the hourly timer is already enabled.
The latest service run succeeded with exit status 0; status reports eight
observations, no pending work, and no last error. User lingering is disabled,
so collection runs only while the user session is active. The existing local
credential is in use.

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
not apply a full migration.

Restore into a new or empty Grist document using the Grist internal restore
workflow and verify tables, attachments, and views before reconnecting the
collector. The collector's local state backup alone cannot restore Grist data.
Keep a backup copy away from the host running the collector.

## Rollout boundary

This is a pilot for the named private document. Local configuration, sync,
exports, and a full backup restore to a temporary Grist document have been
verified. Complete the Mac browser walkthrough before any full legacy
migration. Acquisition replacement/copy cases remain flagged until the local
collector acknowledgement workflow is implemented and reviewed; do not clear
machine-generated reconciliation fields by hand. The current migration tool is
limited to explicitly selected pilot rows; it does not enable a full migration.

## Verification

Run the register acceptance suite from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

Automated tests use temporary synthetic directories and fake transports. See
[ACCEPTANCE.md](ACCEPTANCE.md) for the live pilot checks and the remaining Mac
browser walkthrough.
