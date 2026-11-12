# Dummy experiment entries

Set Acquisition status to **Dummy** for an intentional placeholder. The registration dialog offers Planned (the default) and Dummy entry. Observed acquisition outcomes remain collector managed. Dummy entries have no attention reason and are excluded from electrostatics matrix records, identity indexes, observed axes, and protein counts. Collector refreshes retain Dummy while continuing to synchronize other machine metadata. No migration marks existing records Dummy.

**Desktop data** reports **Not linked** when a normal entry lacks either AcquisitionKey or SourcePath. Typing Folder alone does not establish a link. A linked folder can contain acquisition setup without raw FID data and does not imply a completed experiment.

`experiments_ui.upgrade_dummy_schema(client, apply=False)` previews the guarded schema upgrade; `apply=True` applies it. Only exact previous shipped attention and effective status formulas are accepted. Existing choices and widget options remain intact. `editable_experiment_status` changes existing grid field references in place, preserving field IDs, widths, positions, and options. `apply_experiment_styles` appends native grey Dummy and pale yellow attention row rules, plus subtle blue text for Not linked cells, retaining custom rules. These two helpers can be run independently of the full layout setup workflow.

## Calibration experiments

Exact normalized pulse names zg, zg30, 15N_calib and 13C_calib are classified as Experiment type Calibration. Acquisition status continues to reflect acquisition progress and outcome. Calibration skips only the routine completed-over-24-hours unassessed reminder. Missing registration, missing samples, reconciliation issues, scan errors and other warnings remain active. Calibration has no special row colour.
