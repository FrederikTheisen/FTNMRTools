#!/usr/bin/env python3
"""
Consolidate 3-column files (residue, rate, error) into a single table.

Usage:
  python relax_consolidate.py <directory> <pattern>

Output goes to stdout with columns:
  residue f1_rate f1_error f2_rate f2_error ...

Notes:
- Lines starting with '#' and blank lines are ignored.
- Residues are sorted numerically when possible, otherwise lexicographically.
- If a file contains duplicate residues, the first occurrence is kept.
- Requires pandas (and its dependency numpy).
- Only files matching <pattern>*results_filtered.out in <directory> are included.
"""
import sys
import os
import glob
import pandas as pd

def read_one(path, idx):
    # Extract base name for column naming
    fname = os.path.basename(path)
    # Remove 'results_filtered.out' from end
    colbase = fname.replace('results_filtered.out', '').rstrip('_')
    try:
        df = pd.read_csv(
            path,
            delim_whitespace=True,
            comment="#",
            header=None,
            usecols=[0, 1, 2],
            names=["residue", f"{colbase}_rate", f"{colbase}_error"],
            dtype={"residue": str}  # keep as string; we sort numerically later if possible
        )
    except FileNotFoundError:
        print(f"[error] file not found: {path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[error] reading {path}: {e}", file=sys.stderr)
        sys.exit(1)

    # Drop duplicate residues within the same file (keep first)
    before = len(df)
    df = df.drop_duplicates(subset=["residue"], keep="first")
    if len(df) < before:
        print(f"[warn] {path}: dropped {before - len(df)} duplicate residue rows", file=sys.stderr)
    return df

def weighted_mean_and_error(rates, errors):
    import numpy as np
    rates = np.array(rates, dtype=float)
    errors = np.array(errors, dtype=float)
    # Only use valid (non-nan, nonzero error) entries
    mask = (~np.isnan(rates)) & (~np.isnan(errors)) & (errors > 0)
    rates = rates[mask]
    errors = errors[mask]
    if len(rates) == 0:
        return (np.nan, np.nan)
    if len(rates) == 1:
        return (rates[0], errors[0])
    weights = 1 / errors**2
    mean = np.sum(weights * rates) / np.sum(weights)
    sem = np.sqrt(1 / np.sum(weights))
    # If repeats differ a lot, inflate error by SD
    sd = np.std(rates, ddof=1)
    error = max(sem, sd)
    return (mean, error)

def main(argv):
    if len(argv) != 3:
        print("Usage: python relax_consolidate.py <directory> <pattern>", file=sys.stderr)
        sys.exit(2)
    directory = argv[1]
    pattern = argv[2]
    search_pattern = os.path.join(directory, f"{pattern}*results_filtered.out")
    files = sorted(glob.glob(search_pattern))
    if not files:
        print(f"No files found matching pattern: {search_pattern}", file=sys.stderr)
        sys.exit(1)
    merged = None
    for i, path in enumerate(files, start=1):
        df = read_one(path, i)
        merged = df if merged is None else merged.merge(df, on="residue", how="outer")

    # Sort residues numerically when possible, else by string
    key = pd.to_numeric(merged["residue"], errors="coerce")
    merged = merged.assign(__key=key).sort_values(by=["__key", "residue"], ascending=[True, True]).drop(columns="__key")

    # Calculate consolidated weighted mean and error for each residue
    rate_cols = [c for c in merged.columns if c.endswith('_rate')]
    error_cols = [c for c in merged.columns if c.endswith('_error')]
    consolidated = []
    for idx, row in merged.iterrows():
        rates = [row[c] for c in rate_cols if not pd.isna(row[c])]
        errors = [row[c.replace('_rate', '_error')] for c in rate_cols if not pd.isna(row[c])]
        mean, err = weighted_mean_and_error(rates, errors)
        consolidated.append((row['residue'], mean, err))
    # Build output DataFrame
    outdf = pd.DataFrame(consolidated, columns=['residue', 'consolidated_rate', 'consolidated_error'])
    # Optionally, merge with original columns for reference
    outdf = outdf.merge(merged, on='residue', how='left')
    outdf.to_csv(sys.stdout, sep="\t", index=False, na_rep="")

if __name__ == "__main__":
    main(sys.argv)
