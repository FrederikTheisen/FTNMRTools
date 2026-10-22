#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from collections import defaultdict

def read_peak_file(path):
    """Return unchanged table bytes and optional metadata from the comment footer."""
    data = Path(path).read_bytes()
    begin = b'# CCPN_METADATA_BEGIN\n'
    if begin not in data:
        if b'CCPN_METADATA_BEGIN' in data or b'CCPN_METADATA_END' in data:
            raise ValueError(f'{path}: malformed metadata footer')
        return data, None
    if data.count(begin) != 1:
        raise ValueError(f'{path}: multiple metadata blocks')
    table, footer = data.split(begin, 1)
    lines = footer.decode('utf-8').splitlines()
    if not table.endswith(b'\n') or not lines or lines[-1] != '# CCPN_METADATA_END':
        raise ValueError(f'{path}: incomplete metadata footer')
    if not all(line.startswith('# ') for line in lines[:-1]):
        raise ValueError(f'{path}: invalid metadata comment line')
    metadata = json.loads('\n'.join(line[2:] for line in lines[:-1]))
    return table, metadata


def read_embedded_noise_sigma(path):
    table, metadata = read_peak_file(path)
    if metadata is None:
        return None
    if metadata.get('schema') != 'ccpn_peak_export_metadata' or metadata.get('schema_version') != 2:
        raise ValueError(f'{path}: unsupported embedded metadata schema')
    if hashlib.sha256(table).hexdigest() != metadata.get('peak_table_sha256'):
        raise ValueError(f'{path}: peak-table hash mismatch; re-export before using noise metadata')
    sigma = float(metadata['noise']['noise_sigma'])
    if metadata['noise'].get('status') != 'estimated' or not math.isfinite(sigma) or sigma <= 0:
        raise ValueError(f'{path}: invalid embedded noise SD')
    return sigma

def to_bool(s):
    if s is None:
        return False
    return s.strip() == "True"

int_re = re.compile(r'(\d+)')

def residue_from_assignment(assign):
    m = int_re.search(assign)
    if not m:
        raise ValueError(f"Could not extract residue number from assignment '{assign}'")
    return int(m.group(1))

def parse_sparky_list(path):
    """
    Expect columns:
      assignment  w1  w2  intensity  overlap_bool  overlap_info...
    Lines starting with '#' are ignored. Whitespace-delimited.
    Returns dict: residue -> dict(fields)
    In case of duplicates per residue, keep the row with largest |intensity|.
    """
    data = {}
    with open(path, 'r', encoding='utf-8') as f:
        print(f"# Parsing Sparky list: {path}")
        idx_int = 3
        idx_overlap = 4
        for ln in f:
            line = ln.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if parts[0].lower() == 'assignment':
                if 'intensity' in parts:
                    idx_int = parts.index('intensity')
                if 'overlap' in parts:
                    idx_overlap = parts.index('overlap')
                continue
            if len(parts) < 4:
                # skip malformed lines
                continue
            assignment = parts[0]
            try:
                w1 = float(parts[1])
                w2 = float(parts[2])
                intensity = float(parts[idx_int])
            except ValueError:
                print(f'LINE ERROR: {line}')
                continue
            overlap_bool = False
            overlap_info = ""
            if len(parts) > idx_overlap:
                overlap_bool = to_bool(parts[idx_overlap])
            if len(parts) > idx_overlap:
                overlap_info = " ".join(parts[idx_overlap:])

            try:
                resid = residue_from_assignment(assignment)
            except ValueError:
                # skip entries that don't encode a residue number
                continue

            # Keep the entry with largest absolute intensity if duplicates
            prev = data.get(resid)
            if (prev is None) or (abs(intensity) > abs(prev["intensity"])):
                data[resid] = {
                    "assignment": assignment,
                    "w1": w1,
                    "w2": w2,
                    "intensity": intensity,
                    "overlap": bool(overlap_bool),
                    "overlap_info": overlap_info
                }
    return data

def build_output_prefix_from_ref(ref_file):
    """
    From ref filename like '500mM_25C_600MHz_HetNOE1.1.list',
    return '500mM_25C_600MHz_HetNOE'.
    Strategy: take basename, strip directory, then cut at the first occurrence of '_HetNOE'
              (case-insensitive), remove trailing underscores, then append '_HetNOE'.
    If pattern not found, use stem + '_HetNOE'.
    """
    name = Path(ref_file).name
    m = re.search(r'(_+)HetNOE\d*', name, flags=re.IGNORECASE)
    if m:
        prefix = name[:m.start()].rstrip('_')
    else:
        prefix = Path(ref_file).stem
    return f"{prefix}_HetNOE"

def compute_ratio_and_error(A_noe, B_ref, sigma, allow_zero_ref=False, sigma_ref=None):
    """
    R = A/B
    sigma_R = sqrt((sigma_NOE/B)^2 + (A*sigma_REF/B^2)^2).
    sigma is the NOE SD; sigma_ref defaults to the same value for legacy calls.
    If B ~ 0, return (nan, nan) unless allow_zero_ref=True (then tiny epsilon guard).
    """
    if B_ref == 0.0:
        if allow_zero_ref:
            B_ref = 1e-12
        else:
            return (float('nan'), float('nan'))
    R = A_noe / B_ref
    ref_sd = sigma if sigma_ref is None else sigma_ref
    sigma_R = math.hypot(sigma / B_ref, R * ref_sd / B_ref)
    return (R, sigma_R)

def write_out(path, rows, simple=True):
    """
    rows: list of dicts with keys:
      resid, I_ref, I_noe, ratio, ratio_sigma, overlap_ref, overlap_noe
    """
    with open(path, 'w', encoding='utf-8') as out:
        if not simple: out.write("resid I_ref I_noe ratio error overlap_ref overlap_noe\n")
        for r in sorted(rows, key=lambda x: x["resid"]):
            if not simple: out.write(f'{r["resid"]:d} {r["I_ref"]:.5f} {r["I_noe"]: .6g} {r["ratio"]: .6g} {r["error"]: .6g} {int(r["overlap_ref"])} {int(r["overlap_noe"])}\n')
            else: out.write(f'{r["resid"]:d} {r["ratio"]:.5f} {r["error"]:.5f}\n')


def main():
    ap = argparse.ArgumentParser(description="Compute HetNOE intensity ratios (NOE/REF) from two Sparky peak lists.")
    ap.add_argument("-ref", help="Reference (noe OFF) Sparky .list file (used for output naming).", required=True)
    ap.add_argument("-noe", help="NOE (noe ON) Sparky .list file.", required=True)
    ap.add_argument("-noise", type=float, required=False,
                    help="Override noise SD for both spectra; otherwise read embedded noise metadata when present.", default=None)
    ap.add_argument("--allow-zero-ref", action="store_true",
                    help="If REF intensity is zero, use a tiny epsilon instead of returning NaN.")
    args = ap.parse_args()

    ref = parse_sparky_list(args.ref)
    noe = parse_sparky_list(args.noe)

    if args.noise is not None:
        if not math.isfinite(args.noise) or args.noise < 0:
            ap.error('-noise must be finite and nonnegative')
        sigma_ref = sigma_noe = args.noise
    else:
        sigma_ref = read_embedded_noise_sigma(args.ref)
        sigma_noe = read_embedded_noise_sigma(args.noe)
        if (sigma_ref is None) != (sigma_noe is None):
            ap.error('Both lists need embedded noise metadata, or supply -noise explicitly')
        if sigma_ref is None:
            sigma_ref = sigma_noe = 0.0  # Preserve legacy plain-list behaviour.
    print(f"# Noise SD: REF={sigma_ref:g}; NOE={sigma_noe:g} (independent noise-only errors)")

    # Intersect residues (they're expected to match; we enforce intersection)
    common = sorted(set(ref.keys()) & set(noe.keys()))
    missing_ref = sorted(set(noe.keys()) - set(ref.keys()))
    missing_noe = sorted(set(ref.keys()) - set(noe.keys()))
    if missing_ref or missing_noe:
        print(f"# Warning: Residue mismatch. In NOE but not REF: {missing_ref}; in REF but not NOE: {missing_noe}")

    rows_all = []
    rows_filtered = []
    for resid in common:
        I_ref = ref[resid]["intensity"]
        I_noe = noe[resid]["intensity"]

        ratio, sigma_ratio = compute_ratio_and_error(I_noe, I_ref, sigma_noe, allow_zero_ref=args.allow_zero_ref, sigma_ref=sigma_ref)

        overlap_note = "overlap" if ref[resid]["overlap"] or noe[resid]["overlap"] else ""
        print(f"""{resid}\t{ratio:+.2f} ± {sigma_ratio:.2f} {overlap_note}""")

        row = {
            "resid": resid,
            "I_ref": I_ref,
            "I_noe": I_noe,
            "ratio": ratio,
            "error": sigma_ratio,
            "overlap_ref": ref[resid]["overlap"],
            "overlap_noe": noe[resid]["overlap"]
        }
        rows_all.append(row)
        if (not ref[resid]["overlap"]) and (not noe[resid]["overlap"]):
            rows_filtered.append(row)

    out_prefix = build_output_prefix_from_ref(args.ref)
    out_all = Path(f"{out_prefix}_results.out")
    out_filtered = Path(f"{out_prefix}_results_filtered.out")

    write_out(out_all, rows_all)
    write_out(out_filtered, rows_filtered)

    print(f"# Wrote: {out_all}")
    print(f"# Wrote: {out_filtered}")
    print(f"# Residues total (common): {len(common)}  |  filtered (non-overlap in both): {len(rows_filtered)}")

if __name__ == "__main__":
    main()
