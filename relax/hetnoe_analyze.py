#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path
from collections import defaultdict

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
            if 'intensity' in line:
                idx_int = line.split().index('intensity')
            if 'overlap' in line:
                idx_overlap = line.split().index('overlap')
            if not line or line.startswith('#'):
                continue
            parts = line.split()
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

def compute_ratio_and_error(A_noe, B_ref, sigma, allow_zero_ref=False):
    """
    R = A/B
    sigma_R = sqrt((sigma/B)^2 + (A*sigma/B^2)^2)
    If B ~ 0, return (nan, nan) unless allow_zero_ref=True (then tiny epsilon guard).
    """
    if B_ref == 0.0:
        if allow_zero_ref:
            B_ref = 1e-12
        else:
            return (float('nan'), float('nan'))
    R = A_noe / B_ref
    sigma_R = abs(sigma / B_ref) * math.sqrt(1.0 + R*R)
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
                    help="Noise SD to use for both spectra (same sigma for REF and NOE).", default=0.0)
    ap.add_argument("--allow-zero-ref", action="store_true",
                    help="If REF intensity is zero, use a tiny epsilon instead of returning NaN.")
    args = ap.parse_args()

    ref = parse_sparky_list(args.ref)
    noe = parse_sparky_list(args.noe)

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

        ratio, sigma_ratio = compute_ratio_and_error(I_noe, I_ref, args.noise, allow_zero_ref=args.allow_zero_ref)

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
