"""CSPanalysis — compute CSPs and fit a simple 1:1 binding model.

Minimal, readable utilities to parse Sparky HSQC peak lists (handles
pseudo-dimensions and nucleus-labelled w# columns), extract concentrations
from filenames (converted to µM), compute CSPs with nucleus weighting, and
fit binding curves. See the CLI at the bottom of the file for usage.
"""

import argparse
import glob
import os
import re
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares

# Maximum allowed CSP value when fitting – prevents runaway fits
FITTING_CSP_MAX = 10.0

# Weighting factors for nitrogen and hydrogen chemical shifts (1H and 15N)
SHIFT_WEIGHTS = {'H': 1.0, 'N': 0.154}


def extract_variable_from_filename(filename: str) -> float:
    pattern = re.compile(r"(?i)(\d+(?:\.\d+)?)(mM|uM|nM|M|µM)")
    for num_str, unit in pattern.findall(filename):
        try:
            value = float(num_str)
        except ValueError:
            continue
        unit = unit.lower().replace('µ', 'u')
        if unit == 'mm':
            return value * 1e3
        elif unit == 'um':
            return value
        elif unit == 'nm':
            return value / 1e3
        elif unit == 'm':
            return value * 1e6
    
    print("No concentration found in filename:", filename)
    quit()


def extract_residue_number(assignment: str) -> int:
    match = re.search(r"(\d+)", str(assignment))
    if match:
        return int(match.group(1))
    
    quit()
    return 10 ** 9  # large number to send unknowns to the end


def read_peak_list(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep='\t')
    if df.empty:
        return df
    # Inspect the assignment of the first row to determine how to rename w#
    assignments = [tok.strip() for tok in str(df['assignment'].iloc[0]).split('-')]
    print(assignments)
    # Build a mapping from original w# columns to new names; skip pseudo dimensions.
    # To ensure subsequent merging works, collapse the dimension index and only
    # encode the nucleus (H or N).  For HSQC peak lists there should be at most
    # one H and one N dimension.
    col_mapping: Dict[str, str] = {}
    columns_to_drop: List[str] = []
    used_nuclei: Dict[str, int] = {}
    for idx, assign in enumerate(assignments):
        w_col = f"w{idx + 1}"
        if assign == '?' or assign == '':
            if w_col in df.columns:
                columns_to_drop.append(w_col)
            continue
        atom_type = assign[-1]  # 'H' or 'N'
        if w_col in df.columns:
            count = used_nuclei.get(atom_type, 0) + 1
            used_nuclei[atom_type] = count
            # If multiple tokens of the same atom type appear, append an index
            suffix = '' if count == 1 else str(count)
            col_mapping[w_col] = f"shift_{atom_type}{suffix}"
    df = df.rename(columns=col_mapping)
    if columns_to_drop:
        df = df.drop(columns=columns_to_drop)
    df['residue'] = df['assignment'].apply(extract_residue_number)
    df = df.sort_values('residue').reset_index(drop=True)
    return df


def compute_intensity_ratio(df_ref: pd.DataFrame, df_target: pd.DataFrame) -> pd.DataFrame:
    merged = pd.merge(df_ref[['assignment', 'intensity']],
                     df_target[['assignment', 'intensity']],
                     on='assignment', suffixes=('_ref', '_target'))
    result = merged[['assignment']].copy()
    result['int_ratio'] = merged['intensity_target'] / merged['intensity_ref']
    return result


def compute_csp(ref, target, minimize: bool = False, ref_file: str = None, target_file: str = None) -> pd.DataFrame:
    print(f"Computing CSP between peak lists:")
    print(f"  Reference: {ref['file']}")
    print(f"  Target:    {target['file']}")
    merged = pd.merge(ref['df'], target['df'], on='assignment', suffixes=('_ref', '_target'))
    # Identify nitrogen and proton columns for reference and target separately.
    n_col_ref = None
    h_col_ref = None
    for col in merged.columns:
        if col.endswith('_ref'):
            base = col[:-4]  # strip '_ref'
            # recognise names like shift_N, shift_N2, shift_H, shift_H1, etc.
            if base.startswith('shift_N'):
                n_col_ref = col
            elif base.startswith('shift_H'):
                h_col_ref = col
    # Identify target columns: columns with suffix '_target'
    n_col_target = None
    h_col_target = None
    for col in merged.columns:
        if col.endswith('_target'):
            base = col[:-7]  # strip '_target'
            if base.startswith('shift_N'):
                n_col_target = col
            elif base.startswith('shift_H'):
                h_col_target = col
    if None in (n_col_ref, h_col_ref, n_col_target, h_col_target):
        raise RuntimeError(f"Could not identify N/H chemical shift columns in merged data: {list(merged.columns)}")
    dn = merged[n_col_ref] - merged[n_col_target]
    dh = merged[h_col_ref] - merged[h_col_target]
    if minimize:
        drift_n = np.median(dn)
        drift_h = np.median(dh)
        dn = dn - drift_n
        dh = dh - drift_h
        print(f"\t-drift correction: N={drift_n}, H={drift_h}")
    total_csp = np.sqrt((SHIFT_WEIGHTS['N'] * dn) ** 2 + (SHIFT_WEIGHTS['H'] * dh) ** 2)
    result = merged[['assignment']].copy()
    result['csp'] = total_csp
    return result


def binding_model(L: np.ndarray, Delta_max: float, Kd: float, P: float) -> np.ndarray:
    return Delta_max * (((L + P + Kd) - np.sqrt((L + P + Kd) ** 2 - 4 * L * P)) / (2 * P))


def fit_binding_models(result_df: pd.DataFrame, prot_conc: float, global_fit: bool = False) -> pd.DataFrame:
    ligand_cols = [col for col in result_df.columns if col.endswith('_csp')]
    ligand_concs = np.array([float(col[:-4]) for col in ligand_cols])
    fit_results: List[Dict[str, float]] = []
    if not global_fit:
        for _, row in result_df.iterrows():
            y_data = row[ligand_cols].values.astype(float)
            p0 = [np.max(y_data), 10.0]
            print(p0)
            try:
                popt, _ = curve_fit(
                    lambda L, Delta_max, Kd: binding_model(L, Delta_max, Kd, prot_conc),
                    ligand_concs,
                    y_data,
                    p0=p0,
                    bounds=([0, 0.001], [FITTING_CSP_MAX, 100000]),
                )
                Delta_max_fit, Kd_fit = popt
            except RuntimeError:
                Delta_max_fit, Kd_fit = np.nan, np.nan
            fit_results.append({
                'assignment': row['assignment'],
                'residue': row['residue'],
                'Delta_max': Delta_max_fit,
                'Kd': Kd_fit,
                'prot_conc': prot_conc,
            })
        return pd.DataFrame(fit_results)
    else:
        y_all = [row[ligand_cols].values.astype(float) for _, row in result_df.iterrows()]
        Delta_max_inits = [np.max(y) for y in y_all]
        p0 = np.array(Delta_max_inits + [10.0])
        N = len(result_df)
        print("Starting global fitting...")
        print("Number of assignments:", N, "parameters:", len(p0))
        print("Initial parameters:", p0)
        def global_residuals(params: np.ndarray) -> np.ndarray:
            Kd_global = params[-1]
            residuals: List[np.ndarray] = []
            for i, y in enumerate(y_all):
                Delta_max_i = params[i]
                pred = binding_model(ligand_concs, Delta_max_i, Kd_global, prot_conc)
                residuals.append(y - pred)
            return np.concatenate(residuals)
        res = least_squares(global_residuals, p0)
        params_opt = res.x
        Kd_global_opt = params_opt[-1]
        for i, (_, row) in enumerate(result_df.iterrows()):
            fit_results.append({
                'assignment': row['assignment'],
                'residue': row['residue'],
                'Delta_max': params_opt[i],
                'Kd': Kd_global_opt,
                'prot_conc': prot_conc,
            })
        return pd.DataFrame(fit_results)


def plot_binding_curves(result_df: pd.DataFrame, fit_results: Optional[pd.DataFrame], pdf_filename: str) -> None:
    ligand_cols_csp = [c for c in result_df.columns if c.endswith('_csp')]
    ligand_cols_int = [c for c in result_df.columns if c.endswith('_int') and c[:-4] != '0']
    ligand_concs = np.array([float(c[:-4]) for c in ligand_cols_csp])
    prot_conc = fit_results['prot_conc'][0] if fit_results is not None else None
    with PdfPages(pdf_filename) as pdf:
        # Max CSP summary
        fig, ax = plt.subplots(figsize=(6, 4))
        residues = result_df['residue'].values
        csp_max = result_df[ligand_cols_csp].max(axis=1).values
        ax.bar(residues, csp_max)
        ax.set_xlabel('Residue Number')
        ax.set_ylabel('Measured CSP (max)')
        pdf.savefig(fig)
        plt.close(fig)
        # Intensity ratio summary
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hlines(1, residues.min(), residues.max(), colors='k', linestyles='dashed')
        for label in ligand_cols_int:
            heights = result_df[label].values
            ax.bar(residues, heights, label=label[:-4])
        ax.set_xlabel('Residue Number')
        ax.set_ylabel('Intensity ratio (I/I0)')
        pdf.savefig(fig)
        plt.close(fig)
        # Fit summaries
        if fit_results is not None:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.bar(fit_results['residue'], fit_results['Delta_max'], color='green')
            ax.set_xlabel('Residue Number')
            ax.set_ylabel('Fitted CSP')
            pdf.savefig(fig)
            plt.close(fig)
            fig, ax = plt.subplots(figsize=(6, 4))
            kd_vals = fit_results['Kd'].astype(float).values
            kd_plot = np.where((kd_vals > 0) & np.isfinite(kd_vals), kd_vals, np.nan)
            ax.bar(fit_results['residue'], kd_plot, color='purple')
            ax.set_xlabel('Residue Number')
            ax.set_ylabel('Fitted Kd')
            ax.set_yscale('log')
            ax.set_ylim(1e-3, 1e4)
            pdf.savefig(fig)
            plt.close(fig)
        for _, row in result_df.iterrows():
            fig, ax = plt.subplots(figsize=(6, 4))
            y_data = row[ligand_cols_csp].values.astype(float)
            ax.scatter(ligand_concs, y_data, color='blue')
            assignment = row['assignment']
            residue = row['residue']
            ax.set_title(f'Assignment: {assignment} (Residue: {residue})')
            ax.set_xlabel('Ligand Concentration (µM)')
            ax.set_ylabel('CSP')
            if fit_results is not None:
                fit_row = fit_results[fit_results['assignment'] == assignment]
                if not fit_row.empty:
                    Delta_max = fit_row['Delta_max'].values[0]
                    Kd = fit_row['Kd'].values[0]
                    L_fit = np.linspace(ligand_concs.min(), ligand_concs.max(), 100)
                    y_fit = binding_model(L_fit, Delta_max, Kd, prot_conc)
                    ax.plot(L_fit, y_fit, color='red')
                    textstr = f'CSP = {Delta_max:.3f}\nKd = {Kd:.3f}'
                    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
                            fontsize=10, verticalalignment='top',
                            bbox=dict(boxstyle='round', fc='w'))
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Binding curves and summary plots saved to {pdf_filename}")


def print_significant_residues(fit_results: Optional[pd.DataFrame], csp_threshold: float = 0.01, kd_threshold: float = 10000.0) -> None:
    if fit_results is None:
        print("No fitting results available.")
        return
    significant = fit_results[(fit_results['Delta_max'] >= csp_threshold) & (fit_results['Kd'] <= kd_threshold)]
    if significant.empty:
        print(f"No residues found with fitted CSP >= {csp_threshold}")
        return
    significant = significant.sort_values('Delta_max', ascending=False)
    print(f"\nResidues with fitted CSP >= {csp_threshold}:")
    print("Residue\tCSP\tKd")
    print("-" * 30)
    for _, row in significant.iterrows():
        print(f"{int(row['residue'])}\t{row['Delta_max']:.3f}\t{row['Kd']:.3f}")
    residue_list = ','.join(str(int(r)) for r in significant['residue'])
    print(f"\nFor --residues argument: {residue_list}")


def process_folder(folder: str, minimize: bool = False, fit_mode: Optional[str] = None,
                   protein_conc: Optional[float] = None, filter_overlap: bool = False,
                   min_csp: float = 0.0001, residues: Optional[str] = None) -> None:
    file_pattern = os.path.join(folder, '*')
    all_files = [f for f in glob.glob(file_pattern) if f.lower().endswith(('.list'))]
    if not all_files:
        print("No peak list files (.list or .txt) found in folder:", folder)
        return
    print("\nPeak list ligand concentrations (uM):")
    peak_lists: List[Dict[str, object]] = []
    for fp in all_files:
        df = read_peak_list(fp)
        base = os.path.basename(fp)
        var = extract_variable_from_filename(base)
        label = str(var) if var is not None else base
        peak_lists.append({'file': fp, 'df': df, 'var': var, 'label': label})
        print(f"{base}\t{var}")
    
    if filter_overlap:
        print("\nFiltering overlapping peaks...")
        for entry in peak_lists:
            if 'overlap' in entry['df'].columns:
                before = len(entry['df'])
                entry['df'] = entry['df'][~entry['df']['overlap']].reset_index(drop=True)
                after = len(entry['df'])
                print(f"{os.path.basename(entry['file'])}: removed {before - after} overlapping rows")
    
    numeric_entries = [p for p in peak_lists if p['var'] is not None]
    ref_entry = min(numeric_entries, key=lambda x: x['var']) if numeric_entries else peak_lists[0]
    print("\nReference file:", os.path.basename(ref_entry['file']))
    
    df_ref = ref_entry['df']
    ref_label = ref_entry['label']
    csp_dict: Dict[str, Dict[str, Dict[str, float]]] = {}
    for _, row in df_ref.iterrows():
        ass = row['assignment']
        csp_dict.setdefault(ass, {})[ref_label] = {'csp': 0.0, 'int': 1.0}
    for entry in peak_lists:
        if entry['file'] == ref_entry['file']:
            continue
        print('\nProcessing Peak List:\t', os.path.basename(entry['file']))
        csp_df = compute_csp(ref_entry, entry, minimize=minimize)
        int_df = compute_intensity_ratio(df_ref, entry['df'])
        combined = pd.merge(csp_df, int_df, on='assignment')
        variable = entry['label']
        for _, row in combined.iterrows():
            ass = row['assignment']
            csp_dict.setdefault(ass, {})[variable] = {'csp': row['csp'], 'int': row['int_ratio']}
    df = pd.DataFrame.from_dict(csp_dict, orient='index')
    df.index.name = 'assignment'
    df.reset_index(inplace=True)
    flat: Dict[str, List[float]] = {'assignment': df['assignment']}
    for label in df.columns.drop('assignment'):
        flat[f'{label}_csp'] = df[label].apply(lambda d: d['csp'])
        flat[f'{label}_int'] = df[label].apply(lambda d: d['int'])
    result_df = pd.DataFrame(flat)
    result_df['residue'] = result_df['assignment'].apply(extract_residue_number)
    cols = ['assignment', 'residue'] + sorted([c for c in result_df.columns if c not in ('assignment', 'residue')], key=lambda x: float(x[:-4]) if x[:-4].replace('.', '', 1).isdigit() else -1)
    result_df = result_df[cols]
    ligand_cols_csp = [c for c in result_df.columns if c.endswith('_csp')]
    result_df['csp_max'] = result_df[ligand_cols_csp].max(axis=1) if ligand_cols_csp else 0.0
    csv_name = f"{folder}/csp_results{'_minimised' if minimize else ''}.csv"
    result_df.to_csv(csv_name, index=False)
    print(f"CSP results written to {csv_name}")
    fit_results: Optional[pd.DataFrame] = None
    if fit_mode:
        prot_conc = protein_conc if protein_conc is not None else float(input("Enter the constant protein concentration (µM): "))
        to_fit = result_df.copy()
        if residues:
            try:
                tokens = re.split(r"[,;\s]+", residues.strip())
                include = {int(t) for t in tokens if t}
            except Exception:
                print(f"Invalid residue list: '{residues}'. No residues will be fitted.")
                include = set()
            to_fit = result_df[result_df['residue'].isin(include)].copy() if include else result_df.iloc[0:0].copy()
        else:
            to_fit = result_df[result_df['csp_max'] >= min_csp].copy()
        if to_fit.empty:
            print(f"No residues meet the minimum CSP threshold ({min_csp}). Skipping fitting.")
        else:
            fit_results = fit_binding_models(to_fit, prot_conc, global_fit=(fit_mode == 'global'))
    pdf_name = f"{folder}/graphs{'_minimised' if minimize else ''}.pdf"
    plot_binding_curves(result_df, fit_results, pdf_name)
    if fit_results is not None and fit_mode == 'individual':
        print_significant_residues(fit_results, csp_threshold=min_csp)


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate CSPs and fit binding models from HSQC peak lists.")
    parser.add_argument('folder', help="Folder containing peak list files (\*.list or \*.txt).")
    parser.add_argument('-m', '--minimize', action='store_true', help="Apply drift minimisation when computing CSPs.")
    parser.add_argument('--fit', nargs='?', const='individual', choices=['global', 'individual'], help="Perform binding model fitting: 'individual' (default) or 'global'.")
    parser.add_argument('-conc', type=float, nargs='?', help="Protein concentration in µM for fitting.  If omitted and fitting is requested, you will be prompted.")
    parser.add_argument('-fo', '--filter_overlap', action='store_true', help="Remove rows flagged as overlapping (requires 'overlap' column).")
    parser.add_argument('--mincsp', type=float, default=0.0001, help="Minimum measured CSP to include a residue in fitting.")
    parser.add_argument('--residues', type=str, help="Comma/space/semicolon separated list of residue numbers for fitting; overrides --mincsp.")
    args = parser.parse_args()
    process_folder(args.folder, minimize=args.minimize, fit_mode=args.fit, protein_conc=args.conc,
                   filter_overlap=args.filter_overlap, min_csp=args.mincsp, residues=args.residues)


if __name__ == '__main__':
    main()