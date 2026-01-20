import pandas as pd
import numpy as np
import re
import argparse
import os
import glob
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

FITTING_CSP_MAX = 10
SHIFT_WEIGHTS = {'H':1, 'N':0.154}

def split_filename_properties(filename):
    """
    Break a filename into property-like tokens.
    Splits on '_' first (user-declared properties) and then on '.' to catch HSQC.# chunks.
    """
    base = os.path.basename(filename)
    stem = re.sub(r'\.list$', '', base, flags=re.IGNORECASE)
    tokens = []
    for chunk in stem.split('_'):
        tokens.extend(chunk.split('.'))
    return [t for t in tokens if t]

def extract_variable_from_filename(filename):
    # Look for numeric value followed by a concentration unit
    # Pattern matches: number + unit (e.g., "500uM", "50mM", "5nM")
    pattern = re.compile(r"(?i)(\d+(?:\.\d+)?)(mM|uM|nM|M|µM)")
    matches = []
    for token in split_filename_properties(filename):
        matches.extend(pattern.findall(token))
    if not matches:
        matches = pattern.findall(filename)
    if not matches:
        return None

    # Prefer the first occurrence (typically the ligand concentration at start of filename)
    for num_str, unit in matches:
        try:
            val = float(num_str)
        except ValueError:
            continue

        u = unit.lower().replace('µ', 'u')
        if u == 'mm':
            return val * 1000.0
        elif u == 'um':
            return val
        elif u == 'nm':
            return val / 1000.0
        elif u == 'm':
            return val * 1e6

    return None

def extract_temperature_from_filename(filename):
    """
    Extract temperature from filename tokens. Returns Celsius if found, else None.
    Accepts tokens like '25C' or '298K'.
    """
    for token in split_filename_properties(filename):
        m = re.search(r'(?i)(\d+(?:\.\d+)?)([ck])', token)
        if not m:
            continue
        try:
            temp_val = float(m.group(1))
        except ValueError:
            continue
        unit = m.group(2).lower()
        if unit == 'k':
            return temp_val - 273.15
        return temp_val
    return None

def get_variable_concentration(filename):
    var = extract_variable_from_filename(filename)
    if var is None:
        while True:
            user_input = input(f"Could not extract variable concentration from '{filename}'. Please provide a value: ")
            try:
                var = float(user_input)
                break
            except ValueError:
                print("Invalid input. Please enter a numeric value.")
    return var

def extract_dimensions(assignment):
    parts = [p.strip() for p in str(assignment).split('-') if p.strip() and p.strip() != '?']
    dims = [p[-1] for p in parts if p]
    
    if len(dims) >= 2:
        return dims[:2]
    else:
        print(f"Warning: couldn't parse two dimensions from assignment '{assignment}'; defaulting to ['H','N']")
        return ['H', 'N']

def extract_residue(assignment):
    m = re.search(r'(\d+)', assignment)
    if m:
        return int(m.group(1))
    return float('inf')

def col_sort_key(col):
    col = col.split('_')[0]
    try:
        return float(col)
    except ValueError:
        return -1

def read_peak_list(file_path):
    print(f'\nReading Peak List:\t', os.path.basename(file_path))
    df = pd.read_csv(file_path, sep='\t')
    
    if len(df) == 0:
        return df
    
    # Derive column mapping from the first assignment token string
    assign_column = str(df['assignment'].iloc[0])
    assignments = [tok.strip() for tok in assign_column.split('-')]

    print(f'\tDetected assignment tokens: {assignments}')

    col_mapping = {}
    drop_cols = []
    for idx, assignment in enumerate(assignments, start=1):
        w_col = f'w{idx}'
        if assignment == '?' or assignment == '':
            drop_cols.append(w_col)
            continue
        nucleus = assignment[-1].upper()
        if nucleus in ('N', 'H'):
            # Map this w# to a standardized chemical shift column
            target = f'CS_{nucleus}'
            if target in col_mapping.values():
                print(f"Warning: duplicate nucleus '{nucleus}' found in assignment '{assign_column}'.")
            col_mapping[w_col] = target
        else:
            print(f"Warning: unrecognized nucleus in token '{assignment}' from assignment '{assign_column}'.")

        print(f'\tMapping column {w_col} to {col_mapping.get(w_col, "DROP")} (1st entry = {df[w_col].iloc[0]}ppm)')

    # Drop pseudo-dimension columns if they exist
    for col in drop_cols:
        if col in df.columns:
            df = df.drop(columns=[col])

    # Apply renaming for real dimensions
    if col_mapping:
        df = df.rename(columns=col_mapping)

    # Remove the pseudo-dimension token from all assignments if present at the front
    if assignments and assignments[0] == '?':
        df['assignment'] = df['assignment'].astype(str).str.replace(r'^\?-', '', regex=True)
    
    df['residue'] = df['assignment'].apply(extract_residue)
    df.sort_values(by='residue', inplace=True)

    return df

def compute_intensity_ratio(ref, target):
    """Compute I/I0 using an outer merge so missing peaks become NaN (create gaps)."""
    print(f'\tComputing intensity ratios...')
    merged = pd.merge(
        ref['df'][['assignment', 'intensity']],
        target['df'][['assignment', 'intensity']],
        on='assignment', suffixes=('_ref', '_target')
    )

    result = merged[['assignment']].copy()
    ref_int = merged['intensity_ref']
    tgt_int = merged['intensity_target']
    result['int_ratio'] = np.where(ref_int == 0, np.nan, tgt_int / ref_int)
    
    return result

def compute_csp(ref, target, minimize=False):
    print(f'\tComputing CSPs...')

    # Outer merge keeps peaks that are missing in one list so they plot as gaps.
    merged = pd.merge(ref['df'], target['df'], on='assignment', suffixes=('_ref', '_target'))

    # Locate standardized chemical shift columns (CS_N / CS_H) on the reference side.
    cs_cols_ref = [c for c in merged.columns if c.startswith('CS_') and c.endswith('_ref')]
    n_col_ref = 'CS_N_ref' if 'CS_N_ref' in cs_cols_ref else None
    h_col_ref = 'CS_H_ref' if 'CS_H_ref' in cs_cols_ref else None

    if not n_col_ref or not h_col_ref:
        raise ValueError(f"Could not identify CS_N/CS_H columns. Found: {cs_cols_ref}")
    
    n_col_target = n_col_ref.replace('_ref', '_target')
    h_col_target = h_col_ref.replace('_ref', '_target')

    if n_col_target not in merged.columns or h_col_target not in merged.columns:
        raise ValueError("Missing CS_N/CS_H target columns after merge.")
    
    # Calculate chemical shift differences
    dn = merged[n_col_ref] - merged[n_col_target]
    dh = merged[h_col_ref] - merged[h_col_target]

    if minimize:
        # Calculate constant offsets for the target relative to the reference
        drift_n = np.median(dn)
        drift_h = np.median(dh)

        dn -= drift_n
        dh -= drift_h

        print(f'\t-drift correction: N={drift_n}, H={drift_h}')
    
    # Calculate total CSP with nitrogen weighting
    total_csp = np.sqrt((SHIFT_WEIGHTS['N']*dn)**2 + (SHIFT_WEIGHTS['H']*dh)**2)
    result = merged[['assignment']].copy()
    result['csp'] = total_csp
    
    return result

def binding_model(L, Delta_max, Kd, P):
    """
    Binding model function.
    L: ligand concentration (array)
    Delta_max: maximum CSP
    Kd: dissociation constant
    P: constant protein concentration
    """
    return Delta_max * (((L + P + Kd) - np.sqrt((L + P + Kd)**2 - 4 * L * P)) / (2 * P))

def plot_binding_curves(result_df, fit_results, pdf_filename):
    """
    Create a multi-page PDF of binding curves and summary plots.
    
    Parameters:
      result_df: DataFrame with columns 'assignment', 'residue', and CSP columns (with ligand concentration labels).
      ligand_cols: list of column names corresponding to ligand concentrations (convertible to float).
      prot_conc: constant protein concentration (needed for plotting the fitted curves).
      ligand_concs: array of ligand concentration values (floats) corresponding to ligand_cols.
      fit_results: Optional DataFrame with columns 'assignment', 'residue', 'Delta_max', 'Kd'.
    """
    figsize=(6,4)

    ligand_cols_csp = [col for col in result_df.columns if col.endswith('_csp')]
    ligand_cols_int = [col for col in result_df.columns if col.endswith('_int') and col_sort_key(col) != 0]

    print(ligand_cols_int)

    ligand_concs = np.array([col_sort_key(col) for col in ligand_cols_csp])
    if fit_results is not None: 
        prot_conc = fit_results['prot_conc'][0]

    # Build a plotting frame with all residue numbers present so gaps show explicitly.
    plot_df = result_df.copy()
    if not plot_df.empty:
        residues_full = np.arange(plot_df['residue'].min(), plot_df['residue'].max() + 1)
        plot_df = plot_df.set_index('residue').reindex(residues_full).reset_index().rename(columns={'index': 'residue'})

    with PdfPages(pdf_filename) as pdf:

        # Summary plot: CSPs vs residue for all ligand points (line + markers per series).
        fig, ax = plt.subplots(figsize=figsize)
        residues = plot_df['residue'].values
        for label in ligand_cols_csp:
            y_vals = plot_df[label].values.astype(float)
            y_vals = np.where(np.isfinite(y_vals), y_vals, np.nan)  # break lines across missing points
            ax.plot(residues, y_vals, marker='o', label=col_sort_key(label), markersize=2)
        ax.set_xlabel("Residue Number")
        ax.set_ylabel("CSP")
        ax.legend(title="Ligand")
        pdf.savefig(fig)
        plt.close(fig)

        # Summary plot: Measured I/I0 vs residue number (line + markers per series).
        fig, ax = plt.subplots(figsize=figsize)
        residues = plot_df['residue'].values
        ax.hlines(1, residues.min(), residues.max(), colors='k', linestyles='dashed')

        for label in ligand_cols_int:
            heights = plot_df[label].values.astype(float)
            heights = np.where(np.isfinite(heights), heights, np.nan)  # break lines across missing points
            ax.plot(residues, heights, marker='o', label=col_sort_key(label), markersize=2)

        ax.set_xlabel("Residue Number")
        ax.set_ylabel("Intensity ratio (I/I0)")
        ax.legend(title="Ligand")
        pdf.savefig(fig)
        plt.close(fig)

        if fit_results is not None:

            # Fitted CSP vs residue.
            fig, ax = plt.subplots(figsize=figsize)
            ax.bar(fit_results['residue'], fit_results['Delta_max'], color='green')
            ax.set_xlabel("Residue Number")
            ax.set_ylabel("Fitted CSP")
            ax.set_title("Fitted CSP vs Residue")
            pdf.savefig(fig)
            plt.close(fig)
            
            # Fitted Kd vs residue (log scale).
            fig, ax = plt.subplots(figsize=figsize)
            kd = fit_results['Kd'].values.astype(float)
            # Avoid log-scale errors by masking non-finite or non-positive values
            kd_plot = kd.copy()
            kd_plot[~np.isfinite(kd_plot)] = np.nan
            kd_plot[kd_plot <= 0] = np.nan
            ax.bar(fit_results['residue'], kd_plot, color='purple')
            ax.set_xlabel("Residue Number")
            ax.set_ylabel("Fitted Kd")
            ax.set_yscale('log')
            ax.set_ylim(0.001, 1e4)
            ax.set_title("Fitted Kd vs Residue (log scale)")
            pdf.savefig(fig)
            plt.close(fig)

            # Individual binding curves for each assignment.
            for idx, row in result_df.iterrows():
                fig, ax = plt.subplots(figsize=figsize)
                # Get the measured data points.
                y_data = row[ligand_cols_csp].values.astype(float)
                ax.scatter(ligand_concs, y_data, color='blue', label='Data')
                
                assignment = row['assignment']
                residue = row['residue']
                title = f"Assignment: {assignment} (Residue: {residue})"
                ax.set_title(title)
                ax.set_xlabel("Ligand Concentration")
                ax.set_ylabel("CSP")
                
                # If fitting was performed, overlay the fitted curve and annotate parameters.
                
                fit_row = fit_results[fit_results['assignment'] == assignment]
                if not fit_row.empty:
                    Delta_max = fit_row['Delta_max'].values[0]
                    Kd = fit_row['Kd'].values[0]
                    # Generate a smooth curve.
                    L_fit = np.linspace(min(ligand_concs), max(ligand_concs), 100)
                    y_fit = binding_model(L_fit, Delta_max, Kd, prot_conc)
                    ax.plot(L_fit, y_fit, color='red', label='Fit')
                    # Add annotation with the fitted parameters.
                    textstr = f"CSP = {Delta_max:.3f}\nKd = {Kd:.3f}"
                    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                            verticalalignment='top', bbox=dict(boxstyle="round", fc="w"))
                    #ax.legend()
                
                pdf.savefig(fig)
                plt.close(fig)
                
    print(f"Binding curves and summary plots saved to {pdf_filename}")

def fit_binding_models(result_df, prot_conc, global_fit=False):
    # Identify ligand concentration columns (those not 'assignment' or 'residue')
    ligand_cols = [col for col in result_df.columns if col.endswith('_csp')]
    # Convert ligand column names to floats (they represent the ligand concentrations)
    ligand_concs = np.array([float(col[0:-4]) for col in ligand_cols])

    # Build a reusable default Kd guess from the ligand range (fallback to 10 uM if unknown).
    finite_ligand = ligand_concs[np.isfinite(ligand_concs) & (ligand_concs > 0)]
    kd_guess_default = float(np.nanmedian(finite_ligand)) if finite_ligand.size else 10.0
    kd_guess_default = kd_guess_default if np.isfinite(kd_guess_default) and kd_guess_default > 0 else 10.0

    def make_initial_guess(y_vals):
        """Return (Delta_max, Kd) guesses that stay positive and ignore NaNs."""
        finite_y = y_vals[np.isfinite(y_vals)]
        delta_guess = max(finite_y)
        if not np.isfinite(delta_guess) or delta_guess <= 0:
            delta_guess = 0.05
        kd_guess = kd_guess_default
        return delta_guess, kd_guess
    
    fit_results = []

    if not global_fit:
        # Individual fit for each assignment
        print(f'\nPerforming individual fits for each assignment...')
        for idx, row in result_df.iterrows():
            # y data are the CSP values for this assignment
            y_data = row[ligand_cols].values.astype(float)
            print(f"{row['assignment']} (residue {row['residue']})")
            finite_mask = np.isfinite(y_data)
            if finite_mask.sum() < 2:
                # Not enough finite points to perform a fit.
                Delta_max_fit, Kd_fit = np.nan, np.nan
                fit_results.append({
                    'assignment': row['assignment'],
                    'residue': row['residue'],
                    'Delta_max': Delta_max_fit,
                    'Kd': Kd_fit,
                    'prot_conc': prot_conc
                })
                continue

            delta_guess, kd_guess = make_initial_guess(y_data[finite_mask])

            try:
                popt, _ = curve_fit(lambda L, Delta_max, Kd: binding_model(L, Delta_max, Kd, prot_conc),
                                    ligand_concs[finite_mask], 
                                    y_data[finite_mask], 
                                    p0=[delta_guess, kd_guess],
                                    bounds=([0,0.001],[FITTING_CSP_MAX,100000]))
                Delta_max_fit, Kd_fit = popt
            except RuntimeError:
                Delta_max_fit, Kd_fit = np.nan, np.nan
            
            fit_results.append({
                'assignment': row['assignment'],
                'residue': row['residue'],
                'Delta_max': Delta_max_fit,
                'Kd': Kd_fit,
                'prot_conc': prot_conc
            })
        return pd.DataFrame(fit_results)
    else:
        # Global fit: one global Kd and individual Delta_max for each assignment.
        N = len(result_df)
        # Initial guesses: per-assignment Delta_max uses finite max; global Kd uses median ligand concentration.
        def finite_max(vals):
            finite_vals = vals[np.isfinite(vals)]
            if finite_vals.size == 0:
                return 0.05
            guess = float(np.nanmax(finite_vals))
            return guess if np.isfinite(guess) and guess > 0 else 0.05

        Delta_max_inits = [finite_max(row[ligand_cols].values.astype(float)) for _, row in result_df.iterrows()]
        initial_Kd = kd_guess_default
        # Parameters vector: [Delta_max_1, Delta_max_2, ..., Delta_max_N, Kd]
        p0 = np.array(Delta_max_inits + [initial_Kd])
        
        # Concatenate observed data for all assignments into one vector.
        y_all = [row[ligand_cols].values.astype(float) for _, row in result_df.iterrows()]

        if not any(np.isfinite(arr).any() for arr in y_all):
            print("No finite CSP data available for global fitting.")
            return pd.DataFrame()

        print("Starting global fitting...")
        print("Number of assignments:", N, "parameters:", len(p0))
        print("Initial parameters:", p0)

        def global_residuals(params):
            Kd_global = params[-1]
            residuals = []
            for i, row in enumerate(y_all):
                Delta_max_i = params[i]
                pred = binding_model(ligand_concs, Delta_max_i, Kd_global, prot_conc)
                mask = np.isfinite(row)
                if not np.any(mask):
                    continue
                residuals.append(row[mask] - pred[mask])
            return np.concatenate(residuals)
        
        lower_bounds = [0.0]*N + [0.001]
        upper_bounds = [FITTING_CSP_MAX]*N + [100000.0]
        res = least_squares(global_residuals, p0, bounds=(lower_bounds, upper_bounds))
        params_opt = res.x
        Kd_global_opt = params_opt[-1]
        for i in range(len(result_df)):
            fit_results.append({
                'assignment': result_df.iloc[i]['assignment'],
                'residue': result_df.iloc[i]['residue'],
                'Delta_max': params_opt[i],
                'Kd': Kd_global_opt,
                'prot_conc': prot_conc
            })
        return pd.DataFrame(fit_results)

def print_significant_residues(fit_results, csp_threshold=0.01, kd_threshold=10000):
    """Print residue numbers with fitted CSP (Delta_max) values above the threshold.
    Useful for identifying residues to include in global fitting."""
    if fit_results is None:
        print("No fitting results available.")
        return
    
    # Get residues with significant CSPs
    significant = fit_results[fit_results['Delta_max'] >= csp_threshold]
    significant = significant[significant['Kd'] <= kd_threshold]
    if significant.empty:
        print(f"No residues found with fitted CSP >= {csp_threshold}")
        return
    
    # Sort by CSP value descending
    significant = significant.sort_values('Delta_max', ascending=False)
    
    print(f"\nResidues with fitted CSP >= {csp_threshold}:")
    print("Residue\tCSP\tKd")
    print("-" * 30)
    for _, row in significant.iterrows():
        print(f"{int(row['residue'])}\t{row['Delta_max']:.3f}\t{row['Kd']:.3f}")
    
    # Print a comma-separated list for easy --residues input
    residue_list = ','.join(map(str, map(int, significant['residue'].tolist())))
    print(f"\nFor --residues argument: {residue_list}")

def main(folder, args):
    # Find all .list files in the folder.
    file_pattern = os.path.join(folder, "*.list")
    file_paths = glob.glob(file_pattern)
    if not file_paths:
        print("No .list files found in the folder:", folder)
        return
    
    # Process each file: read the data and extract the variable.
    peak_lists = []
    for file_path in file_paths:
        base_name = os.path.basename(file_path)
        temp_c = extract_temperature_from_filename(base_name)

        if args.temp is not None:
            if temp_c is None:
                print(f"Skipping {base_name}: no temperature found in filename while --temp was requested.")
                continue
            if not np.isclose(temp_c, args.temp, atol=0.5):
                print(f"Skipping {base_name}: temperature {temp_c:.2f}C does not match requested {args.temp}C.")
                continue

        df = read_peak_list(file_path)
        var = get_variable_concentration(base_name)
        col_label = str(var)
        peak_lists.append({'file': file_path, 'df': df, 'var': var, 'label': col_label, 'temp_c': temp_c})

    if not peak_lists:
        if args.temp is not None:
            print(f"No .list files matching temperature {args.temp}C found in the folder: {folder}")
        return

    print("\nIdentified peak lists and variable concentrations:")
    for entry in peak_lists:
        temp_info = f"{entry['temp_c']:.2f}C" if entry.get('temp_c') is not None else "temp=?"
        print(f"\t{os.path.basename(entry['file'])}\t{entry['var']}\t{temp_info}")

    # Choose the reference peak list as the one with the lowest numeric variable value.
    numeric_lists = [p for p in peak_lists if p['var'] is not None]
    if numeric_lists:
        ref_entry = min(numeric_lists, key=lambda x: x['var'])
    else:
        ref_entry = peak_lists[0]
    
    print("\nReference file:", ref_entry['file'])
    df_ref = ref_entry['df']
    ref_label = ref_entry['label']
    
    # Create a dictionary to hold CSP results per assignment.
    csp_dict = {}
    
    # The reference file has a CSP of 0 and int ratio of 1 for every assignment.
    for _, row in df_ref.iterrows():
        ass = row['assignment']
        csp_dict.setdefault(ass, {})[ref_label] = {'csp':0.0, 'int': 1}

    # Compute CSP for each non-reference peak list.
    for entry in peak_lists:
        print('\nProcessing Peak List:\t',os.path.basename(entry['file']))
        if entry['file'] == ref_entry['file']:
            print('\tReference peak list')
            continue
        
        csp_df = compute_csp(ref_entry, entry, minimize=args.minimize)
        int_df = compute_intensity_ratio(ref_entry, entry)
        
        data_df = pd.merge(csp_df, int_df, on='assignment')

        variable = entry['label']

        for _, row in data_df.iterrows():
            ass = row['assignment']
            csp_val = row['csp']
            int_val = row['int_ratio']
            csp_dict.setdefault(ass, {})[variable] = {'csp':csp_val, 'int':int_val}

    # Convert the results into a DataFrame.
    df = pd.DataFrame.from_dict(csp_dict, orient='index')
    df.index.name = 'assignment'
    df.reset_index(inplace=True)
    
    flat_data = {'assignment': df['assignment']}
    def safe_lookup(cell, key):
        """Return nested dict value or NaN when the cell is missing."""
        return cell.get(key, np.nan) if isinstance(cell, dict) else np.nan
    for label in df.columns.drop('assignment'):
        flat_data[f'{label}_csp'] = df[label].apply(lambda d: safe_lookup(d, 'csp'))
        flat_data[f'{label}_int'] = df[label].apply(lambda d: safe_lookup(d, 'int'))

    result_df = pd.DataFrame(flat_data)
    result_df['residue'] = result_df['assignment'].apply(extract_residue)
    
    # Grab CSP columns to separate df (we ignore int ratios here [for now])
    cols = ['assignment', 'residue'] + sorted([col for col in result_df if col not in ['assignment','residue']], key=col_sort_key)
    result_df = result_df[cols]
    # Compute measured max CSP per residue for filtering/inspection
    ligand_cols_csp = [col for col in result_df.columns if col.endswith('_csp')]
    if ligand_cols_csp:
        result_df['csp_max'] = result_df[ligand_cols_csp].max(axis=1)
    else:
        result_df['csp_max'] = 0.0

    result_df.to_csv(f"{folder}/csp_results{'_minimised' if args.minimize else ''}.csv", index=False)

    fit_results=None

    if args.fit is not None:
        prot_conc = args.conc if args.conc is not None else float(input("Enter the constant protein concentration (uM): "))
        # Determine which residues to fit.
        # If the user provided an explicit residue list via --residues, prefer that (overrides --mincsp).
        to_fit = result_df.copy()
        if args.residues:
            # Accept comma, space or semicolon separated lists. Ignore empty tokens.
            try:
                tokens = re.split(r"[,;\s]+", args.residues.strip())
                print(tokens)
                include = set(int(t) for t in tokens if t != "")
            except Exception:
                print(f"Invalid residue list provided: '{args.residues}'. Expected comma/space/semicolon-separated integers.")
                include = set()

            if not include:
                print("No valid residues parsed from --residues; no residues will be fitted.")
                to_fit = result_df.iloc[0:0].copy()
            else:
                to_fit = result_df[result_df['residue'].isin(include)].copy()
        else:
            # Fall back to mincsp threshold filtering
            ligand_cols_csp = [col for col in result_df.columns if col.endswith('_csp')]
            if ligand_cols_csp:
                to_fit = result_df[result_df['csp_max'] >= args.mincsp].copy()
            else:
                to_fit = result_df.copy()

        if to_fit.empty:
            print(f"No residues meet the minimum CSP threshold ({args.mincsp}). Skipping fitting.")
            fit_results = None
        else:
            fit_results = fit_binding_models(to_fit, prot_conc, global_fit=(args.fit == 'global'))  # or global_fit=True
    
    plot_binding_curves(result_df, fit_results=fit_results, pdf_filename=f"{folder}/graphs{'_minimised' if args.minimize else ''}.pdf")
    
    # If we did individual fits, print residues that might be good candidates for global fitting
    if fit_results is not None and args.fit == 'individual':
        print_significant_residues(fit_results, args.mincsp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate CSP from sparky HSQC peak lists in a folder."
    )
    parser.add_argument("folder", help="Folder containing .list files to process.")
    parser.add_argument("-m","--minimize", action="store_true",
                        help="Enable drift compensation using least squares minimization of offsets.")
    parser.add_argument("--temp", type=float, default=None,
                        help="Filter input files to a specific temperature in Celsius. Temperature is read from filename tokens containing a value with C or K.")
    parser.add_argument("--fit", nargs="?", const="individual", choices=["global", "individual"],
                    help="Activate binding model fitting. Optionally specify 'global' for a global Kd, or 'individual' (default) for per-assignment fits.")
    parser.add_argument("-conc", type=float, nargs="?", help="Provide the concentration of labeled protein for the fitting algorithm. ")
    parser.add_argument("--mincsp", type=float, default=0.0001,
                        help="Minimum measured maximum CSP (across ligand conditions) required for a residue to be included in fitting.")
    parser.add_argument("--residues", type=str, default=None,
                        help="Comma-, space- or semicolon-separated list of residue numbers to include in fitting. If provided, this overrides --mincsp.")

    args = parser.parse_args()
    main(args.folder, args)
