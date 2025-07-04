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

def extract_variable_from_filename(filename):
    """
    Look for numbers followed by mM, uM, or C in the filename.
    Returns the last match as a float if possible, otherwise None.
    """
    pattern = r'(\d+(?:\.\d+)?)(?:mM|uM|C|K)'
    matches = re.findall(pattern, filename)
    if matches:
        try:
            return float(matches[-1])
        except ValueError:
            return None
    return 

def get_variable_concentration(filename):
    """
    Try to extract the variable concentration from the filename.
    If extraction fails, ask the user to input a numeric value.
    """
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
    dims = assignment.split('-')

    return [dim[-1] for dim in dims]

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
    """
    Read the sparky HSQC peak list file.
    Expected tab-delimited format with columns:
    assignment    w1    w2    intensity
    """
    df = pd.read_csv(file_path, sep='\t')

    df['residue'] = df['assignment'].apply(extract_residue)
    df.sort_values(by='residue', inplace=True)

    return df

def compute_intensity_ratio(df_ref, df_target):
    """
    Compute intensity ratio I/I0 (target/reference) for each assignment.
    Returns a Series indexed by residue with the ratio.
    """
    merged = pd.merge(
        df_ref[['assignment', 'intensity']],
        df_target[['assignment', 'intensity']],
        on='assignment', suffixes=('_ref', '_target')
    )

    result = merged[['assignment']].copy()
    result['int_ratio'] = merged['intensity_target'] / merged['intensity_ref']
    
    return result

def compute_csp(df_ref, df_target, minimize=False):
    """
    For common assignments between df_ref and df_target, compute the differences 
    and then the total CSP using nitrogen weighting.
    
    If minimize is True, compute a constant offset (dx, dy) for the target 
    such that the mean difference is compensated before calculating the CSP.
    """
    merged = pd.merge(df_ref, df_target, on='assignment', suffixes=('_ref', '_target'))

    dims = extract_dimensions(merged['assignment'].iloc[0])

    print("\t-dim","atom","weight")
    print("\t-dw1:",dims[0],SHIFT_WEIGHTS[dims[0]])
    print("\t-dw2:",dims[1],SHIFT_WEIGHTS[dims[1]])

    dw1 = merged['w1_ref'] - merged['w1_target']
    dw2 = merged['w2_ref'] - merged['w2_target']

    if minimize:
        # Calculate constant offsets for the target relative to the reference.
        drift_dw1 = np.median(merged['w1_ref'] - merged['w1_target'])
        drift_dw2 = np.median(merged['w2_ref'] - merged['w2_target'])

        dw1 -= drift_dw1
        dw2 -= drift_dw2

        print(f'\t-drift correction: {dims[0]}={drift_dw1}, {dims[1]}={drift_dw2}')
    
    # Calculate total CSP with nitrogen weighting
    total_csp = np.sqrt((SHIFT_WEIGHTS[dims[0]]*dw1)**2 + (SHIFT_WEIGHTS[dims[1]] * dw2)**2)
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

    with PdfPages(pdf_filename) as pdf:

        # Summary plot: Measured maximum CSP vs residue number.
        fig, ax = plt.subplots(figsize=figsize)
        residues = result_df['residue'].values
        csp_max = result_df[ligand_cols_csp].max(axis=1).values
        ax.bar(residues, csp_max)
        ax.set_xlabel("Residue Number")
        ax.set_ylabel("Measured CSP (max)")
        #ax.set_title("Measured CSP vs Residue")
        pdf.savefig(fig)
        plt.close(fig)

        # Summary plot: Measured I/I0 vs residue number.
        fig, ax = plt.subplots(figsize=figsize)
        residues = result_df['residue'].values
        ax.hlines(1,residues.min(),residues.max(), colors='k', linestyles='dashed')

        for i, label in enumerate(ligand_cols_int):
            heights = result_df[label].values
            ax.bar(residues, heights, label=col_sort_key(label))

        ax.set_xlabel("Residue Number")
        ax.set_ylabel("Intensity ratio (I/I0)")
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
            
            # Fitted Kd vs residue.
            fig, ax = plt.subplots(figsize=figsize)
            ax.bar(fit_results['residue'], fit_results['Kd'], color='purple')
            ax.set_xlabel("Residue Number")
            ax.set_ylabel("Fitted Kd")
            ax.set_title("Fitted Kd vs Residue")
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
    ligand_concs = np.array([float(col) for col in ligand_cols])
    
    fit_results = []

    if not global_fit:
        # Individual fit for each assignment
        for idx, row in result_df.iterrows():
            # y data are the CSP values for this assignment
            y_data = row[ligand_cols].values.astype(float)
            # Initial guess: Delta_max as max observed CSP and Kd arbitrarily as 10
            p0 = [np.max(y_data), 10.0]
            
            try:
                popt, _ = curve_fit(lambda L, Delta_max, Kd: binding_model(L, Delta_max, Kd, prot_conc),
                                    ligand_concs, 
                                    y_data, 
                                    p0=p0,
                                    bounds=([0,0],[FITTING_CSP_MAX,np.inf]))
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
        # Initial guesses: for each assignment, Delta_max = max(y) and global Kd = 10.
        Delta_max_inits = [np.max(row[ligand_cols].values.astype(float)) for _, row in result_df.iterrows()]
        initial_Kd = 10.0
        # Parameters vector: [Delta_max_1, Delta_max_2, ..., Delta_max_N, Kd]
        p0 = np.array(Delta_max_inits + [initial_Kd])
        
        # Concatenate observed data for all assignments into one vector.
        y_all = np.concatenate([row[ligand_cols].values.astype(float) for _, row in result_df.iterrows()])
        
        def global_residuals(params):
            Kd_global = params[-1]
            residuals = []
            for i, row in enumerate(result_df.itertuples(index=False)):
                Delta_max_i = params[i]
                pred = binding_model(ligand_concs, Delta_max_i, Kd_global, prot_conc)
                # Assuming the order of ligand columns is preserved
                y_obs = np.array(row[2:])  
                residuals.append(y_obs - pred)
            return np.concatenate(residuals)
        
        res = least_squares(global_residuals, p0)
        params_opt = res.x
        Kd_global_opt = params_opt[-1]
        for i, row in result_df.iterrows():
            fit_results.append({
                'assignment': row['assignment'],
                'residue': row['residue'],
                'Delta_max': params_opt[i],
                'Kd': Kd_global_opt,
                'prot_conc': prot_conc
            })
        return pd.DataFrame(fit_results)

def main(folder, args):
    # Find all .list files in the folder.
    file_pattern = os.path.join(folder, "*.list")
    file_paths = glob.glob(file_pattern)
    if not file_paths:
        print("No .list files found in the folder:", folder)
        return
    
    # Process each file: read the data and extract the variable.
    print("\nPeak List Variable Values:")
    peak_lists = []
    for file_path in file_paths:
        df = read_peak_list(file_path)
        base_name = os.path.basename(file_path)
        var = get_variable_concentration(base_name)
        col_label = str(var)
        peak_lists.append({'file': file_path, 'df': df, 'var': var, 'label': col_label})
        
        print(file_path.split('/')[-1],'\t',var)

    # Filter overlap?
    if args.filter_overlap:
        print("\nFiltering overlapping peaks...")
        for idx, plist in enumerate(peak_lists):
            print(os.path.basename(plist['file']))
            print(df[df.overlap]['residue'].tolist())
            df = plist['df']
            peak_lists[idx]['df'] = df.drop(df[df.overlap].index)
    
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
        print('Processing Peak List:\t',entry['file'])
        if entry['file'] == ref_entry['file']:
            print('\t-reference peak list')
            continue
        
        csp_df = compute_csp(df_ref, entry['df'], minimize=args.minimize)
        int_df = compute_intensity_ratio(df_ref, entry['df'])
        
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
    for label in df.columns.drop('assignment'):
        flat_data[f'{label}_csp'] = df[label].apply(lambda d: d['csp'])
        flat_data[f'{label}_int'] = df[label].apply(lambda d: d['int'])

    result_df = pd.DataFrame(flat_data)
    result_df['residue'] = result_df['assignment'].apply(extract_residue)
    
    # Grab CSP columns to separate df (we ignore int ratios here [for now])
    cols = ['assignment', 'residue'] + sorted([col for col in result_df if col not in ['assignment','residue']], key=col_sort_key)
    result_df = result_df[cols]
    result_df.to_csv(f"{folder}/csp_results{'_minimised' if args.minimize else ''}.csv", index=False)

    fit_results=None

    if args.fit is not None:
        prot_conc = args.conc if args.conc is not None else float(input("Enter the constant protein concentration: "))
        fit_results = fit_binding_models(result_df, prot_conc, global_fit=(args.fit == 'global'))  # or global_fit=True
    
    plot_binding_curves(result_df, fit_results=fit_results, pdf_filename=f"{folder}/graphs{'_minimised' if args.minimize else ''}.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate CSP from sparky HSQC peak lists in a folder."
    )
    parser.add_argument("folder", help="Folder containing .list files to process.")
    parser.add_argument("-m","--minimize", action="store_true",
                        help="Enable drift compensation using least squares minimization of offsets.")
    parser.add_argument("--fit", nargs="?", const="individual", choices=["global", "individual"],
                    help="Activate binding model fitting. Optionally specify 'global' for a global Kd, or 'individual' (default) for per-assignment fits.")
    parser.add_argument("-conc", type=float, nargs="?", help="Provide the concentration of labeled protein for the fitting algorithm. ")
    parser.add_argument("-fo", "--filter_overlap", action="store_true", help="Remove data from peaks that overlap too much to determine CSPs. Determined by 'overlap' column in peak list file.")

    args = parser.parse_args()
    main(args.folder, args)
