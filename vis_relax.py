import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Patch

# Pattern to parse filenames: e.g. 50mM_5C_600MHz_T1_results_filtered.out
FNAME_RE = re.compile(
    r"(?P<salt>\d+mM)_(?P<temp>\d+)C_(?P<field>\d+)MHz_(?P<exp>\w+)_results_filtered\.out"
)

def collect_data(directory):
    """
    Walks through `directory`, finds all `*_results_filtered.out` files,
    parses metadata from filename, and loads the data into a single DataFrame.
    """
    rows = []
    for path in glob.glob(os.path.join(directory, "*_results_filtered.out")):
        fname = os.path.basename(path)
        print(f"Processing file: {fname}")
        m = FNAME_RE.match(fname)
        print(f"Matched: {m.groupdict() if m else 'No match'}")
        
        # Skip files that don't match the expected pattern
        if not m:
            continue
        info = m.groupdict()
        df = pd.read_csv(path, sep=" ", comment="#")
        
        # Name columns: residue, rate, error
        ncols = df.shape[1]
        col_names = ['residue', 'rate', 'error']
        df.columns = col_names
        
        df['salt'] = info['salt']
        df['temp_C'] = int(info['temp'])
        df['field_MHz'] = int(info['field'])
        df['experiment'] = info['exp']
        rows.append(df)
    data = pd.concat(rows, ignore_index=True)
    return data


def get_protein_name(directory):
    """Extracts the protein name from the directory path."""
    return os.path.basename(os.path.abspath(directory))

def plot_relaxation(data, experiment='T1', protein_name=None, plot_type='bar'):
    """
    Creates a multi-panel figure of relaxation rates vs field strength,
    with one panel per temperature, and lines colored by salt concentration.
    plot_type: 'bar' (default) or 'line'
    """
    subset = data[data['experiment'] == experiment]
    temps = sorted(subset['temp_C'].unique())
    fields = sorted(subset['field_MHz'].unique())
    salts = sorted(subset['salt'].unique(), key=lambda x: int(x.replace('mM','')))

    # Create a color map for salts
    cmap = plt.get_cmap('viridis')
    norm = mpl.colors.Normalize(vmin=min([int(s.replace('mM','')) for s in salts]), vmax=max([int(s.replace('mM','')) for s in salts]))
    salt_to_color = {salt: cmap(norm(int(salt.replace('mM','')))) for salt in salts}

    nrows = len(temps)
    ncols = len(fields)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), sharex=True, sharey=True, num=experiment)
    # If only one row or column, axes may not be 2D, so ensure it's always 2D
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]

    for i, temp in enumerate(temps):
        for j, field in enumerate(fields):
            ax = axes[i, j]
            df_t = subset[(subset['temp_C'] == temp) & (subset['field_MHz'] == field)]
            for salt in salts:
                df_ts = df_t[df_t['salt'] == salt]
                if not df_ts.empty:
                    if plot_type == 'bar':
                        ax.bar(
                            df_ts['residue'], df_ts['rate'], yerr=df_ts.get('error'),
                            color=salt_to_color[salt], label=f"{salt}"
                        )
                    elif plot_type == 'line':
                        # Sort by residue for line plot
                        df_ts_sorted = df_ts.sort_values('residue')
                        ax.errorbar(
                            df_ts_sorted['residue'], df_ts_sorted['rate'], yerr=df_ts_sorted.get('error'),
                            color=salt_to_color[salt], label=f"{salt}", marker='o', markersize=1, linestyle='-'
                        )
            if i == 0:
                ax.set_title(f"{field} MHz")
            if j == 0:
                ax.set_ylabel(f"{temp} °C\n{experiment} (s$^{{-1}}$)")
            # Clean up x-axis: show as 'Residue (ProteinName)'
            if i == len(temps) - 1:
                if protein_name:
                    ax.set_xlabel(f"{protein_name} Residue Number")
                else:
                    ax.set_xlabel('Residue Number')

            

    # Hide any unused axes (if any)
    for i in range(nrows):
        for j in range(ncols):
            if i >= len(temps) or j >= len(fields):
                axes[i, j].axis('off')

    # Custom legend for salts (only once, outside loop)
    salt_patches = [Patch(color=salt_to_color[salt], label=f"{salt}") for salt in salts]
    fig.legend(handles=salt_patches, title='Salt (mM)', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Plot NMR relaxation rates from filtered output files'
    )
    parser.add_argument(
        'directory', help='Folder containing *_results_filtered.out files'
    )
    parser.add_argument(
        '-exp', default='T1', help='Experiment type to plot, e.g. T1, T1rho'
    )
    parser.add_argument(
        '-style', default='line', help='Plot style: "bar" or "line" (default: "line")'
    )
    args = parser.parse_args()

    data = collect_data(args.directory)
    protein_name = get_protein_name(args.directory)
    plot_relaxation(data, experiment=args.exp, protein_name=protein_name, plot_type=args.style)
