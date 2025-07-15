import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Patch

# Pattern to parse filenames: e.g. 50mM_5C_600MHz_T1_results_filtered.out, 0.5uM_5C_600MHz_T1rho_extra_results_filtered.out, or _5C_600MHz_T1_results_filtered.out
FNAME_RE = re.compile(
    r"(?P<conc>[^_]+)?_?(?P<temp>\d+)C_(?P<field>\d+)MHz_(?P<exp>[^_]+)(?:_(?P<extra>.+))?_results_filtered\.out"
)

def filefield_to_concentration(s):
    """
    Converts a string representation of concentration (e.g. '50mM', '0.5uM') to a float in Molar.
    Handles various units: mM, uM, nM, M.
    """
    if s is None or s == '':
        return 0.0
    match = re.match(r"([\d\.]+)([a-zA-Z]*)", s)
    if match:
        value, unit = match.groups()
        try:
            value = float(value)
        except ValueError:
            value = 0.0
        if unit == 'mM':
            return value / 1e3
        elif unit == 'uM' or unit == 'µM':
            return value / 1e6
        elif unit == 'nM':
            return value / 1e9
        elif unit == 'M':
            return value
    return 0.0

def collect_data(directory):
    """
    Walks through `directory`, finds all `*_results_filtered.out` files,
    parses metadata from filename, and loads the data into a single DataFrame.
    """
    rows = []

    print(f"Collecting data from directory: {directory}")
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

        # Concentration (could be any string, or missing)
        conc = info.get('conc')
        if conc is None or conc == '':
            conc = '0mM'
        df['concentration'] = filefield_to_concentration(conc)
        df['temp_C'] = int(info['temp'])
        df['field_MHz'] = int(info['field'])
        df['experiment'] = info['exp']
        # Store extra field for later labeling
        df['extra_label'] = info.get('extra') if info.get('extra') else ''
        rows.append(df)

    data = pd.concat(rows, ignore_index=True)

    print(f"Collected data: {len(rows)}")
    print()
    print()

    return data


def get_protein_name(directory):
    """Extracts the protein name from the directory path."""
    return os.path.basename(os.path.abspath(directory))

def plot_relaxation(data, experiment='T1', protein_name=None, plot_type='bar'):

    """
    Unified entry point for plotting NMR relaxation data.
    Determines panel/series keys based on experiment selection and delegates to grid plotting helper.
    """
    if experiment == 'all':
        # Plot all experiment types as series, panels by temp/field
        panel_keys = ['temp_C', 'field_MHz']
        series_keys = ['experiment', 'concentration', 'extra_label']
        plot_grid_panels(data, panel_keys, series_keys, protein_name, plot_type)
    else:
        # Plot by variable/field/temperature for a single experiment
        subset = data[data['experiment'] == experiment]
        panel_keys = ['temp_C', 'field_MHz']
        # If only one concentration, color by extra_label
        if len(subset['concentration'].unique()) == 1:
            series_keys = ['extra_label']
        else:
            series_keys = ['concentration', 'extra_label']
        plot_grid_panels(subset, panel_keys, series_keys, protein_name, plot_type, experiment=experiment)


def plot_grid_panels(data, panel_keys, series_keys, protein_name, plot_type, experiment=None):
    """
    Helper function to plot grid panels for NMR relaxation data.
    panel_keys: list of column names to use for grid panels (e.g. ['temp_C', 'field_MHz'])
    series_keys: list of column names to use for series within each panel (e.g. ['experiment', 'concentration', 'extra_label'])
    """
    # Get unique panels
    panels = data.groupby(panel_keys)
    panel_vals = sorted(data[panel_keys[0]].unique())
    panel_vals2 = sorted(data[panel_keys[1]].unique())
    nrows = len(panel_vals)
    ncols = len(panel_vals2)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), sharex=True, sharey=True, num=experiment if experiment else 'all')
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]

    # Determine coloring
    # If coloring by experiment, use tab10; else use viridis for concentration
    color_keys = []
    if 'experiment' in series_keys:
        color_keys = sorted(data['experiment'].unique())
        cmap = plt.get_cmap('tab10')
        key_to_color = {key: cmap(i % 10) for i, key in enumerate(color_keys)}
    elif 'concentration' in series_keys:
        color_keys = sorted(data['concentration'].unique())
        cmap = plt.get_cmap('viridis')
        norm = mpl.colors.Normalize(vmin=min(color_keys), vmax=max(color_keys))
        key_to_color = {var: cmap(norm(var)) for var in color_keys}
    elif 'extra_label' in series_keys:
        color_keys = sorted(data['extra_label'].unique())
        cmap = plt.get_cmap('tab10')
        key_to_color = {key: cmap(i % 10) for i, key in enumerate(color_keys)}
    else:
        key_to_color = {}

    legend_handles = []
    legend_labels = []

    for i, val1 in enumerate(panel_vals):
        for j, val2 in enumerate(panel_vals2):
            ax = axes[i, j]
            panel_df = data[(data[panel_keys[0]] == val1) & (data[panel_keys[1]] == val2)]
            # Group by series_keys
            for keys, df_series in panel_df.groupby(series_keys):
                # Build label and color
                if isinstance(keys, tuple):
                    label = ' '.join([str(k) for k in keys if k != '' and k is not None])
                    # Color by first key in series_keys
                    color_key = keys[0]
                else:
                    label = str(keys)
                    color_key = keys
                color = key_to_color.get(color_key, 'gray')
                if not df_series.empty:
                    if plot_type == 'bar':
                        ax.bar(
                            df_series['residue'], df_series['rate'], yerr=df_series.get('error'),
                            color=color, label=label
                        )
                    elif plot_type == 'line':
                        df_series_sorted = df_series.sort_values('residue')
                        ax.errorbar(
                            df_series_sorted['residue'], df_series_sorted['rate'], yerr=df_series_sorted.get('error'),
                            color=color, label=label, marker='o', markersize=1, linestyle='-'
                        )
                if label not in legend_labels:
                    legend_handles.append(Patch(color=color, label=label))
                    legend_labels.append(label)
            # Panel titles and axis labels
            if i == 0:
                ax.set_title(f"{val2} MHz")
            if j == 0:
                ax.set_ylabel(f"{val1} °C" + (f"\n{experiment} (s$^{{-1}}$)" if experiment else ""))
            if i == len(panel_vals) - 1:
                if protein_name:
                    ax.set_xlabel(f"{protein_name} Residue Number")
                else:
                    ax.set_xlabel('Residue Number')

    # Hide unused axes
    for i in range(nrows):
        for j in range(ncols):
            if i >= len(panel_vals) or j >= len(panel_vals2):
                axes[i, j].axis('off')

    # Legend in top-right
    ax = axes[0, -1]
    ax.legend(handles=legend_handles)
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
        '-exp', default='T1', help='Experiment type to plot, e.g. T1, T1rho, or "all" for all experiments'
    )
    parser.add_argument(
        '-style', default='line', help='Plot style: "bar" or "line" (default: "line")'
    )
    parser.add_argument(
        '-protein', help='Observed protein name (optional, derived from directory name)'
    )
    args = parser.parse_args()

    data = collect_data(args.directory)
    protein_name = get_protein_name(args.directory)
    if args.protein:
        protein_name = args.protein
    plot_relaxation(data, experiment=args.exp, protein_name=protein_name, plot_type=args.style)
