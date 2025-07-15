import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt

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


def plot_relaxation(data, experiment='T1'):
    """
    Creates a multi-panel figure of relaxation rates vs field strength,
    with one panel per temperature, and lines colored by salt concentration.
    """
    subset = data[data['experiment'] == experiment]
    temps = sorted(subset['temp_C'].unique())
    salts = sorted(subset['salt'].unique(), key=lambda x: int(x.replace('mM','')))

    n = len(temps)
    cols = min(3, n)
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(4*cols,3*rows), sharey=True)
    axes = axes.flatten()

    for ax, temp in zip(axes, temps):
        df_t = subset[subset['temp_C'] == temp]
        for salt in salts:
            df_ts = df_t[df_t['salt'] == salt]
            ax.bar(
                df_ts['residue'], df_ts['rate'], yerr=df_ts.get('error'), label=salt
            )
        ax.set_title(f"{temp} °C")
        ax.set_xlabel('Field (MHz)')
    axes[0].set_ylabel('Relaxation rate (s^-1)')
    # hide any unused axes
    for ax in axes[len(temps):]:
        ax.axis('off')

    plt.legend(title='Salt (mM)', bbox_to_anchor=(1.05, 1), loc='upper left')
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
        '--exp', default='T1', help='Experiment type to plot, e.g. T1, T1rho'
    )
    args = parser.parse_args()

    data = collect_data(args.directory)
    plot_relaxation(data, experiment=args.exp)
