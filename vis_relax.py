import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Patch

# Options structure for carrying plot settings
from dataclasses import dataclass

@dataclass
class PlotOptions:
    experiment: str = None
    protein_name: str = None
    plot_type: str = 'bar'
    tempunit: str = 'C'
    tempunit_string: str = '°C'
    concentration_unit: str = 'mM'  # Default unit for concentration
    # Add more fields as needed

    def get_conc_factor(self):
        if self.concentration_unit == 'mM': return 1e3
        elif self.concentration_unit in ['uM', 'µM']: return 1e6
        elif self.concentration_unit == 'nM': return 1e9
        elif self.concentration_unit == 'pM': return 1e12
        else: return 1

    def get_string_conc(self, conc):
        conc = self.get_conc_factor() * conc
        return f"{conc:.0f} {self.concentration_unit}"

# Pattern to parse filenames: e.g. 50mM_5C_600MHz_T1_results_filtered.out, 0.5uM_5C_600MHz_T1rho_extra_results_filtered.out, or _5C_600MHz_T1_results_filtered.out
FNAME_RE = re.compile(
    r"(?P<conc>\d+(?:\.\d+)?(?:mM|uM|µM|nM|M))?.*?_?(?P<temp>\d+[CK])_(?P<field>\d+)MHz_(?P<exp>[^_]+)(?:_(?P<extra>.+))?_results_filtered\.out"
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
        print(f"Parsing concentration: value={value}, unit={unit}")
        try:
            value = float(value)
        except ValueError:
            value = 0.0
        if 'mM' in unit:
            return value / 1e3
        elif 'uM' in unit or 'µM' in unit:
            return value / 1e6
        elif 'nM' in unit:
            return value / 1e9
        elif 'M' in unit:
            return value
    return 0.0

def filefield_to_temperature(s, output_unit='C'):
    """
    Converts a string representation of temperature to degrees Celsius.
    Handles various units: Celcius and Kelvin.
    """
    if s is None or s == '':
        return 0.0
    match = re.match(r"([\d\.]+)([CK]*)", s)
    if match:
        value, unit = match.groups()
        try:
            value = float(value)
        except ValueError:
            value = 0.0
        if unit == 'C': # Celsius
            if output_unit == 'C':
                return value
            elif output_unit == 'K': # Convert Celsius to Kelvin
                return value + 273.15
        elif unit == 'K': # Convert Kelvin to Celsius
            if output_unit == 'C':
                return value - 273.15
            elif output_unit == 'K':
                return value
    # Default case, if no match or conversion fails
    if output_unit == 'C':
        return 25 # Default to 25C
    elif output_unit == 'K':
        return 298.15

def get_best_common_conc_unit(concs):
    average = sum(concs) / len(concs)
    magnitude = 10 ** (np.floor(np.log10(average)))
    if magnitude >= 1e-3 and magnitude < 1:
        return 'mM'
    elif magnitude >= 1e-6 and magnitude < 1e-3:
        return 'uM'
    elif magnitude >= 1e-9 and magnitude < 1e-6:
        return 'nM'
    elif magnitude >= 1e-12 and magnitude < 1e-9:
        return 'pM'
    else:
        return 'M'

def collect_data(directory, options=None):
    """
    Walks through `directory`, finds all `*_results_filtered.out` files,
    parses metadata from filename, and loads the data into a single DataFrame.
    """
    rows = []
    variable_concentrations = []

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
        print(conc)
        if conc is None or conc == '':
            conc = '0mM'
        df['concentration'] = filefield_to_concentration(conc)
        df['temp'] = filefield_to_temperature(info['temp'], output_unit=options.tempunit)
        df['field_MHz'] = int(info['field'])
        df['experiment'] = info['exp']
        # Store extra field for later labeling
        df['extra_label'] = info.get('extra') if info.get('extra') else ''
        rows.append(df)

        variable_concentrations.append(filefield_to_concentration(conc))

    data = pd.concat(rows, ignore_index=True)
    common_conc_unit = get_best_common_conc_unit(variable_concentrations)

    print(f"Collected data: {len(rows)}")
    print(f"Common concentration: {common_conc_unit}")
    print()
    print()

    return data, common_conc_unit

def get_protein_name(directory):
    """Extracts the protein name from the directory path."""
    return os.path.basename(os.path.abspath(directory))

def plot_relaxation(data, options: PlotOptions):
    print(f"plot_relaxation: experiment={options.experiment}, protein_name={options.protein_name}, plot_type={options.plot_type}, tempunit={options.tempunit}")
    print(f"Data shape: {data.shape}, columns: {list(data.columns)}")
    """
    Unified entry point for plotting NMR relaxation data.
    Determines panel/series keys based on experiment selection and delegates to grid plotting helper.
    """
    if options.experiment == None:
        # Plot all experiment types as separate panels (rows), with temp/field as columns
        panel_keys = ['experiment', 'temp', 'field_MHz']
        if len(data['concentration'].unique()) == 1:
            series_keys = ['extra_label']
        else:
            series_keys = ['concentration', 'extra_label']
        plot_grid_panels(data, panel_keys, series_keys, options)
    else:
        # Plot by variable/field/temperature for a single experiment
        subset = data[data['experiment'] == options.experiment]
        panel_keys = ['temp', 'field_MHz']
        # If only one concentration, color by extra_label
        if len(subset['concentration'].unique()) == 1:
            series_keys = ['extra_label']
        else:
            series_keys = ['concentration', 'extra_label']
        plot_grid_panels(subset, panel_keys, series_keys, options)

def plot_grid_panels(data, panel_keys, series_keys, options: PlotOptions):
    """
    Helper function to plot grid panels for NMR relaxation data.
    panel_keys: list of column names to use for grid panels (e.g. ['temp', 'field_MHz'])
    series_keys: list of column names to use for series within each panel (e.g. ['experiment', 'concentration', 'extra_label'])
    """
    print(f"plot_grid_panels: panel_keys={panel_keys}, series_keys={series_keys}, protein_name={options.protein_name}, plot_type={options.plot_type}, experiment={options.experiment}")
    print(f"Unique values for panel_keys:")
    fig_title = options.protein_name or ''
    fig_title += "_" + options.experiment if options.experiment else ''
    
    for k in panel_keys:
        print(f"  {k}: {sorted(data[k].unique())}")
        if len(data[k].unique()) == 1:
            if k == 'temp':
                fig_title += f"_{data[k].unique()[0]}{options.tempunit}"
            elif k == 'field_MHz':
                fig_title += f"_{data[k].unique()[0]}MHz"

    # Get unique panels
    panels = data.groupby(panel_keys)
    panel_vals = [sorted(data[k].unique()) for k in panel_keys]
    nrows = len(panel_vals[0])
    ncols = len(panel_vals[1]) if len(panel_keys) > 1 else 1
    fig, axes = plt.subplots(
        nrows, 
        ncols, 
        figsize=(4*ncols, 3*nrows), 
        sharex=False, 
        sharey=False, 
        num=fig_title)
    
    # Ensure axes is always 2D
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]

    axes_ranges = [None] * nrows

    for i, val0 in enumerate(panel_vals[0]):
        print(f"Processing panel row {i+1}/{nrows} for {panel_keys[0]}={val0}")

        if panel_keys[0] == 'experiment':
            df = data[data[panel_keys[0]] == val0]
            delta = df['rate'].max() - df['rate'].min()
            axes_ranges[i] = round(10*(df['rate'].min() - 0.1*delta))/10, round(10*(df['rate'].max() + 0.1*delta))/10

    print(f'Axes ranges: {axes_ranges}')

    # Determine coloring: combine all series_keys into a tuple key for each series
    print(f"Coloring method: based on series_keys: {series_keys}")
    # Get all unique combinations of series_keys present in the data
    unique_keys = sorted(
        set(tuple(row[k] for k in series_keys) for _, row in data.iterrows())
    )
    print(f"Color keys (series_keys tuples): {unique_keys}")
    # Use tab10 or tab20 colormap depending on number of series
    if len(unique_keys) <= 10:
        cmap = plt.get_cmap('tab10')
    elif len(unique_keys) <= 20:
        cmap = plt.get_cmap('tab20')
    else:
        cmap = plt.get_cmap('hsv')
    key_to_color = {key: cmap(i % cmap.N) for i, key in enumerate(unique_keys)}
    print(f"Key to color mapping: {key_to_color}")

    print(f"Panel grid: {[(v0, v1) for v0 in panel_vals[0] for v1 in panel_vals[1]]}")
    print(f"Panel loop: nrows={nrows}, ncols={ncols}")

    legend_handles = []
    legend_labels = []

    print()
    print("Starting panel processing...")
    for i, val0 in enumerate(panel_vals[0]):
        print(f"Processing panel row {i+1}/{nrows} for {panel_keys[0]}={val0}")
        for j in range(ncols):
            print(f"Processing column {j+1}/{ncols} for {panel_keys[1]}={panel_vals[1][j] if ncols > 1 else 'N/A'}")
            val1 = panel_vals[1][j] if ncols > 1 else None
            print(f"Panel [{i},{j}]: {panel_keys[0]}={val0}, {panel_keys[1] if ncols > 1 else ''}={val1}")
            ax = axes[i, j]
            # Build panel filter
            panel_filter = (data[panel_keys[0]] == val0)
            if ncols > 1:
                panel_filter &= (data[panel_keys[1]] == val1)
            if len(panel_keys) > 2:
                val2 = None
                if len(panel_vals) > 2:
                    val2 = panel_vals[2][0] if len(panel_vals[2]) == 1 else None
                if val2 is not None:
                    print(f"  Filtering for {panel_keys[2]}={val2}")
                    panel_filter &= (data[panel_keys[2]] == val2)
            panel_df = data[panel_filter]
            print(f"  Panel data shape: {panel_df.shape}")
            # Group by series_keys
            for keys, df_series in panel_df.groupby(series_keys):
                print(f"      n={len(df_series)}")
                # Build label and color
                label = ''
                for key,val in zip(series_keys, keys):
                    print(f"      {key}: {val}")
                    if key == 'concentration': label += f"{options.get_string_conc(val)} "
                    else: label += f"{val} "

                color = key_to_color.get(keys, 'gray')
                print(f"      Color: {color}, Label: {label}")
                if not df_series.empty:
                    if options.plot_type == 'bar':
                        ax.bar(
                            df_series['residue'], df_series['rate'], yerr=df_series.get('error'),
                            color=color, label=label
                        )
                    elif options.plot_type == 'line':
                        df_series_sorted = df_series.sort_values('residue')
                        ax.errorbar(
                            df_series_sorted['residue'], df_series_sorted['rate'], yerr=df_series_sorted.get('error'),
                            color=color, label=label, marker='o', markersize=1, linestyle='-'
                        )
                if label not in legend_labels:
                    legend_handles.append(Patch(color=color, label=label))
                    legend_labels.append(label)

                if axes_ranges[i]:
                    print(f"  [ROW] Setting axes limits for row {i+1}: {axes_ranges[i][0]} - {axes_ranges[i][1]}")
                    ax.set_ylim(axes_ranges[i])
            
            # Panel titles and axis labels
            title = ''
            if i == 0: # Only set title for the first row
                if options.experiment == None:
                    # Add temperature and field with units if present
                    key = panel_keys[1]
                    if key == 'temp' and (val1 is not None):
                        title += f"{val1} {options.tempunit_string}"
                    elif key == 'field_MHz' and (val1 is not None):
                        title += f"{val1} MHz"
                else:
                    if val1 is not None:
                        # Add units for temp/field if present
                        key = panel_keys[1]
                        if key == 'temp':
                            title = f"{val1} {options.tempunit_string}"
                        elif key == 'field_MHz':
                            title = f"{val1} MHz"
                        else:
                            title = f"{val1}"
                print(f"  [COLUMN] Setting title: {title}")
            
            ax.set_title(title)

            if j == 0:
                print(f"  [ROW] Setting ylabel for row [{i},{j}]: {val0} {options.experiment}")
                if options.experiment == None:
                    # val0 should be experiment name, but if tuple, take first element
                    exp_label = val0[0] if isinstance(val0, tuple) else val0
                    ax.set_ylabel(f"{exp_label} (s$^{{-1}}$)")
                else:
                    ax.set_ylabel(f"{val0} {options.tempunit_string}" + (f"\n{options.experiment} (s$^{{-1}}$)" if options.experiment else ""))
            
            if i == nrows - 1:
                if protein_name:
                    ax.set_xlabel(f"{protein_name} Residue Number")
                else:
                    ax.set_xlabel('Residue Number')
        
        print(f"Finished processing panel row {i+1}/{nrows}")
        print()
        print()

    # Hide unused axes
    for i in range(nrows):
        for j in range(ncols):
            if i >= len(panel_vals[0]) or (ncols > 1 and j >= len(panel_vals[1])):
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
        '-exp', default=None, help='Specific experiment type to plot, e.g. T1, T1rho'
    )
    parser.add_argument(
        '-style', default='line', help='Plot style: "bar" or "line" (default: "line")'
    )
    parser.add_argument(
        '-protein', help='Observed protein name (optional, derived from directory name)'
    )
    parser.add_argument(
        '-temp', type=int, help='Temperature (Celsius) to plot (optional)'
    )
    parser.add_argument(
        '-field', type=int, help='Field (MHz) to plot (optional)'
    )
    parser.add_argument(
        '-tempunit', choices=['C', 'K'], default='C',
    )
    args = parser.parse_args()

    # Determine protein name from directory or command line argument
    protein_name = get_protein_name(args.directory)
    if args.protein:
        protein_name = args.protein

    # Build options object
    options = PlotOptions(
        experiment=args.exp,
        protein_name=protein_name,
        plot_type=args.style,
        tempunit=args.tempunit,
        tempunit_string=f"°{args.tempunit}" if args.tempunit == 'C' else args.tempunit
    )

    # Collect data from the specified directory
    data, conc_unit = collect_data(args.directory, options)

    options.concentration_unit = conc_unit

    # Filter by temperature and field if specified
    if args.temp is not None:
        data = data[data['temp'] == args.temp]
    if args.field is not None:
        data = data[data['field_MHz'] == args.field]

    # If plotting all experiments, check number of variables
    if args.exp == 'all':
        n_temps = len(data['temp'].unique())
        n_fields = len(data['field_MHz'].unique())
        n_conc = len(data['concentration'].unique())
        n_extra = len(data['extra_label'].unique())
        # Only allow if at least one variable is fixed (max 3 vary)
        n_varying = sum([n_temps > 1, n_fields > 1, n_conc > 1, n_extra > 1])
        if n_varying > 3:
            raise ValueError("Too many variables vary for 'all' experiment plotting. Please fix temperature, field, or concentration using -temp or -field.")

    plot_relaxation(data, options)
