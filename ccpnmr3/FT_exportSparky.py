"""
CCPNMR3 Macro: Export Peak List as Sparky Format Peak List
Frederik Theisen 2025

This macro exports a given peak list as a Sparky formatted text file.
The output path must be designated and writeable.

Peaks can be excluded on basis of Euclidean distance.
A pseudo 3D options exists to produce third dimension necessary for ancient relaxation data analyses.
The carbon dimension is curently not implemented.

"""

PEAKLIST = 'PL:T1rho_plane.2'


OUTPUTPATH = '/home/ftheisen/Documents/' + PEAKLIST.split(':')[1] + '.list'
INCLUDE_VOLUME = False
RECALCULATE_HEIGHTS = True
WRITE_NOISE_METADATA = True
# None estimates noise automatically across the spectrum. To select a reviewed
# background region, use exact axis codes, e.g. {'H': (9.0, 10.0), 'N': (110.0, 130.0)}.
NOISE_REGION = None
NOISE_MAX_SAMPLES = 65536
NOISE_RANDOM_SEED = 0
AUTO_FILTER_OVERLAPS = True
COMMENT_OVERLAP = True
PSEUDO3D = True
DIMORDER = {'H':2,'N':1,'PSEUDO':0, 'C':1}
SEQUENCE_CODE = 'oneLetter' #'oneLetter' or 'threeLetter'
OVERLAP_DISTANCE_THRESHOLD = 0.02
ATOM_CS_WEIGHTS = {'H':1, 'N':0.154, 'C':0.3}

from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import io
import json
import os
import tempfile
from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar, notificationEchoBlocking
import numpy as np

# Dictionary to convert three-letter residue codes to one-letter.
chemCompCodesDict = {
    'Ala': ('A', 'ALA', 'ALANINE', 'C3H7N1O2'),
    'Cys': ('C', 'CYS', 'CYSTEINE', 'C3H7N1O2S1'),
    'Asp': ('D', 'ASP', 'ASPARTIC ACID', 'C4H6N1O4'),
    'Glu': ('E', 'GLU', 'GLUTAMIC ACID', 'C5H8N1O4'),
    'Phe': ('F', 'PHE', 'PHENYLALANINE', 'C9H11N1O2'),
    'Gly': ('G', 'GLY', 'GLYCINE', 'C2H5N1O2'),
    'His': ('H', 'HIS', 'L-Histidine', 'C6H10N3O2'),
    'Ile': ('I', 'ILE', 'ISOLEUCINE', 'C6H13N1O2'),
    'Lys': ('K', 'LYS', 'LYSINE', 'C6H15N2O2'),
    'Leu': ('L', 'LEU', 'LEUCINE', 'C6H13N1O2'),
    'Met': ('M', 'MET', 'METHIONINE', 'C5H11N1O2S1'),
    'Asn': ('N', 'ASN', 'ASPARAGINE', 'C4H8N2O3'),
    'Pro': ('P', 'PRO', 'PROLINE', 'C5H9N1O2'),
    'Gln': ('Q', 'GLN', 'GLUTAMINE', 'C5H10N2O3'),
    'Arg': ('R', 'ARG', 'ARGININE', 'C6H15N4O2'),
    'Ser': ('S', 'SER', 'SERINE', 'C3H7N1O3'),
    'Thr': ('T', 'THR', 'THREONINE', 'C4H9N1O3'),
    'Val': ('V', 'VAL', 'VALINE', 'C5H11N1O2'),
    'Trp': ('W', 'TRP', 'TRYPTOPHAN', 'C11H12N2O2'),
    'Tyr': ('Y', 'TYR', 'TYROSINE', 'C9H11N1O3'),
}

def get_sparky_header(peakList):
    header = 'assignment'
    dims = len(peakList.peaks[0].assignmentsByDimensions)

    if PSEUDO3D: dims += 1

    for idx in range(dims):
        dim = '\tw' + str(1+idx)
        header += dim

    header += '\tintensity'

    if INCLUDE_VOLUME: 
        header += '\tvolume'

    if AUTO_FILTER_OVERLAPS:
        header += '\toverlap'
        if COMMENT_OVERLAP:
            header += '\tcomment'

    return header, dims

def sort_by_key(data):
    '''
        This function sorts the assignment based on the orderkey defined further up. This order information is trown away afterwards.
    '''
    return [[y,z,k] for x, y, z, k in sorted(data, key=lambda item: item[0])]

def convertResidueCode(residueName, inputCodeType='threeLetter', outputCodeType='oneLetter'):
    modes = ['oneLetter', 'threeLetter', 'synonym', 'molFormula']
    if inputCodeType not in modes or outputCodeType not in modes:
        return residueName
    for key, values in chemCompCodesDict.items():
        dd = dict(zip(modes, values))
        if residueName == dd.get(inputCodeType):
            return dd.get(outputCodeType)
    return residueName

def flag_overlapping_peaks(peaklist):
    num_peaks = len(peaklist)

    # Initialize overlap flag
    for peak in peaklist:
        peak['overlap'] = False
        peak['comment'] += "|" # Divider between user comment and overlap auto comment

    # Compute distances and flag overlaps
    for i in range(num_peaks):
        for j in range(num_peaks):  # To avoid redundant calculations, set start to i + 1
            if i == j: continue # don't check against self
            peak_i = peaklist[i]
            peak_j = peaklist[j]

            # Calculate weighted Euclidean distance
            distance_squared = 0
            for (atom_i, cs_i), (atom_j, cs_j) in zip(peak_i['assign'], peak_j['assign']):
                if atom_i != atom_j or atom_i == '?':
                    continue  # Ensure same atom type and that we are not looking at a pseudo atom
                weight = ATOM_CS_WEIGHTS[atom_i]
                distance_squared += (weight * (cs_i - cs_j)) ** 2

            distance = np.sqrt(distance_squared)

            if distance < OVERLAP_DISTANCE_THRESHOLD:
                peaklist[i]['overlap'] = True
                peaklist[j]['overlap'] = True  # Mark both as overlapping (unnecessary)

                peaklist[i]['comment'] += str(peaklist[j]['residuenumber']) + ','

                print(f"Peak overlap detected: {peaklist[i]['label']} {peaklist[j]['label']} {round(distance,5)}")

        peaklist[i]['comment'] = peaklist[i]['comment'].strip('|').strip(',') # Remove delimiters if not used for anything

    return peaklist  # Returns modified list with overlap flags

def recalculate_heights_at_positions(peakList):
    """Refresh signed, interpolated heights without fitting or moving peaks."""
    spectrum = peakList.spectrum
    if not spectrum.hasValidPath():
        raise RuntimeError(f"Cannot recalculate heights: {spectrum.pid} has no valid data path")

    # Read every value before updating peaks, so a read failure leaves them alone.
    heights = []
    for peak in peakList.peaks:
        try:
            height = spectrum.getHeight(peak.position)
            if height is None or not np.isfinite(height):
                raise ValueError("height is missing or non-finite")
            heights.append((peak, float(height)))
        except Exception as error:
            raise RuntimeError(f"Cannot recalculate height for {peak.pid}: {error}") from error

    with undoBlockWithoutSideBar(), notificationEchoBlocking():
        for peak, height in heights:
            peak.height = height

def estimate_spectrum_noise(spectrum):
    """Estimate Gaussian-equivalent background SD from signed grid intensities.

    Median absolute deviation limits contamination from sparse peaks of either
    sign. This is a global background estimate, not a fitted-peak uncertainty.
    """
    if not spectrum.hasValidPath():
        raise RuntimeError(f"Cannot estimate noise: {spectrum.pid} has no valid data path")
    if NOISE_MAX_SAMPLES < 128:
        raise ValueError("NOISE_MAX_SAMPLES must be at least 128")
    rng = np.random.default_rng(NOISE_RANDOM_SEED)
    point_counts = tuple(int(n) for n in spectrum.pointCounts)
    total_points = int(np.prod(point_counts))
    region = dict(zip(spectrum.axisCodes, spectrum.spectrumLimits))

    if NOISE_REGION is not None:
        if not NOISE_REGION or set(NOISE_REGION) - set(spectrum.axisCodes):
            raise ValueError("NOISE_REGION must use this spectrum's exact axis codes")
        for axis, bounds in NOISE_REGION.items():
            limits = sorted(region[axis])
            if (len(bounds) != 2 or not np.all(np.isfinite(bounds)) or
                    not limits[0] <= bounds[0] < bounds[1] <= limits[1]):
                raise ValueError(f"NOISE_REGION bounds must be inside the spectrum: {axis}")
        region.update(NOISE_REGION)
        # CCPN 3.3.4 getRegion returns data-source units; apply Spectrum.scale.
        values = np.asarray(spectrum.getRegion(**region), dtype=np.float64).ravel()
        total_points = values.size
        selection = "configured_ppm_region"
    elif spectrum.dimensionCount == 1:
        values = np.asarray(spectrum.dataSource.getSliceData(), dtype=np.float64).ravel()
        selection = "whole_spectrum"
    elif spectrum.dimensionCount == 2:
        values = np.asarray(spectrum.dataSource.getPlaneData(), dtype=np.float64).ravel()
        selection = "whole_spectrum"
    else:
        # Bound data access for large nD spectra instead of loading an entire cube.
        indices = rng.choice(total_points, size=min(total_points, 10000, NOISE_MAX_SAMPLES),
                             replace=False)
        positions = np.column_stack(np.unravel_index(indices, point_counts)) + 1
        values = np.array([spectrum.dataSource.getPointData(position=tuple(map(int, p)))
                           for p in positions], dtype=np.float64)
        selection = "whole_spectrum_random_grid_points"

    if values.size > NOISE_MAX_SAMPLES:
        values = values[rng.choice(values.size, size=NOISE_MAX_SAMPLES, replace=False)]
    # All data-source reads above already include their format's dataScale.
    # Spectrum.scale brings noise into the same units as exported peak heights.
    values = values * float(spectrum.scale)
    nonfinite = int(np.count_nonzero(~np.isfinite(values)))
    values = values[np.isfinite(values)]
    if values.size < 128:
        raise RuntimeError("Too few finite data points for a spectrum noise estimate")
    baseline = float(np.median(values))
    mad = float(np.median(np.abs(values - baseline)))
    sigma = mad / 0.6744897501960817
    if not np.isfinite(sigma) or sigma <= 0:
        raise RuntimeError("Spectrum noise estimate is zero or invalid; review the noise region")
    return {
        "status": "estimated", "noise_sigma": sigma,
        "method": "median_absolute_deviation_gaussian_scaled",
        "definition": "median(abs(x - median(x))) / 0.6744897501960817",
        "baseline_median": baseline, "mad": mad,
        "units": "same scaled intensity units as the exported peak heights",
        "selection": selection, "region_ppm": region,
        "region_user_selected": NOISE_REGION is not None,
        "available_grid_points": total_points, "finite_sample_count": int(values.size),
        "discarded_nonfinite_samples": nonfinite, "random_seed": NOISE_RANDOM_SEED,
        "limitations": ["Automatic selection assumes most sampled points are background.",
                        "Grid noise is a proxy for interpolated peak-height uncertainty; processing correlates neighbouring points.",
                        "Overlap, baseline artefacts, peak placement and systematic errors are not quantified.",
                        "For nD data this is pooled spectrum noise, not a separate estimate per plane."],
    }

def write_peak_export(peakList, noise, table_text):
    """Write the peak table and a commented JSON footer as one atomic file."""
    output = Path(OUTPUTPATH).expanduser().absolute()
    content = table_text
    if noise is not None:
        spectrum = peakList.spectrum
        metadata = {
            "schema": "ccpn_peak_export_metadata", "schema_version": 2,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "export_file": output.name,
            "peak_table_sha256": hashlib.sha256(table_text.encode('utf-8')).hexdigest(),
            "peak_list_pid": peakList.pid, "peak_count": len(peakList.peaks),
            "spectrum_pid": spectrum.pid, "spectrum_path": str(spectrum.filePath),
            "spectrum_scale": float(spectrum.scale),
            "axis_codes": list(spectrum.axisCodes), "point_counts": list(spectrum.pointCounts),
            "heights_recalculated_at_positions": bool(RECALCULATE_HEIGHTS),
            "height_sampling": "CCPN spectrum.getHeight(peak.position)" if RECALCULATE_HEIGHTS else "stored peak.height",
            "noise": noise,
        }
        comment = '\n'.join('# ' + line for line in json.dumps(metadata, indent=2, allow_nan=False).splitlines())
        content += '# CCPN_METADATA_BEGIN\n' + comment + '\n# CCPN_METADATA_END\n'
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=output.parent,
                                         prefix='.' + output.name + '.', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        if output.exists():
            os.chmod(temporary, output.stat().st_mode & 0o777)
        else:
            os.chmod(temporary, 0o644)
        os.replace(temporary, output)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    if noise is not None:
        print(f"Noise SD: {noise['noise_sigma']:.6g}; metadata embedded in: {output}")

def runMacro():
    # Retrieve the currently active peak list.
    # (getCurrentPeakList() is provided by the CCPNMR3 macro environment.)
    peakList = get(PEAKLIST)
    if not peakList:
        print("No active peak list found!")
        return

    if not peakList.peaks:
        print(f"No peaks found in {PEAKLIST}; nothing exported.")
        return
    noise = estimate_spectrum_noise(peakList.spectrum) if WRITE_NOISE_METADATA else None
    if RECALCULATE_HEIGHTS:
        recalculate_heights_at_positions(peakList)

    # Get and write the header line in Sparky format.
    sparky_header,dims = get_sparky_header(peakList)
    try:
        sortedPeakList = sorted(peakList.peaks, key=lambda item: item.assignments[0][0].nmrResidue.sequenceCode)
    except:
        sortedPeakList = peakList.peaks

    data=[]

    # Loop over each peak in the current peak list.
    for pk in sortedPeakList:
        ass = [] # [[priority, assignment, cs]]
        dat = {}
        dat['comment'] = pk.comment.replace(' ','_')

        if PSEUDO3D: ass.append([DIMORDER['PSEUDO'],'?','?',10]) # First dim will be pseudo 10ppm, 

        # For each dimension, take the first available assignment.
        for dimIdx, assignments in enumerate(pk.assignmentsByDimensions):
            cs = pk.position[dimIdx] # Get CS for peak and dimension

            if len(assignments) == 0 or '@' in assignments[0].nmrResidue.sequenceCode:
                pkNum = int(pk.serial) + 10000 #.split('.')[-1]
                atomGuess = 'H' if dimIdx == 0 else 'N'
                ass.append([DIMORDER[atomGuess],f"X{pkNum}{atomGuess}", atomGuess, cs])
                dat['residuenumber'] = pkNum
            else:
                residue = assignments[0] # Use the first assignment for this dimension (not sure how there could be more)
                residueNumber = assignments[0].nmrResidue.sequenceCode
                # Convert the residue type from three-letter. Get residue number. Get atom type.
                resType = convertResidueCode(residue.nmrResidue.residueType, inputCodeType='threeLetter', outputCodeType=SEQUENCE_CODE)
                atomName = residue.name.replace('%','')
                
                # Save assignment
                ass.append([DIMORDER[atomName],f"{resType}{residueNumber}{atomName}", atomName, cs])
                print(residueNumber, resType, atomName)
                dat['residuenumber'] = int(residueNumber)
                

        # Sort the assignment according to priority
        sorted_ass = sort_by_key(ass)

        # Construct the assignment label X-Y-Z
        label = '-'.join([ass[0] for ass in sorted_ass])

        dat['label'] = label
        dat['assign'] = [[atm,cs] for lbl, atm, cs in sorted_ass]

        # Get the intensity
        dat['height'] = pk.height

        # If selected, export the volume
        if INCLUDE_VOLUME:
            if pk.volume != None: dat['volume'] = pk.volume
            else: dat['volume'] = None

        data.append(dat)
        
    sorted_data = sorted(data, key=lambda item: item['residuenumber'])
    
    with io.StringIO() as f:
        f.write(sparky_header + '\n')

        data = flag_overlapping_peaks(sorted_data)

        for peak in data:
            f.write(f"{peak['label']}")

            for atm,cs in peak['assign']:
                f.write(f"\t{round(cs,4)}")
            
            f.write(f"\t{peak['height']}")
            
            if INCLUDE_VOLUME: f.write(f"\t{peak['volume']}")
            
            if AUTO_FILTER_OVERLAPS:
                f.write(f"\t{peak['overlap']}")
                if COMMENT_OVERLAP: f.write(f"\t{peak['comment']}")

            f.write(f'\n')

        table_text = f.getvalue()

    write_peak_export(peakList, noise, table_text)

    print("Sparky peak list exported successfully to:", OUTPUTPATH)

if __name__ == '__main__':
    runMacro()
