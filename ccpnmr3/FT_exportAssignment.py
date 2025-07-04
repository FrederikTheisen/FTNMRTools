#!/usr/bin/env python
"""
CCPNMR3 Macro: Export an assignment based on a collection of assignment spectra
Frederik Theisen 2025

This script groups chemical shifts by residue in an nmrChain. For each residue (identified
by chain and sequence code) an assignment label is generated.
For each atom type (e.g. H, N, CA, CB, etc.) present in the chemical shift list,
three possible columns are created: <Atom>_shift, <Atom>_error, <Atom>_peaks
shift is the chemical shift
error is the standard deviation of the chemical shift if the resonance has been assigned in multiple spectra
peaks is the number of peak which are included in the given chemical shift

Peak count information can be reported for each shift or for the residue as a 
whole as minimum number of assignments for the given shift

If a residue does not have data for a particular atom, those fields are left blank.

Usage (within CCPNMR3 macro environment):
    - Edit macro to select chemical shift list
    - Set valid output path
    - Change output parameters if desired
"""

### SETTING ###

NMRCHAIN = 'NC:heTDR1'
ASSIGNMENT_SPECTRA_COLLECTION = 'CO:AssignmentSpectra'
PERFORM_TROSY_CORRECTION = True
CHEMICALSHIFT_WEIGHTING = True

OUTPUTPATH = '/home/ftheisen/Documents/' + NMRCHAIN.split(':')[1]
OUTPUT_FORMAT = "TALOS" #CS or TALOS

# Settings for non-talos output
THREELETTER = True
INCLUDE_ERRORS = True
INCLUDE_PEAK_COUNT = True
INCLUDE_CS_SPECIFIC_PEAK_COUNT = False
DELIMITER = '\t'

### SCRIPT ###

DEBUG = True

# Adjust settings
if not INCLUDE_PEAK_COUNT: SHOW_CS_SPECIFIC_PEAK_COUNT = False

if OUTPUT_FORMAT == 'TALOS': OUTPUTPATH += '.tab'
else: OUTPUTPATH += '.list'

# Fetch spectra from collection
SPECTRA = [item.name for item in project.getByPid(ASSIGNMENT_SPECTRA_COLLECTION).items if 'SP' in item.pid]

from collections import defaultdict
import csv
import math

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

ccp2talos = {
    'CA':'CA',
    'CB':'CB',
    'N':'N',
    'H':'HN',
    'C':'C'
}

def debug_print(s1, s2 = "", s3 = "", s4 = ""):
    if DEBUG:
        print("DEBUG: ",s1,s2,s3,s4)

def convertResidueCode(residueName, inputCodeType='threeLetter', outputCodeType='oneLetter'):
    """
    Convert a residue name/code from one format to another.
    For example, converting "Asn" (three-letter) to "N" (one-letter).
    """
    modes = ['oneLetter', 'threeLetter', 'synonym', 'molFormula']
    if inputCodeType not in modes or outputCodeType not in modes:
        return residueName
    for key, values in chemCompCodesDict.items():
        dd = dict(zip(modes, values))
        if residueName == dd.get(inputCodeType):
            return dd.get(outputCodeType)
    return residueName

def calcChemicalShift(peaks, weigting = True):
    cs = 0
    cs_n = 0
    cs_e = 0

    shifts = [f'{p["val"]:.2f}+{p["err"]:.2f}' for p in peaks]
    debug_print(' '.join(shifts))

    if not weigting:
        for peak in peaks:
            cs += peak['val']
            cs_n += 1
        cs = cs / cs_n

        ss = 0
        for peak in peaks:
            ss += (peak['val'] - cs)**2
        cs_e = math.sqrt(ss) / cs_n

    else:
        w = 0
        for peak in peaks:
            cs += peak['val'] / peak['err']
            w += 1 / peak['err']
            cs_n += 1
        cs = cs / w
        cs_e = w / cs_n

    return cs, cs_n, cs_e

def correctTrosyShift(observed_shift, atom_name, field_strength_mhz):
    """
    Corrects chemical shifts from TROSY-selected peaks to match the multiplet center.
    
    Args:
        atom_name (str): Atom name (e.g., 'H', 'HN', 'N').
        observed_shift (float): The shift in ppm from the TROSY spectrum.
        field_strength_mhz (float): The proton frequency in MHz (e.g., 600, 800).

    Returns:
        float: Corrected shift in ppm.
    """
    j_coupling_hz = 95.0
    half_j = j_coupling_hz / 2.0
    ppm_offset = half_j / field_strength_mhz

    if atom_name == 'H':
        return observed_shift + ppm_offset  # TROSY H is lower
    elif atom_name == 'N':
        return observed_shift - ppm_offset  # TROSY N is higher
    else:
        return observed_shift  # no correction for other nuclei

def getChemicalShifts():
    """
    Export a CCPNMR3 chemical shift list to a file.
    
    Groups chemical shifts by residue (using chain and sequence code) and then,
    for each residue, outputs:
        - Assignment: one-letter code + sequence code,
        - For each atom type present in the overall csList:
              * <Atom>_shift, <Atom>_error, <Atom>_nPeaks
  
    If a residue does not have data for a given atom type, blank entries are written.
    """
    # Dictionary for grouping by residue.
    # Key: string "ChainSequence" (e.g. "A53")
    # Value: {'residue': residueObject, 'atoms': {atomType: (value, error, peaks)}}
    data = {}
    # Set of all atom types encountered.
    atomTypes = set()

    chain = project.getByPid(NMRCHAIN)

    for res in chain.nmrResidues:

        seqCode = res.sequenceCode
        res_code = int(seqCode)
        rdat = {'residue': res, 'atoms': {}}

        debug_print(res_code, res.residueType)

        for atm in res.nmrAtoms:
            debug_print(atm.name)
            pdat = []
            for cs in atm.chemicalShifts:
                for peak in cs.assignedPeaks:
                    debug_print(peak,'| include =',peak.spectrum.name in SPECTRA)
                    if peak.spectrum.name in SPECTRA:
                        dim_idx = None
                        for i in range(len(peak.assignments[0])):
                            ass = peak.assignments[0][i]
                            if ass and ass.atom.id == atm.id:
                                dim_idx = i

                                # Save the atom type as an assigned type that needs to be exported
                                atmType = atm.name.strip().replace('%','')
                                if atmType not in atomTypes: atomTypes.add(atmType)

                                # Found the correct dimension index, stop searching
                                break

                        if dim_idx is not None:
                            # Get the assignment linewidth 
                            err = peak.lineWidths[dim_idx]
                            frq = peak.spectrum.spectrometerFrequencies[dim_idx]
                            if err: err = err / frq
                            else: err = peak.spectrum.ppmPerPoints[dim_idx] * 3 # Error estimation attempt based on spectral resolution

                            shift = peak.position[dim_idx]

                            if PERFORM_TROSY_CORRECTION: shift = correctTrosyShift(shift, atmType, frq)

                            pdat.append({'val':shift, 'err':err})

            if len(pdat) > 0:
                # Calculate (linewidth weighted) average chemical shift
                cs, cs_n, cs_e = calcChemicalShift(pdat, weigting = CHEMICALSHIFT_WEIGHTING)
                debug_print(cs,cs_e,n)

                rdat['atoms'][atm.name] = {'cs':cs,'n':cs_n,'err':cs_e}
            else: print(res_code,res.residueType,atm.name,'NO PEAKS FOUND')

        # Determine minimum peaks for residue
        tmp = []
        for k in rdat['atoms']:
            tmp.append(rdat['atoms'][k]['n'])
        rdat['minAssignedPeak'] = min(tmp) if len(tmp) > 0 else 0

        data[res_code] = rdat

    # Sort atom types for consistent ordering.
    sorted_atoms = sorted(list(atomTypes))

    return data, sorted_atoms

def exportStandardFormat(data, sorted_atoms, outFilePath):

    # Build header row: start with Assignment, then for each atom type add three columns.
    header = ["assign"]
    for at in sorted_atoms:
        header.append(f"{at}_cs")
        if INCLUDE_ERRORS: header.append(f"{at}_err")
        if INCLUDE_CS_SPECIFIC_PEAK_COUNT: header.append(f"{at}_n")

    if INCLUDE_PEAK_COUNT and not INCLUDE_CS_SPECIFIC_PEAK_COUNT: header.append(f"n_peaks")
    
    # Write CSV file.
    with open(outFilePath, "w", newline='') as f:
        writer = csv.writer(f, delimiter=DELIMITER)
        writer.writerow(header)
        
        # Process residues in sorted order.
        for key in sorted(data.keys()):
            rdat = data[key]
            # Build assignment label: convert the three-letter code to one-letter.
            code = convertResidueCode(rdat['residue'].residueType,
                inputCodeType='threeLetter',
                outputCodeType= 'threeLetter' if THREELETTER else 'oneLetter')
            assignment = f"{code.capitalize()}{rdat['residue'].sequenceCode}"

            row = [assignment]
            # For each atom type, output the shift, error, and nPeaks if available.
            for at in sorted_atoms:
                if at in rdat['atoms']:
                    val = rdat['atoms'][at]['cs']
                    err = rdat['atoms'][at]['err']
                    npeaks = rdat['atoms'][at]['n']

                    val_str = f"{round(val, 3)}" if val is not None else ""
                    err_str = f"{round(err, 3)}" if err is not None else ""
                    npeaks_str = f"{npeaks}"

                    row.append(val_str)
                    if INCLUDE_ERRORS: row.append(err_str)
                    if INCLUDE_CS_SPECIFIC_PEAK_COUNT: row.append(npeaks_str)
                else: # Atom not found, add empty string
                    row.append("")
                    if INCLUDE_ERRORS: row.append("")
                    if INCLUDE_CS_SPECIFIC_PEAK_COUNT: row.append("")

                if INCLUDE_PEAK_COUNT and not INCLUDE_CS_SPECIFIC_PEAK_COUNT: row.append(rdat['minAssignedPeak'])

                writer.writerow(row)

def exportTALOSFormat(data, outFilePath):
    sequence = []
    residues_with_shifts = []

    for resnum in sorted(data.keys()):
        rdat = data[resnum]

        one_letter = convertResidueCode(rdat['residue'].residueType,
            inputCodeType='threeLetter',
            outputCodeType= 'oneLetter')
        sequence.append(one_letter)
        
        for atm in rdat['atoms']:
            residues_with_shifts.append((resnum, one_letter, ccp2talos[atm], rdat['atoms'][atm]['cs']))

    with open(outFilePath, 'w') as f:
        f.write(f'REMARK CCPNMR {NMRCHAIN}:{ASSIGNMENT_SPECTRA_COLLECTION} Chemical Shift Table Export\n\n')
        f.write("DATA FIRST_RESID 1\n\n")

        # Write the sequence in 10-letter blocks
        f.write("DATA SEQUENCE ")
        for i in range(0, len(sequence), 10):
            f.write(''.join(sequence[i:i+10]) + ' ')
        f.write("\n\n")

        f.write("VARS   RESID RESNAME ATOMNAME SHIFT\n")
        f.write("FORMAT %4d   %1s     %4s      %8.3f\n\n")

        for resid, resname, atomname, shift in residues_with_shifts:
            f.write(f"{resid:4d} {resname:>2}   {atomname:>4}  {shift:8.3f}\n")

# Main execution block.
if __name__ == "__main__":

    data, sorted_atoms = getChemicalShifts()

    if OUTPUT_FORMAT == 'TALOS': exportTALOSFormat(data, OUTPUTPATH)
    else: exportChemicalShiftList(data, sorted_atoms, OUTPUTPATH)

    print(f'Exported {OUTPUT_FORMAT} chemical shift list to: {OUTPUTPATH}')
 
