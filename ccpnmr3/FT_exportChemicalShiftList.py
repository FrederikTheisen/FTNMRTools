#!/usr/bin/env python
"""
CCPNMR3 Macro: Export a Chemical Shift List
Frederik Theisen 2025

This script groups chemical shifts by residue. For each residue (identified
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

CHEMICALSHIFTLIST = 'CL:default'
OUTPUTPATH = '/home/ftheisen/Documents/' + CHEMICALSHIFTLIST + '.list'
THREELETTER = True
INCLUDE_ERRORS = True
INCLUDE_PEAK_COUNT = True
INCLUDE_CS_SPECIFIC_PEAK_COUNT = False
DELIMITER = '\t'

from collections import defaultdict
import csv

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

if not INCLUDE_PEAK_COUNT: SHOW_CS_SPECIFIC_PEAK_COUNT = False

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

def exportChemicalShiftList(csList, outFilePath):
    """
    Export a CCPNMR3 chemical shift list to a CSV file.
    
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
    res_data = {}
    # Set of all atom types encountered.
    atomTypes = set()
    
    for cs in csList.chemicalShifts:
        atom = cs.nmrAtom
        if not atom:
            continue
        residue = atom.nmrResidue
        if not residue:
            continue
        chainCode = residue.chain.code if hasattr(residue, 'chain') and residue.chain else ""
        seqCode = residue.sequenceCode
        res_key = f"{chainCode}{seqCode}"
        
        # Use the atom's name (trimmed) as its type (e.g. "H", "N", "CA", etc.)
        atType = atom.name.strip().replace('%','')
        atomTypes.add(atType)
        
        # Retrieve values.
        shift_val = cs.value
        shift_err = cs.valueError if hasattr(cs, 'valueError') else None
        try:
            nPeaks = len(cs.assignedPeaks)
        except AttributeError:
            nPeaks = 1
        
        if res_key not in res_data:
            res_data[res_key] = {'residue': residue, 'atoms': {}}
        res_data[res_key]['atoms'][atType] = (shift_val, shift_err, nPeaks)

    for res_key in res_data:
        atms = [res_data[res_key]['atoms'][at] for at in res_data[res_key]['atoms']]
        res_data[res_key]['minAssignedPeak'] = min([cs[2] for cs in atms])
    
    # Sort atom types for consistent ordering.
    sorted_atoms = sorted(list(atomTypes))
    
    # Build header row: start with Assignment, then for each atom type add three columns.
    header = ["assign"]
    for at in sorted_atoms:
        header.append(f"{at}_shift")
        if INCLUDE_ERRORS: header.append(f"{at}_error")
        if INCLUDE_CS_SPECIFIC_PEAK_COUNT: header.append(f"{at}_peaks")

    if INCLUDE_PEAK_COUNT and not INCLUDE_CS_SPECIFIC_PEAK_COUNT: header.append(f"peaks")
    
    # Write CSV file.
    with open(outFilePath, "w", newline='') as f:
        writer = csv.writer(f, delimiter=DELIMITER)
        writer.writerow(header)
        
        # Process residues in sorted order.
        for res_key in sorted(res_data.keys()):
            residue = res_data[res_key]['residue']
            # Build assignment label: convert the three-letter code to one-letter.
            code = convertResidueCode(residue.residueType,
                                           inputCodeType='threeLetter',
                                           outputCodeType= 'threeLetter' if THREELETTER else 'oneLetter')
            assignment = f"{code.capitalize()}{residue.sequenceCode}"
            row = [assignment]
            # For each atom type, output the shift, error, and nPeaks if available.
            for at in sorted_atoms:
                if at in res_data[res_key]['atoms']:
                    val, err, npeaks = res_data[res_key]['atoms'][at]
                    val_str = f"{round(val, 3)}" if val is not None else ""
                    err_str = f"{round(err, 3)}" if err is not None else ""
                    npeaks_str = f"{npeaks}"
                    #row.extend([val_str, err_str, npeaks_str]) #Old code

                    row.append(val_str)
                    if INCLUDE_ERRORS: row.append(err_str)
                    if INCLUDE_CS_SPECIFIC_PEAK_COUNT: row.append(npeaks_str)
                else:
                    row.append("")
                    if INCLUDE_ERRORS: row.append("")
                    if INCLUDE_CS_SPECIFIC_PEAK_COUNT: row.append("")

            if INCLUDE_PEAK_COUNT and not INCLUDE_CS_SPECIFIC_PEAK_COUNT: row.append(res_data[res_key]['minAssignedPeak'])

            writer.writerow(row)
    print("Exported chemical shift list to:", outFilePath)

# Main execution block.
if __name__ == "__main__":
    try:
        csList = get(CHEMICALSHIFTLIST)
    except Exception as e:
        print("Error retrieving chemical shift list:", e)
        csList = None
    
    if csList is not None:
        outFile = OUTPUTPATH

        exportChemicalShiftList(csList, outFile)
    else:
        print("No chemical shift list available.")
