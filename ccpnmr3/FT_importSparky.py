"""
CCPNMR3 Macro: Assign Peaks from a Sparky HSQC Peak List

This script reads a Sparky-format peak list file and for each peak:
  - Extracts a residue number from the assignment column (first integer found),
    then applies an integer offset.
  - Reads the chemical shift positions from the specified columns.
  - Uses an NMR chain (with id NMR_CHAIN_ID) to obtain the NMR residue by its
    sequence code and then looks up the corresponding atoms (using TRANSLATION_TABLE)
    for each dimension.
  - Adds a new peak to the target peak list (with id PEAK_LIST_ID) at the sparky
    positions and assigns the NMR atoms.
  
User-adjustable parameters are defined near the top.
"""

import os
import re
import sys
import ccpn

# --- User-adjustable parameters ---

# Sparky peak list file
SPARKY_PEAKLIST = "/home/ftheisen/Documents/hANP32a_FL_500mMNaCl_20C.list"

# IDs of objects in Analysis (set these to your actual project IDs)
PEAK_LIST_ID = "PL:500mMNaCl_20C_hANP32aFL.1"      # The (empty) peak list in Analysis to which peaks will be added
NMR_CHAIN_ID = "NC:FL"         # The NMR chain from which residues/atoms are obtained

# Translation table: map sparky dimension names to NMR atom names.
# For example, if sparky column "w2" corresponds to the H nucleus and "w1" to N.
TRANSLATION_TABLE = {
    "w2": "H",
    "w1": "N"
}

# Integer offset for residue numbering (applied to the number extracted from the sparky assignment)
RESIDUE_NUMBER_OFFSET = 0

# --- End of user-adjustable parameters ---

# Built-in functions from the CCPNMR3 macro environment are assumed available.
# For example, project.getByPid() to obtain objects by ID.

def parse_residue_number(assignment_str):
    """
    Extract the first integer found in the assignment string.
    Return None if no integer is found.
    """
    m = re.search(r'(\d+)', assignment_str)
    if m:
        return int(m.group(1))
    return None

def read_sparky_file(filepath):
    """
    Reads a Sparky peak list file in a generic manner.
    
    Expected format:
      - The first column is "assignment", which contains a string with at least one integer (the residue number).
      - The following columns are named "wX" (e.g., w1, w2, w3, ...) representing chemical shift dimensions.
      - Other columns (e.g. intensity) are ignored.
    
    The function extracts:
      - The residue number from the assignment column.
      - For each column starting with 'w', converts the column name to an atom name using the TRANSLATION_TABLE.
        If no translation exists for a given column, that column's header is used as the atom name.
      - The chemical shift values (converted to float) for those dimensions.
    
    Returns:
      A list of dictionaries, one per peak, with structure:
          {
              "resnum": <int>, 
              "dims": {
                  <atom_name>: <chemical_shift>, 
                  ...
              }
          }
    """
    peaks = []
    with open(filepath, 'r') as f:
        print("Reading file: filepath")
        # Read header line
        header_line = f.readline().strip()
        if not header_line:
            print(f"Error: File {filepath} is empty or missing header.")
            return []
        headers = header_line.split()
        # Normalize headers: first column should be assignment.
        headers = [h.strip().lower() for h in headers]
        if len(headers) < 2:
            print(f"Error: File {filepath} must have at least two columns.")
            return []
        
        # Identify all dimension columns (those that start with "w")
        # Also create a mapping from file column index -> atom name via the translation table.
        dim_indices = {}
        for i, col in enumerate(headers[1:], start=1):
            if col.startswith("w"):
                # Look up the translation; if not found, use the original column name.
                atom_name = TRANSLATION_TABLE.get(col, col)
                dim_indices[i] = atom_name
            else:
                # Ignore non-w columns (e.g. intensity)
                continue
        
        # Process each line
        for line in f:
            if not line.strip():
                continue
            
            parts = line.split()
            d = {}
            # First column: assignment. Extract residue number.
            assignment_str = parts[0]
            # Extract integers from the assignment field
            numbers = re.findall(r'\d+', assignment_str)
            if not numbers:
                d["resnum_candidates"] = [None]
                print(f"Warning: No residue number found in assignment '{assignment_str}'")
            elif len(numbers) == 1:
                d["resnum_candidates"] = [int(numbers[0])]
            else:
                d["resnum_candidates"] = [int(n) for n in numbers]

            d["resnum"] = d["resnum_candidates"][0]
            dims = {}
            # Process each dimension column
            for idx, atom_name in dim_indices.items():
                try:
                    shift_val = float(parts[idx])
                    dims[atom_name] = shift_val
                except ValueError:
                    print(f"Warning: Could not convert value '{parts[idx]}' to float in file {filepath}")
                    dims[atom_name] = None
            d["dims"] = dims
            peaks.append(d)
            print(f'{parts[0]}:\t{d["resnum_candidates"]}\t{dims}')
    
    return peaks

def get_nmr_chain(nmr_chain_id):
    chain = project.getByPid(nmr_chain_id)
    if chain is None:
        print(f"Error: NMR chain with id '{nmr_chain_id}' not found.")
        sys.exit(1)
    return chain

def get_peak_list(peak_list_id):
    pl = project.getByPid(peak_list_id)
    if pl is None:
        print(f"Error: Peak list with id '{peak_list_id}' not found.")
        sys.exit(1)
    return pl

def assign_peak_to_chain(peak_data, chain):
    peak = {}
    assignments = []

    resnum = peak_data.get("resnum")
    resnum_corrected = resnum + RESIDUE_NUMBER_OFFSET

    # Fetch the corresponding NMR residue from the chain.
    # TODO: This should be moved to within the dimension loop to account for peaks originating from more than one residue
    # TODO: Some dimension offset table could be used to determine the actual residue. Alternatively, a not shit Sparky peak list should contain this info
    residue = chain.fetchNmrResidue(resnum_corrected)
    if residue is None:
        print(f"Warning: Residue with sequence code {resnum_corrected} not found in chain {chain.pid}")
        return assignments
    
    # Iterate over dimensions in the order defined by the translation_table.
    for dim in peak_data['dims']:
        assigment = {}
        # Retrieve the chemical shift for the given atom
        assigment['shift'] = peak_data["dims"][dim]
        # Fetch the corresponding NMR atom from the residue.
        assigment['atom'] = residue.fetchNmrAtom(dim)

        assignments.append(assigment)

    if len(peak_data["resnum_candidates"]) > 1:
        peak['comment'] = f'{",".join(peak_data["resnum_candidates"])}'

    peak['assignments'] = assignments
    
    return peak

def sort_assignment_object(assignments, spectrum_axes):
    """
    Reorder the assignment object so that the list of assignments matches
    the order of the spectrum_axes.
    """
    # Make a copy of the assignment list that we can remove from.
    available = assignments.copy()
    sorted_assignment = []
    
    for axis in spectrum_axes:
        found = None
        for candidate in available:
            # Assume axisCodeMatch returns the best matching axis for candidate's nmr_atom.
            candidate_axis = ccpn.core.lib.AxisCodeLib.axisCodeMatch(candidate['atom'].name, spectrum_axes)
            if candidate_axis == axis:
                found = candidate
                break
        if found is not None:
            sorted_assignment.append(found)
            available.remove(found)
        else:
            print(f"Warning: No assignment found for axis '{axis}'")
    
    # Optionally, you might want to append any remaining assignments that didn't match any axis.
    if available:
        print("Note: Some assignments were not matched to a spectrum axis and will be appended at the end.")
        sorted_assignment.extend(available)
    
    return sorted_assignment

def add_peak_to_peaklist(peak_list_obj, peak):
    """
    Adds a new peak to the provided peak list using built-in methods and assigns NMR atoms.
    """
    sorted_assignment = sort_assignment_object(peak['assignments'], peak_list_obj.spectrum.axisCodes)

    ppmPositions = [a['shift'] for a in sorted_assignment]
    assigned_atoms = [a['atom'] for a in sorted_assignment]

    new_peak = peak_list_obj.newPeak(ppmPositions=ppmPositions, comment=peak['comment'])

    for i in range(len(peak_list_obj.spectrum.axisCodes)):
        try:
            new_peak.assignDimension(peak_list_obj.spectrum.axisCodes[i], assigned_atoms[i])
        except Exception as e:
            print(f"Warning: Could not assign atoms for peak at {ppmPositions}: {e}")

def main():
    print("Peak list file:",SPARKY_PEAKLIST)
    print("Target peak list:", PEAK_LIST_ID)
    print("NMR chain:", NMR_CHAIN_ID)
    print("Residue number offset:", RESIDUE_NUMBER_OFFSET)
    print("Translation table:", TRANSLATION_TABLE)

    # Get the target peak list and NMR chain objects.
    pl_obj = get_peak_list(PEAK_LIST_ID)
    nmr_chain = get_nmr_chain(NMR_CHAIN_ID)
    
    peaks_data = read_sparky_file(SPARKY_PEAKLIST)
    print(f"Found {len(peaks_data)} peaks in file.")

    for peak_data in peaks_data:
        # Get ppm positions and assigned NMR atoms for this peak.
        assignments = assign_peak_to_chain(peak_data, nmr_chain)
        if len(assignments) == 0:
            continue  # Skip peaks with no residue number
        # Add the peak to the peak list.
        add_peak_to_peaklist(pl_obj, assignments)
    
    print("Peak import complete.")

if __name__ == "__main__":
    main()
