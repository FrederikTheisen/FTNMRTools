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
AUTO_FILTER_OVERLAPS = True
COMMENT_OVERLAP = True
PSEUDO3D = True
DIMORDER = {'H':2,'N':1,'PSEUDO':0, 'C':1}
SEQUENCE_CODE = 'oneLetter' #'oneLetter' or 'threeLetter'
OVERLAP_DISTANCE_THRESHOLD = 0.02
ATOM_CS_WEIGHTS = {'H':1, 'N':0.154, 'C':0.3}

from collections import defaultdict
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

def runMacro():
    # Retrieve the currently active peak list.
    # (getCurrentPeakList() is provided by the CCPNMR3 macro environment.)
    peakList = get(PEAKLIST)
    if not peakList:
        print("No active peak list found!")
        return

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
    
    with open(OUTPUTPATH, "w") as f:
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

    print("Sparky peak list exported successfully to:", OUTPUTPATH)

if __name__ == '__main__':
    runMacro()
