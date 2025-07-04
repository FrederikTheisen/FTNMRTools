# CCPNMRv3 macro script to synchronize a Molecular Chain (MC) and NMR Chain (NC)
# Edit these PIDs to match your project setup

pid = 'heTDR1'

mc_pid = 'MC:' + pid    # Molecular Chain PID
nc_pid = 'NC:' + pid    # NMR Chain PID

mc = project.getByPid(mc_pid)
nc = project.getByPid(nc_pid)

# Step 1: Add missing NmrResidues based on MolecularChain residues
print("Adding missing NmrResidues...")
for residue in mc.residues:
    seq_code = str(residue.sequenceCode)
    res_type = str(residue.residueType)
    existing = nc.fetchNmrResidue(sequenceCode=seq_code, residueType=res_type)
    if not existing:
        print(f"Added NmrResidue for seqCode {seq_code}")
        #nc.fetchNmrResidue(sequenceCode=seq_code)

# Step 2: Merge offset residues with base residues
print("Merging offset residues...")
for nmr_res in list(nc.nmrResidues):
    if nmr_res.relativeOffset == -1:
        main_res = nmr_res.mainNmrResidue
        if main_res and main_res.previousNmrResidue:
            merge_target = main_res.previousNmrResidue
            print(f"Merging {nmr_res.pid} into {merge_target.pid}")
            for nmr_atom in nmr_res.nmrAtoms:
                try:
                    nmr_atom.assignTo(sequenceCode=merge_target.sequenceCode, mergeToExisting=True)
                except Exception as e:
                    print(f"  Failed to assign {nmr_atom.name}: {e}")
            project.deleteObjects(nmr_res)

    elif nmr_res.relativeOffset == +1:
        main_res = nmr_res.mainNmrResidue
        if main_res and main_res.nextNmrResidue:
            merge_target = main_res.nextNmrResidue
            print(f"Merging {nmr_res.pid} into {merge_target.pid}")
            for nmr_atom in nmr_res.nmrAtoms:
                try:
                    nmr_atom.assignTo(sequenceCode=merge_target.sequenceCode, mergeToExisting=True)
                except Exception as e:
                    print(f"  Failed to assign {nmr_atom.name}: {e}")
            project.deleteObjects(nmr_res)

print("Finished updating NmrChain.")
