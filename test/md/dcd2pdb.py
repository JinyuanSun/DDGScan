def dcd2pdb(dcd_file, topol_file, out_file, stride=1, noWater=True, superimpose=True):

    top = mt.load(topol_file)
    if noWater:
        indices = top.topology.select("protein")
    else:
        indices = top.topology.select("all")

    traj = mt.load_dcd(dcd_file, top=topol_file, stride=stride, atom_indices=indices)

    if superimpose:
        print("INFO: Superimposing to topology ......")
        CA_indices = top.topology.select("protein and name CA")
        traj.superpose(top, ref_atom_indices=CA_indices, atom_indices=CA_indices)

    traj.save_pdb(out_file)

    return None

if __name__ == '__main__':
    import sys
    import mdtraj as mt
    dcd_file, topol_file, out_file = sys.argv[1:]
    dcd2pdb(dcd_file, topol_file, out_file)
