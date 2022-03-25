#!/usr/bin/env python

# https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html
# http://docs.openmm.org/latest/userguide/application/02_running_sims.html
# Thanks Dr. Zheng Liangzhen for help with MDs.

from __future__ import print_function

import os
os.environ["OPENMM_CPU_THREADS"] = "4"

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
import mdtraj as mt


def fix(pdbfile):
    fixed = pdbfile.replace(".pdb", "_fixed.pdb")
    fixer = PDBFixer(pdbfile)
    numChains = len(list(fixer.topology.chains()))
    fixer.removeChains(range(1, numChains))
    fixer.findMissingResidues()

    # only add missing residues in the middle of the chain, do not add terminal ones
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    missingResidues = dict()
    for key in keys:
        chain = chains[key[0]]
        if not (key[1] == 0 or key[1] == len(list(chain.residues()))):
            missingResidues[key] = fixer.missingResidues[key]
    fixer.missingResidues = missingResidues

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    PDBFile.writeFile(
        fixer.topology,
        fixer.positions,
        open(pdbfile.replace(".pdb", "_fixed.pdb"), "w"),
    )
    return fixed


def produciton(pdbfilename, platform="CUDA"):

    # load in input PDB file and force field XML files
    pdb = app.PDBFile(pdbfilename)
    forcefield = app.ForceField("amber99sbildn.xml", "tip3p.xml")

    # use app.Modeller to add hydrogens and solvent
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, model="tip3p", padding=1.0 * unit.nanometers)
    topname = pdbfilename.replace("fixed", "modeller_tip3p")
    app.PDBFile.writeFile(modeller.topology, modeller.positions, open(topname, "w"))

    # prepare system and integrator
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    integrator = mm.LangevinIntegrator(
        300 * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds
    )
    integrator.setConstraintTolerance(0.00001)

    # prepare simulation
    if platform == "CUDA":
        platform = mm.Platform.getPlatformByName(platform)
        properties = {"CudaPrecision": "mixed"}
        simulation = app.Simulation(
            modeller.topology, system, integrator, platform, properties
        )
    else:
        platform = mm.Platform.getPlatformByName("CPU")
        simulation = app.Simulation(modeller.topology, system, integrator, platform)

    simulation.context.setPositions(modeller.positions)

    # minimize
    print("Minimizing...")
    simulation.minimizeEnergy()

    # equilibrate for 100 steps
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
    print("Equilibrating...")
    simulation.step(100)

    # append reporters
    dcdname = pdbfilename.replace("fixed.pdb", "_tip3p.dcd")
    simulation.reporters.append(app.DCDReporter(dcdname, 1000))
    simulation.reporters.append(
        app.StateDataReporter(
            stdout,
            1000,
            step=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=500000,
            separator="\t",
        )
    )

    # run 1 ns of production simulation
    print("Running Production...")
    simulation.step(500000)
    print("Done!")
    return topname, dcdname


def dcd2pdb(dcd_file, topol_file, out_file, stride=100, noWater=True, superimpose=True):

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


def main(pdbfile, out_file, platform):
    fixed = fix(pdbfile)
    topname, dcdname = produciton(fixed, platform)
    dcd2pdb(dcdname, topname, out_file)


if __name__ == "__main__":
    import sys

    pdbfile, out_file, platform = sys.argv[1:]
    main(pdbfile, out_file, platform)
