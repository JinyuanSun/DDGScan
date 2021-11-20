from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

# load in input PDB file and force field XML files
pdb = app.PDBFile('1o9s_fixed.pdb')
forcefield = app.ForceField('amber99sbildn.xml', 'tip4pew.xml')

# use app.Modeller to add hydrogens, extra particles, and solvent
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens()
modeller.addExtraParticles(forcefield)
modeller.addSolvent(forcefield, model='tip4pew', padding=1.0*unit.nanometers)
app.PDBFile.writeFile(modeller.topology, modeller.positions, open('1o9s_modeller_tip4pew.pdb', 'w'))

# prepare system and integrator
system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)

# prepare simulation
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(modeller.topology, system, integrator, platform, 
    properties)
simulation.context.setPositions(modeller.positions)

# minimize
print('Minimizing...')
simulation.minimizeEnergy()

# equilibrate for 100 steps
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulation.step(100)

# append reporters
simulation.reporters.append(app.DCDReporter('trajectory_tip4pew.dcd', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=25000000, separator='\t'))

# run 50 ns of production simulation
print('Running Production...')
simulation.step(25000000)
print('Done!')
