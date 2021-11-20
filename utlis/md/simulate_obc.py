from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

# load in input PDB file and force field XML files
pdb = app.PDBFile('6JTT_fixed.pdb')
#forcefield = app.ForceField('amber10.xml')
forcefield = app.ForceField('amber99sb.xml', 'amber99_obc.xml')

# use app.Modeller to add hydrogens and solvent
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
#modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometers)
app.PDBFile.writeFile(modeller.topology, modeller.positions, open('6JTT_obc.pdb', 'w'))

# prepare system and integrator
system=forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, soluteDielectric=2.0,
        solventDielectric=80.0)

#system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, 
#    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, 
#    ewaldErrorTolerance=0.0005)
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
simulation.reporters.append(app.DCDReporter('trajectory_obc.dcd', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=500000, separator='\t'))

# run 50 ns of production simulation
print('Running Production...')
simulation.step(500000)
print('Done!')
