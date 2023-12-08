import sys
from simtk.openmm import app, LangevinIntegrator, MonteCarloBarostat
from simtk import openmm
import openmm.unit as u
from openmm import DrudeLangevinIntegrator, MonteCarloMembraneBarostat
from simtk.unit import *
import os

# Load the system state from the previous simulation
with open('min.xml', 'r') as f:
    system = openmm.XmlSerializer.deserialize(f.read())



integrator = DrudeLangevinIntegrator(
    300.0 * u.kelvin,  # temperature
    5.0 / u.picosecond,  # friction coeff
    1.0 * u.kelvin,  # drude temperature
    20.0 / u.picosecond,  # drude friction
    0.001 * u.picosecond  # timestep
)
integrator.setMaxDrudeDistance(0.025) # Drude Hardwall

barostat = MonteCarloMembraneBarostat(
    1.0 * u.atmosphere,  # Default pressure
    0.0,                 # Default surface tension
    integrator.getTemperature(),  # Temperature
    MonteCarloMembraneBarostat.XYIsotropic,  # XY mode
    MonteCarloMembraneBarostat.ZFree,        # Z mode
    10                  # Frequency of Monte Carlo volume changes
    )
system.addForce(barostat)

# Load PSF file (for topology)
psf = app.CharmmPsfFile('step2_drude_.psf')
A = 10
B = 10
C = 20
psf.setBox(A* u.nanometers, B * u.nanometers, C *u.nanometers)


# Create a Simulation object
simulation = app.Simulation(psf.topology, system, integrator)

# Load the checkpoint file to continue from the previous state
checkpoint_file = 'min_system.chk'  # Replace with your checkpoint file
simulation.loadCheckpoint(checkpoint_file)

# Add reporters as needed
output_file = 'out.txt'
simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(app.DCDReporter('bulk/trajectory.dcd', 1000))
simulation.reporters.append(app.StateDataReporter(output_file, 1000, step=True, potentialEnergy=True, temperature=True, volume=True,density=True,speed=True))
simulation.reporters.append(app.CheckpointReporter('bulk/checkpoint.chk', 1000))  # Save checkpoints every 10000 steps

# Continue the simulation
simulation.step(500000)

# Save system and integrator states in XML format at the end
## rename to prefered output file
with open('bulk/continued_system.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.system))
with open('bulk/continued_integrator.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.integrator))

# Save the final coordinates in PDB format
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open('bulk/continued_output.pdb', 'w'))
