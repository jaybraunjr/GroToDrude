import sys
import os
import openmm
import openmm.app as app
import openmm.unit as u
from openmm import DrudeLangevinIntegrator

# Load the system state from the previous simulation
with open('min_system.xml', 'r') as f:
    system = openmm.XmlSerializer.deserialize(f.read())

# Integrator setup
integrator = DrudeLangevinIntegrator(
    400.0 * u.kelvin,      # Temperature for non-Drude particles
    5.0 / u.picosecond,    # Friction coefficient for non-Drude particles
    1.0 * u.kelvin,        # Temperature for Drude particles
    20.0 / u.picosecond,   # Friction coefficient for Drude particles
    0.001 * u.picosecond   # Time step
)
integrator.setMaxDrudeDistance(0.025)  # Drude hardwall constraint

# Load PSF file for topology and set the simulation box dimensions
psf = app.CharmmPsfFile('step2_drude_.psf')
psf.setBox(8 * u.nanometers, 8 * u.nanometers, 17 * u.nanometers)

# Create a Simulation object
simulation = app.Simulation(psf.topology, system, integrator)

# Load the checkpoint file to continue from the previous state
checkpoint_file = 'min_system.chk'
simulation.loadCheckpoint(checkpoint_file)

# Add reporters for output
output_file = 'nvt/out.txt'
simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(app.DCDReporter('nvt/trajectory.dcd', 1000))
simulation.reporters.append(app.StateDataReporter(output_file, 1000, step=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True))
simulation.reporters.append(app.CheckpointReporter('nvt/checkpoint.chk', 1000))  # Save checkpoints every 1000 steps

# Run the simulation
simulation.step(50000)

# Save system and integrator states in XML format at the end
with open('nvt/system.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.system))
with open('nvt/integrator.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.integrator))

# Save the final coordinates in PDB format
positions = simulation.context.getState(getPositions=True).getPositions()
with open('nvt/output.pdb', 'w') as f:
    app.PDBFile.writeFile(simulation.topology, positions, f)
