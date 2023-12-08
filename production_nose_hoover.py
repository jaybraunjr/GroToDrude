import sys
import openmm
from simtk.openmm import app, LangevinIntegrator
from openmm import MonteCarloMembraneBarostat, MonteCarloBarostat

#from openmm.app import MonteCarloMembraneBarostat
from simtk import openmm
import openmm.unit as u
from openmm import DrudeLangevinIntegrator, DrudeNoseHooverIntegrator
from simtk.unit import *
import os

# Load the system state from the previous simulation
with open('npt/system2.xml', 'r') as f:
    system = openmm.XmlSerializer.deserialize(f.read())


temperature = 323.15 * u.kelvin
drude_temperature = 1.0 * u.kelvin
timestep = 0.001 * u.picosecond
friction = 1.0 / u.picosecond # That's how I interpret tau = 0.1 ps from Lamoureux and Roux
drude_friction = 1 #/ u. picosecond  

integrator = DrudeNoseHooverIntegrator(
       temperature, friction, drude_temperature, drude_friction, timestep,
       3, 3, 3
)

integrator.setMaxDrudeDistance(0.025)  # hard wall constraint


pressure = 1.0 * u.atmosphere  # Desired pressure
pressureFrequency = 25  # Frequency of Monte Carlo volume changes

barostat = MonteCarloBarostat(pressure, temperature, pressureFrequency)
system.addForce(barostat)



# Load PSF file (for topology)
psf = app.CharmmPsfFile('step2_drude_.psf')
psf.setBox(65.493126 * nanometers,65.493126 * nanometers, 241.439224 * nanometers)


# Create a Simulation object
simulation = app.Simulation(psf.topology, system, integrator)

# Load the checkpoint file to continue from the previous state
checkpoint_file = 'npt/checkpoint2.chk'  # Replace with your checkpoint file
simulation.loadCheckpoint(checkpoint_file)

# Add reporters as needed
output_file = 'npt.txt'
simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(app.DCDReporter('npt/trajectory3.dcd', 1000))
simulation.reporters.append(app.StateDataReporter(output_file, 1000, step=True, potentialEnergy=True, temperature=True, volume=True,density=True,speed=True))
simulation.reporters.append(app.CheckpointReporter('npt/checkpoint3.chk', 1000))  # Save checkpoints every 10000 steps

# Continue the simulation
simulation.step(5000000)

# Save system and integrator states in XML format at the end
with open('npt/system3.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.system))
with open('npt/integrator3.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.integrator))

# Save the final coordinates in PDB format
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open('npt/output_.pdb', 'w'))
with open('npt/output.pdb', 'w') as f:
    app.PDBFile.writeFile(simulation.topology, positions, f)
    # Write PBC box dimensions
    box = state.getPeriodicBoxVectors()
    f.write("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n".format(
        box[0][0] / u.angstrom, box[1][1] / u.angstrom, box[2][2] / u.angstrom,
        90.0, 90.0, 90.0))  # Assumes an orthogonal box
