import os
import openmm.app as app
import openmm
import openmm.unit as unit

# Load the PSF file
psf_file = 'step2_drude_.psf'  # Replace with the actual path to your PSF file
psf = app.CharmmPsfFile(psf_file)

# Set the box dimensions
A = 8.68698  # Replace with your actual value
B = 8.68698  # Replace with your actual value
C = 18.13459 # Replace with your actual value
psf.setBox(A * unit.nanometers, B * unit.nanometers, C * unit.nanometers)

# Load the system state from the XML file
state_file = 'min_system.xml'
with open(state_file, 'r') as f:
    state_xml = f.read()
state = openmm.XmlSerializer.deserialize(state_xml)

# Load force field parameters
toppath = './toppar_drude'
toppar = [
    os.path.join(toppath, "toppar_drude_master_protein_2022c.str"),
    os.path.join(toppath, "lipid_drud.str")
]
params = app.CharmmParameterSet(*toppar)

# Create the system
system = psf.createSystem(params, nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers)

# Define the DrudeLangevinIntegrator
temperature = 300 * unit.kelvin
frictionCoeff = 1/unit.picosecond
timestep = 1 * unit.femtoseconds
drudeTemperature = 1 * unit.kelvin
drudeFrictionCoeff = 20/unit.picosecond
integrator = openmm.DrudeLangevinIntegrator(temperature, frictionCoeff, drudeTemperature, drudeFrictionCoeff, timestep)

# Create a Simulation object
simulation = app.Simulation(psf.topology, system, integrator)

# Set the initial positions from the state
simulation.context.setPositions(state.getPositions())

# Define the maximum number of minimization steps
max_minimization_steps = 5000  # Adjust this number as needed

# Minimize the energy with a specified number of steps
print("Further minimizing energy with a maximum of {} steps...".format(max_minimization_steps))
simulation.minimizeEnergy(maxIterations=max_minimization_steps)

# Get the state after minimization with positions and energy
final_state = simulation.context.getState(getEnergy=True, getPositions=True)
final_energy = final_state.getPotentialEnergy()

# Print the final potential energy
print("Final potential energy after minimization: ", final_energy)

# Save the further minimized coordinates
positions = final_state.getPositions()
app.PDBFile.writeFile(psf.topology, positions, open('min2.pdb', 'w'))

# Optionally, save the new state in XML and checkpoint formats
with open('min_system2.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.system))
with open('min2.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.integrator))
with open('min2.chk', 'wb') as f:
    f.write(simulation.context.createCheckpoint())

print("Further minimization completed.")
