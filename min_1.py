import os
import openmm.app as app
import openmm
import openmm.unit as u
from openmm.app import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet, LJPME,\
                       HBonds, Simulation, DCDReporter, StateDataReporter
from openmm import Platform, MonteCarloMembraneBarostat, DrudeLangevinIntegrator



# Load the PSF and CRD files
psf = app.CharmmPsfFile('step2_drude_.psf')
crd = app.CharmmCrdFile('step2_drude_.crd')

# Set the box dimensions to a little bigger than what system is 
A = 10
B = 10
C = 20
psf.setBox(A * u.nanometers, B * u.nanometers, C * u.nanometers)

# Load force field parameters
toppath = './toppar_drude'
toppar = [
    os.path.join(toppath, "toppar_drude_master_protein_2022c.str"),
    os.path.join(toppath, "lipid_drud.str")
]
params = app.CharmmParameterSet(*toppar)



# Set platform
#platform = Platform.getPlatformByName('CUDA')
#prop = dict(CudaPrecision='mixed')

# Build simulation context
#simulation = Simulation(psf.topology, system, integrator, platform, prop)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)

# Drude VirtualSites
simulation.context.computeVirtualSites()

# Calculate initial system energy
print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Energy minimization
simulation.minimizeEnergy(maxIterations=1000)

# Get the state after minimization with positions and energy
final_state = simulation.context.getState(getEnergy=True, getPositions=True)
final_energy = final_state.getPotentialEnergy()

# Print the final potential energy
print("Final potential energy after minimization: ", final_energy)

# Save the further minimized coordinates
positions = final_state.getPositions()
app.PDBFile.writeFile(psf.topology, positions, open('min.pdb', 'w'))

# Optionally, save the new state in XML and checkpoint formats
with open('min.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(simulation.system))
# with open('min_integrator.xml', 'w') as f:
#     f.write(openmm.XmlSerializer.serialize(simulation.integrator))
with open('min.chk', 'wb') as f:
    f.write(simulation.context.createCheckpoint())
