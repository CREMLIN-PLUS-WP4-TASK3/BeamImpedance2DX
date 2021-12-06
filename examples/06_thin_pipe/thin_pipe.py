from bi2d.materials import vacuum, beam
from mpi4py import MPI
import bi2d
import logging
import numpy as np

# Assign material indices. When setting index `-1` to SIBC material, boundary will be applied to all boundaries
beam.index = 1
vacuum.index = 2
steel = bi2d.Material(3, sigma=1e6)
outer_vacuum = vacuum.copy()
outer_vacuum.index = 4
wall = steel.copy()
wall.index = -1

# Setup mesh and materials
m_v = bi2d.Mesh("vacuum.xdmf")
m_vm = bi2d.Mesh("vacuum_metal.xdmf")
m_vmv = bi2d.Mesh("vacuum_metal_vacuum.xdmf")
mc_v = bi2d.MaterialMap(m_m, [beam, vacuum])
mc_vm = bi2d.MaterialMap(m_mv, [beam, vacuum, steel])
mc_vmv = bi2d.MaterialMap(m_mvm, [beam, vacuum, steel, outer_vacuum])

# Configure solution
solution_v = bi2d.Solution(mc_v, Hcurl_order=2, H1_order=2)
solution_vm = bi2d.Solution(mc_vm, Hcurl_order=2, H1_order=2)
solution_vmv = bi2d.Solution(mc_vmv, Hcurl_order=2, H1_order=2)
solution_v.logger.setLevel(logging.INFO)
solution_vm.logger.setLevel(logging.INFO)
solution_vmv.logger.setLevel(logging.INFO)

# Calculate impedance in frequency range and save the results
for solution, name, f_range, sibc_boundaries, source_function \
    in [(solution_vm, "metal", np.logspace(3, 8, num=100), [], bi2d.SourceFunction.MONOPOLE),
        (solution_v, "sibc", np.logspace(6, 12, num=100), [wall], bi2d.SourceFunction.MONOPOLE),
        (solution_vmv, "metal", np.logspace(3, 8, num=100), [], bi2d.SourceFunction.DIPOLE),
        (solution_v, "sibc", np.logspace(6, 12, num=100), [wall], bi2d.SourceFunction.DIPOLE)]:
    data = solution.get_z(f_range, beta=0.999999, source_function=source_function)
    if MPI.COMM_WORLD.rank == 0:
        np.savetxt(f"{source_function.name.lower()}_{name}.dat", data)
