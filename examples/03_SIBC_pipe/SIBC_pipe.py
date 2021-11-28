from bi2d.materials import vacuum, beam, steel
from itertools import product
from mpi4py import MPI
import bi2d
import logging
import numpy as np

# Assign material indices. When setting index `-1` to SIBC material, boundary will be applied to all boundaries
beam.index = 1
vacuum.index = 2
steel.index = 3
wall = steel.copy()
wall.index = -1

# Setup mesh and materials
m_metal = bi2d.Mesh("metal.xdmf")
m_sibc = bi2d.Mesh("sibc.xdmf")
mc_metal = bi2d.MaterialMap(m_metal, [beam, vacuum, steel])
mc_sibc = bi2d.MaterialMap(m_sibc, [beam, vacuum])

# Configure solution
solution_metal = bi2d.Solution(mc_metal, Hcurl_order=2, H1_order=2)
solution_metal.logger.setLevel(logging.INFO)
solution_sibc = bi2d.Solution(mc_sibc, Hcurl_order=2, H1_order=2)
solution_sibc.logger.setLevel(logging.INFO)

# Calculate impedance in frequency range and save the results
for source_function in [bi2d.SourceFunction.MONOPOLE, bi2d.SourceFunction.DIPOLE]:
    data = solution_metal.get_z(np.logspace(1, 13, num=20), beta=0.99999, source_function=source_function)
    if MPI.COMM_WORLD.rank == 0:
        np.savetxt(f"{source_function.name.lower()}_metal.dat", data)

for source_function in [bi2d.SourceFunction.MONOPOLE, bi2d.SourceFunction.DIPOLE]:
    data = solution_sibc.get_z(np.logspace(1, 13, num=20), beta=0.99999, source_function=source_function, sibc=[wall])
    if MPI.COMM_WORLD.rank == 0:
        np.savetxt(f"{source_function.name.lower()}_sibc.dat", data)
