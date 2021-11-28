from bi2d.materials import vacuum, beam
from itertools import product
from mpi4py import MPI
import bi2d
import logging
import numpy as np

# Assign material indices
beam.index = 1
vacuum.index = 2

# Setup mesh and materials
m = bi2d.Mesh("mesh.xdmf")
mc = bi2d.MaterialMap(m, [beam, vacuum])

# Configure solution
solution = bi2d.Solution(mc, Hcurl_order=2, H1_order=2)
# Enable info level logging
solution.logger.setLevel(logging.INFO)

# Calculate impedance in frequency range and save the results
for beta, source_function in product([0.1, 0.9],
                                     [bi2d.SourceFunction.MONOPOLE, bi2d.SourceFunction.DIPOLE]):
    data = solution.get_z(np.logspace(5, 12, num=20), beta=beta, source_function=source_function)
    if MPI.COMM_WORLD.rank == 0:
        np.savetxt(f"{source_function.name.lower()}_{beta:.1f}.dat", data)
