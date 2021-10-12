#!/usr/bin/env python3

import logging
import bi2d
from bi2d.esum import Esum
from bi2d.materials import vacuum, beam, steel_isotropic
import numpy as np
from mpi4py import MPI

# Import materials
beam.index = 1
vacuum.index = 2
vacuum2 = vacuum.copy()
vacuum2.index = 3
# steel_isotropic.index = 3


# Setup mesh and materials
m = bi2d.Mesh("mesh.xdmf")
mc = bi2d.MaterialMap(m, [beam, vacuum, vacuum2])
# Configure solution
solution = bi2d.Solution(mc, Hcurl_order=2 H1_order=2)
# Enable logging
solution.logger.setLevel(logging.INFO)

# Solve and visualize the fields for one frequency point
print(solution.get_z([1e5], beta=0.1, source_function=bi2d.SourceFunction.DIPOLE))
esum = Esum(solution)  # Sum rotational and irrotational fields
esum.solve()
solution.save("solution.xdmf")

# Calculate impedance in frequency range and save the results
# data = solution.get_z(np.logspace(5, 12, num=3), beta=0.1, source_function=bi2d.SourceFunction.MONOPOLE_GAUSSIAN)
# if MPI.COMM_WORLD.rank == 0:
#     np.savetxt("out.csv", data)
