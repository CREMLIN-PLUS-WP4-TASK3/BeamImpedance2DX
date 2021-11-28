from bi2d.materials import vacuum, beam, mischung43
from mpi4py import MPI
import bi2d
import logging
import numpy as np

# Assign material indices
beam.index = 1
vacuum.index = 2
mischung43.index = 3

# Setup mesh and materials
m = bi2d.Mesh("mesh.xdmf")
mc = bi2d.MaterialMap(m, [beam, vacuum, mischung43])

# Configure solution
solution = bi2d.Solution(mc, Hcurl_order=2, H1_order=2)
# Enable info level logging
solution.logger.setLevel(logging.INFO)

betas = [0.01, 0.05, 0.1, 0.5, 1.0]
for beta in betas:
    data = solution.get_z(np.logspace(5, 9, num=50), beta=beta, source_function=bi2d.SourceFunction.MONOPOLE)
    data[:, (1, 2)] *= 0.0254
    if MPI.COMM_WORLD.rank == 0:
        np.savetxt(f"monopole_{beta:.2f}.dat", data)
