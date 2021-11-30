from bi2d.materials import vacuum, beam
from itertools import product
from mpi4py import MPI
import bi2d
import logging
import numpy as np

beam.index = 1
vacuum.index = 2
carbon = bi2d.Material(3, sigma=1e4)
wall = carbon.copy()
wall.index = -1

# Setup mesh and materials
m_sibc = bi2d.Mesh("collimator_sibc.xdmf")
mc_sibc = bi2d.MaterialMap(m_sibc, [beam, vacuum])
m_metal = bi2d.Mesh("collimator_metal.xdmf")
mc_metal = bi2d.MaterialMap(m_metal, [beam, vacuum, carbon])

solution_sibc = bi2d.Solution(mc_sibc, Hcurl_order=2, H1_order=2)
solution_sibc.logger.setLevel(logging.INFO)
solution_metal = bi2d.Solution(mc_metal, Hcurl_order=2, H1_order=2)
solution_metal.logger.setLevel(logging.INFO)

## Calculate impedance

# Calculate impedance in frequency range and save the results
for (solution, name, f_range, sibc_boundaries), source_function, rotation \
    in product([(solution_metal, "metal", np.logspace(2, 7, num=100), []),
                (solution_sibc, "sibc", np.logspace(6, 11, num=100), [wall])],
               [bi2d.SourceFunction.MONOPOLE, bi2d.SourceFunction.DIPOLE],
               [0, np.pi/2]):
    data = solution.get_z(f_range, beta=0.999999, rotation=rotation, source_function=source_function)
    if MPI.COMM_WORLD.rank == 0:
        np.savetxt(f"{source_function.name.lower()}_{name}_{rotation*180/np.pi:.0f}.dat", data)
