from bi2d.materials import vacuum, beam
from itertools import product
from mpi4py import MPI
import bi2d
import logging
import numpy as np

# Assign material indices. When setting index `-1` to SIBC material, boundary will be applied to all boundaries
beam.index = 1
vacuum.index = 2
steel = bi2d.Material(3, sigma=1e6)
wall = steel.copy()
wall.index = -1
vacuum_outer = vacuum.copy()
vacuum_outer.index = 4

# Setup mesh and materials
circle_sibc_mesh = bi2d.Mesh("circle_sibc.xdmf")
circle_metal_mesh = bi2d.Mesh("circle_metal.xdmf")
ellipse_sibc_mesh = bi2d.Mesh("ellipse_sibc.xdmf")
ellipse_metal_mesh = bi2d.Mesh("ellipse_metal.xdmf")
circle_sibc_map = bi2d.MaterialMap(circle_sibc_mesh, [beam, vacuum])
circle_metal_map = bi2d.MaterialMap(circle_metal_mesh, [beam, vacuum, steel, vacuum_outer])
ellipse_sibc_map = bi2d.MaterialMap(ellipse_sibc_mesh, [beam, vacuum])
ellipse_metal_map = bi2d.MaterialMap(ellipse_metal_mesh, [beam, vacuum, steel, vacuum_outer])

# Configure solution
circle_sibc_solution = bi2d.Solution(circle_sibc_map, Hcurl_order=2, H1_order=2)
circle_sibc_solution.logger.setLevel(logging.INFO)
circle_metal_solution = bi2d.Solution(circle_metal_map, Hcurl_order=2, H1_order=2)
circle_metal_solution.logger.setLevel(logging.INFO)
ellipse_sibc_solution = bi2d.Solution(ellipse_sibc_map, Hcurl_order=2, H1_order=2)
ellipse_sibc_solution.logger.setLevel(logging.INFO)
ellipse_metal_solution = bi2d.Solution(ellipse_metal_map, Hcurl_order=2, H1_order=2)
ellipse_metal_solution.logger.setLevel(logging.INFO)

# Calculate impedance in frequency range and save the results
for (solution, name, f_range, sibc_boundaries), source_function, rotation \
    in product([(circle_metal_solution, "circle_metal", np.logspace(3, 8, num=100), []),
                (circle_sibc_solution, "circle_sibc", np.logspace(6, 12, num=100), [wall]),
                (ellipse_metal_solution, "ellipse_metal", np.logspace(3, 8, num=100), []),
                (ellipse_sibc_solution, "ellipse_sibc", np.logspace(6, 12, num=100), [wall])],
               [bi2d.SourceFunction.MONOPOLE, bi2d.SourceFunction.DIPOLE],
               [0, np.pi/2]):
    data = solution.get_z(f_range, beta=0.999999, rotation=rotation ,source_function=source_function)
    if MPI.COMM_WORLD.rank == 0:
        np.savetxt(f"{source_function.name.lower()}_{name}_{rotation*180/np.pi:.0f}.dat", data)
