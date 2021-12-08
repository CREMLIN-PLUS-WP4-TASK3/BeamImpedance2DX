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
circle_m_v = bi2d.Mesh("circle_vacuum.xdmf")
circle_m_vm = bi2d.Mesh("circle_vacuum_metal.xdmf")
circle_m_vmv = bi2d.Mesh("circle_vacuum_metal_vacuum.xdmf")
circle_mc_v = bi2d.MaterialMap(circle_m_v, [beam, vacuum])
circle_mc_vm = bi2d.MaterialMap(circle_m_vm, [beam, vacuum, steel])
circle_mc_vmv = bi2d.MaterialMap(circle_m_vmv, [beam, vacuum, steel, outer_vacuum])
ellipse_m_v = bi2d.Mesh("ellipse_vacuum.xdmf")
ellipse_m_vm = bi2d.Mesh("ellipse_vacuum_metal.xdmf")
ellipse_m_vmv = bi2d.Mesh("ellipse_vacuum_metal_vacuum.xdmf")
ellipse_mc_v = bi2d.MaterialMap(ellipse_m_v, [beam, vacuum])
ellipse_mc_vm = bi2d.MaterialMap(ellipse_m_vm, [beam, vacuum, steel])
ellipse_mc_vmv = bi2d.MaterialMap(ellipse_m_vmv, [beam, vacuum, steel, outer_vacuum])

# Configure solution
solution_circle_v = bi2d.Solution(circle_mc_v, Hcurl_order=2, H1_order=2)
solution_circle_vm = bi2d.Solution(circle_mc_vm, Hcurl_order=2, H1_order=2)
solution_circle_vmv = bi2d.Solution(circle_mc_vmv, Hcurl_order=2, H1_order=2)
solution_circle_v.logger.setLevel(logging.INFO)
solution_circle_vm.logger.setLevel(logging.INFO)
solution_circle_vmv.logger.setLevel(logging.INFO)
solution_ellipse_v = bi2d.Solution(ellipse_mc_v, Hcurl_order=2, H1_order=2)
solution_ellipse_vm = bi2d.Solution(ellipse_mc_vm, Hcurl_order=2, H1_order=2)
solution_ellipse_vmv = bi2d.Solution(ellipse_mc_vmv, Hcurl_order=2, H1_order=2)
solution_ellipse_v.logger.setLevel(logging.INFO)
solution_ellipse_vm.logger.setLevel(logging.INFO)
solution_ellipse_vmv.logger.setLevel(logging.INFO)

# Calculate impedance in frequency range and save the results
for solution, name, f_range, sibc_boundaries, source_function, rotation \
    in [(solution_circle_vm, "circle_metal", np.logspace(3, 8, num=100), [], bi2d.SourceFunction.MONOPOLE, 0),
        (solution_circle_v, "circle_sibc", np.logspace(6, 12, num=100), [wall], bi2d.SourceFunction.MONOPOLE, 0),
        (solution_circle_vmv, "circle_metal", np.logspace(3, 8, num=100), [], bi2d.SourceFunction.DIPOLE, 0),
        (solution_circle_v, "circle_sibc", np.logspace(6, 12, num=100), [wall], bi2d.SourceFunction.DIPOLE, 0),
        (solution_circle_vmv, "circle_metal", np.logspace(3, 8, num=100), [], bi2d.SourceFunction.DIPOLE, np.pi/2),
        (solution_circle_v, "circle_sibc", np.logspace(6, 12, num=100), [wall], bi2d.SourceFunction.DIPOLE, np.pi/2),
        (solution_ellipse_vm, "ellipse_metal", np.logspace(3, 8, num=100), [], bi2d.SourceFunction.MONOPOLE, 0),
        (solution_ellipse_v, "ellipse_sibc", np.logspace(6, 12, num=100), [wall], bi2d.SourceFunction.MONOPOLE, 0),
        (solution_ellipse_vmv, "ellipse_metal", np.logspace(3, 8, num=100), [], bi2d.SourceFunction.DIPOLE, 0),
        (solution_ellipse_v, "ellipse_sibc", np.logspace(6, 12, num=100), [wall], bi2d.SourceFunction.DIPOLE, 0),
        (solution_ellipse_vmv, "ellipse_metal", np.logspace(3, 8, num=100), [], bi2d.SourceFunction.DIPOLE, np.pi/2),
        (solution_ellipse_v, "ellipse_sibc", np.logspace(6, 12, num=100), [wall], bi2d.SourceFunction.DIPOLE, np.pi/2)]:
    data = solution.get_z(f_range, beta=0.999999, source_function=source_function)
    if MPI.COMM_WORLD.rank == 0:
        np.savetxt(f"{source_function.name.lower()}_{name}_{rotation*180/np.pi:.0f}.dat", data)
