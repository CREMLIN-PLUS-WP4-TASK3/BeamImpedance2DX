"""Source definition."""

from enum import Enum
import dolfinx
import ufl
from ufl import inner, dx
import numpy as np
from petsc4py import PETSc
from mpi4py import MPI


class SourceFunction(Enum):
    """Source function type."""

    MONOPOLE = 1
    MONOPOLE_CONSTANT = 1
    MONOPOLE_GAUSSIAN = 2
    DIPOLE_SIN = 3


class Js():
    """Beam current."""

    def __source_function_monopole_constant(self, _, test_function):
        self.solution._monopole = True
        (minx, _), (maxx, _) = self.mesh.get_limits(self.material_map.beam_index)
        area = np.pi * (maxx - minx)**2 / 4
        return 1 / area * self.material_map.beam * test_function * dx

    def __source_function_monopole_gaussian(self, function_space, test_function):
        self.solution._monopole = True
        (minx, miny), (maxx, maxy) = self.mesh.get_limits(self.material_map.beam_index)
        centerx = (minx + maxx) / 2
        centery = (miny + maxy) / 2
        sigma = (maxx - minx) / 2 / 2
        x = ufl.SpatialCoordinate(function_space)
        A = 1 / (2 * np.pi * sigma**2)
        return A * ufl.exp(- (x[0] - centerx)**2 / (2 * sigma**2)
                           - (x[1] - centery)**2 / (2 * sigma**2)) * \
                           test_function * dx

    def __source_function_multipole_sin(self, function_space, test_function, N):
        self.solution._monopole = False
        (minx, miny), (maxx, maxy) = self.mesh.get_limits(self.material_map.beam_index)
        centerx = (minx + maxx) / 2
        centery = (miny + maxy) / 2
        R = (maxx - minx) / 2
        self.solution.beam_center = (centerx, centery)
        x = ufl.SpatialCoordinate(function_space)
        r = ufl.sqrt(x[0]**2+x[1]**2) / R
        theta = ufl.atan_2(x[1], x[0])
        # FIXME: Fix the normalized amplitude
        return ufl.sin(np.pi*r)*ufl.sin(N * theta + self.rotation) * \
            self.material_map.beam * test_function * dx

    def __source_function_dipole_sin(self, function_space, test_function):
        return self.__source_function_multipole_sin(function_space, test_function, 2)

    def __init__(self, solution, rotation=0, source_function=SourceFunction.MONOPOLE):
        """Initialize."""
        self.source_functions = {
            SourceFunction.MONOPOLE_CONSTANT: self.__source_function_monopole_constant,
            SourceFunction.MONOPOLE_GAUSSIAN: self.__source_function_monopole_gaussian,
            SourceFunction.DIPOLE_SIN: self.__source_function_dipole_sin,
        }
        self.solution = solution
        self.source_function = source_function
        self.material_map = solution.material_map
        self.mesh = solution.mesh
        self.rotation = rotation
        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Setting source function")
        self._V = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1)
        u, v = ufl.TrialFunction(self._V), ufl.TestFunction(self._V)
        self._a_p = inner(u, v) * dx
        self._L_p = self.source_functions[self.source_function](self._V, v)
        self.solution.Js = dolfinx.Function(self._V)
        self._A = dolfinx.fem.create_matrix(self._a_p)
        self._b = dolfinx.fem.create_vector(self._L_p)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Set source function")

    def solve(self, petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}):
        """Solve equation."""

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solving source function")

        self.solution._solve(self._a_p, self._L_p, self._A, self._b, self.solution.Js,
                             petsc_options=petsc_options)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solved source function")

        self.solution._Js_stale = False
