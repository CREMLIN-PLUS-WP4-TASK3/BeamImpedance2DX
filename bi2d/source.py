"""Source definition."""

from enum import Enum
import dolfinx
import ufl
from ufl import inner, dx
import numpy as np
from mpi4py import MPI


class SourceFunction(Enum):
    """Source function type."""

    MONOPOLE = 1
    MONOPOLE_CONSTANT = 1
    MONOPOLE_GAUSSIAN = 2
    DIPOLE = 3
    DIPOLE_RING_LINEAR = 3


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
        x0 = (minx + maxx) / 2
        y0 = (miny + maxy) / 2
        sigma = (maxx - minx) / 2 / 2
        x = ufl.SpatialCoordinate(function_space)
        A = 1 / (2 * np.pi * sigma**2)
        return A * ufl.exp(- (x[0] - x0)**2 / (2 * sigma**2)
                           - (x[1] - y0)**2 / (2 * sigma**2)) * \
                           test_function * dx

    def __source_function_dipole_ring_linear(self, function_space, test_function):
        self.solution._monopole = False
        (minx, miny), (maxx, maxy) = self.mesh.get_limits(self.material_map.beam_index)
        x0 = (minx + maxx) / 2
        y0 = (miny + maxy) / 2
        R = (maxx - minx) / 2
        # We are assuming DOFs located on triangle vertices. Thus, this function is only valid in linear H1 space.
        V = dolfinx.FunctionSpace(self.mesh.mesh, ("Lagrange", 1))
        circle = dolfinx.Function(V)
        dofs = dolfinx.fem.locate_dofs_geometrical(
            V,
            lambda x: np.isclose((x[0] - x0)**2+(x[1] - y0)**2, R**2)
        )
        if dofs.size == 0:
            raise ValueError("Cannot set source function. No vertices on beam region boundary.")
        with circle.vector.localForm() as loc:
            loc.set(0)
            loc.setValues(dofs, np.full(dofs.size, 1))

        x = ufl.SpatialCoordinate(V)
        r = ufl.sqrt((x[0]-x0)**2 + (x[1]-y0)**2)
        theta = ufl.atan_2(x[1]-x0, x[0]-y0)
        d = r * ufl.sin(theta - self.rotation)
        function = ufl.sin(theta - self.rotation) / R * circle
        A = dolfinx.fem.assemble_scalar(d * function * dx)
        A = MPI.COMM_WORLD.gather(A)
        if MPI.COMM_WORLD.rank == 0:
            A = np.abs(np.sum(A))
        A = MPI.COMM_WORLD.bcast(A)
        return function / A * test_function * dx

    def __integral_function_monopole(self, _):
        q = dolfinx.fem.assemble_scalar(self.solution.Js * ufl.dx)
        q = MPI.COMM_WORLD.gather(q)
        if MPI.COMM_WORLD.rank == 0:
            q = np.sum(q)
        self.solution.q = MPI.COMM_WORLD.bcast(q)

    def __integral_function_dipole(self, function_space):
        (minx, miny), (maxx, maxy) = self.mesh.get_limits(self.material_map.beam_index)
        x0 = (minx + maxx) / 2
        y0 = (miny + maxy) / 2
        x = ufl.SpatialCoordinate(function_space)
        r = ufl.sqrt((x[0]-x0)**2 + (x[1]-y0)**2)
        theta = ufl.atan_2(x[1], x[0])
        d = r * ufl.sin(theta - self.rotation)
        A = dolfinx.fem.assemble_scalar(d * self.solution.Js * dx)
        A = MPI.COMM_WORLD.gather(A)
        if MPI.COMM_WORLD.rank == 0:
            A = np.abs(np.sum(A))
        self.solution.dm = MPI.COMM_WORLD.bcast(A)

    def __init__(self, solution, rotation=0, source_function=SourceFunction.MONOPOLE):
        """Initialize."""
        self.source_functions = {
            SourceFunction.MONOPOLE_CONSTANT: self.__source_function_monopole_constant,
            SourceFunction.MONOPOLE_GAUSSIAN: self.__source_function_monopole_gaussian,
            SourceFunction.DIPOLE_RING_LINEAR: self.__source_function_dipole_ring_linear,
        }
        self.integral_functions = {
            SourceFunction.MONOPOLE_CONSTANT: self.__integral_function_monopole,
            SourceFunction.MONOPOLE_GAUSSIAN: self.__integral_function_monopole,
            SourceFunction.DIPOLE_RING_LINEAR: self.__integral_function_dipole,
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

        self.integral_functions[self.source_function](self._V)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solved source function")

        self.solution._Js_stale = False
