"""Source definition."""

from enum import Enum
import dolfinx
import ufl
from ufl import inner, dx
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


class SourceFunction(Enum):
    """Source function type."""

    MONOPOLE = 1
    MONOPOLE_CONSTANT = 1
    MONOPOLE_GAUSSIAN = 2
    DIPOLE = 3
    DIPOLE_RING_SIN = 3
    QUADRUPOLE = 4
    QUADRUPOLE_RING_SIN = 4


class Js():
    """Beam current."""

    def __source_function_monopole_constant(self, _, test_function):
        self.solution.source_function = SourceFunction.MONOPOLE
        (minx, _), (maxx, _) = self.mesh.get_limits(self.material_map.beam_index)
        area = np.pi * (maxx - minx)**2 / 4
        return inner(1 / area * self.material_map.beam, test_function) * dx

    def __source_function_monopole_gaussian(self, function_space, test_function):
        self.solution.source_function = SourceFunction.MONOPOLE
        (minx, miny), (maxx, maxy) = self.mesh.get_limits(self.material_map.beam_index)
        x0 = (minx + maxx) / 2
        y0 = (miny + maxy) / 2
        sigma = (maxx - minx) / 2 / 2
        x = ufl.SpatialCoordinate(function_space)
        A = 1 / (2 * np.pi * sigma**2)
        return inner(A * ufl.exp(- (x[0] - x0)**2 / (2 * sigma**2)
                                 - (x[1] - y0)**2 / (2 * sigma**2)), test_function) * dx

    def __source_function_dipole_ring_linear(self, function_space, test_function):
        self.solution.source_function = SourceFunction.DIPOLE
        (minx, miny), (maxx, maxy) = self.mesh.get_limits(self.material_map.beam_index)
        x0 = (minx + maxx) / 2
        y0 = (miny + maxy) / 2
        R = (maxx - minx) / 2
        # We are assuming DOFs located on triangle vertices. Thus, this function is only valid in linear H1 space.
        V = dolfinx.FunctionSpace(self.mesh.mesh, ("Lagrange", 1))
        circle = dolfinx.Function(V)
        func = dolfinx.Function(function_space)
        d = dolfinx.Function(function_space)
        dofs = dolfinx.fem.locate_dofs_geometrical(
            V,
            lambda x: np.isclose((x[0] - x0)**2 + (x[1] - y0)**2, R**2)
        )
        if dofs.size == 0:
            raise ValueError("Cannot set source function. No vertices on beam region boundary.")
        with circle.vector.localForm() as loc:
            loc.set(0)
            loc.setValues(dofs, np.full(dofs.size, 1))

        def sinfunc(arg):
            x = np.real(arg[0])
            y = np.real(arg[1])
            theta = np.arctan2(x - x0, y - y0)
            return np.sin(theta + self.rotation) / R

        def dfunc(arg):
            x = np.real(arg[0])
            y = np.real(arg[1])
            r = np.sqrt((x - x0)**2 + (y - y0)**2)
            theta = np.arctan2(x - x0, y - y0)
            return r * np.sin(theta + self.rotation)

        func.interpolate(sinfunc)
        func.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        d.interpolate(dfunc)
        d.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        A = dolfinx.fem.assemble_scalar(func * d * circle * dx)
        A = MPI.COMM_WORLD.gather(A)
        if MPI.COMM_WORLD.rank == 0:
            A = np.sum(A)
        A = MPI.COMM_WORLD.bcast(A)
        return inner(func * circle / A, test_function) * dx

    def __source_function_quadrupole_ring_linear(self, function_space, test_function):
        self.solution.source_function = SourceFunction.QUADRUPOLE
        (minx, miny), (maxx, maxy) = self.mesh.get_limits(self.material_map.beam_index)
        x0 = (minx + maxx) / 2
        y0 = (miny + maxy) / 2
        R = (maxx - minx) / 2
        # We are assuming DOFs located on triangle vertices. Thus, this function is only valid in linear H1 space.
        V = dolfinx.FunctionSpace(self.mesh.mesh, ("Lagrange", 1))
        circle = dolfinx.Function(V)
        func = dolfinx.Function(function_space)
        d = dolfinx.Function(function_space)
        dofs = dolfinx.fem.locate_dofs_geometrical(
            V,
            lambda x: np.isclose((x[0] - x0)**2 + (x[1] - y0)**2, R**2)
        )
        if dofs.size == 0:
            raise ValueError("Cannot set source function. No vertices on beam region boundary.")
        with circle.vector.localForm() as loc:
            loc.set(0)
            loc.setValues(dofs, np.full(dofs.size, 1))

        def sinfunc(arg):
            x = np.real(arg[0])
            y = np.real(arg[1])
            theta = np.arctan2(x - x0, y - y0)
            return np.sin(2 * theta + np.pi / 2 - self.rotation) / R

        def dfunc(arg):
            x = np.real(arg[0])
            y = np.real(arg[1])
            r = np.sqrt((x - x0)**2 + (y - y0)**2)
            theta = np.arctan2(x - x0, y - y0)
            return r * np.sin(2 * theta + np.pi / 2 - self.rotation)

        func.interpolate(sinfunc)
        func.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        d.interpolate(dfunc)
        d.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        A = dolfinx.fem.assemble_scalar(func * d * circle * dx)
        A = MPI.COMM_WORLD.gather(A)
        if MPI.COMM_WORLD.rank == 0:
            A = np.sum(A)
        A = MPI.COMM_WORLD.bcast(A)
        return inner(func * circle / A, test_function) * dx

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

        d = dolfinx.Function(function_space)

        def dfunc(arg):
            x = np.real(arg[0])
            y = np.real(arg[1])
            r = np.sqrt((x - x0)**2 + (y - y0)**2)
            theta = np.arctan2(x - x0, y - y0)
            return r * np.sin(theta + self.rotation)

        d.interpolate(dfunc)
        d.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        A = dolfinx.fem.assemble_scalar(d * self.solution.Js * dx)
        A = MPI.COMM_WORLD.gather(A)
        if MPI.COMM_WORLD.rank == 0:
            A = np.sum(A)
        self.solution.dm = MPI.COMM_WORLD.bcast(A)

    def __integral_function_quadrupole(self, function_space):
        raise NotImplementedError("Not implemented")

    def __init__(self, solution, rotation=0, source_function=SourceFunction.MONOPOLE):
        """Initialize."""
        self.source_functions = {
            SourceFunction.MONOPOLE_CONSTANT: self.__source_function_monopole_constant,
            SourceFunction.MONOPOLE_GAUSSIAN: self.__source_function_monopole_gaussian,
            SourceFunction.DIPOLE_RING_SIN: self.__source_function_dipole_ring_linear,
            SourceFunction.QUADRUPOLE_RING_SIN: self.__source_function_quadrupole_ring_linear,
        }
        self.integral_functions = {
            SourceFunction.MONOPOLE_CONSTANT: self.__integral_function_monopole,
            SourceFunction.MONOPOLE_GAUSSIAN: self.__integral_function_monopole,
            SourceFunction.DIPOLE_RING_SIN: self.__integral_function_dipole,
            SourceFunction.QUADRUPOLE_RING_SIN: self.__integral_function_quadrupole,
        }
        self.solution = solution
        self.source_function = source_function
        self.material_map = solution.material_map
        self.mesh = solution.mesh
        self.rotation = rotation
        self.solution._Js_stale = True
        self.solution._Ediv_stale = True
        self.solution._Ecurl_stale = True
        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Setting source function")
        self._V = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1)
        u, v = ufl.TrialFunction(self._V), ufl.TestFunction(self._V)
        self._a_p = inner(u, v) * dx
        self._L_p = self.source_functions[self.source_function](self._V, v)
        if type(self.solution.Js) != dolfinx.Function:
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
