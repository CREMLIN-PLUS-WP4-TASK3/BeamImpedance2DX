"""Source definition."""

from enum import Enum
import dolfinx
import ufl
import numpy as np


class SourceFunction(Enum):
    """Source function type."""

    MONOPOLE = 1
    MONOPOLE_CONSTANT = 1
    MONOPOLE_GAUSSIAN = 2
    DIPOLE = 3


class Js():
    """Beam current."""

    def __source_function_monopole_constant(self, _, test_function):
        (minx, _), (maxx, _), _ = self.mesh.get_limits(self.material_map.beam_index)
        area = np.pi * (maxx - minx)**2 / 4
        return 1 / area * self.material_map.beam * test_function * ufl.dx

    def __source_function_monopole_gaussian(self, function_space, test_function):
        (minx, miny), (maxx, maxy), _ = self.mesh.get_limits(self.material_map.beam_index)
        centerx = (minx + maxx) / 2
        centery = (miny + maxy) / 2
        sigma = (maxx - minx) / 2 / 2
        x = ufl.SpatialCoordinate(function_space)
        A = 1 / (2 * np.pi * sigma**2)
        return A * ufl.exp(-(x[0]-centerx)**2 / (2 * sigma**2)
                           -(x[1]-centery)**2 / (2 * sigma**2)) * test_function * ufl.dx

    def __source_function_dipole(self, function_space, test_function):
        raise NotImplementedError()

    def __init__(self, solution, source_function=SourceFunction.MONOPOLE):
        """Initialize."""
        self.source_functions = {
            SourceFunction.MONOPOLE_CONSTANT: self.__source_function_monopole_constant,
            SourceFunction.MONOPOLE_GAUSSIAN: self.__source_function_monopole_gaussian,
            SourceFunction.DIPOLE: self.__source_function_dipole,
        }
        self.solution = solution
        self.source_function = source_function
        self.material_map = solution.material_map
        self.mesh = solution.mesh

    def solve(self):
        """Solve equation."""
        V = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1)
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
        a_p = ufl.inner(u, v) * ufl.dx
        L_p = self.source_functions[self.source_function](V, v)
        self.solution.Js = dolfinx.Function(V)
        projection = dolfinx.fem.LinearProblem(a_p, L_p, u=self.solution.Js)
        projection.solve()
        self.solution.q = dolfinx.fem.assemble_scalar(self.solution.Js * ufl.dx)
