"""Curl definition."""

from enum import Enum, auto
import dolfinx
import ufl
from ufl import inner, dx
import numpy as np
from petsc4py import PETSc
from mpi4py import MPI


class Esum():
    """Total field."""

    def __init__(self, solution):
        """Initialize."""
        self.solution = solution
        self.mesh = self.solution.mesh

        for name in ["Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im",
                     "Ecurl_perp_re", "Ecurl_perp_im", "Ecurl_z_re", "Ecurl_z_im"]:
            attr = getattr(self.solution, name)
            if attr is None:
                raise AttributeError(f"{name} function is not available")

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Setting field summator")

        Vperp = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.Hcurl_vis)
        Vz = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1_vis)

        for name, expr, V in [("Eperp_re", self.solution.Ediv_perp_re+self.solution.ecurl_perp_re, Vperp),
                              ("Eperp_im", self.solution.Ediv_perp_im+self.solution.ecurl_perp_im, Vperp),
                              ("Ez_re", self.solution.Ediv_z_re+self.solution.ecurl_z_re, Vz),
                              ("Ez_im", self.solution.Ediv_z_im+self.solution.ecurl_z_im, Vz)]:
            u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
            a_p = inner(v, u) * dx
            L_p = inner(v, expr) * dx
            A = dolfinx.fem.create_matrix(a_p)
            b = dolfinx.fem.create_vector(L_p)
            setattr(self, f"_a_{name}", a_p)
            setattr(self, f"_L_{name}", L_p)
            setattr(self, f"_A_{name}", A)
            setattr(self, f"_b_{name}", b)
            setattr(self.solution, name, dolfinx.Function(V))

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Set field summator")

    def solve(self, petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}):
        """Solve equation."""

        if self.solution._Ediv_stale:
            raise ValueError("Ediv solutions are stale and need to be recalculated")
        if self.solution._Ecurl_stale:
            raise ValueError("Ecurl solutions are stale and need to be recalculated")

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Prepare field summation")

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Summing fields")

        for name in ["Eperp_re", "Eperp_im", "Ez_re", "Ez_im"]:
            a_p = getattr(self, f"_a_{name}")
            L_p = getattr(self, f"_L_{name}")
            A = getattr(self, f"_A_{name}")
            b = getattr(self, f"_b_{name}")
            f = getattr(self.solution, name)

            self.solution._solve(a_p, L_p, A, b, f, petsc_options=petsc_options)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Fields summed")

        self.solution._E_stale = False
