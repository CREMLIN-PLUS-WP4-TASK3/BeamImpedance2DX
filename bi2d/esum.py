"""Curl definition."""

import dolfinx
import ufl
import numpy as np
from ufl import inner, dx
from mpi4py import MPI
from petsc4py import PETSc


class Esum():
    """Total field."""

    def __init__(self, solution, H1_order=2, Hcurl_order=2):
        """Initialize."""
        self.solution = solution
        self.mesh = self.solution.mesh

        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            for name in ["Ediv_perp", "Ediv_z", "Ecurl_perp", "Ecurl_z"]:
                attr = getattr(self.solution, name)
                if attr is None:
                    raise AttributeError(f"{name} function is not available")
        else:
            for name in ["Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im",
                         "Ecurl_perp_re", "Ecurl_perp_im", "Ecurl_z_re", "Ecurl_z_im"]:
                attr = getattr(self.solution, name)
                if attr is None:
                    raise AttributeError(f"{name} function is not available")

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Setting field summator")

        H1 = ufl.FiniteElement("Lagrange", self.mesh.mesh.ufl_cell(), H1_order)
        Hcurl = ufl.VectorElement(ufl.FiniteElement("Lagrange", self.mesh.mesh.ufl_cell(), Hcurl_order), dim=2)
        Vperp = dolfinx.FunctionSpace(self.mesh.mesh, Hcurl)
        Vz = dolfinx.FunctionSpace(self.mesh.mesh, H1)

        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            expr_list = [("Eperp", self.solution.Ediv_perp + self.solution.ecurl_perp, Vperp),
                         ("Ez", self.solution.Ediv_z + self.solution.ecurl_z, Vz)]
        else:
            expr_list = [("Eperp_re", self.solution.Ediv_perp_re + self.solution.ecurl_perp_re, Vperp),
                         ("Eperp_im", self.solution.Ediv_perp_im + self.solution.ecurl_perp_im, Vperp),
                         ("Ez_re", self.solution.Ediv_z_re + self.solution.ecurl_z_re, Vz),
                         ("Ez_im", self.solution.Ediv_z_im + self.solution.ecurl_z_im, Vz)]

        for name, expr, V in expr_list:
            u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
            a_p = inner(u, v) * dx
            L_p = inner(expr, v) * dx
            A = dolfinx.fem.create_matrix(a_p)
            b = dolfinx.fem.create_vector(L_p)
            setattr(self, f"_a_{name}", a_p)
            setattr(self, f"_L_{name}", L_p)
            setattr(self, f"_A_{name}", A)
            setattr(self, f"_b_{name}", b)
            if type(getattr(self.solution, name)) != dolfinx.Function:
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

        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            names = ["Eperp", "Ez"]
        else:
            names = ["Eperp_re", "Eperp_im", "Ez_re", "Ez_im"]

        for name in names:
            a_p = getattr(self, f"_a_{name}")
            L_p = getattr(self, f"_L_{name}")
            A = getattr(self, f"_A_{name}")
            b = getattr(self, f"_b_{name}")
            f = getattr(self.solution, name)

            self.solution._solve(a_p, L_p, A, b, f, petsc_options=petsc_options)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Fields summed")

        self.solution._E_stale = False
