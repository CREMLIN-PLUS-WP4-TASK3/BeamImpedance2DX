"""Poisson definition."""

import dolfinx
import ufl
from ufl import inner, grad, dx
import numpy as np
from petsc4py import PETSc
from mpi4py import MPI


class Ediv():
    """Irrotational electric field solver."""

    def __init__(self, solution):
        """Initialize."""
        self.material_map = solution.material_map
        self.mesh = solution.mesh
        self.solution = solution

        if self.solution.Js is None:
            raise AttributeError("Source function Js is not initialized")

        self.solution._Ediv_stale = True
        self.solution._Ecurl_stale = True

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Setting poisson function")

        self._a_phi = 0
        self._L_phi = 0

        omega = self.solution._omega
        beta = self.solution._beta
        c0 = self.solution.c0
        eps = self.material_map.eps

        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            V = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1)
            phi, vv = ufl.TrialFunction(V), ufl.TestFunction(V)

            # $$\underline{\varepsilon}\beta c_0\int_\Omega{\nabla_\perp v \cdot \nabla_\perp\Phi \;d\Omega}$$
            self._a_phi += eps * beta * c0 * inner(grad(phi), grad(vv)) * dx
            # $$\frac{\omega^2\underline{\varepsilon}}{\beta c_0}\int_\Omega{v\Phi \;d\Omega}$$
            self._a_phi += omega**2 * eps / (beta * c0) * inner(phi, vv) * dx
            # $$\int_\Omega{v J_s \;d\Omega}$$
            self._L_phi += inner(self.solution.Js, vv) * dx

        else:
            V = dolfinx.FunctionSpace(self.mesh.mesh, ufl.MixedElement(self.solution.H1,
                                                                       self.solution.H1))
            phi_re, phi_im = ufl.TrialFunctions(V)
            v_re, v_im = ufl.TestFunctions(V)

            r"""
            $$
            \varepsilon \beta c_0\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Re \;d\Omega}
            $$
            """
            self._a_phi += eps * beta * c0 * inner(grad(v_re), grad(phi_re)) * dx
            r"""
            $$
            \frac{\omega^2\varepsilon}{\beta c_0}\int_\Omega{v^\Re\Phi^\Re \;d\Omega}
            $$
            """
            self._a_phi += omega**2 * eps / (beta * c0) * inner(v_re, phi_re) * dx
            r"""
            $$
            \varepsilon \beta c_0\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Im \;d\Omega}
            $$
            """
            self._a_phi += eps * beta * c0 * inner(grad(v_im), grad(phi_im)) * dx
            r"""
            $$
            \frac{\omega^2\varepsilon}{\beta c_0}\int_\Omega{v^\Im\Phi^\Im \;d\Omega}
            $$
            """
            self._a_phi += omega**2 * eps / (beta * c0) * inner(v_im, phi_im) * dx

            if self.material_map.sigma is not None:
                sigma = self.material_map.sigma

                r"""
                $$
                \frac{\sigma \beta c_0}{\omega}\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Im \;d\Omega}
                $$
                """
                self._a_phi += sigma * beta * c0 / omega * inner(grad(v_re), grad(phi_im)) * dx
                r"""
                $$
                \frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Re\Phi^\Im \;d\Omega}
                $$
                """
                self._a_phi += omega * sigma / (beta * c0) * inner(v_re, phi_im) * dx
                r"""
                $$
                -\frac{\sigma \beta c_0}{\omega}\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Re \;d\Omega}
                $$
                """
                self._a_phi += -sigma * beta * c0 / omega * inner(grad(v_im), grad(phi_re)) * dx
                # $$-\frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Im\Phi^\Re \;d\Omega}$$
                self._a_phi += -omega * sigma / (beta * c0) * inner(v_im, phi_re) * dx

            # $$\int_\Omega{v^\Re J_s \;d\Omega}$$
            self._L_phi += inner(v_re, self.solution.Js) * dx

        self._A_phi = dolfinx.fem.create_matrix(self._a_phi)
        self._b_phi = dolfinx.fem.create_vector(self._L_phi)

        u_bc = dolfinx.Function(V)
        bc_facets = np.where(np.array(dolfinx.cpp.mesh.compute_boundary_facets(self.mesh.mesh.topology)) == 1)[0]
        bc_dofs = dolfinx.fem.locate_dofs_topological(V, self.mesh.mesh.topology.dim - 1, bc_facets)
        with u_bc.vector.localForm() as loc:
            loc.setValues(bc_dofs, np.full(bc_dofs.size, 0))
        u_bc.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        self._bc = dolfinx.DirichletBC(u_bc, bc_dofs)

        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            if type(self.solution._phi) != dolfinx.Function:
                self.solution.Phi = dolfinx.Function(V)
                self.solution._phi = self.solution.Phi

            Vperp = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.Hcurl)
            # Vperp = dolfinx.FunctionSpace(self.mesh.mesh, ufl.VectorElement(self.solution.H1, dim=2))
            Vz = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1)

            loop_list = [("Ediv_perp", -grad(self.solution._phi), Vperp),
                         ("Ediv_z", 1j * omega / (beta * c0) * self.solution._phi, Vz)]
        else:
            if type(self.solution._phi) != dolfinx.Function:
                self.solution._phi = dolfinx.Function(V)
                self.solution.Phi_re, self.solution.Phi_im = self.solution._phi.split()
                # FIXME: workaround for https://github.com/FEniCS/dolfinx/issues/1577
                Phi_re, Phi_im = ufl.split(self.solution._phi)

            Vperp = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.Hcurl)
            # Vperp = dolfinx.FunctionSpace(self.mesh.mesh, ufl.VectorElement(self.solution.H1, dim=2))
            Vz = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1)

            loop_list = [("Ediv_perp_re", -grad(Phi_re), Vperp),
                         ("Ediv_perp_im", -grad(Phi_im), Vperp),
                         ("Ediv_z_re", -omega / (beta * c0) * Phi_im, Vz),
                         ("Ediv_z_im", omega / (beta * c0) * Phi_re, Vz)]

        for name, expr, V in loop_list:
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
            self.solution.logger.debug("Set poisson function")

    def solve(self, petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}):
        """Solve equation."""
        if self.solution._Js_stale:
            raise ValueError("Source function Js solution is stale and needs to be recalculated first")

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solving poisson function")

        self.solution._solve(self._a_phi, self._L_phi, self._A_phi,
                             self._b_phi, self.solution._phi, bcs=[self._bc],
                             petsc_options=petsc_options)

        Efields = (["Ediv_perp", "Ediv_z"] if np.issubdtype(PETSc.ScalarType, np.complexfloating)
                   else ["Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im"])

        for name in Efields:
            a_p = getattr(self, f"_a_{name}")
            L_p = getattr(self, f"_L_{name}")
            A = getattr(self, f"_A_{name}")
            b = getattr(self, f"_b_{name}")
            f = getattr(self.solution, name)

            self.solution._solve(a_p, L_p, A, b, f, petsc_options=petsc_options)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solved poisson function")

        self.solution._Ediv_stale = False
