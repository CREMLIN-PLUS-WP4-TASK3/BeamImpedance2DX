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

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Setting poisson function")

        V = dolfinx.FunctionSpace(self.mesh.mesh, ufl.MixedElement(self.solution.H1,
                                                                   self.solution.H1))
        phi_re, phi_im = ufl.TrialFunctions(V)
        v_re, v_im = ufl.TestFunctions(V)

        omega = self.solution._omega
        v = self.solution._beta * self.solution.c0

        # $$\varepsilon_r \varepsilon_0 \beta c_0$$
        eps_v = self.material_map.eps * v

        # $$\frac{\omega^2 \varepsilon_r \varepsilon_0}{\beta c_0}$$
        omega2_eps_v = omega**2 * self.material_map.eps / v

        self._a_phi = 0

        # $$\varepsilon_0 \varepsilon_r \beta c_0\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Re \;d\Omega}$$
        self._a_phi += eps_v * inner(grad(v_re), grad(phi_re)) * dx
        # $$\frac{\omega^2\varepsilon_0 \varepsilon_r}{\beta c_0}\int_\Omega{v^\Re\Phi^\Re \;d\Omega}$$
        self._a_phi += omega2_eps_v * inner(v_re, phi_re) * dx
        # $$ \varepsilon_0 \varepsilon_r \beta c_0\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Im \;d\Omega}$$
        self._a_phi += eps_v * inner(grad(v_im), grad(phi_im)) * dx
        # $$\frac{\omega^2\varepsilon_0 \varepsilon_r}{\beta c_0}\int_\Omega{v^\Im\Phi^\Im \;d\Omega}$$
        self._a_phi += omega2_eps_v * inner(v_im, phi_im) * dx

        if self.material_map.sigma is not None:
            # $$\frac{\sigma \beta c_0}{\omega}$$
            sigma_v_omega = self.material_map.sigma * v / omega

            # $$\frac{\omega \sigma}{\beta c_0}$$
            omega_sigma_v = omega * self.material_map.sigma / v

            # $$\frac{\sigma \beta c_0}{\omega}\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Im \;d\Omega}$$
            self._a_phi += sigma_v_omega * inner(grad(v_re), grad(phi_im)) * dx
            # $$\frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Re\Phi^\Im \;d\Omega}$$
            self._a_phi += omega_sigma_v * inner(v_re, phi_im) * dx
            # $$-\frac{\sigma \beta c_0}{\omega}\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Re \;d\Omega}$$
            self._a_phi += -sigma_v_omega * inner(grad(v_im), grad(phi_re)) * dx
            # $$-\frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Im\Phi^\Re \;d\Omega}$$
            self._a_phi += -omega_sigma_v * inner(v_im, phi_re) * dx

        # $$\int_\Omega{v^\Re J_s \;d\Omega}$$
        self._L_phi = inner(v_re, self.solution.Js) * dx

        self._A_phi = dolfinx.fem.create_matrix(self._a_phi)
        self._b_phi = dolfinx.fem.create_vector(self._L_phi)

        u_bc = dolfinx.Function(V)
        bc_facets = np.where(np.array(dolfinx.cpp.mesh.compute_boundary_facets(self.mesh.mesh.topology)) == 1)[0]
        bc_dofs = dolfinx.fem.locate_dofs_topological(V, self.mesh.mesh.topology.dim-1, bc_facets)
        with u_bc.vector.localForm() as loc:
            loc.setValues(bc_dofs, np.full(bc_dofs.size, 0))
        u_bc.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        self._bc = dolfinx.DirichletBC(u_bc, bc_dofs)

        self._phi = dolfinx.Function(V)
        self.solution.Phi_re, self.solution.Phi_im = self._phi.split()
        # FIXME: workaround for https://github.com/FEniCS/dolfinx/issues/1577
        Phi_re, Phi_im = ufl.split(self._phi)

        Vperp = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.Hcurl)
        # Vperp = dolfinx.FunctionSpace(self.mesh.mesh, ufl.VectorElement(self.solution.H1, dim=2))
        Vz = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1)

        for name, expr, V in [("Ediv_perp_re", -grad(Phi_re), Vperp),
                              ("Ediv_perp_im", -grad(Phi_im), Vperp),
                              ("Ediv_z_re", -omega / v * Phi_im, Vz),
                              ("Ediv_z_im", omega / v * Phi_re, Vz)]:
            u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
            a_p = inner(v, u) * dx
            L_p = inner(v, expr) * dx
            A = dolfinx.fem.create_matrix(a_p)
            b = dolfinx.fem.create_vector(L_p)
            setattr(self, f"_a_{name}", a_p)
            setattr(self, f"_L_{name}", L_p)
            setattr(self, f"_A_{name}", A)
            setattr(self, f"_b_{name}", b)
            # setattr(self, f"_bc_{name}", self.__set_bc(V, bcs))
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
                             self._b_phi, self._phi, bcs=[self._bc],
                             petsc_options=petsc_options)

        for name in ["Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im"]:
            a_p = getattr(self, f"_a_{name}")
            L_p = getattr(self, f"_L_{name}")
            A = getattr(self, f"_A_{name}")
            b = getattr(self, f"_b_{name}")
            # bc = getattr(self, f"_bc_{name}")
            f = getattr(self.solution, name)

            self.solution._solve(a_p, L_p, A, b, f, petsc_options=petsc_options)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solved poisson function")

        self.solution._Ediv_stale = False
