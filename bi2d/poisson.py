"""Poisson definition."""

import dolfinx
import ufl
from ufl import inner, grad, dx
import numpy as np


class Ediv():
    """Irrotational electric field solver."""

    def __init__(self, solution, boundary_index=None):
        """Initialize."""
        self.material_map = solution.material_map
        self.mesh = solution.mesh
        self.solution = solution
        self.boundary_index = boundary_index

    def __set_bc(self, V, value=0):
        u_bc = dolfinx.Function(V)
        if self.boundary_index is not None:
            bc_dofs = dolfinx.fem.locate_dofs_topological(V,
                                                          self.mesh.mesh.topology.dim-1,
                                                          self.mesh.get_boundary(self.boundary_index))
        else:
            bc_facets = np.where(
                np.array(dolfinx.cpp.mesh.compute_boundary_facets(self.mesh.mesh.topology)) == 1)[0]
            bc_dofs = dolfinx.fem.locate_dofs_topological(V,
                                                          self.mesh.mesh.topology.dim-1,
                                                          bc_facets)

        with u_bc.vector.localForm() as loc:
            loc.setValues(bc_dofs, np.full(len(bc_dofs), value))
        return dolfinx.DirichletBC(u_bc, bc_dofs)

    def solve(self):
        """Solve equation."""
        if self.solution.Js is None:
            raise ValueError("Term Js is not available")

        V = dolfinx.FunctionSpace(self.mesh.mesh, ufl.MixedElement(self.solution.H1,
                                                                   self.solution.H1))

        phi_re, phi_im = ufl.TrialFunctions(V)
        v_re, v_im = ufl.TestFunctions(V)

        omega = 2 * np.pi * self.solution.f
        v = self.solution.beta * self.solution.c0

        # $$\varepsilon_r \varepsilon_0 \beta c_0$$
        eps_v = dolfinx.Function(self.material_map.eps.function_space)
        with self.material_map.eps.vector.localForm() as eps, eps_v.vector.localForm() as vec:
            vec += eps * v

        # $$\frac{\omega^2 \varepsilon_r \varepsilon_0}{\beta c_0}$$
        omega2_eps_v = dolfinx.Function(self.material_map.eps.function_space)
        with self.material_map.eps.vector.localForm() as eps, omega2_eps_v.vector.localForm() as vec:
            vec += omega**2 * eps / v

        a_p = 0

        # $$\varepsilon_0 \varepsilon_r \beta c_0\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Re \;d\Omega}$$
        a_p += eps_v * inner(grad(v_re), grad(phi_re)) * dx
        # $$\frac{\omega^2\varepsilon_0 \varepsilon_r}{\beta c_0}\int_\Omega{v^\Re\Phi^\Re \;d\Omega}$$
        a_p += omega2_eps_v * inner(v_re, phi_re) * dx
        # $$ \varepsilon_0 \varepsilon_r \beta c_0\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Im \;d\Omega}$$
        a_p += eps_v * inner(grad(v_im), grad(phi_im)) * dx
        # $$\frac{\omega^2\varepsilon_0 \varepsilon_r}{\beta c_0}\int_\Omega{v^\Im\Phi^\Im \;d\Omega}$$
        a_p += omega2_eps_v * inner(v_im, phi_im) * dx

        if self.material_map.sigma is not None:
            # $$\frac{\sigma \beta c_0}{\omega}$$
            sigma_v_omega = dolfinx.Function(self.material_map.sigma.function_space)
            with self.material_map.sigma.vector.localForm() as sigma, sigma_v_omega.vector.localForm() as vec:
                vec += sigma * v / omega

            # $$\frac{\omega \sigma}{\beta c_0}$$
            omega_sigma_v = dolfinx.Function(self.material_map.sigma.function_space)
            with self.material_map.sigma.vector.localForm() as sigma, omega_sigma_v.vector.localForm() as vec:
                vec += omega * sigma / v

            # $$\frac{\sigma \beta c_0}{\omega}\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Im \;d\Omega}$$
            a_p += sigma_v_omega * inner(grad(v_re), grad(phi_im)) * dx
            # $$\frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Re\Phi^\Im \;d\Omega}$$
            a_p += omega_sigma_v * inner(v_re, phi_im) * dx
            # $$-\frac{\sigma \beta c_0}{\omega}\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Re \;d\Omega}$$
            a_p += - sigma_v_omega * inner(grad(v_im), grad(phi_re)) * dx
            # $$-\frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Im\Phi^\Im \;d\Omega}$$
            a_p += - omega_sigma_v * inner(v_im, phi_im) * dx

        # $$\int_\Omega{v^\Re J_s \;d\Omega}$$
        L_p = inner(v_re, self.solution.Js) * dx

        bc = self.__set_bc(V)
        system = dolfinx.fem.LinearProblem(a_p, L_p, bcs=[bc],
                                           petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        solution = system.solve()
        # FIXME: workaround for https://github.com/FEniCS/dolfinx/issues/1577
        self.solution.Phi_re, self.solution.Phi_im = solution.split()
        Phi_re, Phi_im = ufl.split(solution)

        Vperp = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.Hcurl)
        Vperp_bc = self.__set_bc(Vperp)
        Vz = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.H1)
        Vz_bc = self.__set_bc(Vz)
        for name, expr, V, bc in [("Ediv_re_perp", -grad(Phi_re), Vperp, Vperp_bc),
                                  ("Ediv_im_perp", -grad(Phi_im), Vperp, Vperp_bc),
                                  ("Ediv_re_z", -omega / v * Phi_im, Vz, Vz_bc),
                                  ("Ediv_im_z", omega / v * Phi_re, Vz, Vz_bc)]:
            u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
            a_p = inner(v, u) * dx
            L_p = inner(v, expr) * dx
            # projection = dolfinx.fem.LinearProblem(a_p, L_p, bcs=[bc],
            #                                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
            projection = dolfinx.fem.LinearProblem(a_p, L_p,
                                                   petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
            solution = projection.solve()
            setattr(self.solution, name, solution)
