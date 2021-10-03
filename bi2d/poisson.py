"""Poisson definition."""

import dolfinx
import ufl
from ufl import inner, grad, dx
import numpy as np
from petsc4py import PETSc

class Ediv():
    """Irrotational electric field solver."""

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

    def __init__(self, solution, boundary_index=None):
        """Initialize."""
        self.material_map = solution.material_map
        self.mesh = solution.mesh
        self.solution = solution
        self.boundary_index = boundary_index
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
            self._a_phi += - sigma_v_omega * inner(grad(v_im), grad(phi_re)) * dx
            # $$-\frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Im\Phi^\Im \;d\Omega}$$
            self._a_phi += - omega_sigma_v * inner(v_im, phi_im) * dx

        # $$\int_\Omega{v^\Re J_s \;d\Omega}$$
        self._L_phi = inner(v_re, self.solution.Js) * dx

        self._A_phi = dolfinx.fem.create_matrix(self._a_phi)
        self._b_phi = dolfinx.fem.create_vector(self._L_phi)
        self._bc = self.__set_bc(V)
        self._phi = dolfinx.Function(V)
        self.solution.Phi_re, self.solution.Phi_im = self._phi.split()
        # FIXME: workaround for https://github.com/FEniCS/dolfinx/issues/1577
        Phi_re, Phi_im = ufl.split(self._phi)

        Vperp = dolfinx.FunctionSpace(self.mesh.mesh, self.solution.Hcurl)
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
            setattr(self.solution, name, dolfinx.Function(V))

    def solve(self, petsc_options={"linear_solver": "preonly",
                                   "preconditioner": "la"}):
        """Solve equation."""
        if self.solution.Js is None:
            raise AttributeError("Term Js is not available")

        self.solution.solver.reset()
        self.solution._set_solver_options(petsc_options)

        self._A_phi.zeroEntries()
        dolfinx.fem.assemble_matrix(self._A_phi, self._a_phi, bcs=[self._bc])
        self._A_phi.assemble()
        with self._b_phi.localForm() as b_loc:
            b_loc.set(0)
        dolfinx.fem.assemble_vector(self._b_phi, self._L_phi)
        self.solution.solver.setOperators(self._A_phi)

        dolfinx.fem.apply_lifting(self._b_phi, [self._a_phi], [[self._bc]])
        self._b_phi.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.set_bc(self._b_phi, [self._bc])

        self.solution.solver.solve(self._b_phi, self._phi.vector)
        self._phi.x.scatter_forward()

        for name in ["Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im"]:
            self.solution.solver.reset()

            a_p = getattr(self, f"_a_{name}")
            L_p = getattr(self, f"_L_{name}")
            A = getattr(self, f"_A_{name}")
            b = getattr(self, f"_b_{name}")
            f = getattr(self.solution, name)

            A.zeroEntries()
            dolfinx.fem.assemble_matrix(A, a_p)
            A.assemble()
            with b.localForm() as b_loc:
                b_loc.set(0)
            dolfinx.fem.assemble_vector(b, L_p)
            self.solution.solver.setOperators(A)

            self.solution.solver.solve(b, f.vector)
            f.x.scatter_forward()
