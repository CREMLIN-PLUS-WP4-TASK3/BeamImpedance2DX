"""Curl definition."""

from enum import Enum, auto
import dolfinx
import ufl
from ufl import inner, dx
import numpy as np
from petsc4py import PETSc


class BoundaryType(Enum):
    """Boundary type."""

    DIRICHLET = auto()
    SIBC = auto()


class Ecurl():
    """Solenoidal electric field solver."""

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

    def A(self, scal):
        """Set A matrix."""
        rx = scal.dx(1)
        ry = -scal.dx(0)
        return ufl.as_vector((rx, ry))

    def B(self, vec):
        """Set B matrix."""
        rz = -vec[0].dx(1) + vec[1].dx(0)
        return rz

    def Z(self, vec):
        """Set Z matrix."""
        rx = self._omega_v * vec[1]
        ry = -self._omega_v * vec[0]
        return ufl.as_vector((rx, ry))

    def __init__(self, solution, boundary_index=None, boundary_type=BoundaryType.DIRICHLET):
        """Initialize."""
        self.material_map = solution.material_map
        self.mesh = solution.mesh
        self.solution = solution
        self.boundary_index = boundary_index
        self.boundary_type = boundary_type

        for name in ["Js", "Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im"]:
            attr = getattr(self.solution, name)
            if attr is None:
                raise AttributeError(f"Term {name} is not available")

        V = dolfinx.FunctionSpace(self.mesh.mesh, ufl.MixedElement(self.solution.Hcurl,
                                                                   self.solution.Hcurl,
                                                                   self.solution.H1,
                                                                   self.solution.H1))

        omega = self.solution._omega
        # $$\frac{\omega}{\beta c_0}$$
        omega_v = omega / (self.solution._beta * self.solution.c0)
        self._omega_v = omega_v
        # $$\frac{\omega^2}{\beta^2 c_0^2}$$
        omega2_v2 = omega_v ** 2
        A = self.A
        B = self.B
        Z = self.Z
        sigma = self.material_map.sigma
        nu_re = self.material_map.nu_re
        nu_im = self.material_map.nu_im
        # $$\omega^2 \varepsilon_0 \varepsilon_r$$
        omega2_eps = omega**2 * self.material_map.eps

        Eperp_re, Eperp_im, Ez_re, Ez_im = ufl.TrialFunctions(V)
        w_re, w_im, v_re, v_im = ufl.TestFunctions(V)

        self._a_p = 0

        """
        1
        $$
        S^{\Re\Re}_{\perp\perp}
        =\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Re\right)
        \left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\;d\Omega}
        +\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Re\nu^\Re\vec{E}^\Re_\perp\;d\Omega}
        $$
        """
        self._a_p += inner(B(w_re), nu_re * B(Eperp_re)) * dx
        self._a_p += omega2_v2 * inner(w_re, nu_re * Eperp_re) * dx

        """
        4
        $$
        S^{\Im\Re}_{z\perp}
        =-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega}
        $$
        """
        self._a_p += -inner(w_re, Z(nu_re * A(Ez_im))) * dx

        """
        6
        $$
        S^{\Im\Im}_{\perp\perp}
        =\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Im\right)
        \left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\;d\Omega}
        +\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Im\nu^\Re\vec{E}^\Im_\perp\;d\Omega}
        $$
        """
        self._a_p += inner(B(w_im), nu_re * B(Eperp_im)) * dx
        self._a_p += omega2_v2 * inner(w_im, nu_re * Eperp_im) * dx

        """
        7
        $$
        S^{\Re\Im}_{z\perp}
        =\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega}
        $$
        """
        self._a_p += inner(w_im, Z(nu_re * A(Ez_re))) * dx

        """
        10-
        $$
        S^{\Im\Re}_{\perp z}
        =\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
        \left(\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\right)\;d\Omega}
        $$
        """
        self._a_p += -inner(A(v_re), nu_re * Z(Eperp_im)) * dx

        """
        11
        $$
        S^{\Re\Re}_{z z}
        =\int_\Omega{\left(\hat{\operatorname{A}} v^\Re\right)
        \left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;d\Omega}
        $$
        """
        self._a_p += inner(A(v_re), nu_re * A(Ez_re)) * dx

        """
        13-
        $$
        S^{\Re\Im}_{\perp z}
        =-\int_\Omega{\left(\hat{\operatorname{A}}v^\Im\right)
        \left(\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\right)\;d\Omega}
        $$
        """
        self._a_p += inner(A(v_im), nu_re * Z(Eperp_re)) * dx

        """
        16
        $$
        S^{\Im\Im}_{z z}
        =\int_\Omega{\left(\hat{\operatorname{A}} v^\Im\right)
        \left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;d\Omega}
        $$
        """
        self._a_p += inner(A(v_im), nu_re * A(Ez_im)) * dx

        """
        1
        $$
        M^{\Re\Re}_{\varepsilon\perp}=
        -\omega^2\varepsilon_0\varepsilon_r\int_{\Omega}{\vec{w}^\Re\vec{E}^\Re_\perp\;d\Omega}
        $$
        """
        self._a_p += -omega2_eps * inner(w_re, Eperp_re) * dx

        """
        4
        $$
        M^{\Im\Im}_{\varepsilon\perp}=
        -\omega^2\varepsilon_0\varepsilon_r\int_{\Omega}{\vec{w}^\Im\vec{E}^\Im_\perp\;d\Omega}
        $$
        """
        self._a_p += -omega2_eps * inner(w_im, Eperp_im) * dx

        """
        5
        $$
        M^{\Re\Re}_{\varepsilon z}=
        -\omega^2\varepsilon_0\varepsilon_r\int_{\Omega}{v^\Re E^\Re_z\;d\Omega}
        $$
        """
        self._a_p += -omega2_eps * inner(v_re, Ez_re) * dx

        """
        8
        $$
        M^{\Im\Im}_{\varepsilon z}=
        -\omega^2\varepsilon_0\varepsilon_r\int_{\Omega}{v^\Im E^\Im_z\;d\Omega}
        $$
        """
        self._a_p += -omega2_eps * inner(v_im, Ez_im) * dx

        if nu_im is not None:
            pass
            """
            2
            $$
            S^{\Im\Re}_{\perp\perp}
            =-\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Re\right)
            \left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\;d\Omega}
            -\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Re\nu^\Im\vec{E}^\Im_\perp\;d\Omega}
            $$
            """
            # self._a_p += -inner(B(w_re), nu_re * B(Eperp_re)) * dx
            # self._a_p += omega2_v2 * inner(w_re, nu_re * Eperp_re) * dx

            """
            3
            $$
            S^{\Re\Re}_{z\perp}
            =-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega}
            $$
            """
            self._a_p += -inner(w_re, nu_re * Eperp_re) * dx

            """
            5
            $$
            S^{\Re\Im}_{\perp\perp}
            =\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Im\right)
            \left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\;d\Omega}
            +\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Im\nu^\Im\vec{E}^\Re_\perp\;d\Omega}
            $$
            """

            """
            8
            $$
            S^{\Im\Im}_{z\perp}
            =-\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega}
            $$
            """

            """
            9
            $$
            S^{\Re\Re}_{\perp z}
            =\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
            \left(\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\right)\;d\Omega}
            $$
            """

            """
            12
            $$
            S^{\Im\Re}_{z z}
            =-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
            \left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;d\Omega}
            $$
            """

            """
            14
            $$
            S^{\Im\Im}_{\perp z}
            =\int_\Omega{\left(\hat{\operatorname{A}}v^\Im\right)
            \left(\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\right)\;d\Omega}
            $$
            """

            """
            15
            $$
            S^{\Re\Im}_{z z}
            =\int_\Omega{\left(\hat{\operatorname{A}} v^\Im\right)
            \left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;d\Omega}
            $$
            """

        if sigma is not None:
            # $$\omega^2 \varepsilon_0 \varepsilon_r$$
            omega_sigma = dolfinx.Function(sigma.function_space)
            with sigma.vector.localForm() as sigma, omega_sigma.vector.localForm() as vec:
                vec += omega * sigma

            """
            2
            $$
            M^{\Im\Re}_{\sigma\perp}=
            \omega\sigma\int_{\Omega}{\vec{w}^\Re\vec{E}^\Im_\perp\;d\Omega}
            $$
            """

            """
            3
            $$
            M^{\Re\Im}_{\sigma\perp}=
            -\omega\sigma\int_{\Omega}{\vec{w}^\Im\vec{E}^\Re_\perp\;d\Omega}
            $$
            """

            """
            6
            $$
            M^{\Im\Re}_{\sigma z}=
            \omega\sigma\int_{\Omega}{v^\Re E^\Im_z\;d\Omega}
            $$
            """

            """
            7
            $$
            M^{\Re\Im}_{\sigma z}=
            -\omega\sigma\int_{\Omega}{v^\Im E^\Re_z\;d\Omega}
            $$
            """

        if self.boundary_type == BoundaryType.SIBC:
            pass
            """
            $$
            B^{\Re\Re}_{\perp\perp}
            =\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)
            \vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Im\Re}_{\perp\perp}
            =-\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Re\Im}_{\perp\perp}
            =\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Im\Im}_{\perp\perp}
            =\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Re\Re}_{\perp z}
            =\int_{\partial\Omega}{v^\Re\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Im\Re}_{\perp z}
            =\int_{\partial\Omega}{v^\Re\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Re\Re}_{z z}
            =-\int_{\partial\Omega}{v^\Re\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Im\Re}_{z z}
            =\int_{\partial\Omega}{v^\Re\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Re\Im}_{\perp z}
            =-\int_{\partial\Omega}{v^\Im\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Im\Im}_{\perp z}
            =\int_{\partial\Omega}{v^\Im\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Re\Im}_{z z}
            =-\int_{\partial\Omega}{v^\Im\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\vec{\tau}\;dS}
            $$
            """

            """
            $$
            B^{\Im\Im}_{z z}
            =-\int_{\partial\Omega}{v^\Im\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\vec{\tau}\;dS}
            $$
            """

        """
        $$
        J_s^\Im = -\omega\int_\Omega{v^\Im J_s \;d\Omega}
        $$
        """
        self._L_p = -omega * inner(v_im, self.solution.Js) * dx

        """
        1
        $$
        N^{\Re\Re}_{\varepsilon\perp}=
        \omega^2\varepsilon_0\varepsilon_r\int_{\Omega}{\vec{w}^\Re\vec{E}^\Re_\perp\;d\Omega}
        $$
        """
        self._L_p += omega2_eps * inner(w_re, self.solution.Ediv_perp_re) * dx

        """
        4
        $$
        N^{\Im\Im}_{\varepsilon\perp}=
        \omega^2\varepsilon_0\varepsilon_r\int_{\Omega}{\vec{w}^\Im\vec{E}^\Im_\perp\;d\Omega}
        $$
        """
        self._L_p += omega2_eps * inner(w_im, self.solution.Ediv_perp_im) * dx

        """
        5
        $$
        N^{\Re\Re}_{\varepsilon z}=
        \omega^2\varepsilon_0\varepsilon_r\int_{\Omega}{v^\Re E^\Re_z\;d\Omega}
        $$
        """
        self._L_p += omega2_eps * inner(v_re, self.solution.Ediv_z_re) * dx

        """
        8
        $$
        N^{\Im\Im}_{\varepsilon z}=
        \omega^2\varepsilon_0\varepsilon_r\int_{\Omega}{v^\Im E^\Im_z\;d\Omega}
        $$
        """
        self._L_p += omega2_eps * inner(v_im, self.solution.Ediv_z_im) * dx

        if sigma is not None:
            pass
            """
            2
            $$
            N^{\Im\Re}_{\sigma\perp}=
            -\omega\sigma\int_{\Omega}{\vec{w}^\Re\vec{E}^\Im_\perp\;d\Omega}
            $$
            """

            """
            3
            $$
            N^{\Re\Im}_{\sigma\perp}=
            \omega\sigma\int_{\Omega}{\vec{w}^\Im\vec{E}^\Re_\perp\;d\Omega}
            $$
            """

            """
            6
            $$
            N^{\Im\Re}_{\sigma z}=
            -\omega\sigma\int_{\Omega}{v^\Re E^\Im_z\;d\Omega}
            $$
            """

            """
            7
            $$
            N^{\Re\Im}_{\sigma z}=
            \omega\sigma\int_{\Omega}{v^\Im E^\Re_z\;d\Omega}
            $$
            """

        self._bc = self.__set_bc(V)

        self._A = dolfinx.fem.create_matrix(self._a_p)
        self._b = dolfinx.fem.create_vector(self._L_p)

        self._Ecurl = dolfinx.Function(V)

        # FIXME: workaround for https://github.com/FEniCS/dolfinx/issues/1577
        (self.solution.Ecurl_perp_re,
         self.solution.Ecurl_perp_im,
         self.solution.Ecurl_z_re,
         self.solution.Ecurl_z_im) = self._Ecurl.split()
        (self.solution.ecurl_perp_re,
         self.solution.ecurl_perp_im,
         self.solution.ecurl_z_re,
         self.solution.ecurl_z_im) = ufl.split(self._Ecurl)

    def solve(self, petsc_options={"linear_solver": "preonly",
                                   "preconditioner": "la"}):
        """Solve equation."""

        self.solution.solver.reset()
        self.solution._set_solver_options(petsc_options)

        self._A.zeroEntries()
        dolfinx.fem.assemble_matrix(self._A, self._a_p, bcs=[self._bc])
        self._A.assemble()
        with self._b.localForm() as b_loc:
            b_loc.set(0)
        dolfinx.fem.assemble_vector(self._b, self._L_p)
        self.solution.solver.setOperators(self._A)

        dolfinx.fem.apply_lifting(self._b, [self._a_p], [[self._bc]])
        self._b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.set_bc(self._b, [self._bc])

        self.solution.solver.solve(self._b, self._Ecurl.vector)
        self._Ecurl.x.scatter_forward()
