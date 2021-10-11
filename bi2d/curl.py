"""Curl definition."""

from enum import Enum, auto
import dolfinx
import ufl
from ufl import inner, dx
import numpy as np
from petsc4py import PETSc
from mpi4py import MPI


class BoundaryType(Enum):
    """Boundary type."""

    DIRICHLET = auto()
    SIBC = auto()


class Ecurl():
    """Solenoidal electric field solver."""

    def __set_bc(self, V, bcs, value=0):
        u_bc = dolfinx.Function(V)
        bc_dofs = []
        for boundary_index, boundary_type in bcs:
            if boundary_type == BoundaryType.DIRICHLET:
                if boundary_index == -1:
                    bc_facets = np.where(
                        np.array(dolfinx.cpp.mesh.compute_boundary_facets(self.mesh.mesh.topology)) == 1)[0]
                    dofs = dolfinx.fem.locate_dofs_topological(V,
                                                               self.mesh.mesh.topology.dim-1,
                                                               bc_facets)
                else:
                    dofs = dolfinx.fem.locate_dofs_topological(V,
                                                               self.mesh.mesh.topology.dim-1,
                                                               self.mesh.get_boundary(boundary_index))
                bc_dofs.append(dofs)

        bc_dofs = np.hstack(bc_dofs)
        with u_bc.vector.localForm() as loc:
            loc.setValues(bc_dofs, np.full(bc_dofs.size, value))
        u_bc.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        return dolfinx.DirichletBC(u_bc, bc_dofs)

    def A(self, scal):
        """Apply A operator."""
        rx = scal.dx(1)
        ry = -scal.dx(0)
        return ufl.as_vector((rx, ry))

    def B(self, vec):
        """Apply B operator."""
        rz = -vec[0].dx(1) + vec[1].dx(0)
        return rz

    def Z(self, vec):
        """Apply Z operator."""
        rx = self.solution._omega / (self.solution._beta * self.solution.c0) * vec[1]
        ry = -self.solution._omega / (self.solution._beta * self.solution.c0) * vec[0]
        return ufl.as_vector((rx, ry))

    def __init__(self, solution, bcs=[(-1, BoundaryType.DIRICHLET)]):
        """Initialize."""
        self.material_map = solution.material_map
        self.mesh = solution.mesh
        self.solution = solution

        for name in ["Js", "Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im"]:
            attr = getattr(self.solution, name)
            if attr is None:
                raise AttributeError(f"{name} function is not available")

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Setting curl function")

        V = dolfinx.FunctionSpace(self.mesh.mesh, ufl.MixedElement(self.solution.Hcurl,
                                                                   self.solution.Hcurl,
                                                                   self.solution.H1,
                                                                   self.solution.H1))

        omega = self.solution._omega
        # $$\frac{\omega}{\beta c_0}$$
        omega_v = omega / (self.solution._beta * self.solution.c0)
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

        r"""
        1
        $$
        S^{\Re\Re}_{\perp\perp}
        =\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Re\right)
        \left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\;d\Omega}
        +\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Re\nu^\Re\vec{E}^\Re_\perp\;d\Omega}
        $$
        """
        self._a_p += inner(B(w_re), nu_re * B(Eperp_re)) * dx
        self._a_p += inner(w_re, omega2_v2 * nu_re * Eperp_re) * dx

        r"""
        4
        $$
        S^{\Im\Re}_{z\perp}
        =-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega}
        $$
        """
        self._a_p += -inner(w_re, Z(nu_re * A(Ez_im))) * dx

        r"""
        6
        $$
        S^{\Im\Im}_{\perp\perp}
        =\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Im\right)
        \left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\;d\Omega}
        +\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Im\nu^\Re\vec{E}^\Im_\perp\;d\Omega}
        $$
        """
        self._a_p += inner(B(w_im), nu_re * B(Eperp_im)) * dx
        self._a_p += inner(w_im, omega2_v2 * nu_re * Eperp_im) * dx

        r"""
        7
        $$
        S^{\Re\Im}_{z\perp}
        =\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega}
        $$
        """
        self._a_p += inner(w_im, Z(nu_re * A(Ez_re))) * dx

        r"""
        10
        $$
        S^{\Im\Re}_{\perp z}
        =-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
        \left(\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\right)\;d\Omega}
        $$
        """
        self._a_p += -inner(A(v_re), nu_re * Z(Eperp_im)) * dx

        r"""
        11
        $$
        S^{\Re\Re}_{z z}
        =\int_\Omega{\left(\hat{\operatorname{A}} v^\Re\right)
        \left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;d\Omega}
        $$
        """
        self._a_p += inner(A(v_re), nu_re * A(Ez_re)) * dx

        r"""
        13
        $$
        S^{\Re\Im}_{\perp z}
        =\int_\Omega{\left(\hat{\operatorname{A}}v^\Im\right)
        \left(\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\right)\;d\Omega}
        $$
        """
        self._a_p += inner(A(v_im), nu_re * Z(Eperp_re)) * dx

        r"""
        16
        $$
        S^{\Im\Im}_{z z}
        =\int_\Omega{\left(\hat{\operatorname{A}} v^\Im\right)
        \left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;d\Omega}
        $$
        """
        self._a_p += inner(A(v_im), nu_re * A(Ez_im)) * dx

        r"""
        1
        $$
        M^{\Re\Re}_{\varepsilon\perp}=
        -\omega^2\int_{\Omega}{\vec{w}^\Re\varepsilon_0\varepsilon_r\vec{E}^\Re_\perp\;d\Omega}
        $$
        """
        self._a_p += -inner(w_re, omega2_eps * Eperp_re) * dx

        r"""
        6
        $$
        M^{\Im\Im}_{\varepsilon\perp}=
        -\omega^2\int_{\Omega}{\vec{w}^\Im\varepsilon_0\varepsilon_r\vec{E}^\Im_\perp\;d\Omega}
        $$
        """
        self._a_p += -inner(w_im, omega2_eps * Eperp_im) * dx

        r"""
        11
        $$
        M^{\Re\Re}_{\varepsilon z}=
        -\omega^2\int_{\Omega}{v^\Re\varepsilon_0\varepsilon_r E^\Re_z\;d\Omega}
        $$
        """
        self._a_p += -inner(v_re, omega2_eps * Ez_re) * dx

        r"""
        16
        $$
        M^{\Im\Im}_{\varepsilon z}=
        -\omega^2\int_{\Omega}{v^\Im\varepsilon_0\varepsilon_r E^\Im_z\;d\Omega}
        $$
        """
        self._a_p += -inner(v_im, omega2_eps * Ez_im) * dx

        if nu_im is not None:
            r"""
            2
            $$
            S^{\Im\Re}_{\perp\perp}
            =-\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Re\right)
            \left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\;d\Omega}
            -\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Re\nu^\Im\vec{E}^\Im_\perp\;d\Omega}
            $$
            """
            self._a_p += -inner(B(w_re), nu_im * B(Eperp_im)) * dx
            self._a_p += -omega2_v2 * inner(w_re, nu_im * Eperp_im) * dx

            r"""
            3
            $$
            S^{\Re\Re}_{z\perp}
            =-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega}
            $$
            """
            self._a_p += -inner(w_re, Z(nu_im * A(Ez_re))) * dx

            r"""
            5
            $$
            S^{\Re\Im}_{\perp\perp}
            =\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Im\right)
            \left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\;d\Omega}
            +\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Im\nu^\Im\vec{E}^\Re_\perp\;d\Omega}
            $$
            """
            self._a_p += inner(B(w_im), nu_im * B(Eperp_re)) * dx
            self._a_p += omega2_v2 * inner(w_im, nu_im * Eperp_re) * dx

            r"""
            8
            $$
            S^{\Im\Im}_{z\perp}
            =-\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega}
            $$
            """
            self._a_p += -inner(w_im, Z(nu_im * A(Ez_im))) * dx

            r"""
            9
            $$
            S^{\Re\Re}_{\perp z}
            =-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
            \left(\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\right)\;d\Omega}
            $$
            """
            self._a_p += -inner(A(v_re), nu_im * Z(Eperp_re)) * dx

            r"""
            12
            $$
            S^{\Im\Re}_{z z}
            =-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
            \left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;d\Omega}
            $$
            """
            self._a_p += -inner(A(v_re), nu_im * A(Ez_im)) * dx

            r"""
            14
            $$
            S^{\Im\Im}_{\perp z}
            =-\int_\Omega{\left(\hat{\operatorname{A}}v^\Im\right)
            \left(\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\right)\;d\Omega}
            $$
            """
            self._a_p += -inner(A(v_im), nu_im * Z(Eperp_im)) * dx

            r"""
            15
            $$
            S^{\Re\Im}_{z z}
            =\int_\Omega{\left(\hat{\operatorname{A}} v^\Im\right)
            \left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;d\Omega}
            $$
            """
            self._a_p += inner(A(v_im), nu_im * A(Ez_re)) * dx

        if sigma is not None:
            # $$\omega\sigma$$
            omega_sigma = omega * sigma

            r"""
            2
            $$
            M^{\Im\Re}_{\sigma\perp}=
            -\omega\int_{\Omega}{\vec{w}^\Re\sigma\vec{E}^\Im_\perp\;d\Omega}
            $$
            """
            self._a_p += -inner(w_re, omega_sigma * Eperp_im) * dx

            r"""
            5
            $$
            M^{\Re\Im}_{\sigma\perp}=
            \omega\int_{\Omega}{\vec{w}^\Im\sigma\vec{E}^\Re_\perp\;d\Omega}
            $$
            """
            self._a_p += inner(w_im, omega_sigma * Eperp_re) * dx

            r"""
            12
            $$
            M^{\Im\Re}_{\sigma z}=
            -\omega\int_{\Omega}{v^\Re\sigma E^\Im_z\;d\Omega}
            $$
            """
            self._a_p += -inner(v_re, omega_sigma * Ez_im) * dx

            r"""
            15
            $$
            M^{\Re\Im}_{\sigma z}=
            \omega\int_{\Omega}{v^\Im\sigma E^\Re_z\;d\Omega}
            $$
            """
            self._a_p += inner(v_im, omega_sigma * Ez_re) * dx

        # if self.boundary_type == BoundaryType.SIBC:
        #     pass
            r"""
            $$
            B^{\Re\Re}_{\perp\perp}
            =\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)
            \vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Im\Re}_{\perp\perp}
            =-\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Re\Im}_{\perp\perp}
            =\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Im\Im}_{\perp\perp}
            =\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Re\Re}_{\perp z}
            =\int_{\partial\Omega}{v^\Re\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Im\Re}_{\perp z}
            =\int_{\partial\Omega}{v^\Re\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Re\Re}_{z z}
            =-\int_{\partial\Omega}{v^\Re\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Im\Re}_{z z}
            =\int_{\partial\Omega}{v^\Re\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Re\Im}_{\perp z}
            =-\int_{\partial\Omega}{v^\Im\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Im\Im}_{\perp z}
            =\int_{\partial\Omega}{v^\Im\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Re\Im}_{z z}
            =-\int_{\partial\Omega}{v^\Im\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\vec{\tau}\;dS}
            $$
            """

            r"""
            $$
            B^{\Im\Im}_{z z}
            =-\int_{\partial\Omega}{v^\Im\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\vec{\tau}\;dS}
            $$
            """

        r"""
        $$
        J_s^\Im = -\omega\int_\Omega{v^\Im J_s \;d\Omega}
        $$
        """
        self._L_p = -inner(v_im, omega * self.solution.Js) * dx

        r"""
        1
        $$
        N^{\Re\Re}_{\varepsilon\perp}=
        \omega^2\int_{\Omega}{\vec{w}^\Re\varepsilon_0\varepsilon_r\vec{E}^\Re_\perp\;d\Omega}
        $$
        """
        self._L_p += inner(w_re, omega2_eps * self.solution.Ediv_perp_re) * dx

        r"""
        6
        $$
        N^{\Im\Im}_{\varepsilon\perp}=
        \omega^2\int_{\Omega}{\vec{w}^\Im\varepsilon_0\varepsilon_r\vec{E}^\Im_\perp\;d\Omega}
        $$
        """
        self._L_p += inner(w_im, omega2_eps * self.solution.Ediv_perp_im) * dx

        r"""
        11
        $$
        N^{\Re\Re}_{\varepsilon z}=
        \omega^2\int_{\Omega}{v^\Re\varepsilon_0\varepsilon_r E^\Re_z\;d\Omega}
        $$
        """
        self._L_p += inner(v_re, omega2_eps * self.solution.Ediv_z_re) * dx

        r"""
        16
        $$
        N^{\Im\Im}_{\varepsilon z}=
        \omega^2\int_{\Omega}{v^\Im\varepsilon_0\varepsilon_r E^\Im_z\;d\Omega}
        $$
        """
        self._L_p += inner(v_im, omega2_eps * self.solution.Ediv_z_im) * dx

        if sigma is not None:
            r"""
            2
            $$
            N^{\Im\Re}_{\sigma\perp}=
            \omega\int_{\Omega}{\vec{w}^\Re\sigma\vec{E}^\Im_\perp\;d\Omega}
            $$
            """
            self._L_p += inner(w_re, omega_sigma * self.solution.Ediv_perp_im) * dx

            r"""
            5
            $$
            N^{\Re\Im}_{\sigma\perp}=
            -\omega\int_{\Omega}{\vec{w}^\Im\sigma\vec{E}^\Re_\perp\;d\Omega}
            $$
            """
            self._L_p += -inner(w_im, omega_sigma * self.solution.Ediv_perp_re) * dx

            r"""
            12
            $$
            N^{\Im\Re}_{\sigma z}=
            \omega\int_{\Omega}{v^\Re\sigma E^\Im_z\;d\Omega}
            $$
            """
            self._L_p += inner(v_re, omega_sigma * self.solution.Ediv_z_im) * dx

            r"""
            15
            $$
            N^{\Re\Im}_{\sigma z}=
            -\omega\int_{\Omega}{v^\Im\sigma E^\Re_z\;d\Omega}
            $$
            """
            self._L_p += -inner(v_im, omega_sigma * self.solution.Ediv_z_re) * dx

        self._bc = self.__set_bc(V, bcs)

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

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Set curl function")

    def solve(self, petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}):
        """Solve equation."""

        if self.solution._Js_stale:
            raise ValueError("Source function Js solution is stale and needs to be recalculated first")
        if self.solution._Ediv_stale:
            raise ValueError("Ediv solutions are stale and need to be recalculated")

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solving curl function")

        self.solution._solve(self._a_p, self._L_p, self._A,
                             self._b, self._Ecurl, bcs=[self._bc],
                             petsc_options=petsc_options)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solved curl function")

        self.solution._Ecurl_stale = False
