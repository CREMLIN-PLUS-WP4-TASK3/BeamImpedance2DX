"""Curl definition."""

from inspect import signature
from numbers import Number
import dolfinx
import ufl
from ufl import inner, dx, dot
import numpy as np
from petsc4py import PETSc
from mpi4py import MPI


def simplify(expr):
    if issubclass(type(expr), ufl.core.operator.Operator):
        return PETSc.ScalarType(expr.__class__(*[simplify(e) for e in expr.ufl_operands]))
    elif issubclass(type(expr), dolfinx.Constant):
        return PETSc.ScalarType(expr.value)
    elif issubclass(type(expr), ufl.constantvalue.ConstantValue):
        return PETSc.ScalarType(expr)
    elif type(expr) == PETSc.ScalarType:
        return expr
    else:
        raise TypeError(f"Unsupported type {type(expr)}")


class Ecurl():
    """Solenoidal electric field solver."""

    def __set_dirichlet(self, V, sibc):
        if -1 in [bc.index for bc in sibc]:
            # assigned SIBC to all boundaries
            return []

        bc_dofs = np.array([])
        bc_facets = np.where(np.array(dolfinx.cpp.mesh.compute_boundary_facets(self.mesh.mesh.topology)) == 1)[0]
        bc_dofs = dolfinx.fem.locate_dofs_topological(V, self.mesh.mesh.topology.dim - 1, bc_facets)

        # remove SIBC regions
        for bc in sibc:
            dofs = dolfinx.fem.locate_dofs_topological(V, self.mesh.mesh.topology.dim - 1,
                                                       self.mesh.get_boundary(bc.index))
            bc_dofs = np.setdiff1d(bc_dofs, dofs)

        if bc_dofs.size == 0:
            return []
        else:
            u_bc = dolfinx.Function(V)
            with u_bc.vector.localForm() as loc:
                loc.setValues(bc_dofs, np.full(bc_dofs.size, 0))
            u_bc.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
            return [dolfinx.DirichletBC(u_bc, bc_dofs)]

    def __create_sibc_measures(self, sibc):
        if len(sibc) == 0:
            return [], []
        if -1 in [material.index for material in sibc]:
            # assign SIBC to all boundaries
            measures = [ufl.ds]
            materials = [material for material in sibc if material.index == -1]
            deltas = [dolfinx.Constant(self.mesh.mesh, PETSc.ScalarType(0))]
        else:
            if self.material_map.mesh.boundaries is None:
                raise ValueError("Boundary data was not loaded")
            for material in sibc:
                if material.index not in self.material_map.mesh.boundaries.values:
                    raise ValueError(f"No boundary data with index {material.index}")

            ds = ufl.Measure('ds', domain=self.material_map.mesh.mesh,
                             subdomain_data=self.material_map.mesh.boundaries)
            measures = [ds(material.index) for material in sibc]
            materials = [material for material in sibc]
            deltas = [dolfinx.Constant(self.mesh.mesh, PETSc.ScalarType(0)) for _ in sibc]
        return measures, list(zip(deltas, materials))

    def __set_skin_depth_from_material(self, k, material, omega):
        if isinstance(material.sigma, Number):
            sigma = material.sigma
        elif ((type(material.sigma).__name__ == "function" or
               type(material.sigma).__name__ == "method") and
              (len(signature(material.sigma).parameters) == 1)):
            sigma = material.sigma(2 * np.pi * omega)
        else:
            raise AttributeError("Only constant or dispersive conductivity is allowed for SIBC material")

        r"""
        $$
        k=\frac{1}{\mu_0\delta}
        =\frac{1}{\mu_0\sqrt{\frac{2}{\mu_0\omega\sigma}}}
        =\sqrt{\frac{\omega\sigma}{2\mu_0}}
        $$
        """
        k.value = PETSc.ScalarType(np.sqrt((omega * sigma) / (2 * self.material_map.mu0)))

    def A(self, scal):
        """Apply A operator."""
        rx = scal.dx(1)
        ry = -scal.dx(0)
        return ufl.as_vector((rx, ry))

    def B(self, vec):
        """Apply B operator."""
        rz = -vec[0].dx(1) + vec[1].dx(0)
        return rz

    def Z_real(self, vec):
        """Apply Z operator."""
        rx = self.omega / (self.beta * self.c0) * vec[1]
        ry = -self.omega / (self.beta * self.c0) * vec[0]
        return ufl.as_vector((rx, ry))

    def Z_complex(self, vec):
        """Apply Z operator."""
        rx = PETSc.ScalarType(1j) * self.omega / (self.beta * self.c0) * vec[1]
        ry = PETSc.ScalarType(-1j) * self.omega / (self.beta * self.c0) * vec[0]
        return ufl.as_vector((rx, ry))

    def __init__(self, solution, sibc=[]):
        """Initialize."""
        self.material_map = solution.material_map
        self.mesh = solution.mesh
        self.solution = solution

        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            for name in ["Js", "Ediv_perp", "Ediv_z"]:
                attr = getattr(self.solution, name)
                if attr is None:
                    raise AttributeError(f"{name} function is not available")
        else:
            for name in ["Js", "Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im"]:
                attr = getattr(self.solution, name)
                if attr is None:
                    raise AttributeError(f"{name} function is not available")

        for i in np.unique([bc.index for bc in sibc]):
            if np.sum([bc.index for bc in sibc] == 1) > 1:
                raise(IndexError(f"Duplicate use of boundary index {i}"))

        self.solution._Ecurl_stale = True

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Setting curl function")

        omega = self.solution._omega
        self.omega = omega
        beta = self.solution._beta
        self.beta = beta
        c0 = self.solution.c0
        self.c0 = c0
        eps = self.material_map.eps
        A = self.A
        B = self.B
        # TODO: experiment with different values and probably remove
        RHS_mult = PETSc.ScalarType(1)
        self.RHS_mult = RHS_mult
        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            Z = self.Z_complex
            mu = self.material_map.mu
        else:
            Z = self.Z_real
            sigma = self.material_map.sigma
            nu_re = self.material_map.nu_re
            nu_im = self.material_map.nu_im

        self._a_p = 0
        self._L_p = 0

        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            V = dolfinx.FunctionSpace(self.mesh.mesh, ufl.MixedElement(self.solution.Hcurl,
                                                                       self.solution.H1))

            Eperp, Ez = ufl.TrialFunctions(V)
            w, v = ufl.TestFunctions(V)

            r"""
            1
            $$\hat{S}_{\perp\perp}
            =\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}\right)\left(\frac{1}{\mu}
            \hat{\operatorname{B}}\underline{\vec{E}}_\perp\right)\;d\Omega}
            +\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}\frac{1}{\mu}\underline{\vec{E}}_\perp\;d\Omega}$$
            """
            self._a_p += inner(1 / mu * B(Eperp), B(w)) * dx
            self._a_p += inner(omega**2 / (beta * c0)**2 * 1 / mu * Eperp, w) * dx

            r"""
            2
            $$\hat{S}_{\perp z}
            =\int_{\Omega}{\vec{w}\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}
            \hat{\operatorname{A}}\underline{E}_z\;d\Omega}$$
            """
            self._a_p += inner(Z(1 / mu * A(Ez)), w) * dx

            r"""
            3
            $$\hat{S}_{z\perp}
            =\int_\Omega{\left(\hat{\operatorname{A}}v\right)
            \left(\frac{1}{\mu}\hat{\operatorname{Z}}\underline{\vec{E}}_\perp\right)\;d\Omega}$$
            """
            self._a_p += inner(1 / mu * Z(Eperp), A(v)) * dx

            r"""
            4
            $$\hat{S}_{zz}\underline{E}_z
            =\int_\Omega{\left(\hat{\operatorname{A}} v\right)
            \left(\frac{1}{\mu}\hat{\operatorname{A}}\underline{E}_z\right)\;d\Omega}$$
            """
            self._a_p += inner(1 / mu * A(Ez), A(v)) * dx

            r"""
            1
            $$
            M_{\perp}=
            -\omega^2\int_{\Omega}{\vec{w}\underline{\varepsilon}\underline{\vec{E}}_\perp\;d\Omega}
            $$
            """
            self._a_p += -inner(omega**2 * eps * Eperp, w) * dx

            r"""
            4
            $$
            M_{z}=
            -\omega^2\int_{\Omega}{v\underline{\varepsilon}\underline{E}_z\;d\Omega}
            $$
            """
            self._a_p += -inner(omega**2 * eps * Ez, v) * dx

            r"""
            $$
            J_s = -j\omega\int_\Omega{v J_s \;d\Omega}
            $$
            """
            self._L_p += RHS_mult * PETSc.ScalarType(-1j) * inner(omega * self.solution.Js, v) * dx

            r"""
            1
            $$
            N_{\perp}=
            \omega^2\int_{\Omega}{\vec{w}\underline{\varepsilon}\underline{\vec{E}}_\perp\;d\Omega}
            $$
            """
            self._L_p += RHS_mult * inner(omega**2 * eps * self.solution.Ediv_perp, w) * dx

            r"""
            4
            $$
            N_{z}=
            \omega^2\int_{\Omega}{v\underline{\varepsilon}\underline{E}_z\;d\Omega}
            $$
            """
            self._L_p += RHS_mult * inner(omega**2 * eps * self.solution.Ediv_z, v) * dx

        else:
            V = dolfinx.FunctionSpace(self.mesh.mesh, ufl.MixedElement(self.solution.Hcurl,
                                                                       self.solution.Hcurl,
                                                                       self.solution.H1,
                                                                       self.solution.H1))

            Eperp_re, Eperp_im, Ez_re, Ez_im = ufl.TrialFunctions(V)
            w_re, w_im, v_re, v_im = ufl.TestFunctions(V)

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
            self._a_p += omega**2 / (beta * c0)**2 * inner(w_re, nu_re * Eperp_re) * dx

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
            self._a_p += omega**2 / (beta * c0)**2 * inner(w_im, nu_re * Eperp_im) * dx

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
            self._a_p += -omega**2 * inner(w_re, eps * Eperp_re) * dx

            r"""
            6
            $$
            M^{\Im\Im}_{\varepsilon\perp}=
            -\omega^2\int_{\Omega}{\vec{w}^\Im\varepsilon_0\varepsilon_r\vec{E}^\Im_\perp\;d\Omega}
            $$
            """
            self._a_p += -omega**2 * inner(w_im, eps * Eperp_im) * dx

            r"""
            11
            $$
            M^{\Re\Re}_{\varepsilon z}=
            -\omega^2\int_{\Omega}{v^\Re\varepsilon_0\varepsilon_r E^\Re_z\;d\Omega}
            $$
            """
            self._a_p += -omega**2 * inner(v_re, eps * Ez_re) * dx

            r"""
            16
            $$
            M^{\Im\Im}_{\varepsilon z}=
            -\omega^2\int_{\Omega}{v^\Im\varepsilon_0\varepsilon_r E^\Im_z\;d\Omega}
            $$
            """
            self._a_p += -omega**2 * inner(v_im, eps * Ez_im) * dx

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
                self._a_p += -omega**2 / (beta * c0)**2 * inner(w_re, nu_im * Eperp_im) * dx

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
                self._a_p += omega**2 / (beta * c0)**2 * inner(w_im, nu_im * Eperp_re) * dx

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

                r"""
                2
                $$
                M^{\Im\Re}_{\sigma\perp}=
                -\omega\int_{\Omega}{\vec{w}^\Re\sigma\vec{E}^\Im_\perp\;d\Omega}
                $$
                """
                self._a_p += -omega * inner(w_re, sigma * Eperp_im) * dx

                r"""
                5
                $$
                M^{\Re\Im}_{\sigma\perp}=
                \omega\int_{\Omega}{\vec{w}^\Im\sigma\vec{E}^\Re_\perp\;d\Omega}
                $$
                """
                self._a_p += omega * inner(w_im, sigma * Eperp_re) * dx

                r"""
                12
                $$
                M^{\Im\Re}_{\sigma z}=
                -\omega\int_{\Omega}{v^\Re\sigma E^\Im_z\;d\Omega}
                $$
                """
                self._a_p += -omega * inner(v_re, sigma * Ez_im) * dx

                r"""
                15
                $$
                M^{\Re\Im}_{\sigma z}=
                \omega\int_{\Omega}{v^\Im\sigma E^\Re_z\;d\Omega}
                $$
                """
                self._a_p += omega * inner(v_im, sigma * Ez_re) * dx

            r"""
            $$
            J_s^\Im = -\omega\int_\Omega{v^\Im J_s \;d\Omega}
            $$
            """
            self._L_p += -RHS_mult * omega * inner(v_im, self.solution.Js) * dx

            r"""
            1
            $$
            N^{\Re\Re}_{\varepsilon\perp}=
            \omega^2\int_{\Omega}{\vec{w}^\Re\varepsilon_0\varepsilon_r\vec{E}^\Re_\perp\;d\Omega}
            $$
            """
            self._L_p += RHS_mult * omega**2 * inner(w_re, eps * self.solution.Ediv_perp_re) * dx

            r"""
            6
            $$
            N^{\Im\Im}_{\varepsilon\perp}=
            \omega^2\int_{\Omega}{\vec{w}^\Im\varepsilon_0\varepsilon_r\vec{E}^\Im_\perp\;d\Omega}
            $$
            """
            self._L_p += RHS_mult * inner(w_im, omega**2 * eps * self.solution.Ediv_perp_im) * dx

            r"""
            11
            $$
            N^{\Re\Re}_{\varepsilon z}=
            \omega^2\int_{\Omega}{v^\Re\varepsilon_0\varepsilon_r E^\Re_z\;d\Omega}
            $$
            """
            self._L_p += RHS_mult * omega**2 * inner(v_re, eps * self.solution.Ediv_z_re) * dx

            r"""
            16
            $$
            N^{\Im\Im}_{\varepsilon z}=
            \omega^2\int_{\Omega}{v^\Im\varepsilon_0\varepsilon_r E^\Im_z\;d\Omega}
            $$
            """
            self._L_p += RHS_mult * omega**2 * inner(v_im, eps * self.solution.Ediv_z_im) * dx

            if sigma is not None:
                r"""
                2
                $$
                N^{\Im\Re}_{\sigma\perp}=
                \omega\int_{\Omega}{\vec{w}^\Re\sigma\vec{E}^\Im_\perp\;d\Omega}
                $$
                """
                self._L_p += RHS_mult * omega * inner(w_re, sigma * self.solution.Ediv_perp_im) * dx

                r"""
                5
                $$
                N^{\Re\Im}_{\sigma\perp}=
                -\omega\int_{\Omega}{\vec{w}^\Im\sigma\vec{E}^\Re_\perp\;d\Omega}
                $$
                """
                self._L_p += -RHS_mult * omega * inner(w_im, sigma * self.solution.Ediv_perp_re) * dx

                r"""
                12
                $$
                N^{\Im\Re}_{\sigma z}=
                \omega\int_{\Omega}{v^\Re\sigma E^\Im_z\;d\Omega}
                $$
                """
                self._L_p += RHS_mult * omega * inner(v_re, sigma * self.solution.Ediv_z_im) * dx

                r"""
                15
                $$
                N^{\Re\Im}_{\sigma z}=
                -\omega\int_{\Omega}{v^\Im\sigma E^\Re_z\;d\Omega}
                $$
                """
                self._L_p += -RHS_mult * omega * inner(v_im, sigma * self.solution.Ediv_z_re) * dx

        sibc_measures, self._sibc_deltas = self.__create_sibc_measures(sibc)
        if len(sibc_measures) > 0:
            n = -ufl.FacetNormal(self.material_map.mesh.mesh)
            tau = ufl.as_vector((-n[1], n[0]))
            for ds, (k, _) in zip(sibc_measures, self._sibc_deltas):
                if np.issubdtype(PETSc.ScalarType, np.complexfloating):
                    r"""
                    1
                    $$B_{\perp\perp}
                    =-\int_{\partial\Omega}{\left(\vec{w}\vec{\tau}\right)\left(\frac{1+j}{\delta\mu}\vec{E}_\perp\vec{\tau}\right)\;dS}$$
                    """
                    self._a_p += -k * PETSc.ScalarType(1 + 1j) * inner(dot(Eperp, tau), dot(w, tau)) * ds

                    r"""
                    4
                    $$
                    B_{zz}
                    =\int_{\partial\Omega}{v\frac{1+j}{\delta\mu}\underline{E}_z\;dS}
                    $$
                    """
                    self._a_p += k * PETSc.ScalarType(1 + 1j) * inner(Ez, v) * ds

                else:
                    r"""
                    1-1
                    $$
                    -\int_{\partial\Omega}{\left(\vec{w}^\Re\vec{\tau}\right)
                    \left(\frac{1}{\mu_0\delta}\vec{E}_\perp^\Re\vec{\tau}\right)\;dS}
                    $$
                    """
                    self._a_p += -k * inner(dot(w_re, tau), dot(Eperp_re, tau)) * ds

                    r"""
                    1-2
                    $$
                    -\int_{\partial\Omega}{\left(\vec{w}^\Im\vec{\tau}\right)
                    \left(\frac{1}{\mu_0\delta}\vec{E}_\perp^\Re\vec{\tau}\right)\;dS}
                    $$
                    """
                    self._a_p += -k * inner(dot(w_im, tau), dot(Eperp_re, tau)) * ds

                    r"""
                    6-1
                    $$
                    -\int_{\partial\Omega}{\left(\vec{w}^\Im\vec{\tau}\right)
                    \left(\frac{1}{\mu_0\delta}\vec{E}_\perp^\Im\vec{\tau}\right)\;dS}
                    $$
                    """
                    self._a_p += -k * inner(dot(w_im, tau), dot(Eperp_im, tau)) * ds

                    r"""
                    6-2
                    $$
                    \int_{\partial\Omega}{\left(\vec{w}^\Re\vec{\tau}\right)
                    \left(\frac{1}{\mu_0\delta}\vec{E}_\perp^\Im\vec{\tau}\right)\;dS}
                    $$
                    """
                    self._a_p += k * inner(dot(w_re, tau), dot(Eperp_im, tau)) * ds

                    r"""
                    11-1
                    $$
                    \int_{\partial\Omega}{v^\Re\left(\frac{1}{\mu_0\delta}\vec{E}_z^\Re\right)\;dS}
                    $$
                    """
                    self._a_p += k * inner(v_re, Ez_re) * ds

                    r"""
                    11-2
                    $$
                    \int_{\partial\Omega}{v^\Im\left(\frac{1}{\mu_0\delta}\vec{E}_z^\Re\right)\;dS}
                    $$
                    """
                    self._a_p += k * inner(v_im, Ez_re) * ds

                    r"""
                    16-1
                    $$
                    \int_{\partial\Omega}{v^\Im\left(\frac{1}{\mu_0\delta}\vec{E}_z^\Im\right)\;dS}
                    $$
                    """
                    self._a_p += k * inner(v_im, Ez_im) * ds

                    r"""
                    16-2
                    $$
                    -\int_{\partial\Omega}{v^\Re\left(\frac{1}{\mu_0\delta}\vec{E}_z^\Im\right)\;dS}
                    $$
                    """
                    self._a_p += -k * inner(v_re, Ez_im) * ds

        self._bcs = self.__set_dirichlet(V, sibc)

        self._A = dolfinx.fem.create_matrix(self._a_p)
        self._b = dolfinx.fem.create_vector(self._L_p)

        if type(self.solution.Ecurl_perp) != dolfinx.Function:
            self.solution._Ecurl = dolfinx.Function(V)

            if np.issubdtype(PETSc.ScalarType, np.complexfloating):
                # FIXME: workaround for https://github.com/FEniCS/dolfinx/issues/1577
                (self.solution.Ecurl_perp,
                 self.solution.Ecurl_z) = self.solution._Ecurl.split()
                (self.solution.ecurl_perp,
                 self.solution.ecurl_z) = ufl.split(self.solution._Ecurl)
            else:
                # FIXME: workaround for https://github.com/FEniCS/dolfinx/issues/1577
                self.solution._Ecurl = dolfinx.Function(V)
                (self.solution.Ecurl_perp_re,
                 self.solution.Ecurl_perp_im,
                 self.solution.Ecurl_z_re,
                 self.solution.Ecurl_z_im) = self.solution._Ecurl.split()
                (self.solution.ecurl_perp_re,
                 self.solution.ecurl_perp_im,
                 self.solution.ecurl_z_re,
                 self.solution.ecurl_z_im) = ufl.split(self.solution._Ecurl)

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

        for k, material in self._sibc_deltas:
            self.__set_skin_depth_from_material(k, material, self.solution._omega.value)

        self.solution._solve(self._a_p, self._L_p, self._A,
                             self._b, self.solution._Ecurl, bcs=self._bcs,
                             petsc_options=petsc_options)
        self.solution._Ecurl.x.array[:] /= simplify(self.RHS_mult)

        if MPI.COMM_WORLD.rank == 0:
            self.solution.logger.debug("Solved curl function")

        self.solution._Ecurl_stale = False
