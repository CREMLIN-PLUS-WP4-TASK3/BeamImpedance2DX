"""Solution definition."""

from binascii import crc32
import logging
import dolfinx
import ufl
from ufl import dx
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from .source import Js, SourceFunction
from .poisson import Ediv
from .curl import Ecurl
from .esum import Esum


class Solution():
    """This class holds solutions to Maxwell's equations."""

    c0 = PETSc.ScalarType(299792458)

    def __init__(self, material_map, H1_order=1, Hcurl_order=1, logger_name=None):
        """Initialize."""
        self.logger = logging.getLogger(logger_name if logger_name is not None
                                        else f"{type(self).__name__}_{crc32(str(id(self)).encode()):x}")
        # Mesh and materials
        self.material_map = material_map
        self.mesh = material_map.mesh
        # Scalar variables
        self._omega = dolfinx.Constant(self.mesh.mesh, PETSc.ScalarType(2 * np.pi * self.material_map.f))
        self._beta = dolfinx.Constant(self.mesh.mesh, PETSc.ScalarType(1))
        # FEM elements
        self.H1 = ufl.FiniteElement("Lagrange", material_map.mesh.mesh.ufl_cell(), H1_order)
        self.Hcurl = ufl.FiniteElement("Nedelec 1st kind H(curl)", material_map.mesh.mesh.ufl_cell(), Hcurl_order)
        # Solver setup
        self.solver = PETSc.KSP().create(MPI.COMM_WORLD)
        self._solver_prefix = "bi2d_solver"
        self.solver.setOptionsPrefix(self._solver_prefix)
        self.solver_options = PETSc.Options()
        self.solver_viewer = PETSc.Viewer()
        # Functions and function flags
        self.Js = None
        self._phi = None
        self.Phi = None
        self.Phi_re = None
        self.Phi_im = None
        self.Ediv_perp = None
        self.Ediv_perp_re = None
        self.Ediv_perp_im = None
        self.Ediv_z = None
        self.Ediv_z_re = None
        self.Ediv_z_im = None
        self._Ecurl = None
        self.Ecurl_perp = None
        self.Ecurl_perp_re = None
        self.Ecurl_perp_im = None
        self.Ecurl_z = None
        self.Ecurl_z_re = None
        self.Ecurl_z_im = None
        self.Eperp = None
        self.Eperp_re = None
        self.Eperp_im = None
        self.Ez = None
        self.Ez_re = None
        self.Ez_im = None
        self._Js_stale = True
        self._Ediv_stale = True
        self._Ecurl_stale = True
        self._E_stale = True
        # Solution type flag
        self.source_function = None
        # Monopole charge
        self.q = 0.0
        # Dipole moment
        self.dm = 0.0
        # Quadrupole moment
        self.qm = 0.0
        # get_z attributes
        self.__sibc_hash = None
        self.__source_hash = None
        self.__Js_solver = None
        self.__Ediv_solver = None
        self.__Ecurl_solver = None
        # sum_fields attributes
        self.__Esum_solver = None
        self.__sum_order_hash = None

    def _set_solver_options(self, petsc_options={}):
        self.solver_options.prefixPush(self._solver_prefix)
        for k, v in petsc_options.items():
            self.solver_options[k] = v
        self.solver_options.prefixPop()
        self.solver.setFromOptions()

    def _solve(self, a_p, L_p, A, b, f, bcs=[],
               petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}):
        self.solver.reset()
        self._set_solver_options(petsc_options)
        self.solver.setOperators(A)

        # Assemble lhs
        A.zeroEntries()
        dolfinx.fem.assemble_matrix(A, a_p, bcs=bcs)
        A.assemble()

        # Assemble rhs
        with b.localForm() as b_loc:
            b_loc.set(0)
        dolfinx.fem.assemble_vector(b, L_p)

        # Apply boundary conditions to the rhs
        dolfinx.fem.apply_lifting(b, [a_p], [bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.set_bc(b, bcs)

        # Solve linear system and update ghost values in the solution
        self.solver.solve(b, f.vector)
        f.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

    @property
    def f(self):
        """Get frequency point."""
        return self._omega.value / 2 / np.pi

    @f.setter
    def f(self, f):
        """Set frequency point."""
        if not np.isclose(self._omega.value, 2 * np.pi * f):
            self._omega.value = PETSc.ScalarType(2 * np.pi * f)
            self.material_map.f = f
            self.material_map.update()
            for name in ["Ediv", "Ecurl", "E"]:
                setattr(self, f"_{name}_stale", True)

    @property
    def beta(self):
        """Get beam particle speed."""
        return self._beta.value

    @beta.setter
    def beta(self, beta):
        """Set beam particle speed."""
        if not np.isclose(self._beta.value, beta):
            self._beta.value = PETSc.ScalarType(beta)
            for name in ["Ediv", "Ecurl", "E"]:
                setattr(self, f"_{name}_stale", True)

    @property
    def z(self):
        """Get impedance."""
        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            for name in ["Js", "Ediv_perp", "Ediv_z", "Ecurl_perp", "Ecurl_z"]:
                attr = getattr(self, name)
                if attr is None:
                    raise AttributeError(f"{name} function is not available")
            for name in ["Js", "Ediv", "Ecurl"]:
                attr = getattr(self, f"_{name}_stale")
                if attr:
                    raise ValueError(f"{name} solutions are stale and need to be recalculated")
        else:
            for name in ["Js", "Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im",
                         "Ecurl_perp_re", "Ecurl_perp_im", "Ecurl_z_re", "Ecurl_z_im"]:
                attr = getattr(self, name)
                if attr is None:
                    raise AttributeError(f"{name} function is not available")
            for name in ["Js", "Ediv", "Ecurl"]:
                attr = getattr(self, f"_{name}_stale")
                if attr:
                    raise ValueError(f"{name} solutions are stale and need to be recalculated")

        if np.issubdtype(PETSc.ScalarType, np.complexfloating):
            E_z = self.Ediv_z + self.ecurl_z
            if self.source_function == SourceFunction.MONOPOLE:
                r"""
                $$
                \underline{Z}_\parallel(\omega)
                =-\frac{1}{q^2}\int_{\Omega}{\vec{E}\vec{J}_\parallel^* \delta \;d\Omega}
                $$
                """
                val = dolfinx.fem.assemble_scalar(-1 / self.q ** 2 *
                                                  ufl.inner(E_z, ufl.conj(self.Js)) * dx)
                Zre = np.real(val)
                Zim = np.imag(val)
            elif self.source_function == SourceFunction.DIPOLE:
                r"""
                $$
                \underline{Z}_\perp(\omega)
                =-\frac{\beta c}{d_r^2\omega}
                \int_{\Omega}{\vec{E} \vec{J}_\perp^* \delta \;d\Omega}
                $$
                """
                val = dolfinx.fem.assemble_scalar(-self._beta * self.c0 / self.dm**2 / self._omega *
                                                  ufl.inner(E_z, ufl.conj(self.Js)) * dx)
                Zre = np.real(val)
                Zim = np.imag(val)
            elif self.source_function == SourceFunction.QUADRUPOLE:
                r"""
                $$
                \underline{Z}_\perp(\omega)
                =-\frac{\beta c}{d_r^2\omega}
                \int_{\Omega}{\vec{E} \vec{J}_\perp^* \delta \;d\Omega}
                $$
                """
                val = dolfinx.fem.assemble_scalar(-self._beta * self.c0 / self.qm**2 / self._omega *
                                                  ufl.inner(E_z, ufl.conj(self.Js)) * dx)
                Zre = np.real(val)
                Zim = np.imag(val)
            else:
                raise AttributeError("Unknown source function")
        else:
            E_z_re = self.Ediv_z_re + self.ecurl_z_re
            E_z_im = self.Ediv_z_im + self.ecurl_z_im
            if self.source_function == SourceFunction.MONOPOLE:
                r"""
                $$
                \underline{Z}_\parallel(\omega)
                =-\frac{1}{q^2}\left(
                \int_{\Omega}{\vec{E}^\Re \vec{J}_\parallel^\Re \delta \;d\Omega}
                +j\int_{\Omega}{\vec{E}^\Im \vec{J}_\parallel^\Re \delta \;d\Omega}
                \right)
                $$
                """
                Zre = dolfinx.fem.assemble_scalar(-1 / self.q ** 2 *
                                                  ufl.inner(E_z_re, self.Js) * dx)
                Zim = dolfinx.fem.assemble_scalar(-1 / self.q ** 2 *
                                                  ufl.inner(E_z_im, self.Js) * dx)
            elif self.source_function == SourceFunction.DIPOLE:
                r"""
                $$
                \underline{Z}_\perp(\omega)
                =-\frac{\beta c}{d_r^2\omega}\left(
                \int_{\Omega}{\vec{E}^\Re \vec{J}_\perp^\Re \delta \;d\Omega}
                +j\int_{\Omega}{\vec{E}^\Im \vec{J}_\perp^\Re \delta \;d\Omega}
                \right)
                $$
                """
                Zre = dolfinx.fem.assemble_scalar(-self._beta * self.c0 / self.dm2**2 / self._omega *
                                                  ufl.inner(E_z_re, self.Js) * dx)
                Zim = dolfinx.fem.assemble_scalar(-self._beta * self.c0 / self.dm2**2 / self._omega *
                                                  ufl.inner(E_z_im, self.Js) * dx)
            elif self.source_function == SourceFunction.QUADRUPOLE:
                r"""
                $$
                \underline{Z}_\perp(\omega)
                =-\frac{\beta c}{d_r^2\omega}\left(
                \int_{\Omega}{\vec{E}^\Re \vec{J}_\perp^\Re \delta \;d\Omega}
                +j\int_{\Omega}{\vec{E}^\Im \vec{J}_\perp^\Re \delta \;d\Omega}
                \right)
                $$
                """
                Zre = dolfinx.fem.assemble_scalar(-self._beta * self.c0 / self.qm**2 / self._omega *
                                                  ufl.inner(E_z_re, self.Js) * dx)
                Zim = dolfinx.fem.assemble_scalar(-self._beta * self.c0 / self.qm**2 / self._omega *
                                                  ufl.inner(E_z_im, self.Js) * dx)
            else:
                raise AttributeError("Unknown source function")
        Zre = MPI.COMM_WORLD.allreduce(Zre)
        Zim = MPI.COMM_WORLD.allreduce(Zim)

        return (Zre, Zim)

    def get_z(self, f, beta=1.0, gamma=None, rotation=0, sibc=[],
              source_function=SourceFunction.MONOPOLE,
              petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}):
        """Get impedance. High level interface to the library functionality."""
        if (beta is None and gamma is None) or (beta is not None and gamma is not None):
            raise AttributeError("Must supply either beta or gamma value")
        elif gamma is not None:
            # $\beta = \sqrt{1-\gamma^{-2}}$
            beta = np.sqrt(1 - 1 / (gamma**2))
        if self.__Js_solver is None or self.__source_hash != hash((source_function, rotation)):
            self.__Js_solver = Js(self, rotation=rotation,
                                  source_function=source_function)
            self.__source_hash = hash((source_function, rotation))
            self.__Ediv_solver = None
        if self.__Ediv_solver is None:
            self.__Ediv_solver = Ediv(self)
            self.__Ecurl_solver = None
        if self.__Ecurl_solver is None or hash(tuple(sibc)) != self.__sibc_hash:
            self.__Ecurl_solver = Ecurl(self, sibc=sibc)
            self.__sibc_hash = hash(tuple(sibc))
        self.beta = beta
        f = np.array(f)

        output = np.zeros((f.size, 3))
        for i, f in enumerate(f):
            if MPI.COMM_WORLD.rank == 0:
                if beta == 1:
                    gamma = "Inf"
                else:
                    # $\gamma = \frac{1}{\sqrt{1-\beta^2}}$
                    gamma = 1/np.sqrt(1-beta**2)
                self.logger.info(f"Solving for f={f:.2e}, ??={beta:.2f}, ??={gamma:.2f}")
            self.f = f
            if self._Js_stale:
                self.__Js_solver.solve(petsc_options=petsc_options)
            if self._Ediv_stale:
                self.__Ediv_solver.solve(petsc_options=petsc_options)
            if self._Ecurl_stale:
                self.__Ecurl_solver.solve(petsc_options=petsc_options)
            z_re, z_im = self.z
            output[i][0] = f
            output[i][1] = z_re
            output[i][2] = z_im

        return output

    def sum_fields(self, H1_order=2, Hcurl_order=2,
                   petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}):
        """Sum irrotational and solenoidal fields. High level interface to the library functionality."""
        if self.__Esum_solver is None or self.__sum_order_hash != hash((H1_order, Hcurl_order)):
            self.Eperp = None
            self.Eperp_re = None
            self.Eperp_im = None
            self.Ez = None
            self.Ez_re = None
            self.Ez_im = None
            self.__Esum_solver = Esum(self, H1_order=H1_order, Hcurl_order=Hcurl_order)
            self.__sum_order_hash = hash((H1_order, Hcurl_order))

        if self._E_stale:
            self.__Esum_solver.solve(petsc_options=petsc_options)

    def save(self, field_file: str):
        """Save solution to XDMF file."""
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, field_file, "w") as xdmf:
            self.material_map.mesh.xdmf_write_mesh(xdmf)
            self.material_map.xdmf_write_all(xdmf)
            self.xdmf_write_all(xdmf)

    def xdmf_write_field(self, xdmf: dolfinx.io.XDMFFile, field_name: str):
        """Write field to XDMF file."""
        field = getattr(self, field_name)
        field.name = field_name
        xdmf.write_function(field)

    def xdmf_write_all(self, xdmf: dolfinx.io.XDMFFile):
        """Save material fields to XDMF file."""
        if not self._Js_stale:
            for n in ["Js"]:
                if getattr(self, n) is not None:
                    self.xdmf_write_field(xdmf, n)
        if not self._Ediv_stale:
            for n in ["Phi_re", "Phi_im", "Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im",
                      "Phi", "Ediv_perp", "Ediv_z"]:
                if getattr(self, n) is not None:
                    self.xdmf_write_field(xdmf, n)
        if not self._Ecurl_stale:
            for n in ["Ecurl_perp_re", "Ecurl_perp_im", "Ecurl_z_re", "Ecurl_z_im",
                      "Ecurl_perp", "Ecurl_z"]:
                if getattr(self, n) is not None:
                    self.xdmf_write_field(xdmf, n)
        if not self._E_stale:
            for n in ["Eperp_re", "Eperp_im", "Ez_re", "Ez_im",
                      "Eperp", "Ez"]:
                if getattr(self, n) is not None:
                    self.xdmf_write_field(xdmf, n)
