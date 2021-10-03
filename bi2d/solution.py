"""Solution definition."""

import dolfinx
import ufl
from numpy import pi
from mpi4py import MPI
from petsc4py import PETSc


class Solution():
    """This class holds solutions to Maxwell's equations."""

    c0 = 299792458

    def __init__(self, material_map, Hcurl_order=1, H1_order=1):
        """Initialize."""
        self.material_map = material_map
        self.mesh = material_map.mesh
        self._omega = dolfinx.Constant(self.mesh.mesh, 2 * pi * self.material_map.f)
        self._beta = dolfinx.Constant(self.mesh.mesh, 1)
        self.Hcurl = ufl.FiniteElement("Nedelec 1st kind H(curl)", material_map.mesh.mesh.ufl_cell(), Hcurl_order)
        self.H1 = ufl.FiniteElement("Lagrange", material_map.mesh.mesh.ufl_cell(), H1_order)
        self.Js = None
        self.solver = PETSc.KSP().create(MPI.COMM_WORLD)
        self._solver_prefix = "dolfinx_solve_{}".format(id(self.solver))
        self.solver.setOptionsPrefix(self._solver_prefix)
        self.solver_options = PETSc.Options()
        # self.Phi_re = None
        # self.Phi_im = None
        # self.Ediv_perp_re = None
        # self.Ediv_perp_im = None
        # self.Ediv_z_re = None
        # self.Ediv_z_im = None
        # self.Ecurl_perp_re = None
        # self.Ecurl_perp_im = None
        # self.Ecurl_z_re = None
        # self.Ecurl_z_im = None
        self.q = 0.0
        self._reset_fields()

    def _set_solver_options(self, petsc_options={}):
        self.solver_options.prefixPush(self._solver_prefix)
        for k, v in petsc_options.items():
            self.solver_options[k] = v
        self.solver_options.prefixPop()
        self.solver.setFromOptions()


    def _reset_fields(self):
        for n in ["Phi_re", "Phi_im",
                  "Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im",
                  "Ecurl_perp_re", "Ecurl_perp_im", "Ecurl_z_re", "Ecurl_z_im"]:
            setattr(self, n, None)

    @property
    def f(self):
        """Get frequency point."""
        return self._omega.value / 2 / pi

    @f.setter
    def f(self, f):
        """Set frequency point."""
        self._omega.value = 2 * pi * f
        self.material_map.f = f
        self.material_map.update()
        # self._reset_fields()

    @property
    def beta(self):
        return self._beta.value

    @beta.setter
    def beta(self, beta):
        """Set beam particle speed."""
        self._beta.value = beta
        # self._reset_fields()

    def get_z(self):
        """Get impedance."""
        for name in ["Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im",
                     "Ecurl_perp_re", "Ecurl_perp_im", "Ecurl_z_re", "Ecurl_z_im"]:
            attr = getattr(self, name)
            if attr is None:
                raise AttributeError(f"Term {name} is not available")

        r"""
        $$
        \underline{Z}_\parallel(\omega)
        =-\frac{1}{q^2}\left(
        \int_{\Omega}{\vec{E}^\Re \vec{J}_\parallel^\Re \delta \;d\Omega}
        +j\int_{\Omega}{\vec{E}^\Im \vec{J}_\parallel^\Re \delta \;d\Omega}
        \right)
        $$
        $$
        \underline{Z}_\perp(\omega)
        =-\frac{\beta c}{(q d_r)^2\omega}\left(
        \int_{\Omega}{\vec{E}^\Re \vec{J}_\perp^\Re \delta \;d\Omega}
        +j\int_{\Omega}{\vec{E}^\Im \vec{J}_\perp^\Re \delta \;d\Omega}
        \right)
        $$
        """
        E_z_re = self.Ediv_z_re + self.ecurl_z_re
        E_z_im = self.Ediv_z_im + self.ecurl_z_im
        Zre = -1 / self.q ** 2 * \
            dolfinx.fem.assemble_scalar(self.material_map.beam * ufl.inner(E_z_re, self.Js) * ufl.dx)
        Zim = -1 / self.q ** 2 * \
            dolfinx.fem.assemble_scalar(self.material_map.beam * ufl.inner(E_z_im, self.Js) * ufl.dx)
        return (Zre, Zim)

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
        for n in ["Js", "Phi_re", "Phi_im",
                  "Ediv_perp_re", "Ediv_perp_im", "Ediv_z_re", "Ediv_z_im",
                  "Ecurl_perp_re", "Ecurl_perp_im", "Ecurl_z_re", "Ecurl_z_im"]:
            if getattr(self, n) is not None:
                self.xdmf_write_field(xdmf, n)
