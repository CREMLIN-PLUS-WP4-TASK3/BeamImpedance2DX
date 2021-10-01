"""Solution definition."""

import dolfinx
from mpi4py import MPI
import ufl


class Solution():
    """This class holds solutions to Maxwell's equations."""

    c0 = 299792458

    def __init__(self, material_map, beta=1, Hcurl_order=1, H1_order=1):
        """Initialize."""
        self.material_map = material_map
        self.mesh = material_map.mesh
        self.f = self.material_map.f
        self.beta = beta
        self.Hcurl = ufl.FiniteElement("Nedelec 1st kind H(curl)", material_map.mesh.mesh.ufl_cell(), Hcurl_order)
        self.H1 = ufl.FiniteElement("Lagrange", material_map.mesh.mesh.ufl_cell(), H1_order)
        self.Js = None
        self.q = 1.0
        self._reset_fields()

    def _reset_fields(self):
        for n in ["Phi_re", "Phi_im",
                  "Ediv_re_perp", "Ediv_im_perp", "Ediv_re_z", "Ediv_im_z",
                  "Ecurl_re_perp", "Ecurl_im_perp", "Ecurl_re_z", "Ecurl_im_z"]:
            setattr(self, n, None)

    def set_frequency(self, f):
        """Set frequency point."""
        self.f = float(f)
        self.material_map.f = self.f
        self.material_map.update()
        self._reset_fields()

    def set_beta(self, beta):
        """Set beam particle speed."""
        self.beta = float(beta)
        self._reset_fields()

    def get_z(self):
        """Get impedance."""
        for name in ["Ediv_re_perp", "Ediv_im_perp", "Ediv_re_z", "Ediv_im_z",
                     "Ecurl_re_perp", "Ecurl_im_perp", "Ecurl_re_z", "Ecurl_im_z"]:
            attr = getattr(self, name)
            if attr is None:
                raise ValueError(f"Term {name} is not available")

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
        E_z_re = self.Ediv_re_z + self.ecurl_re_z
        E_z_im = self.Ediv_im_z + self.ecurl_im_z
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
                  "Ediv_re_perp", "Ediv_im_perp", "Ediv_re_z", "Ediv_im_z",
                  "Ecurl_re_perp", "Ecurl_im_perp", "Ecurl_re_z", "Ecurl_im_z"]:
            if getattr(self, n) is not None:
                self.xdmf_write_field(xdmf, n)
