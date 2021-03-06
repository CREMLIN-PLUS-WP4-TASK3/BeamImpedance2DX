"""Material definitions."""

from copy import deepcopy
from inspect import signature
from numbers import Number
import dolfinx
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


class Material():
    """Material properties."""

    def __init__(self, index, name=None, eps_r=1, sigma=0, mu_r_re=1, mu_r_im=0):
        """Initialize."""
        self.__setattr("index", index)
        self.__setattr("name", name if name is not None else f"material_{index}")
        for n in ["eps_r", "sigma", "mu_r_re", "mu_r_im"]:
            value = locals()[n]
            if not isinstance(value, Number) and \
               not ((type(value).__name__ == "function" or
                     type(value).__name__ == "method") and
                    (len(signature(value).parameters) >= 1 or
                     len(signature(value).parameters) <= 3)):
                raise(AttributeError(f"""Invalid value of argument {n}. \
Should be a number or a function with signature [freq:float] -> val:float or \
[x:float, y:float] -> val:float or [x:float, y:float, freq:float] -> val:float"""))
            self.__setattr(n, value)
            self.__setattr(n + "_dispersive",
                           (type(value).__name__ == "function" or
                            type(value).__name__ == "method") and
                           (len(signature(value).parameters) == 1 or
                            len(signature(value).parameters) == 3))
        self.__setattr("real_mu", mu_r_im == 0)
        self.__setattr("real_eps", sigma == 0)

    def __setattr__(self, name, value):
        if name == "index":
            self.__setattr(name, value)
        else:
            raise AttributeError(f"Cannot set `{name}` attribute. All attributes except for `index` are read only.")

    def __setattr(self, name, value):
        super().__setattr__(name, value)

    def copy(self):
        """Create a copy of a material"""
        return deepcopy(self)

    def __str__(self):
        """String representation."""
        return f"""Material({self.name})[index={self.index},\
eps_r={self.eps_r:.2e},\
sigma={self.sigma:.2e},\
mu=({self.mu_r_re:.2e}, {self.mu_r_im:.2e})]"""


class MaterialMapBase():
    """Map material properties to mesh."""

    eps0 = PETSc.ScalarType(8.8541878128e-12)
    mu0 = PETSc.ScalarType(4 * np.pi * 1e-7)

    def __init__(self, mesh, materials, f=1e5, beam_subdomain_index=1, round_check=True):
        """Initialize."""
        self.f = f
        self._f_previous = f
        self.beam_index = beam_subdomain_index
        self.mesh = mesh
        self.mesh_subdomains = np.unique(mesh.subdomains.values)
        material_indices = np.array([m.index for m in materials])

        for i in np.unique(material_indices):
            if np.sum(material_indices == i) > 1:
                names = ", ".join([m.name for m in materials if m.index == i])
                raise(IndexError(f"Duplicate use of material index {i} ({names})"))

        for i in self.mesh_subdomains:
            if i not in material_indices:
                raise(IndexError(f"No material is set for mesh subdomain index {i}"))

        self.materials = [m for m in materials if m.index in self.mesh_subdomains]

        self.V = dolfinx.FunctionSpace(mesh.mesh, ("Discontinuous Lagrange", 0))
        self.beam = dolfinx.Function(self.V)
        self.beam.name = "beam"
        with self.beam.vector.localForm() as loc:
            loc.set(0)
            cells = mesh.subdomains.indices[mesh.subdomains.values == self.beam_index]
            loc.setValues(cells, np.full(cells.size, 1))

        self.mpiproc = dolfinx.Function(self.V)
        self.mpiproc.name = "MPI_processor"
        with self.mpiproc.vector.localForm() as loc:
            loc.set(MPI.COMM_WORLD.rank)

        (minx, miny), (maxx, maxy) = self.mesh.get_limits(self.beam_index)
        if round_check:
            if not np.isclose(maxx - minx, maxy - miny, atol=(maxx - minx) / 200) or\
               not self.mesh.check_round((minx + maxx) / 2, (miny + maxy) / 2, (maxx - minx) / 2, self.beam_index):
                raise ValueError(f"Beam subdomain {beam_subdomain_index} does not appear to be round")

    def update_field(self, name):
        """Update material field map."""
        field = getattr(self, name)
        subdomains = self.mesh.subdomains
        for m in self.materials:
            value = getattr(m, name)
            index = m.index
            if isinstance(value, Number):
                with field.vector.localForm() as loc:
                    cells = subdomains.indices[subdomains.values == index]
                    loc.setValues(cells, np.full(cells.size, value))
            elif type(value).__name__ == "function" or type(value).__name__ == "method":
                # We have to use temporary function to avoid overwriting data from other subdomains
                f = dolfinx.Function(field.function_space)
                if len(signature(value).parameters) == 1:
                    f.interpolate(lambda x: 0 * x[0] + 0 * x[1] + value(self.f))
                elif len(signature(value).parameters) == 2:
                    f.interpolate(lambda x: 0 * x[0] + 0 * x[1] + value(x[0], x[1]))
                elif len(signature(value).parameters) == 3:
                    f.interpolate(lambda x: 0 * x[0] + 0 * x[1] + value(x[0], x[1], self.f))
                else:
                    raise ValueError(f"Invalid {name} function signature for material {m.name}")
                f.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
                cells = subdomains.indices[subdomains.values == index]
                with f.vector.localForm() as loc_f, field.vector.localForm() as loc_field:
                    val = loc_f.getValues(cells)
                    loc_field.setValues(cells, val)
            else:
                raise ValueError(f"Invalid {name} type for material {m.name}")

    def update(self):
        """Update material map lazily."""
        if self.f != self._f_previous:
            self._update()
            self._f_previous = self.f

    def save(self, field_file: str):
        """Save material fields to XDMF file."""
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, field_file, "w") as xdmf:
            self.mesh.xdmf_write_mesh(xdmf)
            self.xdmf_write_all(xdmf)

    def xdmf_write_field(self, xdmf: dolfinx.io.XDMFFile, field_name: str):
        """Write field to XDMF file."""
        xdmf.write_function(getattr(self, field_name))

    def xdmf_write_all(self, xdmf: dolfinx.io.XDMFFile):
        """Save material fields to XDMF file."""
        self.update()
        for n in self._xdmf_write_all_list:
            if getattr(self, n) is not None:
                self.xdmf_write_field(xdmf, n)


class MaterialMapScalar(MaterialMapBase):
    """Map material properties to mesh."""

    def __init__(self, mesh, materials, f=1e5, beam_subdomain_index=1, round_check=True):
        """Initialize."""
        super().__init__(mesh, materials, f, beam_subdomain_index, round_check=round_check)
        self._xdmf_write_all_list = ["beam", "mpiproc", "eps_r", "eps", "sigma", "mu_r_re",
                                     "mu_r_im", "mu_re", "mu_im", "nu_re", "nu_im"]

        for n in ["eps_r", "mu_r_re"]:
            setattr(self, n, dolfinx.Function(self.V))
            getattr(self, n).name = n
            self.update_field(n)

        for (n, flag) in [("sigma", "real_eps"), ("mu_r_im", "real_mu")]:
            if not np.all([getattr(m, flag) for m in self.materials]):
                setattr(self, n, dolfinx.Function(self.V))
                getattr(self, n).name = n
                self.update_field(n)
            else:
                setattr(self, n, None)

        for n in ["eps", "mu_re", "nu_re"]:
            setattr(self, n, dolfinx.Function(self.V))
            getattr(self, n).name = n
        if self.mu_r_im is not None:
            for n in ["mu_im", "nu_im"]:
                setattr(self, n, dolfinx.Function(self.V))
                getattr(self, n).name = n
        else:
            self.mu_im = None
            self.nu_im = None
        self.calculate_eps()
        self.calculate_mu()
        self.calculate_nu()

    def calculate_eps(self):
        """Calculate epsilon function."""
        with self.eps_r.vector.localForm() as eps_r, self.eps.vector.localForm() as eps:
            val = eps_r * self.eps0
            val.copy(eps)

    def calculate_mu(self):
        """Calculate mu function."""
        with self.mu_r_re.vector.localForm() as mu_r_re, self.mu_re.vector.localForm() as mu_re:
            val = mu_r_re * self.mu0
            val.copy(mu_re)
        if self.mu_r_im is not None:
            with self.mu_r_im.vector.localForm() as mu_r_im, self.mu_im.vector.localForm() as mu_im:
                val = mu_r_im * self.mu0
                val.copy(mu_im)

    def calculate_nu(self):
        """Calculate nu functions."""
        if self.mu_im is None:
            with self.mu_re.vector.localForm() as mu_re, \
                 self.nu_re.vector.localForm() as nu_re:
                val = 1 / mu_re
                val.copy(nu_re)
        else:
            with self.mu_re.vector.localForm() as mu_re, \
                 self.mu_im.vector.localForm() as mu_im, \
                 self.nu_re.vector.localForm() as nu_re, \
                 self.nu_im.vector.localForm() as nu_im:
                val = mu_re / (mu_re * mu_re + mu_im * mu_im)
                val.copy(nu_re)
                val = mu_im / (mu_re * mu_re + mu_im * mu_im)
                val.copy(nu_im)

    def _update(self):
        """Update material map."""
        recalculate_nu = False
        for n in ["eps_r", "sigma", "mu_r_re", "mu_r_im"]:
            if getattr(self, n) is not None and \
               np.any([getattr(m, n + "_dispersive") for m in self.materials]):
                self.update_field(n)
                if n == "eps_r":
                    self.calculate_eps()
                elif n == "mu_r_re" or n == "mu_r_im":
                    self.calculate_mu()
                    recalculate_nu = True
        if recalculate_nu:
            self.calculate_nu()


class MaterialMapComplex(MaterialMapBase):
    """Map material properties to mesh."""

    def __init__(self, mesh, materials, f=1e5, beam_subdomain_index=1, round_check=True):
        """Initialize."""
        super().__init__(mesh, materials, f, beam_subdomain_index, round_check=round_check)
        self._xdmf_write_all_list = ["beam", "mpiproc", "eps_r", "eps", "sigma", "mu_r_re", "mu_r_im", "mu"]

        for n in ["eps_r", "sigma", "mu_r_re", "mu_r_im"]:
            setattr(self, n, dolfinx.Function(self.V))
            getattr(self, n).name = n
            self.update_field(n)

        for n in ["eps", "mu"]:
            setattr(self, n, dolfinx.Function(self.V))
            getattr(self, n).name = n
        self.calculate_eps()
        self.calculate_mu()

    def calculate_eps(self):
        """Calculate epsilon function."""
        with self.eps_r.vector.localForm() as eps_r, \
             self.eps.vector.localForm() as eps, \
             self.sigma.vector.localForm() as sigma:
            val = eps_r * self.eps0 - 1j * sigma / (2 * np.pi * self.f)
            val.copy(eps)

    def calculate_mu(self):
        """Calculate mu function."""
        with self.mu_r_re.vector.localForm() as mu_r_re, \
             self.mu_r_im.vector.localForm() as mu_r_im, \
             self.mu.vector.localForm() as mu:
            val = (mu_r_re - 1j * mu_r_im) * self.mu0
            val.copy(mu)

    def _update(self):
        """Update material map."""
        for n in ["eps_r", "sigma", "mu_r_re", "mu_r_im"]:
            if getattr(self, n) is not None and \
               np.any([getattr(m, n + "_dispersive") for m in self.materials]):
                self.update_field(n)
                if n == "eps_r" or n == "sigma":
                    self.calculate_eps()
                elif n == "mu_r_re" or n == "mu_r_im":
                    self.calculate_mu()
        if self.f != self._f_previous:
            self.calculate_eps()


class ArrayInterpolate():
    """Interpolate data from array."""

    def __init__(self, x, y):
        """Initialize."""
        self.data_x = x
        self.data_y = y

    @classmethod
    def from_numpy(cls, array, xi=0, yi=1):
        """Interpolate data from numpy array."""
        return cls(array[:, 0], array[:, 1])

    @classmethod
    def from_file(cls, f, xi=0, yi=1, delimiter=None):
        """Interpolate data from text file."""
        data = np.genfromtxt(f, delimiter=delimiter)
        return cls(data[:, xi], data[:, yi])

    def interp(self, x):
        """Interpolate point."""
        return np.interp([x], self.data_x, self.data_y)[0]


if np.issubdtype(PETSc.ScalarType, np.complexfloating):
    MaterialMap = MaterialMapComplex
else:
    MaterialMap = MaterialMapScalar
