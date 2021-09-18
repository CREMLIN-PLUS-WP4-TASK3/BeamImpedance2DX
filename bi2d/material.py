from inspect import signature
from numbers import Number
import dolfinx
import numpy as np
from mpi4py import MPI


class Material():
    """Material properties"""
    def __init__(self, index, name=None, eps=1, kappa=1, nu=(1, 0)):
        if not isinstance(index, int) or index < 1:
            raise(ValueError(f"""index must be a positive integer"""))
        self.index = index
        self.name = name if name is not None else f"material_{index}"
        if type(nu).__name__ == "tuple":
            nu_re, nu_im = nu
        else:
            nu_re = nu
            nu_im = 0
        for n in ["eps", "kappa", "nu_re", "nu_im"]:
            value = locals()[n]
            if not isinstance(value, Number) and \
               not (type(value).__name__ == "function" and \
                    (len(signature(value).parameters) == 2 or \
                     len(signature(value).parameters) == 3)):
                raise(ValueError(f"""Invalid value of argument {n}. \
Should be a number or a function with signature [float, float] -> float or \
[float, float, float] -> float"""))
            setattr(self, n, value)
            setattr(self, n + "_dispersive",
                    type(value).__name__ == "function" and
                    len(signature(value).parameters) == 3)
        self.real_nu = nu_im == 0

    def __str__(self):
        return f"""Material({self.name})[index={self.index},\
eps={self.eps},\
kappa={self.kappa},\
nu=({self.nu_re}, {self.nu_im})]"""


class MaterialMap():
    """Map material properties to mesh"""

    def __init__(self, mesh, materials, f=1e5, beam_subdomain_index=1):
        self.f = f
        self._f_previous = f
        self.mesh = mesh
        self.mesh_subdomains = np.unique(mesh.subdomains.values)
        material_indices = np.array([m.index for m in materials])

        for i in np.unique(material_indices):
            if np.sum(material_indices == i) > 1:
                raise(IndexError(f"Duplicate use of material index {i}"))

        for i in self.mesh_subdomains:
            if i not in material_indices:
                raise(IndexError(f"No material is set for mesh subdomain index {i}"))

        self.materials = [m for m in materials if m.index in material_indices]

        functionSpace = dolfinx.FunctionSpace(mesh.mesh, ("Discontinuous Lagrange", 0))
        self.beam = dolfinx.Function(functionSpace)
        self.beam.name = "beam"
        with self.beam.vector.localForm() as loc:
            cells = mesh.subdomains.indices[mesh.subdomains.values == beam_subdomain_index]
            loc.setValues(cells, np.full(len(cells), 1))
            cells = mesh.subdomains.indices[mesh.subdomains.values != beam_subdomain_index]
            loc.setValues(cells, np.full(len(cells), 0))

        for n in ["eps", "kappa", "nu_re"]:
            setattr(self, n, dolfinx.Function(functionSpace))
            getattr(self, n).name = n
            self.update_field(n)

        if not np.all([m.real_nu for m in self.materials]):
            self.nu_im = dolfinx.Function(functionSpace)
            self.nu_im.name = "nu_im"
            self.update_field("nu_im")
        else:
            self.nu_im = None

    def update_field(self, name):
        """Update material field map"""
        field = getattr(self, name)
        subdomains = self.mesh.subdomains
        for m in self.materials:
            value = getattr(m, name)
            index = m.index
            if isinstance(value, Number):
                with field.vector.localForm() as loc:
                    cells = subdomains.indices[subdomains.values == index]
                    loc.setValues(cells, np.full(len(cells), value))
            elif type(value).__name__ == "function":
                f = dolfinx.Function(field.ufl_function_space())
                if len(signature(value).parameters) == 2:
                    f.interpolate(lambda x: 0 * x[0] + 0 * x[1] + value(x[0], x[1]))
                elif len(signature(value).parameters) == 3:
                    f.interpolate(lambda x: 0 * x[0] + 0 * x[1] + value(x[0], x[1], self.f))
                else:
                    raise ValueError(f"Invalid {name} function signature for material {m.name}")
                dolfinx.cpp.la.scatter_forward(f.x)
                cells = subdomains.indices[subdomains.values == index]
                with f.vector.localForm() as loc_f:
                    val = loc_f.getValues(cells)
                with field.vector.localForm() as loc_field:
                    loc_field.setValues(cells, val)
            else:
                raise ValueError(f"Invalid {name} type for material {m.name}")

    def update(self):
        """Update material map lazily"""
        if self.f != self._f_previous:
            for n in ["eps", "kappa", "nu_re"]:
                if np.any([getattr(m, n + "_dispersive")
                           for m in self.materials]):
                    self.update_field(n)
            if self.nu_im != None and np.any([getattr(m, "nu_im_dispersive")
                                       for m in self.materials]):
                self.update_field("nu_im")
        self._f_previous = self.f

    def save(self, field_file: str):
        """Save material fields to XDMF file"""
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, field_file, "w") as xdmf:
            self.mesh.xdmf_write_mesh(xdmf)
            for n in ["eps", "kappa", "nu_re", "beam"]:
                self.xdmf_write_field(xdmf, n)
            if self.nu_im is not None:
                self.xdmf_write_field(xdmf, "nu_im")

    def xdmf_write_field(self, xdmf, field_name):
        """Set mesh in XDMF file"""
        xdmf.write_function(getattr(self, field_name))


        # for m in self.mesh_subdomains:
        #     material = materials[np.where(materials.index == m)[0]]
        #     for n in ["eps", "kappa", "nu_re", "nu_im"]:
        #         value = getattr(material, n)
        #         if type(value).__name__ == "function":
        # # f.interpolate(lambda x: (x[0]**2, x[1]*2))
        # f.interpolate(lambda x: x[0] + 5*x[1])
        # dolfinx.cpp.la.scatter_forward(f.x)
