from mpi4py import MPI
import dolfinx
import ufl

class Ediv():
    """Irrotational electric field solver"""

    def __init__(self, solution):
        self.material_map = solution.material_map
        self.H1 = solution.H1
        self.E_re = dolfinx.Function(self.H1)
        self.E_im = dolfinx.Function(self.H1)
        self.E_re.name = "Ediv_re"
        self.E_im.name = "Ediv_im"

    def solve(self):
        H1= dolfinx.FunctionSpace(self.material_map.mesh.mesh, ("Lagrange", 1))
        u_, v_ = ufl.TrialFunction(H1), ufl.TestFunction(H1)
        a_p = ufl.inner(u_, v_) * ufl.dx
        L_p = self.material_map.beam * v_ * ufl.dx
        projection = dolfinx.fem.LinearProblem(a_p, L_p)
        self.f = projection.solve()
        self.f.name = "current"

    def save(self, field_file: str):
        """Save current field to XDMF file"""
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, field_file, "w") as xdmf:
            self.material_map.mesh.xdmf_write_mesh(xdmf)
            self.xdmf_write_field(xdmf)

    def xdmf_write_field(self, xdmf):
        """Set beam current field in XDMF file"""
        xdmf.write_function(self.f)
