import dolfinx
# from mpi4py import MPI
# import numpy as np
# import ufl


class Solution():
    """Solution class"""

    def __init__(self, material_map, Hcurl_order = 1, H1_order = 1):
        self.material_map = material_map
        self.Hcurl= dolfinx.FunctionSpace(material_map.mesh.mesh,
                                          ("Nedelec 1st kind H(curl)", Hcurl_order))
        self.H1= dolfinx.FunctionSpace(material_map.mesh.mesh,
                                       ("Lagrange", H1_order))

        # help(material_map.mesh.mesh.ufl_domain())
        # f = dolfinx.Function(material_map.eps.ufl_function_space())
        # indices = dolfinx.fem.locate_dofs_geometrical(
        #     material_map.eps.ufl_function_space(),
        #     lambda x: np.isclose((x[0]-0.01)**2+(x[1]-0.01)**2, 0.1) +
        #     np.sqrt((x[0]-0.01)**2+(x[1]-0.01)**2) <= 0.1
        # )
        # f.name = "aaa"
        # with f.vector.localForm() as loc_field:
        #     loc_field.setValues(indices, np.full(indices.size, 2))

        # # material_map.mesh.subdomains.ufl_function_space()
        # # f = dolfinx.Function(material_map.mesh.subdomains.ufl_function_space())
        # u_, v_ = ufl.TrialFunction(self.H1), ufl.TestFunction(self.H1)
        # a_p = ufl.inner(u_, v_) * ufl.dx
        # L_p = ufl.inner(f, v_) * ufl.dx
        # projection = dolfinx.fem.LinearProblem(a_p, L_p)
        # ff = projection.solve()
        # with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "field.xdmf", "w") as xdmf:
        #     xdmf.write_mesh(material_map.mesh.mesh)
        #     xdmf.write_function(ff)
