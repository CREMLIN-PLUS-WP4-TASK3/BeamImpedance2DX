import numpy as np
from mpi4py import MPI
import dolfinx.io

# https://jorgensd.github.io/dolfinx-tutorial/chapter3/subdomains.html
def convert_msh(input_file: str, output_file: str) -> None:
    """Convert mesh to XDMF format"""
    import meshio
    import numpy
    msh = meshio.read(input_file)
    msh.prune_z_0()
    msh.remove_lower_dimensional_cells()
    msh.remove_orphaned_nodes()
    cells = msh.get_cells_type("triangle")
    try:
        cell_data = msh.get_cell_data("gmsh:physical", "triangle")
    except KeyError:
        cell_data = numpy.full(cells.shape[0], 1)
    out_mesh = meshio.Mesh(points=msh.points,
                           cells={"triangle": cells},
                           cell_data={"regions": [cell_data]})
    meshio.write(output_file, out_mesh, file_format="xdmf")

class Mesh (object):
    """Class holding mesh and subdomains information"""

    def __init__(self, mesh_file: str):
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_file, "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="Grid")
            self.subdomains = xdmf.read_meshtags(self.mesh, name="Grid")
            self.mesh.topology.create_connectivity(self.mesh.topology.dim, self.mesh.topology.dim-1)

    def save(self, mesh_file: str):
        """Save mesh to XDMF file"""
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_file, "w") as xdmf:
            self.xdmf_write_mesh(xdmf)

    def xdmf_write_mesh(self, xdmf):
        """Set mesh in XDMF file"""
        xdmf.write_mesh(self.mesh)
        xdmf.write_meshtags(self.subdomains)

    def get_bounds(self, index):
        """Get mesh subdomain bounds"""
        subdomain_indices = self.subdomains.indices[self.subdomains.values == index]
        geometry_indices = np.unique(np.array([self.mesh.geometry.dofmap.links(i)
                                               for i in subdomain_indices]).flatten())
        geometry_points = self.mesh.geometry.x[geometry_indices]
        minx = np.min(geometry_points[:,0])
        maxx = np.max(geometry_points[:,0])
        meanx = np.sum(geometry_points[:,0]) / geometry_points[:,0].size
        miny = np.min(geometry_points[:,1])
        maxy = np.max(geometry_points[:,1])
        meany = np.sum(geometry_points[:,1]) / geometry_points[:,1].size
        return ((minx, miny), (maxx, maxy), (meanx, meany))
