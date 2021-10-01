"""Mesh definitions."""

import numpy as np
from mpi4py import MPI
import dolfinx.io


def convert_msh(input_file: str, output_file: str,
                cell_type="triangle", cell_data_tag="gmsh:physical"):
    """Convert mesh to XDMF format.

    Based on https://jorgensd.github.io/dolfinx-tutorial/chapter3/subdomains.html
    """
    import meshio

    msh = meshio.read(input_file)
    msh.prune_z_0()
    msh.remove_orphaned_nodes()
    cells = msh.get_cells_type(cell_type)
    if cells.size == 0:
        raise AttributeError(f"There's no {cell_type} data in {input_file}")
    regions = msh.get_cell_data(cell_data_tag, cell_type)
    out_mesh = meshio.Mesh(points=msh.points,
                           cells={cell_type: cells},
                           cell_data={"regions": [regions]})
    meshio.write(output_file, out_mesh, file_format="xdmf")


class Mesh:
    """Class holding mesh, subdomain and boundary information."""

    def __init__(self, mesh_file, boundary_file=None):
        """Initialize class."""
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_file, "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="Grid")
            self.subdomains = xdmf.read_meshtags(self.mesh, name="Grid")
            # self.mesh.topology.create_connectivity(self.mesh.topology.dim,
            #                                        self.mesh.topology.dim-1)
            self.mesh.topology.create_connectivity(self.mesh.topology.dim-1,
                                                   self.mesh.topology.dim)
        self.boundaries = None
        if boundary_file is not None:
            with dolfinx.io.XDMFFile(MPI.COMM_WORLD, boundary_file, "r") as xdmf:
                self.boundaries = xdmf.read_meshtags(self.mesh, name="Grid")

    def save(self, mesh_file: str):
        """Save mesh to XDMF file."""
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_file, "w") as xdmf:
            self.xdmf_write_mesh(xdmf)

    def xdmf_write_mesh(self, xdmf):
        """Set mesh in XDMF file."""
        xdmf.write_mesh(self.mesh)
        xdmf.write_meshtags(self.subdomains)

    def get_boundary(self, index):
        """Get boundary with index."""
        return self.boundaries.indices[self.boundaries.values == index] if self.boundaries is not None else None

    def get_limits(self, index):
        """Get mesh subdomain limits."""
        subdomain_indices = self.subdomains.indices[self.subdomains.values == index]
        geometry_indices = np.unique(np.array([self.mesh.geometry.dofmap.links(i) for i in subdomain_indices])
                                     .flatten())
        geometry_points = self.mesh.geometry.x[geometry_indices]
        minx = np.min(geometry_points[:, 0])
        maxx = np.max(geometry_points[:, 0])
        meanx = np.sum(geometry_points[:, 0]) / geometry_points[:, 0].size
        miny = np.min(geometry_points[:, 1])
        maxy = np.max(geometry_points[:, 1])
        meany = np.sum(geometry_points[:, 1]) / geometry_points[:, 1].size
        return ((minx, miny), (maxx, maxy), (meanx, meany))
