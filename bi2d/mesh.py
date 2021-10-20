"""Mesh definitions."""

import numpy as np
from mpi4py import MPI
import dolfinx.io


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
        if self.boundaries is None:
            raise AttributeError("Boundary data was not loaded")
        return self.boundaries.indices[self.boundaries.values == index]

    def get_limits(self, index):
        """Get mesh subdomain limits."""
        subdomain_indices = self.subdomains.indices[self.subdomains.values == index]
        geometry_indices = np.unique(np.array([self.mesh.geometry.dofmap.links(i) for i in subdomain_indices])
                                     .flatten())
        if geometry_indices.size > 0:
            geometry_points = self.mesh.geometry.x[geometry_indices]
            minx = np.min(geometry_points[:, 0])
            maxx = np.max(geometry_points[:, 0])
            miny = np.min(geometry_points[:, 1])
            maxy = np.max(geometry_points[:, 1])
        else:
            minx = np.nan
            maxx = np.nan
            miny = np.nan
            maxy = np.nan
        values = MPI.COMM_WORLD.gather(np.array([minx, miny, maxx, maxy]))
        if MPI.COMM_WORLD.rank == 0:
            values = [v for v in values if np.all(~np.isnan(v))]
            values = np.vstack(values)
            minx = np.min(values[:, 0])
            miny = np.min(values[:, 1])
            maxx = np.max(values[:, 2])
            maxy = np.max(values[:, 3])
        return MPI.COMM_WORLD.bcast(((minx, miny), (maxx, maxy)), root=0)

    def check_round(self, x0, y0, R, index):
        """Check all points are inside circle."""
        subdomain_indices = self.subdomains.indices[self.subdomains.values == index]
        geometry_indices = np.unique(np.array([self.mesh.geometry.dofmap.links(i) for i in subdomain_indices])
                                     .flatten())
        if geometry_indices.size > 0:
            geometry_points = self.mesh.geometry.x[geometry_indices]
            x = geometry_points[:, 0]
            y = geometry_points[:, 1]
            return np.all(np.logical_or((x-x0)**2 + (y-y0)**2 < R**2,
                                        np.isclose((x-x0)**2 + (y-y0)**2, R**2)))
        else:
            return True
