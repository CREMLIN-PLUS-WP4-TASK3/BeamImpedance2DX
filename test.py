#!/usr/bin/env python3

import bi2d
# from multipledispatch import dispatch
# import dolfinx
# from mpi4py import MPI

# from inspect import signature


m1 = bi2d.Material(1, eps=lambda x, y: 1+(x**2 - 2*y))
# print(m1)
m2 = bi2d.Material(2, eps=0.1, kappa=lambda x, y, f: f/1e5)

bi2d.convert_msh("test.msh", "test.xdmf")
m = bi2d.Mesh("test.xdmf")
# print(m.get_bounds(1))
# print(m.get_bounds(2))

mc = bi2d.MaterialMap(m, [m1, m2])
mc.save("material.xdmf")
mc.update()
mc.save("material2.xdmf")
mc.f=2e5
mc.update()
mc.save("material3.xdmf")
mc.update()
mc.save("material4.xdmf")


# Q = dolfinx.FunctionSpace(m.mesh, ("Discontinuous Lagrange", 0))
# f = dolfinx.Function(Q)
# with f.vector.localForm() as loc:
#     help(loc)

# 1

# import inspect

# def a(b,c=1,d=1):
#     print(inspect.getargvalues(inspect.currentframe()))


# inspect.getargspec(a)

# a(1,2,3)

# 1

# m.mesh.ufl_cell().num_vertices()
# m.subdomains.values.shape
# m.mesh
# help(m.subdomains)
# help(m.mesh.geometry)
# m.mesh.geometry.dofmap
# m.mesh.geometry.x.shape
# help(m.mesh.geometry.dofmap)

# import numpy as np

# np.max(m.mesh.geometry.dofmap.array)
# np.min(m.mesh.geometry.dofmap.array)
# m.mesh.geometry.dofmap.num_nodes
# m.mesh.geometry.dofmap.links(1)
# m.mesh.geometry.dofmap
# m.mesh.geometry.dofmap.links(3801)

# [m.mesh.geometry.dofmap.links(i) for i in m.subdomains.indices[m.subdomains.values == 1]]
# np.array([m.mesh.geometry.dofmap.links(i) for i in m.subdomains.indices[m.subdomains.values == 1]]).flatten()
# np.unique(np.array([m.mesh.geometry.dofmap.links(i) for i in m.subdomains.indices[m.subdomains.values == 1]]).flatten())
# (np.array([m.mesh.geometry.dofmap.links(i) for i in m.subdomains.indices[m.subdomains.values == 1]]).flatten()).shape
# ii = np.unique(np.array([m.mesh.geometry.dofmap.links(i) for i in m.subdomains.indices[m.subdomains.values == 1]]).flatten())

# m.mesh.geometry.x[ii]
# xx = m.mesh.geometry.x[ii]

# minx = np.min(xx[:,0])
# maxx = np.max(xx[:,0])
# sumx = np.sum(xx[:,0])/xx[:,0].size
# miny = np.min(xx[:,1])
# maxy = np.max(xx[:,1])
# sumy = np.sum(xx[:,1])/xx[:,1].size

# minx
# maxx
# sumx
# miny
# maxy
# sumy

# m.subdomains.indices[m.subdomains.values == 1]
# m.subdomains.values.shape

# m.mesh.geometry.x.shape

# dolfinx.cpp.mesh.midpoints(m.mesh, 0, np.array([0]))

# help(m.mesh.geometry.dofmap)

 # |  array
 # |  
 # |  num_nodes
 # |  
 # |  offsets



# m1 = bi2d.Material(1, eps=1)
# m2 = bi2d.Material(2, eps=1)

# b = bi2d.MaterialCreator(m, m1, m2)


# import numpy as np

# help(m.mesh.geometry)
# print(np.array(m.mesh.geometry.index_map().global_indices()).shape)
# help(m.mesh.topology)
# print(m.subdomains.values.shape)
# print(m.mesh.geometry.x.shape)
# print(m.mesh.geometry.x[m.subdomains.indices[m.subdomains.values == 1]])
# print(m.mesh.geometry.x)

# Q = dolfinx.VectorFunctionSpace(m.mesh, ("Discontinuous Lagrange", 1))
# Q = dolfinx.FunctionSpace(m.mesh, ("Discontinuous Lagrange", 0))
# f = dolfinx.Function(Q)

# f.interpolate(lambda x: (x[0]**2, x[1]*2))
# f.interpolate(lambda x: x[0] + 5*x[1])
# dolfinx.cpp.la.scatter_forward(f.x)
# with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "out.xdmf", "w") as xdmf:
#     m.xdmf_write_mesh(xdmf)
#     xdmf.write_function(f)


# def a(a):
#     return 2*a

# len(signature(a).parameters)

# def a(a, b):
#     return 2*a + b

# class A():
#     def __init__(self, f):
#         self.f = f

#     def do(self, *args):
#         return self.f(*args)

# b = A(a)

# print(b.do(1))
# print(b.do(1,2))

# print(a(1))
# print(a(1,2))
