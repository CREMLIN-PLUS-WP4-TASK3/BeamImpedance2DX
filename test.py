#!/usr/bin/env python3

import bi2d
import dolfinx
import ufl
from mpi4py import MPI


m1 = bi2d.Material(1, eps=lambda x, y: 1+(x**2 - 2*y))
m2 = bi2d.Material(2, eps=0.1, kappa=lambda x, y, f: f/1e5)

bi2d.convert_msh("test.msh", "test.xdmf")
m = bi2d.Mesh("test.xdmf")

mc = bi2d.MaterialMap(m, [m1, m2])
mc.save("material.xdmf")
s = bi2d.Solution(mc)
j = bi2d.Js(s)
j.solve()
j.save("current.xdmf")


H1= dolfinx.FunctionSpace(mc.mesh.mesh, ("Lagrange", 1))
u_, v_ = ufl.TrialFunction(H1), ufl.TestFunction(H1)
dx = ufl.Measure('dx', domain=mc.mesh.mesh, subdomain_id=0)
a_p = ufl.inner(u_, v_) * dx
L_p = v_ * dx
projection = dolfinx.fem.LinearProblem(a_p, L_p)
f = projection.solve()
self.f.name = "current"


dx = ufl.dx(domain=mesh.m)


H1 = dolfinx.FunctionSpace(m.mesh, ("Lagrange",1))    #Space for the Potential

import numpy
boundary_facets = numpy.where(numpy.array(dolfinx.cpp.mesh.compute_boundary_facets(m.mesh.topology)))[0]
boundary_dofs = dolfinx.fem.locate_dofs_topological(V, fdim, boundary_facets)

u_re, v_re = ufl.TrialFunction(H1), ufl.TestFunction(H1)
u_im, v_im = ufl.TrialFunction(H1), ufl.TestFunction(H1)

# Phi=Function(Mix)
# Phir=Function(H1)
# Phii=Function(H1)

from ufl import inner, grad, dx


a_re = inner(grad(v_re),mc.eps*grad(u_re))*dx +inner(grad(v_re),grad(u_im))*dx

ufl.grad(v_re)

1

# #The stiffness matrices
# ar=inner(grad(vr),epsilon*grad(ur))*dx +inner(grad(vr),(kappa/omega)*grad(ui))*dx
# ai=inner(grad(vi),epsilon*grad(ui))*dx -inner(grad(vi),(kappa/omega)*grad(ur))*dx

# #The mass matrices
# br=(omega/(beta*c0))**2 * (inner(vr,epsilon*ur)*dx +inner(vr,(kappa/omega)*ui)*dx )
# bi=(omega/(beta*c0))**2 * (inner(vi,epsilon*ui)*dx -inner(vi,(kappa/omega)*ur)*dx )

# #The right hand side
# RHSr=1/(beta*c0) *vr*Jszr*dx
# RHSi=1/(beta*c0) *vi*Jszi*dx

# eq= ar+ai +br+bi==RHSr+RHSi


# Zero = Expression(('0.0','0.0'))
# def u0_boundary(x, on_boundary):    # returns boolean if x on boundary
#     return on_boundary

# BC=DirichletBC(Mix, Zero, u0_boundary)


# set_log_level(PROGRESS)
# solve(eq, Phi,BC,solver_parameters={"linear_solver": "mumps","preconditioner": "none"})
# #solve(eq, Phi,BC,solver_parameters={"linear_solver": "gmres","preconditioner": "sor"})
# #solve(eq, Phi,BC,solver_parameters={"linear_solver": "lu","preconditioner": "none"})
# (Phir,Phii)=Phi.split(deepcopy=False)

# if(plot3Dflag):
#     plot(Phir,title='Phir')
#     plot(Phii,title='Phii')
#     interactive()

# return [Phir,Phii]

# H1= dolfinx.FunctionSpace(mc.mesh.mesh, ("Lagrange", 1))
# f = dolfinx.Function(H1)
# u_, v_ = ufl.TrialFunction(H1), ufl.TestFunction(H1)
# c = dolfinx.Constant(m.mesh, 2)
# a_p = ufl.inner(u_, v_) * ufl.dx
# L_p = mc.beam * v_ * ufl.dx
# projection = dolfinx.fem.LinearProblem(a_p, L_p)
# ff = projection.solve()
# with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "field.xdmf", "w") as xdmf:
#     xdmf.write_mesh(mc.mesh.mesh)
#     xdmf.write_function(ff)
