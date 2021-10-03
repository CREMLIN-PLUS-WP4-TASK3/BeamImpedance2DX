#!/usr/bin/env python3

import bi2d
import numpy as np

m1 = bi2d.Material(1)
m2 = bi2d.Material(2)

# bi2d.convert_msh("test.msh", "test.xdmf")
# bi2d.convert_msh("test.msh", "test_boundary.xdmf", cell_type="line")
m = bi2d.Mesh("test.xdmf")

mc = bi2d.MaterialMap(m, [m1, m2])
mc.save("material.xdmf")
solution = bi2d.Solution(mc)
j = bi2d.Js(solution, source_function=bi2d.SourceFunction.MONOPOLE_GAUSSIAN)
j.solve()

ediv = bi2d.Ediv(solution)
ecurl = bi2d.Ecurl(solution)


# solution.f = 1e5
# ediv.solve()
# solution.f = 1e6
# ediv.solve()
# solution.f = 1e7
# ediv.solve()
# solution.f = 1e8
# ediv.solve()
# solution.f = 1e9
# ediv.solve()
# solution.f = 1e10
# ediv.solve()
# solution.f = 1e11
# ediv.solve()


solution.beta = 0.5
# solution.f = 1e11

# ediv.solve()
# ecurl.solve()
# solution.save("solution.xdmf")

# print(solution.get_z())

# solution.f = 2e5

# ediv.solve()
# ecurl.solve()
# solution.save("solution.xdmf")

freq = np.logspace(5, 12, num=20)
z_re = np.zeros(freq.size)
z_im = np.zeros(freq.size)

for i, f in enumerate(freq):
    print(f"Point f={f:.0e} ({i+1}/{freq.size})")
    solution.f = f
    # solution.save("solution.xdmf")
    ediv.solve(petsc_options={"linear_solver": "preonly",
                              "preconditioner": "la"})
    ecurl.solve()
    z_re[i], z_im[i] = solution.get_z()

np.savetxt("out.csv", np.column_stack([freq, z_re, z_im]))
