This program solves Maxwell's equations for the 2D device sections and calculates monopole and dipole beam impedance.
This is a port of <https://bitbucket.org/uniederm/beamimpedance2d/src/master/> program to [FEniCSx](https://fenicsproject.org/).

# Installation

It's recommended to use [Docker](https://www.docker.com/) DolfinX image.
If you're new to docker, run a jupyter server command listed below and follow the link. It will take you to the JupyterLab IDE interface.
Check out `examples` folder for calculation examples.

In order to generate calculation meshes and visualize the simulated electromagnetic fields use
[Gmsh](https://gmsh.info/) and [ParaView](https://www.paraview.org/).

## Jupyter server with real/complex numbers
Build docker image by issuing the command
```bash
docker build --tag bi2d -f bi2d_lab.dockerfile  .
```

Jupyter server may be started with the following line
```bash
docker run --init --rm -p 8888:8888 -v "$(pwd)":/root/shared -w /root/shared bi2d
```

Real and complex numbers can be selected in kernel selection box.

## Shell with real numbers
To run script named `script.py`
```bash
docker run --init --rm -v "$(pwd)":/root/shared -w /root/shared dolfinx/dolfinx python3 script.py
```

## Shell with complex numbers
To run script named `script.py`
```bash
docker run -v $(pwd):/root/shared -w "/root/shared" --rm --env LD_LIBRARY_PATH=/usr/local/dolfinx-complex/lib --env PATH=/usr/local/dolfinx-complex/bin:/usr/local/gmsh-4.6.0-Linux64-sdk/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin --env PKG_CONFIG_PATH=/usr/local/dolfinx-complex/lib/pkgconfig --env PETSC_ARCH=linux-gnu-complex-32 --env PYTHONPATH=/usr/local/dolfinx-complex/lib/python3.8/dist-packages dolfinx/dolfinx python3 script.py
```

# Parallel processing

Calculations may be done in parallel using MPI. Assuming proper Python environment configuration to run script in
parallel use command
```bash
mpirun --host localhost:$(nproc) -n 4 python3 script.py
```
In this example number of used CPUs is `4`.

Note that parallel read/write of mesh and fields is supported. However, if you wish to write, for example, numpy arrays
to file, you'll need to guard file access command from parallel writes like so:
```python
data = solution.get_z(np.logspace(5, 12, num=20), beta=0.999, source_function=bi2d.SourceFunction.DIPOLE)
if MPI.COMM_WORLD.rank == 0:
    np.savetxt("output.dat", data)
```

# Workflow

Basic workflow is this:
1. Create geometry using [Gmsh](https://gmsh.info/)
1. Assign [physical region](https://gmsh.info/doc/texinfo/gmsh.html#t1) numbers to beam and material regions and surface
   impedance boundaries.
1. Generate mesh
1. Convert mesh to XDMF format
1. In Python assign materials to material regions and boundaries
1. Calculate the problem

Now, we will talk about these steps in more detail.

## Creating geometry

Creating geometries in Gmsh takes a bit getting used to but it's not that complicated for the 2D structures.
Go through the [Gmsh tutorial](https://gmsh.info/doc/texinfo/gmsh.html#Tutorial) pages (at the very least through the
first one) to get aquatinted with Gmsh syntax.

## Assign physical region numbers

Assigning physical regions is done using commands
```gmsh
Physical Surface(1) = {s1, s2};
Physical Curve(2) = {l1, l2, l3, l4};
```
Numbers in round braces are the unique physical region identifiers that are used later by BeamImpedance2DX library.
Variables or numbers in curly braces are the internal Gmsh surface and line numbers.

## Mesh generation

Generation of the mesh file `mesh.msh` from Gmsh script `mesh.geo` is done by using command
```bash
gmsh mesh.geo -2 -o mesh.msh
```

## Converting mesh to XDMF format

Support script `./tools/convert_msh.py` may be used to convert `.msh` mesh to `.xdfm`.
To convert `mesh.msh` to `mesh.xdmf` issue the command
```bash
./convert_msh.py mesh.msh mesh.xdmf
```
If calculation problem requires boundary data, convert it using command
```bash
./convert_msh.py mesh.msh mesh_boundary.xdmf line
```

## BeamImpedance2DX syntax

### Create materials

There's a number of materials defined in the library
```python
from bi2d.materials import vacuum, beam, steel, copper, titanium
```

By default the physical region number `1` is a beam region. It's used by a source functions and for the integration.
It should be round. Assign it using the line
```python
beam.index = 1
```
Other imported materials should be assigned a physical region number as well
```python
vacuum.index = 2
copper.index = 3
copper_copy = copper.copy()
copper_copy.index = 4
```
Alternatively, materials may be defined from scratch
```python
gold = Material(5, name="Gold", sigma=4.52e7)
alumina = Material(6, name="Alumina", eps_r=8)
unicorn_horn = Material(6, name="Unicorn horn", sigma=1e8, eps_r=2, mu_r_re=0.1, mu_r_im=0.2)
```
Dispersive properties are defined by supplying material property arguments with function that takes and returns one argument
```python
def dispersive_material_sigma(f):
    if f > 1e6:
        return 1e7
    else:
        return 1e6

dispersive_material = Material(7, name="Dispersive", eps=lambda f: 1.1 * (f/1e8), sigma=dispersive_material_sigma)
```
One would probably want to use a table to interpolate material properties from it. Library provides such functionality
in a form of `ArrayInterpolate` class.
```python
my_table = ArrayInterpolate.from_file("my_table.csv", 0, 2, delimiter=",") # x=column 0, y=column 2
my_dispersive_material = Material(8, name="My dispersive_ material", sigma=my_table.interp)
```

### Loading mesh and creating material map

Creating mesh and material map is just a matter of providing a mesh file name and a list of materials with assigned
material numbers.
```python
mesh = bi2d.Mesh("mesh.xdmf")
mesh_with_boundary = bi2d.Mesh("mesh.xdmf", "mesh_boundary.xdmf")
material_map = bi2d.MaterialMap(mesh, [beam, vacuum, steel])
```

### Create a solution object

`Solution` object holds all the calculated fields. It's instantiated using the created material map
```python
solution = Solution(material_map, H1_order=2, Hcurl_order=2)
```

`H1_order` and `Hcurl_order` arguments control a degree of field interpolation. `1` means linear interpolation, `2` ⸺
parabolic interpolation etc.
Higher the polynomial degree, higher the solution precision, bigger the problem size, longer the solution.

After creating a solution object you might want to enable solution logging
```python
import logging

solution.logger.setLevel(logging.INFO)
```
This will enable information messages about current calculation point.

### Calculate impedance

There's a high level function for calculating a beam impedance in a frequency range
```python
z = solution.get_z([1e5, 1e7], beta=0.1, rotation=0,
                  source_function=bi2d.SourceFunction.MONOPOLE, sibc=[steel, copper, titanium])
```
The first argument is a list of frequency points for the calculation. It's convenient to use
[`numpy.logspace`](https://numpy.org/doc/stable/reference/generated/numpy.logspace.html?highlight=logspace#numpy.logspace)
function to generate this list.
`beta` is a relative beam speed.
`source_function` may be `bi2d.SourceFunction.MONOPOLE` for monopole and `bi2d.SourceFunction.DIPOLE` for
dipole calculations.
If using `bi2d.SourceFunction.DIPOLE` source function, argument `rotation` may be used to specify dipole rotation.

Surface impedance boundary condition is applied to the physical lines specified in `sibc` list.
Note that for these materials only conductivity $\sigma$ is used.

There's a shortcut for application of SIBC to all boundaries.
If you set SIBC material index to `-1`, it will be applied to all the boundaries. And in this case you don't need to
load boundary data to your model.

### Visualize the fields

Calculated fields may be saved to file
```python
solution.save("solution.xdmf")
```
Use `ParaView` to actually view the fields
```bash
paraview solution.xdmf
```
This file contains split irrotational and solenoidal field parts. You may create a sum of these parts using `sum_fields`
function
```python
solution.sum_fields()
solution.save("solution.xdmf)"
```
