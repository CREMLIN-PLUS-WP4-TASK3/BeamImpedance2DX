#!/usr/bin/env python3

"""Mesh conversion tools"""

from sys import argv
import meshio


def convert_msh(input_file: str, output_file: str,
                cell_type="triangle", cell_data_tag="gmsh:physical"):
    """Convert mesh to XDMF format.

    Based on https://jorgensd.github.io/dolfinx-tutorial/chapter3/subdomains.html
    """

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


if __name__ == "__main__":
    if len(argv) < 3 or len(argv) > 4:
        exit("Usage:\t./plot mesh.msh mesh.xdmf [[line|triangle]]")

    region = "triangle"
    if len(argv) == 4:
        if argv[3].lower() == "line":
            region = "line"

    convert_msh(argv[1], argv[2],
                cell_type=region)
