import meshio

mesh = meshio.read("vocal_tract.inp")

meshio.write_points_cells("vocal_tract.msh", mesh.points, mesh.cells)