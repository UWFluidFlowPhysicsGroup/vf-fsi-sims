import meshio

mesh = meshio.read("vocal_tract.inp")

meshio.write_points_cells("vocal_tract.msh", mesh.points, mesh.cells)

# meshio convert --input-format gmsh --output-format abaqus inputfile.msh outputfile.inp