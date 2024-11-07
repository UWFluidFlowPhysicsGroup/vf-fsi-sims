//import dealII libraries
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

//import OpenIFEM libraries
//solid linear elastic solver
#include "linear_elasticity.h"
//fluid incompressible navier stokes solver
#include "insim.h"
//fluid-solid interface solver
#include "fsi.h"

#include "parameters.h"
#include "utilities.h"


//import c++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>
#include <filesystem>

//create solid objects
extern template class Solid::LinearElasticity<2>;
extern template class Solid::LinearElasticity<3>;
//create fluid objects
extern template class Fluid::InsIM<2>;
extern template class Fluid::InsIM<3>;
//fluid-solid interface objects
extern template class FSI<2>;
extern template class FSI<3>;

using namespace dealii;

namespace {
const std::string simMeshSolid = "FSIChannelSolid";
//Ability to set multiple fluid meshes to simplify fluid mesh refinement studies
const std::string simMeshFluid[] = {"FSIChannelFluid3"};
const std::string meshPath = "meshes/";
//TODO simplify parameters strings existing - leave to only 2d form for now?
const std::string paramsPath2d = "parameters2d.prm";
const std::string paramsPath3d = "parameters3d.prm";

//define objects for both 2d and 3d mesh manipulation
Triangulation<2> tria2dFluid;
Triangulation<3,3> tria3dFluid;

Triangulation<2> tria2dSolid;
Triangulation<3,3> tria3dSolid;

GridIn<2> gridIn2d;
GridIn<3> gridIn3d;

GridOut gridOut;
}

//imports a mesh and outputs svg file in the XY plane
int loadMesh2d(std::string meshNameSolid, std::string meshNameFluid){
  //identifies mesh to be imported from meshes folder
  std::ifstream solidPath(meshPath + meshNameSolid + ".msh");
  std::ifstream fluidPath(meshPath + meshNameFluid + ".msh");
  //checks if desired mesh can be read
  if (!solidPath){
    //Display error handler that the solid mesh cannot be found
    std::cerr << "----------------------------------------------------"
              << "ERROR FINDING SOLID MESH FILE " << meshNameSolid 
              << "----------------------------------------------------";
    //return to kill the class
    return -1;
  } else if (!fluidPath){
    //Display error handler that file cannot be found
    std::cerr << "----------------------------------------------------"
              << "ERROR FINDING FLUID MESH FILE " << meshNameFluid
              << "----------------------------------------------------";
    //return to kill the class
    return -1;
  }



  //TODO write if statement to check if 2d or 3d mesh being imported, maybe try catch as 3d and have 2d in the catch segment and merge with extrude class?
  //define 2D GridIn object to receive 2d mesh
  gridIn2d.attach_triangulation(tria2dSolid);
  //imports mesh from selected area
  gridIn2d.read_msh(solidPath);
  
  //repeat same for fluid mesh
  gridIn2d.attach_triangulation(tria2dFluid);
  gridIn2d.read_msh(fluidPath);
  
  //std::cout << std::filesystem::current_path();

  //prepares squareMesh.svg file
  //std::ofstream out(meshName + ".svg");
  //writes refined mesh to svg in XY plane
  //gridOut.write_svg(tria2d, out);
  return 0;
}

int importParams2d(std::string paramName){
    //import params from .prm file and assign to 2d square
    Parameters::AllParameters params(paramName);
    //import params for both solid and fluid meshes separately
    Solid::LinearElasticity<2> solid(tria2dSolid, params);
    Fluid::InsIM<2> fluid(tria2dFluid, params);
    //combine solid and fluid meshes to make FSI simulation
    FSI<2> fsi(fluid, solid, params, true);
    fsi.run();
    
    return 0;
}

/*
Commenting out because focusing on 2d mesh for now
//imports a 3D mesh and outputs svg file in the XY plane
int loadMesh3d(std::string meshName){
  //identifies mesh to be imported from meshes folder
  std::ifstream f(meshPath + meshName + ".msh");
  //checks if desired mesh can be read
  if (!f){
    //Display error handler that file cannot be found
    std::cerr << "----------------------------------------------------"
              << "ERROR FINDING MESH FILE " << meshName
              << "----------------------------------------------------";
    //return to kill the class
    return -1;
  }

  //TODO write if statement to check if 2d or 3d mesh being imported, maybe try catch as 3d and have 2d in the catch segment and merge with extrude class?
  //define 2D GridIn object to receive 2d mesh
  gridIn3d.attach_triangulation(tria3d);
  //imports mesh from selected area
  gridIn3d.read_msh(f);
  
  return 1;
}

int importParams3d(std::string paramName){
    //import params from .prm file
    Parameters::AllParameters params(paramName);
    Solid::LinearElasticity<3> solid(tria3dSolid, params);
    Fluid::InsIM<3> fluid(tria3dFluid, params);
    FSI<3> fsi(fluid, solid, params, true);
    fsi.run();
    
    return 0;
}
*/

/*
Commented out extrude and refine functions because focusing on 2d shape first and refining in gmsh instead of c++ dealii
//takes input 2d mesh from before and extrudes to a 3d shape, exports shape to .geo file
int extrude(){  
  //2d input, number of slices, height, output height, output triangulation
  GridGenerator::extrude_triangulation(tria2d, 7, 12.0, tria3d);
  std::ofstream out(meshPath + "vocalFold3d.msh");
  gridOut.write_msh(tria3d, out);
  return 0;
}

int refine(int i){
  //refine_global is set to 1 subdivision because mesh is subdivided from previous loop 
  //one further step into refinement
  tria3d.refine_global(1);
  //output the refined mesh with a different name based on refinement levevl
  std::ofstream out(meshPath + "vocalFold3d" + std::to_string(i) + ".msh");
  gridOut.write_msh(tria3d, out);
  return 0;
}
*/


int main(){
  //iterate through each fluid mesh that was given
  for(const string &meshFluid : simMeshFluid){
    //load meshes through loadMesh class
    loadMesh2d(simMeshSolid, meshFluid);
    //import parameters through importParams class
    importParams2d(paramsPath2d);
    
    //define path to current file location
    std::filesystem::path p = std::filesystem::current_path();
    //create folder with a title corresponding to the current fluid mesh name
    std::filesystem::create_directory(p / meshFluid);

    //iterate through each file in the main directory
    for(const auto& dirEntry : std::filesystem::directory_iterator(p)){
      //checks if each file is a .vtu or .pvd file
      //since these are main outputs for each test case, want to move them somewhere safe before starting another simulation
      if (dirEntry.path().extension() == ".vtu" || dirEntry.path().extension() == ".pvd"){
        //moves the "selected" outputs to the new folder corresponding to the fluid mesh name
        std::filesystem::rename(p / dirEntry.path().filename(), p / meshFluid / dirEntry.path().filename());
      }
    }
  }
}
