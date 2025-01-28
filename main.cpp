//import dealII libraries
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

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

//Vars in unnamed namespace to avoid reading from other files
template <int dim>
class Sim{
public:
  //extern template class Solid::LinearElasticity<dim>;
  Sim();
  int loadMesh(std::string meshNameFluid);
  int setParams(Parameters::AllParameters params);
  Triangulation<3> extrude();
  int refine(int refinement);
private:
  Triangulation<dim> triaSolid, triaFluid;
  DoFHandler<dim> dof_handler;
  GridIn<dim> gridIn;
};

namespace {
//const std::string simMeshSolid = "leafletSolid";
//Ability to set multiple fluid meshes to simplify fluid mesh refinement studies
const std::string simMeshFluid[] = {"ThroatOpen1mm","ThroatOpen3mm","ThroatOpen5mm"};
const std::string meshPath = "meshes/";
const std::string paramsPath = "parameters.prm";
GridOut gridOut;
}

//only need to define dof_handler once for the sim dimensions, i guess this is how it reads what dim to use?
template <int dim>
Sim<dim>::Sim()
  : dof_handler(triaSolid) 
{}

//imports a mesh and outputs svg file in the XY plane
template <int dim>
int Sim<dim>::loadMesh(std::string meshNameFluid){
  //identifies mesh to be imported from meshes folder
  //std::ifstream solidPath(meshPath + meshNameSolid + ".msh");
  std::ifstream fluidPath(meshPath + meshNameFluid + ".msh");
  //checks if desired mesh can be read
  if (!fluidPath){
    //Display error handler that file cannot be found
    std::cerr << "----------------------------------------------------"
              << "ERROR FINDING MESH FILES " << " OR " << meshNameFluid
              << "----------------------------------------------------";
    //return to kill the class
    return 1;
  }

  //define GridIn object to receive 2d mesh
  //gridIn.attach_triangulation(triaSolid);
  //imports mesh from selected area
  //gridIn.read_msh(solidPath);
  
  //repeat same for fluid mesh
  gridIn.attach_triangulation(triaFluid);
  gridIn.read_msh(fluidPath);

  //Exports meshes to .msh file for debugging
  /*
  std::ofstream out(meshNameSolid + ".msh");
  std::ofstream out(meshNameFluid + ".msh");
  gridOut.write_msh(triaSolid, out);
  gridOut.write_msh(triaFluid, out);
  */
  return 0;
}

template <int dim>
int Sim<dim>::setParams(Parameters::AllParameters params){
  //import params for both solid and fluid meshes separately
  //Solid::LinearElasticity<dim> solid(triaSolid, params);
  Fluid::InsIM<dim> fluid(triaFluid, params);
  //combine solid and fluid meshes to make FSI simulation
  fluid.run();
  //FSI<dim> fsi(fluid, solid, params, true);
  //fsi.run();
  
  return 0;
}


int main(){
  //read parameters file to determine the dimensions present
  Parameters::AllParameters params(paramsPath);
  
  //iterate through each fluid mesh that was given
  for(const std::string &meshFluid : simMeshFluid){
    //This section has to be hard coded, since the creation of the Sim object requires a constant variable input
    //the value of ‘dims’ is not usable in a constant expression
    if (params.dimension == 2){
      Sim<2> sim;
      sim.loadMesh(meshFluid);
      sim.setParams(params);

      /*Keeping extrude and refine functions commented out for future reference
      sim.extrude();
      for(int i = 1; i <= 2; i++){
        sim.refine(i);
      } */ 
    } else if (params.dimension == 3){
      Sim<3> sim;
      sim.loadMesh(meshFluid);
      sim.setParams(params);
      
      /*
      for(int i = 1; i <= 2; i++){
        sim.refine(i);
      }*/
    } else {
      std::cerr << "Cannot find dimension from parameters file" << std::endl
                << "Check if " << paramsPath << "exists";
      return 1;
    }

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

/*
Commented out extrude and refine functions because focusing on 2d shape first and refining in gmsh instead of c++ dealii
//takes input 2d mesh from before and extrudes to a 3d shape, exports shape to .geo file
template <int dim>
Triangulation<3> Sim<dim>::extrude(){
  Triangulation<3> tria3d;  
  //2d input, number of slices, height, output height, output triangulation
  GridGenerator::extrude_triangulation(tria2d, 7, 12.0, tria3d);
  std::ofstream out(meshPath + "vocalFold3d.msh");
  gridOut.write_msh(tria3d, out);
  return 0;
}

template <int dim>
int Sim<dim>::refine(int refinement){
  //refine_global is set to 1 subdivision because mesh is subdivided from previous loop 
  //one further step into refinement
  tria.refine_global(1);
  //output the refined mesh with a different name based on refinement level
  std::ofstream out(meshPath + "vocalFold3d" + std::to_string(refinement) + ".msh");
  gridOut.write_msh(tria, out);
  return 0;
}
*/