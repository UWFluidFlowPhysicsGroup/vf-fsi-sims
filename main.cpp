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
//needed for vector relations for anisotropy
#include <deal.II/physics/vector_relations.h>

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
  int loadMesh(std::string meshNameSolid, std::string meshNameFluid);
  int setParams(Parameters::AllParameters params);
  Triangulation<3> extrude();
  int refine(int refinement);
private:
  Triangulation<dim> triaSolid, triaFluid;
  DoFHandler<dim> dof_handler;
  GridIn<dim> gridIn;
};

namespace {
const std::string simMeshSolid = "leafletSolid";
//Ability to set multiple fluid meshes to simplify fluid mesh refinement studies
const std::string simMeshFluid[] = {"leafletFluid_1799"};
const std::string meshPath = "meshes/";
const std::string paramsPath = "fsi_leaflet.prm";
GridOut gridOut;
}

//only need to define dof_handler once for the sim dimensions, i guess this is how it reads what dim to use?
template <int dim>
Sim<dim>::Sim()
  : dof_handler(triaSolid) 
{}

//imports a mesh and outputs svg file in the XY plane
template <int dim>
int Sim<dim>::loadMesh(std::string meshNameSolid, std::string meshNameFluid){
  //identifies mesh to be imported from meshes folder
  std::ifstream solidPath(meshPath + meshNameSolid + ".msh");
  std::ifstream fluidPath(meshPath + meshNameFluid + ".msh");
  //checks if desired mesh can be read
  if (!solidPath || !fluidPath){
    //Display error handler that file cannot be found
    std::cerr << "----------------------------------------------------"
              << "ERROR FINDING MESH FILES " << meshNameSolid << " OR " << meshNameFluid
              << "----------------------------------------------------";
    //return to kill the class
    return 1;
  }

  //define GridIn object to receive 2d mesh
  gridIn.attach_triangulation(triaSolid);
  //imports mesh from selected area
  gridIn.read_msh(solidPath);
  
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
  Solid::LinearElasticity<dim> solid(triaSolid, params);
  Fluid::InsIM<dim> fluid(triaFluid, params);
  //combine solid and fluid meshes to make FSI simulation
  FSI<dim> fsi(fluid, solid, params, true);
  fsi.run();
  
  return 0;
}


int main(){
  //Variables and principal matrix creation
  double E1 = 0, E2 = 0, G12 = 0, mu12 = 0, mu32 = 0;
  //TODO import variable values from parameters, hard coding for now
  E1 = 10;
  E2 = 1;
  
  //TODO replace this dim var with dim obtained from template
  const int dim = 2;
  //creating array to get fiber coordinates before creating tensor (trying to simulate reading from parameter file)
  double fiberCoords[dim];
  fiberCoords[0] = 1;
  fiberCoords[1] = 1;
  //define z axis only for 3D case (otherwise out of bounds)
  if (dim == 3){
    fiberCoords[2] = 0;
  }

  //linear_elastic_material uses symmetric tensor named "elasticity", keeping same for convenience with integration
  SymmetricTensor<4, dim> elasticity;
  
  /*
  TODO create JSON files that map the results of each if statement condition in 4d space, comprised of 1 and 0s (pseudo identity matrix)
  Then just run through each material type, multiply by that mapping matrix and add to elasticity matrix
  */
  
  for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
          {
            for (unsigned int k = 0; k < dim; ++k)
              {
                for (unsigned int l = 0; l < dim; ++l)
                  {
                    //TODO create if else statement for each case for the equivalent Voigt notation for ij and kl
                    //S11
                    if (i==0 && j==0 && k==0 && l==0)
                    {
                      elasticity[i][j][k][l] = 1/E1;
                    }
                    //S22, S33
                    else if ((i==1 && j==1 && k==1 && l==1)||(i==2 && j==2 && k==2 && l==2))
                    {
                      elasticity[i][j][k][l] = 1/E2;
                    }
                    //S12=S13
                    else if (((i==0 && j==0) && ((k==1 && l==1)||(k==2 && l==2)))||(((i==1 && j==1)||(i==2 && j==2)) && (k==0 && l==0)))
                    {
                      elasticity[i][j][k][l] = -mu12/E1;
                    }
                    //S23
                    else if (((i==1 && j==1) && (k==2 && l==2)) || ((i==2 && j==2) && (k==1 && l==1)))
                    {
                      elasticity[i][j][k][l] = -mu32/E2;
                    }
                    //S44
                    else if (((i==1 && j==2) || (i==2 && j==1)) && ((k==1 && l==2) || (k==2 && l==1)))
                    {
                      elasticity[i][j][k][l] = E2/(2*(1-mu12));
                    }
                    //S55=S66
                    else if (((i==0 && j==2) || (i==2 && j==0)) && ((k==0 && l==2) || (k==2 && l==0))||((i==0 && j==1) || (i==1 && j==0)) && ((k==0 && l==1) || (k==1 && l==0)))
                    {
                      elasticity[i][j][k][l] = 1/G12;
                    }
                  }
              }
          }
      }
      
    //create rotation tensor R, size depends on dimension
    Tensor<2, dim> R;
    //create fiber direction vector and import assigned values from parameters list

    if (dim == 2){
      Tensor<1, 2> fiber;
      for(int i = 0; i < dim; i++){
        fiber[i] = fiberCoords[i];
      }

      //defining x axis to find rotation matrix
      Tensor<1, 2> xaxis({1, 0});

      //can obtain 2d rotation tensor directly using sin and cos wrt x axis
      R[0][0] = cos(Physics::VectorRelations::angle(fiber, xaxis));
      R[0][1] = sin(Physics::VectorRelations::angle(fiber, xaxis));
      R[1][0] = -sin(Physics::VectorRelations::angle(fiber, xaxis));
      R[1][1] = cos(Physics::VectorRelations::angle(fiber, xaxis));

      /* int angleTest = Physics::VectorRelations::angle(fiber, xaxis);
      std::cout << angleTest << "\n";
      
      //display rotation tensor for debugging
      std::cout << R[0][0] << " " << R[1][0] <<  "\n"
      << R[0][1] << " " << R[1][1] <<"\n"; */
    } else {
      //brute forced, cannot call subset of identity matrix easily (ex, identity(:,1) as in MATLAB)
      Tensor<1, 3> fiber;
      for(int i = 0; i < dim; i++){
        fiber[i] = fiberCoords[i];
      }
      
      Tensor<1, 3> xaxis({1, 0, 0});
      Tensor<1, 3> yaxis({0, 1, 0});
      Tensor<1, 3> zaxis({0, 0, 1});

      //declare 1st perpendicular fiber vector
      Tensor<1, 3> fiberP1;
      /*both inputs for angle need to be Tensor<1, dim>, cannot call a portion of 3x3 tensor so each vector has to be declared separately
      https://www.dealii.org/current/doxygen/deal.II/namespacePhysics_1_1VectorRelations.html#a9f05135611d90ad209c97bd73a4c4d20
      */
      
      //use z axis as reference axis unless fiber is parallel to z, then use y
      //For 2D case fiber direciton will always be perpendicular to z axis
      Tensor<1, 3> ref;
      if (cos(Physics::VectorRelations::angle(fiber, xaxis)) == 1){
        ref = yaxis;
      }else{
        ref = zaxis;   
      }

      fiberP1[0] = fiber[1]*ref[2]-fiber[2]*ref[1];
      fiberP1[1] = -(fiber[0]*ref[2]-fiber[2]*ref[0]);
      fiberP1[2] = fiber[0]*ref[1]-fiber[1]*ref[0];

      //take cross product of fiber and fiberP1 to find fiberP2
      //TODO not needed for dim = 2 since this will be z axis?
      Tensor<1, 3> fiberP2;
      fiberP2[0] = fiber[1]*fiberP1[2]-fiber[2]*fiberP1[1];
      fiberP2[1] = -(fiber[0]*fiberP1[2]-fiber[2]*fiberP1[0]);
      fiberP2[2] = fiber[0]*fiberP1[1]-fiber[1]*fiberP1[0];
      
      //create array for fiber and global coordinate systems
      Tensor<1, 3> fiberC[3] = {fiber, fiberP1, fiberP2};
      Tensor<1, 3> globalC[3] = {xaxis, yaxis, zaxis};
      for (unsigned int i = 0; i < dim; i++){
        for (unsigned int j = 0; j < dim; j++){
          R[i][j] = cos(Physics::VectorRelations::angle(fiberC[i], globalC[j]));
        }
      }

      /* int angleTest = Physics::VectorRelations::angle(fiber, xaxis);
      std::cout << angleTest << "\n";

      std::cout << R[0][0] << " " << R[1][0] << " " << R[2][0] << "\n"
      << R[0][1] << " " << R[1][1] << " " << R[2][1] << "\n"
      << R[0][2] << " " << R[1][2] << " " << R[2][2] << "\n"; */
    }
  Tensor<4, dim> elasticityCartesian;
  //Rotate elasticity tensor to cartesian global coordinates
  elasticityCartesian = R*R*elasticity*transpose(R)*transpose(R);

  //read parameters file to determine the dimensions present
  Parameters::AllParameters params(paramsPath);
  
  //iterate through each fluid mesh that was given
  for(const std::string &meshFluid : simMeshFluid){
    //This section has to be hard coded, since the creation of the Sim object requires a constant variable input
    //the value of ‘dims’ is not usable in a constant expression
    if (params.dimension == 2){
      Sim<2> sim;
      sim.loadMesh(simMeshSolid, meshFluid);
      //sim.setParams(params);

      /*Keeping extrude and refine functions commented out for future reference
      sim.extrude();
      for(int i = 1; i <= 2; i++){
        sim.refine(i);
      } */ 
    } else if (params.dimension == 3){
      Sim<3> sim;
      sim.loadMesh(simMeshSolid, meshFluid);
      //sim.setParams(params);
      
      /*
      for(int i = 1; i <= 2; i++){
        sim.refine(i);
      }
      */
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