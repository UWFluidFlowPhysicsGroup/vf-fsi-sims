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
  double E1, E2, G12, mu12, mu32;
  //dim not really needed, but would be good to add for 2D case
  unsigned int dim = 3;
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
                    if (i==1 && j==1 && k==1 && l==1)
                    {
                      elasticity[i][j][k][l] = 1/E1;
                    }
                    //S22, S33
                    else if ((i==2 && j==2 && k==2 && l==2)||(i==3 && j==3 && k==3 && l==3))
                    {
                      elasticity[i][j][k][l] = 1/E2;
                    }
                    //S12=S13
                    else if (((i==1 && j==1) && ((k==2 && l==2)||(k==3 && l==3)))||(((i==2 && j==2)||(i==3 && j==3)) && (k==1 && l==1)))
                    {
                      elasticity[i][j][k][l] = -mu12/E1;
                    }
                    //S23
                    else if (((i==2 && j==2) && (k==3 && l==3)) || ((i==3 && j==3) && (k==2 && l==2)))
                    {
                      elasticity[i][j][k][l] = -mu32/E2;
                    }
                    //S44
                    else if (((i==2 && j==3) || (i==3 && j==2)) && ((k==2 && l==3) || (k==3 && l==2)))
                    {
                      elasticity[i][j][k][l] = E2/(2*(1-mu12));
                    }
                    //S55=S66
                    else if (((i==1 && j==3) || (i==3 && j==1)) && ((k==1 && l==3) || (k==3 && l==1))||((i==1 && j==2) || (i==2 && j==1)) && ((k==1 && l==2) || (k==2 && l==1)))
                    {
                      elasticity[i][j][k][l] = 1/G12;
                    }
                  }
              }
          }
      }
  

      //Rotation matrix (not symmetric), contains coordinate system for fiber in global coordinates
      Tensor<2, dim> R;
      
      //Fiber coordinate system tensor, 1st column is the fiber direction, other two dimensions are undefined for now, will be found shortly
      Tensor<1, dim> fiber({1, 1, 0});
      /*
      fiber[0][0] = 1;
      fiber[1][0] = 1;
      fiber[2][0] = 0;
      */

      //brute forced, cannot call subset of identity matrix easily (ex, identity(:,1) as in MATLAB)
      Tensor<1, dim> xaxis({1, 0, 0});
      Tensor<1, dim> yaxis({0, 1, 0});
      Tensor<1, dim> zaxis({0, 0, 1});
      /*
      Tensor<1, 3> xaxis;
      xaxis[0] = 1;
      xaxis[1] = 0;
      xaxis[2] = 0;
      Tensor<1, 3> yaxis;
      yaxis[0] = 0;
      yaxis[1] = 1;
      yaxis[2] = 0;
      Tensor<1, 3> zaxis;
      zaxis[0] = 0;
      zaxis[1] = 0;
      zaxis[2] = 1;
      */

      //checks if fiber is parallel to z axis
      //Tensor<1, 3> ref;

      //declare 1st perpendicular fiber vector
      Tensor<1, dim> fiberP1;
      /*both inputs for angle need to be Tensor<1, dim>
      https://www.dealii.org/current/doxygen/deal.II/namespacePhysics_1_1VectorRelations.html#a9f05135611d90ad209c97bd73a4c4d20
      */
      //Need to fix calling of identity matrix, need to make a sub tensor
      //use y axis if fiber is along x axis direction, otherwise use x
      //TODO need to fix or add something for 2D case?
      if (Physics::VectorRelations::angle(fiber, xaxis) == 0){
        fiberP1[0] = fiber[1]*yaxis[2]-fiber[2]*yaxis[1];
        fiberP1[1] = -(fiber[0]*yaxis[2]-fiber[2]*yaxis[0]);
        fiberP1[2] = fiber[0]*yaxis[1]-fiber[1]*yaxis[0];
      }else{
        fiberP1[0] = fiber[1]*xaxis[2]-fiber[2]*xaxis[1];
        fiberP1[1] = -(fiber[0]*xaxis[2]-fiber[2]*xaxis[0]);
        fiberP1[2] = fiber[0]*xaxis[1]-fiber[1]*xaxis[0];
      }

      //take cross product of fiber and fiberP1 to find fiberP2
      Tensor<1, dim> fiberP2;
      fiberP2[0] = fiber[1]*fiberP1[2]-fiber[2]*fiberP1[1];
      fiberP2[1] = -(fiber[0]*fiberP1[2]-fiber[2]*fiberP1[0]);
      fiberP2[2] = fiber[0]*fiberP1[1]-fiber[1]*fiberP1[0];
      
      //cannot use for loop here, since calling on different tensors
      //counters for rotation matrix positions, need to find better way to do this
      //not working for now, need to look at documentation
      /*int i = 0, j = 0;
      for (Tensor<1, 3> u = {fiber, fiberP1, fiberP2};){
        //reset value of j for next loop
        j=0;
        for (Tensor<1, 3> v = {xaxis, yaxis, zaxis};){
          //Need to find how to identify back to rotation matrix location
          R[i][j] = cos(Physics::VectorRelations::angle(u, v));
          j++;
        }
        i++;
      }
      */


      R[0][0] = cos(Physics::VectorRelations::angle(fiber, xaxis));
      R[0][1] = cos(Physics::VectorRelations::angle(fiber, yaxis));
      R[0][2] = cos(Physics::VectorRelations::angle(fiber, zaxis));
      R[1][0] = cos(Physics::VectorRelations::angle(fiberP1, xaxis));
      R[1][1] = cos(Physics::VectorRelations::angle(fiberP1, yaxis));
      R[1][2] = cos(Physics::VectorRelations::angle(fiberP1, zaxis));
      R[2][0] = cos(Physics::VectorRelations::angle(fiberP2, xaxis));
      R[2][1] = cos(Physics::VectorRelations::angle(fiberP2, yaxis));
      R[2][2] = cos(Physics::VectorRelations::angle(fiberP2, zaxis));

      std::cout << R[0][0] << " " << R[1][0] << " " << R[2][0] << "\n"
      << R[0][1] << " " << R[1][1] << " " << R[2][1] << "\n"
      << R[0][2] << " " << R[1][2] << " " << R[2][2] << "\n";
      
      /*
      for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
          {
            R[i][j] = cos(Physics::VectorRelations::angle(fiber[i], identity[j]));
          }
      }
      */
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
      sim.setParams(params);

      /*Keeping extrude and refine functions commented out for future reference
      sim.extrude();
      for(int i = 1; i <= 2; i++){
        sim.refine(i);
      } */ 
    } else if (params.dimension == 3){
      Sim<3> sim;
      sim.loadMesh(simMeshSolid, meshFluid);
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