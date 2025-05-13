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
  //TODO remove dim when importing to OpenIFEM code
  const int dim = 2;

  //Variables and principal matrix creation
  double E1 = 0, E2 = 0, G12 = 0, mu12 = 0, mu32 = 0;
  //TODO import variable values from parameters, hard coding for now
  E1 = 10;
  E2 = 1;
  
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
  //creating array to get fiber coordinates before creating tensor (trying to simulate reading from parameter file)
  dealii::Tensor<1, dim> fiber;
  fiber[0] = 1;
  fiber[1] = 1;
  //define z axis only for 3D case (otherwise out of bounds)
  if (dim == 3){
    fiber[2] = 0;
  }

  //TODO: Create new rotation tensor using Euler rotations
  //create fiber projection in xy and x axis in 2D
  dealii::Tensor<1, dim> fiberxy({fiber[0], fiber[1]});
  dealii::Tensor<1, dim> xaxis({1, 0});
  //find theta angle from x axis
  double theta = dealii::Physics::VectorRelations::angle(fiberxy, xaxis);
  //TODO find how to prepopulate R, Rz and Ry from identity tensor using SymmetricTensor::unit_symmetric_tensor() function
  dealii::Tensor<2, dim> R;
  dealii::Tensor<2, dim> Rz;
  for (int i = 0; i < dim; i++){
    R[i][i] = 1;
    Rz[i][i] = 1;
  }

  //Rotate about z axis, same process for both 2d and 3d 
  //TODO find way to condense this part?
  Rz[0][0] = cos(theta);
  Rz[0][1] = sin(theta);
  Rz[1][0] = -sin(theta);
  Rz[1][1] = cos(theta);
  //defining R initially as identity matrix creates R for both 2d and 3d
  if(dim == 2){
    R = Rz;
  }else if(dim == 3){
    //if statement for y axis rotation for 3d
    //find phi angle between xy and full fiber vector
    //find rotation tensor directly or multiply together? (Could verify both are correct, would be easier to get program to do all the work)
    double phi = dealii::Physics::VectorRelations::angle(fiber, fiberxy);

    //same issues with Rz, identity tensor and condense assigning sin and cos values
    dealii::Tensor<2, dim> Ry;
    for (int i = 0; i < dim; i++){
      Ry[i][i] = 1;
    }
    Ry[0][0] = cos(phi);
    Ry[0][2] = sin(phi);
    Ry[2][0] = -sin(phi);
    Ry[2][2] = cos(phi);

    R=Rz*Ry;
  }
  
  //test edge cases (0 deg, 90, 180, -90 and 45s between)
}
