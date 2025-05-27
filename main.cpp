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
  const int dim = 3;

  //Variables and principal matrix creation
  double E1 = 0, E2 = 0, G12 = 0, nu12 = 0, nu23 = 0;
  //TODO import variable values from parameters, hard coding for now
  E1 = 1e9; //MPa
  E2 = E1; //MPa
  nu12 = 0.25;
  nu23 = nu12;
  
  //create elasticity tensor in principal coordinates
  dealii::SymmetricTensor<4, dim> elasticityPrincipal;
  
  /*
  TODO create JSON files that map the results of each if statement condition in 4d space, comprised of 1 and 0s (pseudo identity matrix)
  Then just run through each material type, multiply by that mapping matrix and add to elasticity matrix
  */
  //find k before filling elasticity tensor
  const double constk = 1-2*(E2*(1+nu23)*pow(nu12,2))/E1-pow(nu23,2);
  //std::cout << constk << "\n";
  //defining G12 same as later on in the for loop, isotropic test case rn
  G12 = (E1*(1-pow(nu23,2))-E2*nu12*(1+nu23))/(2*constk);

  //SymmetricTensor object automatically applies symmetries in ijkl=jikl=ijlk, but still need to manually input ijkl=klij
  int m=0, n=0;
  for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
          {
            for (unsigned int k = 0; k < dim; ++k)
              {
                for (unsigned int l = 0; l < dim; ++l)
                  {
                    //Create temporary indices m, n to write equivalent Voigt notation for ijkl, makes it easier to process and visualize each case
                    //TODO make another function to convert Einstein to Voigt? (input a, b and output c to avoid duplicate loops for m, n?)
                    m=0;
                    n=0;

                    if (i==0 && j==0) {
                      m=1;
                    } else if (i==1 && j==1) {
                      m=2;
                    } else if (i==2 && j==2) {
                      m=3;
                    } else if ((i==1 && j==2)||((i==2 && j==1))) {
                      m=4;
                    } else if ((i==0 && j==2)||((i==2 && j==0))) {
                      m=5;
                    } else if ((i==0 && j==1)||((i==1 && j==0))) {
                      m=6;
                    }

                    if (k==0 && l==0) {
                      n=1;
                    } else if (k==1 && l==1) {
                      n=2;
                    } else if (k==2 && l==2) {
                      n=3;
                    } else if ((k==1 && l==2)||((k==2 && l==1))) {
                      n=4;
                    } else if ((k==0 && l==2)||((k==2 && l==0))) {
                      n=5;
                    } else if ((k==0 && l==1)||((k==1 && l==0))) {
                      n=6;
                    }

                    if(m==1 && n==1)
                    {
                      elasticityPrincipal[i][j][k][l] = E1*(1-pow(nu23,2))/constk;
                    }
                    else if((m==2 && n==2)||(m==3 && n==3))
                    {
                      elasticityPrincipal[i][j][k][l] = E2*(1-E2/E1*pow(nu12,2))/constk;
                    }
                    else if(m==1 && (n==2 || n==3))
                    {
                      elasticityPrincipal[i][j][k][l] = E2*nu12*(1+nu23)/constk;
                      //applying symmetry along main diagonal (not automatically done for SymmetricTensor)
                      elasticityPrincipal[k][l][i][j] = elasticityPrincipal[i][j][k][l];
                    }
                    else if(m==2 && n==3)
                    {
                      elasticityPrincipal[i][j][k][l] = E2*(E2/E1*pow(nu12,2)+nu23)/constk;
                      //applying symmetry along main diagonal (not automatically done for SymmetricTensor)
                      elasticityPrincipal[k][l][i][j] = elasticityPrincipal[i][j][k][l];
                    }
                    //
                    else if(m==4 && n==4)
                    {
                      //C11 and C12 are already known based on order of for loops (all of i=1 is done first), but writing explicitly just to be safe
                      elasticityPrincipal[i][j][k][l] = (E1*(1-pow(nu23,2))-E2*nu12*(1+nu23))/(2*constk);
                    }
                    //
                    else if((m==5 && n==5)||(m==6 && n==6))
                    {
                      elasticityPrincipal[i][j][k][l] = G12;
                    }
                  }
              }
          }
      }

  //creating array to get fiber coordinates before creating tensor
  dealii::Tensor<1, dim> fiber;
  fiber[0] = 1;
  fiber[1] = 2;
  //define z axis only for 3D case (otherwise out of bounds)
  if (dim == 3){
    fiber[2] = 1;
  }
  
  //Create projection of fiber onto xy plane to find angles
  //fiberxy will be the same as input fiber for 2d case, creating this way to avoid duplicate code for theta
  dealii::Tensor<1, dim> fiberxy;
  fiberxy[0] = fiber[0];
  fiberxy[1] = fiber[1];
  //Create vector of x axis direction
  dealii::Tensor<1, dim> xaxis;
  xaxis[0] = 1;

  //find theta from angle between x axis and xy projection
  double theta = dealii::Physics::VectorRelations::angle(fiberxy, xaxis);
  //TODO find how to prepopulate R, Rz and Ry from identity tensor using SymmetricTensor::unit_symmetric_tensor() function
  //Prepopulating not really necessary, defaults to 0 value, so identity only ensures the "1" along main diagonal is filled for 3d
  dealii::Tensor<2, dim> R;
  dealii::Tensor<2, dim> Rz;

  std::cout << std::scientific << std::setprecision(3);
  
  //temporary identity matrix function
  for (int i = 0; i < dim; i++){
    R[i][i] = 1;
    Rz[i][i] = 1;
  }

  //Rotate about z axis, same process for both 2d and 3d 
  //TODO find way to condense this part?
  Rz[0][0] = cos(theta);
  Rz[0][1] = -sin(theta);
  Rz[1][0] = sin(theta);
  Rz[1][1] = cos(theta);

  //Rotation tensor creation depends on dimension of simulations
  if(dim == 2){
    //2d case can be rotated about the z axis, since the whole simulation is within the xy plane
    R = Rz;
    //display rotation tensor for debugging
    std::cout << R[0][0] << " " << R[1][0] <<  "\n"
    << R[0][1] << " " << R[1][1] <<"\n\n";
    
  }else if(dim == 3){
    //3d case needs phi to account for components in the z direction, which is found using xy projection and full fiber direction
    //TODO finding phi and theta will cause issues for xy projection when fiber is purely in z direction
    double phi = dealii::Physics::VectorRelations::angle(fiber, fiberxy);
    std::cout << "theta = " << 57.3*theta << "\nphi = " << 57.3*phi << "\n";

    //same issues with Rz, identity tensor and condense assigning sin and cos values
    dealii::Tensor<2, dim> Ry;
    for (int i = 0; i < dim; i++){
      Ry[i][i] = 1;
    }
    Ry[0][0] = cos(phi);
    Ry[0][2] = -sin(phi);
    Ry[2][0] = sin(phi);
    Ry[2][2] = cos(phi);
    
    R=Rz*Ry;
    
    /*
    Trying out explicit rotation tensor solution
    R[0][0] = cos(theta)*cos(phi);
    R[0][1] = -sin(theta);
    R[0][2] = -cos(theta)*sin(phi);
    R[1][0] = sin(theta)*cos(phi);
    R[1][1] = cos(theta);
    R[1][2] = -sin(theta)*sin(phi);
    R[2][0] = sin(phi);
    R[2][1] = 0;
    R[2][2] = cos(phi);

    std::cout << R.norm() << "\n\n";
    std::cout << R[0][0] << " " << R[0][1] << " " << R[0][2] << "\n"
      << R[1][0] << " " << R[1][1] << " " << R[1][2] << "\n"
      << R[2][0] << " " << R[2][1] << " " << R[2][2] << "\n";
    */

    //Testing rotation of x axis, which should produce the fiber direction, used for testing
    /*dealii::Tensor<1, dim> ans = R*xaxis;
    std::cout << ans[0] << "\n" << ans[1] << "\n" << ans[2] << "\n";
    */
  }
  //Create temporary asymmetric tensor for multiplications then converting to symmetric after
  dealii::Tensor<4, dim> temp, temp2;
  SymmetricTensor<4, dim> elasticityCartesian;
  //TODO double check multiplication is correct
  //https://stackoverflow.com/questions/50178156/efficient-tensor-multiplication
  //temp = R*R*elasticityPrincipal*transpose(R)*transpose(R);

  temp = R*elasticityPrincipal;
  //shuffle j to front
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          temp2[j][i][k][l] = temp[i][j][k][l];
        }
      }
    }
  }
  temp2 = R*temp2;
  //move j back to prepare for next
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          temp[i][j][k][l] = temp2[j][i][k][l];
        }
      }
    }
  }

  //shuffle k to front
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          temp2[k][i][j][l] = temp[i][j][k][l];
        }
      }
    }
  }
  temp2 = R*temp2;
  //move k back to prepare for next
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          temp[i][j][k][l] = temp2[k][i][j][l];
        }
      }
    }
  }

  //shuffle l to front
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          temp2[l][i][j][k] = temp[i][j][k][l];
        }
      }
    }
  }
  temp2 = R*temp2;
  //move l back for final rotation
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          temp[i][j][k][l] = temp2[l][i][j][k];
        }
      }
    }
  }


  //for loops required to move generic tensor object to symmetric object for output
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          //for cases where ijkl=jikl=ijlk will be overwritten with the last entry, assumes it is already symmetric
          elasticityCartesian[i][j][k][l] = temp[i][j][k][l];
        }
      }
    }
  }
  
  /*
  //do in matlab instead? better suited for matrix calcs
  dealii::Tensor<2, 3> T, VoigtCP, VoigtCG;
  T[0][0] = pow(cos(theta),2);
  T[0][1] = pow(sin(theta),2);
  T[0][2] = 2*cos(theta)*sin(theta);
  T[1][0] = pow(sin(theta),2);
  T[1][1] = pow(cos(theta),2);
  T[1][2] = -2*cos(theta)*sin(theta);
  T[2][0] = -cos(theta)*sin(theta);
  T[2][1] = cos(theta)*sin(theta);
  T[2][2] = pow(cos(theta),2)-pow(sin(theta),2);
  
  VoigtCP[0][0] = elasticityPrincipal[0][0][0][0];
  VoigtCP[0][1] = elasticityPrincipal[0][0][1][1];
  VoigtCP[1][0] = elasticityPrincipal[1][1][0][0];
  VoigtCP[1][1] = elasticityPrincipal[1][1][1][1];
  VoigtCP[2][2] = elasticityPrincipal[0][1][0][1];
  
  VoigtCG = invert(T)*VoigtCP*T;
  //VoigtCG = T*VoigtCP*invert(T);
  std::cout << VoigtCP[0][0] << "    " << VoigtCP[0][1] << "    " << VoigtCP[0][2] << "    " << "\n"
  << VoigtCP[1][0] << "    " << VoigtCP[1][1] << "    " << VoigtCP[1][2] << "    " << "\n"
  << VoigtCP[2][0] << "    " << VoigtCP[2][1] << "    " << VoigtCP[2][2] << "    " << "\n\n";  

  std::cout << VoigtCG[0][0] << "    " << VoigtCG[0][1] << "    " << VoigtCG[0][2] << "    " << "\n"
  << VoigtCG[1][0] << "    " << VoigtCG[1][1] << "    " << VoigtCG[1][2] << "    " << "\n"
  << VoigtCG[2][0] << "    " << VoigtCG[2][1] << "    " << VoigtCG[2][2] << "    " << "\n\n";  */
  /*
  dealii::Tensor<2, dim> RT = transpose(R);
*/
/*
  std::cout << R[0][0] << "    " << R[0][1] << "    " << R[0][2] << "    " << "\n"
  << R[1][0] << "    " << R[1][1] << "    " << R[1][2] << "    " << "\n"
  << R[2][0] << "    " << R[2][1] << "    " << R[2][2] << "    " << "\n\n";

  std::cout << RT[0][0] << "    " << RT[0][1] << "    " << RT[0][2] << "    " << "\n"
  << RT[1][0] << "    " << RT[1][1] << "    " << RT[1][2] << "    " << "\n"
  << RT[2][0] << "    " << RT[2][1] << "    " << RT[2][2] << "    " << "\n\n";
*/

/*std::cout << elasticityPrincipal[0][0][0][0] << "    " << elasticityPrincipal[0][0][1][1] << "    " << elasticityPrincipal[0][0][0][1] << "    " << "\n"
  << elasticityPrincipal[1][1][0][0] << "    " << elasticityPrincipal[1][1][1][1] << "    " << elasticityPrincipal[1][1][0][1] << "    " << "\n"
  << elasticityPrincipal[0][1][0][0] << "    " << elasticityPrincipal[0][1][1][1] << "    " << elasticityPrincipal[0][1][0][1] << "    " << "\n\n";

std::cout << elasticityCartesian[0][0][0][0] << "    " << elasticityCartesian[0][0][1][1] << "    " << elasticityCartesian[0][0][0][1] << "    " << "\n"
  << elasticityCartesian[1][1][0][0] << "    " << elasticityCartesian[1][1][1][1] << "    " << elasticityCartesian[1][1][0][1] << "    " << "\n"
  << elasticityCartesian[0][1][0][0] << "    " << elasticityCartesian[0][1][1][1] << "    " << elasticityCartesian[0][1][0][1] << "    " << "\n\n";
*/


  //outputs 1-3 square of Voigt notation components for debugging
/*std::cout << elasticityPrincipal[0][0][0][0] << "    " << elasticityPrincipal[0][0][1][1] << "    " << elasticityPrincipal[0][0][2][2] << "    " << "\n"
  << elasticityPrincipal[1][1][0][0] << "    " << elasticityPrincipal[1][1][1][1] << "    " << elasticityPrincipal[1][1][2][2] << "    " << "\n"
  << elasticityPrincipal[2][2][0][0] << "    " << elasticityPrincipal[2][2][1][1] << "    " << elasticityPrincipal[2][2][2][2] << "    " << "\n\n";

  std::cout << elasticityPrincipal[1][2][1][2] << "    " << elasticityPrincipal[1][2][0][2] << "    " << elasticityPrincipal[1][2][0][1] << "    " << "\n"
  << elasticityPrincipal[0][2][1][2] << "    " << elasticityPrincipal[0][2][0][2] << "    " << elasticityPrincipal[0][2][0][1] << "    " << "\n"
  << elasticityPrincipal[0][1][1][2] << "    " << elasticityPrincipal[0][1][0][2] << "    " << elasticityPrincipal[0][1][0][1] << "    " << "\n\n";

  std::cout << elasticityCartesian[0][0][0][0] << "    " << elasticityCartesian[0][0][1][1] << "    " << elasticityCartesian[0][0][2][2] << "    " << "\n"
  << elasticityCartesian[1][1][0][0] << "    " << elasticityCartesian[1][1][1][1] << "    " << elasticityCartesian[1][1][2][2] << "    " << "\n"
  << elasticityCartesian[2][2][0][0] << "    " << elasticityCartesian[2][2][1][1] << "    " << elasticityCartesian[2][2][2][2] << "    " << "\n\n";
  */
  
  //Using isotropic tensor code for testing matrix creation  
  //Pulling same vals of E and nu from anisotropic "isotropic" case
    double E = E1;
    double nu = nu12;

    double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));
    dealii::SymmetricTensor<4, dim> elasticityIso;
    for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
          {
            for (unsigned int k = 0; k < dim; ++k)
              {
                for (unsigned int l = 0; l < dim; ++l)
                  {
                    elasticityIso[i][j][k][l] =
                      (i == k && j == l ? mu : 0.0) +
                      (i == l && j == k ? mu : 0.0) +
                      (i == j && k == l ? lambda : 0.0);
                  }
              }
          }
      }
  
  /*SymmetricTensor<4, dim> elasticityCartesianIso;
  temp = R*R*elasticityIso*transpose(R)*transpose(R);
  //for loops required to move generic tensor object to symmetric object for output
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          //for cases where ijkl=jikl=ijlk will be overwritten with the last entry, assumes it is already symmetric
          elasticityCartesianIso[i][j][k][l] = temp[i][j][k][l];
        }
      }
    }
  }*/

  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          std::cout << (int)elasticityIso[i][j][k][l] << "    " << (int)elasticityPrincipal[i][j][k][l] << "    " << (int)temp[i][j][k][l] << "\n";
        }
      }
    }
  }

  /*
  std::cout << elasticityIso[0][0][0][0] << "    " << elasticityIso[0][0][1][1] << "    " << elasticityIso[0][0][2][2] << "    " << "\n"
  << elasticityIso[1][1][0][0] << "    " << elasticityIso[1][1][1][1] << "    " << elasticityIso[1][1][2][2] << "    " << "\n"
  << elasticityIso[2][2][0][0] << "    " << elasticityIso[2][2][1][1] << "    " << elasticityIso[2][2][2][2] << "    " << "\n\n";

  std::cout << elasticityIso[1][2][1][2] << "    " << elasticityIso[1][2][0][2] << "    " << elasticityIso[1][2][0][1] << "    " << "\n"
  << elasticityIso[0][2][1][2] << "    " << elasticityIso[0][2][0][2] << "    " << elasticityIso[0][2][0][1] << "    " << "\n"
  << elasticityIso[0][1][1][2] << "    " << elasticityIso[0][1][0][2] << "    " << elasticityIso[0][1][0][1] << "    " << "\n\n";

  std::cout << elasticityCartesianIso[0][0][0][0] << "    " << elasticityCartesianIso[0][0][1][1] << "    " << elasticityCartesianIso[0][0][2][2] << "    " << "\n"
  << elasticityCartesianIso[1][1][0][0] << "    " << elasticityCartesianIso[1][1][1][1] << "    " << elasticityCartesianIso[1][1][2][2] << "    " << "\n"
  << elasticityCartesianIso[2][2][0][0] << "    " << elasticityCartesianIso[2][2][1][1] << "    " << elasticityCartesianIso[2][2][2][2] << "    " << "\n\n";
  */

}
