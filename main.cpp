//import dealII libraries
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
//needed for vector relations for anisotropy
#include <deal.II/physics/vector_relations.h>

//import OpenIFEM libraries
//solid linear elastic solver
#include "linear_elasticity.h"

//import c++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>
#include <filesystem>

using namespace dealii;

//TODO remove dim when importing to OpenIFEM code
//dim and fiber are global variables for writeTable
const int dim = 3;
dealii::Tensor<1, dim> fiber;

int main(){
  //Variables and principal matrix creation
  double E1 = 0, E2 = 0, G12 = 0, nu12 = 0, nu23 = 0;
  //TODO import variable values from parameters, hard coding for now
  E1 = 131e9; //GPa
  E2 = 10.3e9; //GPa
  nu12 = 0.22;
  nu23 = 0.30;
  G12 = 6.9e9;
  
  //creating array to get fiber coordinates before creating tensor
  fiber[0] = 1;
  fiber[1] = 0;
  //define z axis only for 3D case (otherwise out of bounds)
  if (dim == 3){
    fiber[2] = 2;
  }

  //create elasticity tensor in principal coordinates
  dealii::SymmetricTensor<4, dim> elasticityPrincipal;
  
  /*
  TODO create JSON files that map the results of each if statement condition in 4d space, comprised of 1 and 0s (pseudo identity matrix)
  Then just run through each material type, multiply by that mapping matrix and add to elasticity matrix
  */
  //find k before filling elasticity tensor, used several times
  const double constk = 1-2*(E2*(1+nu23)*pow(nu12,2))/E1-pow(nu23,2);

  //defining G12 same as later on in the for loop, isotropic test case rn
  //G12 = (E1*(1-pow(nu23,2))-E2*nu12*(1+nu23))/(2*constk);
  
  //SymmetricTensor object automatically applies symmetries in ijkl=jikl=ijlk, but still need to manually input ijkl=klij symmetry
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
                    else if(m==4 && n==4)
                    {
                      //C11 and C12 are already known based on order of for loops (all of i=1 is done first), but writing explicitly just to be safe
                      elasticityPrincipal[i][j][k][l] = (E1*(1-pow(nu23,2))-E2*nu12*(1+nu23))/(2*constk);
                    }
                    else if((m==5 && n==5)||(m==6 && n==6))
                    {
                      elasticityPrincipal[i][j][k][l] = G12;
                    }
                  }
              }
          }
      }
  
  //Create projection of fiber onto xy plane to find angles
  //fiberxy will be the same as input fiber for 2d case, creating this way to avoid duplicate code for theta
  dealii::Tensor<1, dim> fiberxy;
  fiberxy[0] = fiber[0];
  fiberxy[1] = fiber[1];
  //Create vector of x axis direction
  dealii::Tensor<1, dim> xaxis ;
  xaxis[0] = 1;
  
  double theta = fiberxy.norm() == 0 ? 0 : dealii::Physics::VectorRelations::angle(fiberxy, xaxis);
  // if(fiberxy.norm() != 0){
  //   //find theta from angle between x axis and xy projection only if xy projection is non-zero
  //   theta = dealii::Physics::VectorRelations::angle(fiberxy, xaxis);
  // }
  
  dealii::Tensor<2, dim> R = unit_symmetric_tensor<dim>(), Rz = unit_symmetric_tensor<dim>();
  //dealii::Tensor<2, dim> Rz = unit_symmetric_tensor<dim>();

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
    /*std::cout << R[0][0] << " " << R[1][0] <<  "\n"
    << R[0][1] << " " << R[1][1] <<"\n\n";*/
    
  }else if(dim == 3){
    //3d case needs phi to account for components in the z direction, which is found using xy projection and full fiber direction
    const double pi = 3.14159265358979323846;
    double phi = fiberxy.norm() == 0 ? pi/2 : dealii::Physics::VectorRelations::angle(fiber, fiberxy);
    // double phi = dealii::Physics::VectorRelations::angle(fiber, fiberxy);
    // if(fiberxy.norm() == 0){
    //   //if fiberxy is 0, then vector is purely vertical, which is assumed to be 90 degrees
    //   phi = pi/2;
    // }else{
    //   //3d case needs phi to account for components in the z direction, which is found using xy projection and full fiber direction
    //   phi = dealii::Physics::VectorRelations::angle(fiber, fiberxy);
    // }

    std::cout << "theta = " << 57.3*theta << "\nphi = " << 57.3*phi << "\n";

    dealii::Tensor<2, dim> Ry = unit_symmetric_tensor<dim>();

    Ry[0][0] = cos(phi);
    Ry[0][2] = -sin(phi);
    Ry[2][0] = sin(phi);
    Ry[2][2] = cos(phi);
    
    R=Rz*Ry;

    /*std::cout << R[0][0] << " " << R[0][1] << " " << R[0][2] << "\n"
      << R[1][0] << " " << R[1][1] << " " << R[1][2] << "\n"
      << R[2][0] << " " << R[2][1] << " " << R[2][2] << "\n";
    */
  }
  //Create temporary asymmetric tensor for multiplications then converting to symmetric after
  dealii::Tensor<4, dim> temp, temp2;
  dealii::SymmetricTensor<4, dim> elasticityCartesian1, elasticityCartesian2;

  /*//Working rotation code
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
  }*/

//Attempted more condensed rotation code
  //Rotate about i and l
  temp = R*elasticityPrincipal*transpose(R);
  //Flip so j and k are on outside
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          temp2[j][i][l][k] = temp[i][j][k][l];
        }
      }
    }
  }
  //Rotate about j and k
  temp2 = R*temp2*transpose(R);
  //Rotate back to create original temp tensor
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          temp[i][j][k][l] = temp2[j][i][l][k];
          elasticityCartesian1[i][j][k][l] = temp2[j][i][l][k];
        }
      }
    }
  }

  //if elasticityCartesian1 == elasticityCartesian2, then this loop can be removed, should be the same, but double check just in case
  //for loops required to move generic tensor object to symmetric object for output
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          //for cases where ijkl=jikl=ijlk will be overwritten with the last entry, assumes it is already symmetric
          elasticityCartesian2[i][j][k][l] = temp[i][j][k][l];
        }
      }
    }
  }
  
  //Using isotropic tensor code for testing matrix creation  
  //Pulling same vals of E and nu from anisotropic "isotropic" case
  // double E = E1;
  // double nu = nu12;
  // double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
  // double mu = E / (2 * (1 + nu));
  // dealii::SymmetricTensor<4, dim> elasticityIso;
  // for (unsigned int i = 0; i < dim; ++i){
  //   for (unsigned int j = 0; j < dim; ++j){
  //     for (unsigned int k = 0; k < dim; ++k){
  //       for (unsigned int l = 0; l < dim; ++l){
  //         elasticityIso[i][j][k][l] =
  //           (i == k && j == l ? mu : 0.0) +
  //           (i == l && j == k ? mu : 0.0) +
  //           (i == j && k == l ? lambda : 0.0);
  //       }
  //     }
  //   }
  // }

  std::cout << std::scientific << std::setprecision(3);

  //Create output file
  std::string filename;
  if (dim == 2){
    filename = std::to_string((int)fiber[0]) + std::to_string((int)fiber[1]) + "_cpp.csv";
  } else {
    filename = std::to_string((int)fiber[0]) + std::to_string((int)fiber[1]) + std::to_string((int)fiber[2]) + "_cpp.csv";
  }
  std::ofstream file(filename);
  file << "Principal, Asymmetric Rotated, Symmetric Rotated 1, Symmetric Rotated 2\n";

  //for loop to iterate through all values in elasticity tensors
  for (unsigned int i = 0; i < dim; i++){
    for (unsigned int j = 0; j < dim; j++){
      for (unsigned int k = 0; k < dim; k++){
        for (unsigned int l = 0; l < dim; l++){
          //std::cout << (long)elasticityIso[i][j][k][l] << "    " << (long)elasticityPrincipal[i][j][k][l] << "    " << (long)temp[i][j][k][l] << "\n";
          file << elasticityPrincipal[i][j][k][l] << "," << temp[i][j][k][l] << "," << elasticityCartesian1[i][j][k][l] << "," << elasticityCartesian2[i][j][k][l] << "\n";
        }
      }
    }
  }
  std::cout << "Written to file " << filename << " successfully \n";
    
}