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

#include <chrono>
#include <cstdint>
#include <sstream>
#include <random>

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
    void loadSolid(std::string solidMeshName);
    void loadFluid(std::string fluidMeshName);
    void setParams(Parameters::AllParameters params);
    //Triangulation<3> extrude(); 
    //int refine(int refinement);
    //is this expecting to call a function?
    // Solid::LinearElasticity<dim> solid(Triangulation<dim> triaSolid, Parameters::AllParameters params);
    // Fluid::InsIM<dim> fluid(Triangulation<dim> triaFluid, Parameters::AllParameters params);

  private:
    Triangulation<dim> triaSolid, triaFluid;
    DoFHandler<dim> dof_handler;
    GridIn<dim> gridIn;

};

namespace {
const std::string simMeshSolid[] = {"SquareMesh_4", "SquareMesh_16", "SquareMesh_64", "SquareMesh_256", "SquareMesh_1024", "SquareMesh_4096", "SquareMesh_16384"};
// const std::string simMeshSolid[] = {"SquareMesh_1024"};
//Ability to set multiple fluid meshes to simplify fluid mesh refinement studies
const std::string simMeshFluid[] = {""};
const std::string meshPath = "meshes/";
const std::string paramsPath = "parameters.prm";
GridOut gridOut;
}

//only need to define dof_handler once for the sim dimensions, i guess this is how it reads what dim to use?
template <int dim>
Sim<dim>::Sim()
  : dof_handler(triaSolid) 
{}

//Solid object creation
template <int dim>
void Sim<dim>::loadSolid(std::string solidMeshName){
  // Dynamically define path for solid mesh location
  std::ifstream solidPath(meshPath + solidMeshName + ".msh");
  // Check if given mesh path is valid
  if (!solidPath){
    std::cerr << "----------------------------------------------------" << "\n"
              << "ERROR FINDING SOLID MESH FILE " << solidMeshName  << "\n"
              << "----------------------------------------------------" << "\n";
    //exit the program
    exit(0);
  }

  //define GridIn object to receive 2d mesh
  gridIn.attach_triangulation(triaSolid);
  //imports mesh from selected area
  gridIn.read_msh(solidPath);
  
  return;
}

// Fluid object creation
template <int dim>
void Sim<dim>::loadFluid(std::string fluidMeshName){
  // Dynamically define path for fluid mesh location
  std::ifstream fluidPath(meshPath + fluidMeshName + ".msh");
  // Check if given mesh path is valid
  if (!fluidPath){
    std::cerr << "----------------------------------------------------" << "\n"
              << "ERROR FINDING FLUID MESH FILE " << fluidMeshName << "\n"
              << "----------------------------------------------------" << "\n";
    exit(0);
  }

  //define GridIn object to receive fluid mesh
  gridIn.attach_triangulation(triaFluid);
  //imports the fluid mesh from the valid file path
  gridIn.read_msh(fluidPath);

  return;
  //Example of exporting fluid mesh for debugging/checking what is used
  /*
  std::ofstream out(meshNameFluid + ".msh");
  gridOut.write_msh(triaFluid, out);
  */
}

template <int dim>
void Sim<dim>::setParams(Parameters::AllParameters params){
  //import params for both solid and fluid meshes separately
  //Solid::LinearElasticity<dim> solid(triaSolid, params);
  //Fluid::InsIM<dim> fluid(triaFluid, params);
  
  if (params.simulation_type == "Solid"){
    Solid::LinearElasticity<dim> solid(triaSolid, params);
    solid.run();
  }else if(params.simulation_type == "Fluid"){
    Fluid::InsIM<dim> fluid(triaFluid, params);
    fluid.run();
  }else if(params.simulation_type == "FSI"){
    //combine solid and fluid meshes to make FSI simulation
    Solid::LinearElasticity<dim> solid(triaSolid, params);
    Fluid::InsIM<dim> fluid(triaFluid, params);
    FSI<dim> fsi(fluid, solid, params, true);
    fsi.run();
  }else{
    //error occured
    exit(0);
  }
  return;
}

int main(){

  /*
  // std::random_device rd;
  //can replace rd() to a constant instead for seed
  std::mt19937 gen(1);
  //can use other random number generation methods https://en.cppreference.com/w/cpp/numeric/random.html
  //ex. normal distribution for deviatioin from a main fiber direction
  std::uniform_real_distribution<> dist(-M_PI, M_PI);


  //new code generates random angles in n x n format 
  int NCellX = 32;
  
  // std::cout << "Max random value = " << RAND_MAX << "\n";
  // long int seed = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count();
  // std::cout << seed << "\n";
  // srand(seed);
  std::ofstream outfile;  
  outfile.open("Rand dataset_" + std::to_string(NCellX*NCellX) + ".csv");

  // std::string filename = "Angle dataset " + to_string(dataset) + ".csv";
  // outfile.open(filename);
  // outfile.open("Angle dataset", to_string(dataset), ".csv");
  // int NCells = 16384;
  // for (int i = 0; i < NCells; i++){
  //   // 2D ANGLE GENERATION
  //   double gentheta = dist(gen);
  //   double genphi = 0;
  //   dealii::Tensor<1,2> rnd_fiber;
  //   //generate equivalent fiber directions from the generated theta angle
  //   rnd_fiber[0] = cos(gentheta)*cos(genphi);
  //   rnd_fiber[1] = sin(gentheta)*cos(genphi);
  //   // outfile << rnd_fiber[0] << "," << rnd_fiber[1] << "\n";

  //   dealii::Tensor<1, 2> xaxis;
  //   xaxis[0] = 1;

  //   double theta = rnd_fiber.norm() == 0 ? 0 : dealii::Physics::VectorRelations::angle(rnd_fiber, xaxis);
  //   // double theta = fiberxy.norm() == 0 ? 0 : (fiber[1] > 0 ? dealii::Physics::VectorRelations::angle(fiberxy, xaxis) : -dealii::Physics::VectorRelations::angle(fiberxy, xaxis));
  //   theta = rnd_fiber[1] > 0 ? theta : -theta;
  //   // theta = theta*180/M_PI;

  //   outfile << theta << "\n";
  // }


  //columns
  for (int y = 0; y <= NCellX; y++){
    //rows
    for (int x = 0; x <= NCellX; x++){
      double gentheta = dist(gen);
      double genphi = 0;
      dealii::Tensor<1,2> rnd_fiber;
      //generate equivalent fiber directions from the generated theta angle
      rnd_fiber[0] = cos(gentheta)*cos(genphi);
      rnd_fiber[1] = sin(gentheta)*cos(genphi);
      
      dealii::Tensor<1, 2> xaxis;
      xaxis[0] = 1;

      double theta = rnd_fiber.norm() == 0 ? 0 : dealii::Physics::VectorRelations::angle(rnd_fiber, xaxis);
      // double theta = fiberxy.norm() == 0 ? 0 : (fiber[1] > 0 ? dealii::Physics::VectorRelations::angle(fiberxy, xaxis) : -dealii::Physics::VectorRelations::angle(fiberxy, xaxis));
      theta = rnd_fiber[1] > 0 ? theta : -theta;
      // theta = theta*180/M_PI;
      outfile << theta;
      if (x < NCellX){
        outfile << ",";
      }
    }
    outfile << "\n";
  }
  outfile.close();
  */
  
  //   // // 3D ANGLE GENERATION
  //   // double gentheta = dist(gen);
  //   // //need different generator for phi? seems to be biased towards phi = +- pi/2
  //   // double genphi = dist(gen)-M_PI/2;

  //   // dealii::Tensor<1,3> rnd_fiber, rnd_fiberxy, xaxis;
  //   // //generate equivalent fiber directions from the generated theta angle
  //   // rnd_fiber[0] = cos(gentheta)*cos(genphi);
  //   // rnd_fiber[1] = sin(gentheta)*cos(genphi);
  //   // rnd_fiber[2] = sin(genphi);

  //   // rnd_fiberxy[0] = rnd_fiber[0];
  //   // rnd_fiberxy[1] = rnd_fiber[1];
    
  //   // xaxis[0] = 1;
    
  //   // double theta = rnd_fiberxy.norm() == 0 ? 0 : dealii::Physics::VectorRelations::angle(rnd_fiberxy, xaxis);
  //   // theta = rnd_fiber[1] > 0 ? theta : -theta;
  //   // double phi = rnd_fiberxy.norm() == 0 ? M_PI/2 : dealii::Physics::VectorRelations::angle(rnd_fiber, rnd_fiberxy);
  //   // phi = rnd_fiber[2] > 0 ? phi : -phi;

  //   // outfile << theta << "," << phi << "\n";
  // }
  // outfile.close();

  
  //read parameters file to determine the dimensions present
  Parameters::AllParameters params(paramsPath);

  //iterate through each fluid mesh that was given
  for(const std::string &meshFluid : simMeshFluid){
    //iterate through each solid mesh
    for(const std::string &meshSolid : simMeshSolid){
      //This section has to be hard coded, since the creation of the Sim object requires a constant variable input
      //the value of ‘dims’ is not usable in a constant expression
      if (params.dimension == 2){
        Sim<2> sim;

        if(params.simulation_type == "Solid" || params.simulation_type == "FSI")
          sim.loadSolid(meshSolid);

        if (params.simulation_type == "Fluid" || params.simulation_type == "FSI")
          sim.loadFluid(meshFluid);

        //sim.loadMesh(meshSolid, meshFluid);
        sim.setParams(params);

      } else if (params.dimension == 3){
        Sim<3> sim;

        if(params.simulation_type == "Solid" || params.simulation_type == "FSI")
          sim.loadSolid(meshSolid);

        if (params.simulation_type == "Fluid" || params.simulation_type == "FSI")
          sim.loadFluid(meshFluid);
        
        sim.setParams(params);
        
      } else {
        std::cerr << "Cannot find dimension from parameters file" << std::endl
                  << "Check if " << paramsPath << "exists or has valid dimensions";
        exit(0);
      }

      //define path to current file location
      std::filesystem::path p = std::filesystem::current_path();
      
      std::string outputFolder;
      if (params.simulation_type == "Solid"){
        outputFolder = meshSolid;
      }else if(params.simulation_type == "Fluid"){
        outputFolder = meshFluid;
      }else if(params.simulation_type == "FSI"){
        outputFolder = meshSolid + "_" + meshFluid;
      }  
      //TODO change meshFluid call to a new output file name
      //create folder with a title corresponding to the current fluid mesh name
      std::filesystem::create_directory(p / outputFolder);

      //iterate through each file in the main directory
      for(const auto& dirEntry : std::filesystem::directory_iterator(p)){
        //checks if each file is a .vtu or .pvd file
        //since these are main outputs for each test case, want to move them somewhere safe before starting another simulation
        if (dirEntry.path().extension() == ".vtu" || dirEntry.path().extension() == ".pvd"){
          
          //moves the "selected" outputs to the new folder corresponding to the fluid mesh name
          std::filesystem::rename(p / dirEntry.path().filename(), p / outputFolder / dirEntry.path().filename());
        }
      }
    }
  }

}