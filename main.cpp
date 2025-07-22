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

//import dealII libraries for parallel computing
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

//TODO missing library 'libp4est.so.3'
//#include "../OpenIFEM-dependencies/p4est-2.8.6/src/.libs/*"

//import OpenIFEM libraries
//solid linear elastic solver
#include "mpi_shared_linear_elasticity.h"
//fluid incompressible navier stokes solver
#include "mpi_insim.h"
//fluid-solid interface solver
#include "mpi_fsi.h"
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
extern template class Solid::MPI::SharedLinearElasticity<2>;
extern template class Solid::MPI::SharedLinearElasticity<3>;
//create fluid objects
extern template class Fluid::MPI::InsIM<2>;
extern template class Fluid::MPI::InsIM<3>;
//fluid-solid interface objects
extern template class MPI::FSI<2>;
extern template class MPI::FSI<3>;

using namespace dealii;

const double L = 4, H = 1, a = 0.1, b = 0.4, h = 0.05, U = 1.5;

int main(int argc, char *argv[])
{

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      std::string infile("parameters.prm");
      if (argc > 1)
        {
          infile = argv[1];
        }
      Parameters::AllParameters params(infile);

      auto inflow_bc = [U = U](const Point<2> &p,
                               const unsigned int component,
                               const double time) -> double {
        (void)time;
        if (component == 0 && std::abs(p[0]) < 1e-10)
          {
            return U;
          }
        return 0.0;
      };

      auto inflow_bc_3d = [U = U](const Point<3> &p,
                                  const unsigned int component,
                                  const double time) -> double {
        (void)time;
        if (component == 2 && std::abs(p[2]) < 1e-10 &&
            std::abs(p[0]) > 1e-10 && std::abs(p[1]) > 1e-10)
          {
            return U;
          }
        return 0.0;
      };

      if (params.dimension == 2)
        {
          parallel::distributed::Triangulation<2> fluid_tria(MPI_COMM_WORLD);
          
          //load fluid tria
          dealii::GridGenerator::subdivided_hyper_rectangle(
            fluid_tria,
            {static_cast<unsigned int>(L / h),
             static_cast<unsigned int>(H / h)},
            Point<2>(0, 0),
            Point<2>(L, H),
            true);

          Fluid::MPI::InsIM<2> fluid(fluid_tria, params);
          fluid.add_hard_coded_boundary_condition(0, inflow_bc);

          Triangulation<2> solid_tria;
          //load solid tria
          dealii::GridGenerator::subdivided_hyper_rectangle(
            solid_tria,
            {static_cast<unsigned int>(a / h),
             static_cast<unsigned int>(b / h)},
            Point<2>(L / 4, 0),
            Point<2>(a + L / 4, b),
            true);
          Solid::MPI::SharedLinearElasticity<2> solid(solid_tria, params);

          MPI::FSI<2> fsi(fluid, solid, params, true);
          fsi.run();
        }
      else
        {
          parallel::distributed::Triangulation<3> fluid_tria(MPI_COMM_WORLD);
          dealii::GridGenerator::subdivided_hyper_rectangle(
            fluid_tria,
            {static_cast<unsigned int>(H / (2 * h)),
             static_cast<unsigned int>(H / (2 * h)),
             static_cast<unsigned int>(L / (2 * h))},
            Point<3>(0, 0, 0),
            Point<3>(H, H, L),
            true);
          Fluid::MPI::InsIM<3> fluid(fluid_tria, params);
          fluid.add_hard_coded_boundary_condition(4, inflow_bc_3d);

          Triangulation<3> solid_tria;
          dealii::GridGenerator::subdivided_hyper_rectangle(
            solid_tria,
            {static_cast<unsigned int>(b / (1 * h)),
             static_cast<unsigned int>(a / (1 * h)),
             static_cast<unsigned int>(a / (1 * h))},
            Point<3>(0, (H - a) / 2, L / 4),
            Point<3>(b, (H + a) / 2, a + L / 4),
            true);
          Solid::MPI::SharedLinearElasticity<3> solid(solid_tria, params);

          MPI::FSI<3> fsi(fluid, solid, params, true);
          fsi.run();
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}

// //Vars in unnamed namespace to avoid reading from other files
// template <int dim>
// class Sim{
//   public:
//     //extern template class Solid::LinearElasticity<dim>;
//     Sim();
//     void loadSolid(std::string solidMeshName);
//     void loadFluid(std::string fluidMeshName);
//     void setParams(Parameters::AllParameters params);
//     //Triangulation<3> extrude(); 
//     //int refine(int refinement);
//     //is this expecting to call a function?
//     // Solid::LinearElasticity<dim> solid(Triangulation<dim> triaSolid, Parameters::AllParameters params);
//     // Fluid::InsIM<dim> fluid(Triangulation<dim> triaFluid, Parameters::AllParameters params);

//   private:
//     parallel::distributed::Triangulation<dim> triaSolid, triaFluid;
//     DoFHandler<dim> dof_handler;
//     GridIn<dim> gridIn;

// };

// namespace {
// const std::string simMeshSolid[] = {""};
// //Ability to set multiple fluid meshes to simplify fluid mesh refinement studies
// const std::string simMeshFluid[] = {""};
// const std::string meshPath = "meshes/";
// const std::string paramsPath = "parameters.prm";
// GridOut gridOut;
// }

// //only need to define dof_handler once for the sim dimensions, i guess this is how it reads what dim to use?
// template <int dim>
// Sim<dim>::Sim()
// //error: no matching function for call to ‘dealii::parallel::distributed::Triangulation<2, 2>::Triangulation()
//   : dof_handler(triaSolid) 
// {}

// //Solid object creation
// template <int dim>
// void Sim<dim>::loadSolid(std::string solidMeshName){
//   // Dynamically define path for solid mesh location
//   std::ifstream solidPath(meshPath + solidMeshName + ".msh");
//   // Check if given mesh path is valid
//   if (!solidPath){
//     std::cerr << "----------------------------------------------------" << "\n"
//               << "ERROR FINDING SOLID MESH FILE " << solidMeshName  << "\n"
//               << "----------------------------------------------------" << "\n";
//     //exit the program
//     exit(0);
//   }

//   //define GridIn object to receive 2d mesh
//   gridIn.attach_triangulation(triaSolid);
//   //imports mesh from selected area
//   gridIn.read_msh(solidPath);
  
//   return;
// }

// // Fluid object creation
// template <int dim>
// void Sim<dim>::loadFluid(std::string fluidMeshName){
//   // Dynamically define path for fluid mesh location
//   std::ifstream fluidPath(meshPath + fluidMeshName + ".msh");
//   // Check if given mesh path is valid
//   if (!fluidPath){
//     std::cerr << "----------------------------------------------------" << "\n"
//               << "ERROR FINDING FLUID MESH FILE " << fluidMeshName << "\n"
//               << "----------------------------------------------------" << "\n";
//     exit(0);
//   }

//   //define GridIn object to receive fluid mesh
//   gridIn.attach_triangulation(triaFluid);
//   //imports the fluid mesh from the valid file path
//   gridIn.read_msh(fluidPath);

//   return;
//   //Example of exporting fluid mesh for debugging/checking what is used
//   /*
//   std::ofstream out(meshNameFluid + ".msh");
//   gridOut.write_msh(triaFluid, out);
//   */
// }

// template <int dim>
// void Sim<dim>::setParams(Parameters::AllParameters params){
//   //import params for both solid and fluid meshes separately
//   //Solid::LinearElasticity<dim> solid(triaSolid, params);
//   //Fluid::InsIM<dim> fluid(triaFluid, params);
  
//   if (params.simulation_type == "Solid"){
//     Solid::MPI::SharedLinearElasticity<dim> solid(triaSolid, params);
//     solid.run();
//   }else if(params.simulation_type == "Fluid"){
//     Fluid::MPI::InsIM<dim> fluid(triaFluid, params);
//     fluid.run();
//   }else if(params.simulation_type == "FSI"){
//     //combine solid and fluid meshes to make FSI simulation
//     Solid::MPI::SharedLinearElasticity<dim> solid(triaSolid, params);
//     Fluid::MPI::InsIM<dim> fluid(triaFluid, params);
//     MPI::FSI<dim> fsi(fluid, solid, params, true);
//     fsi.run();
//   }else{
//     //error occured
//     exit(0);
//   }
//   return;
// }

// int main(){
//   //read parameters file to determine the dimensions present
//   Parameters::AllParameters params(paramsPath);

//   //iterate through each fluid mesh that was given
//   for(const std::string &meshFluid : simMeshFluid){
//     //iterate through each solid mesh
//     for(const std::string &meshSolid : simMeshSolid){
//       //This section has to be hard coded, since the creation of the Sim object requires a constant variable input
//       //the value of ‘dims’ is not usable in a constant expression
//       if (params.dimension == 2){
//         Sim<2> sim;

//         if(params.simulation_type == "Solid" || params.simulation_type == "FSI")
//           sim.loadSolid(meshSolid);

//         if (params.simulation_type == "Fluid" || params.simulation_type == "FSI")
//           sim.loadFluid(meshFluid);

//         //sim.loadMesh(meshSolid, meshFluid);
//         sim.setParams(params);

//       } else if (params.dimension == 3){
//         Sim<3> sim;

//         if(params.simulation_type == "Solid" || params.simulation_type == "FSI")
//           sim.loadSolid(meshSolid);

//         if (params.simulation_type == "Fluid" || params.simulation_type == "FSI")
//           sim.loadFluid(meshFluid);
        
//         sim.setParams(params);
        
//       } else {
//         std::cerr << "Cannot find dimension from parameters file" << std::endl
//                   << "Check if " << paramsPath << "exists or has valid dimensions";
//         exit(0);
//       }

//       //define path to current file location
//       std::filesystem::path p = std::filesystem::current_path();
      
//       std::string outputFolder;
//       if (params.simulation_type == "Solid"){
//         outputFolder = meshSolid;
//       }else if(params.simulation_type == "Fluid"){
//         outputFolder = meshFluid;
//       }else if(params.simulation_type == "FSI"){
//         outputFolder = meshSolid + "_" + meshFluid;
//       }  
//       //TODO change meshFluid call to a new output file name
//       //create folder with a title corresponding to the current fluid mesh name
//       std::filesystem::create_directory(p / outputFolder);

//       //iterate through each file in the main directory
//       for(const auto& dirEntry : std::filesystem::directory_iterator(p)){
//         //checks if each file is a .vtu or .pvd file
//         //since these are main outputs for each test case, want to move them somewhere safe before starting another simulation
//         if (dirEntry.path().extension() == ".vtu" || dirEntry.path().extension() == ".pvd"){
          
//           //moves the "selected" outputs to the new folder corresponding to the fluid mesh name
//           std::filesystem::rename(p / dirEntry.path().filename(), p / outputFolder / dirEntry.path().filename());
//         }
//       }
    
//     }
//   }
// }