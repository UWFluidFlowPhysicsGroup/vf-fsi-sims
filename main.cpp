// Build commands to set up libraries, need to export p4est library path each time Ubuntu is launched:
// not needed if added to .bashrc file in home directory
// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/p4est_build/local/lib

// modified CMakeLists to not need -I and -L commands anymore, still need -l commands
// sudo make main -I $HOME/p4est_build/local/include -L $HOME/p4est_build/local/lib -lp4est -lsc -lz -lm
// mpiexec -n 4 main

// Replace 4 with number of cores you wish to run

//https://education.molssi.org/parallel-programming/04-distributed-examples.html
//https://stackoverflow.com/questions/23163075/how-to-compile-an-mpi-included-c-program-using-cmake
//https://hpc-discourse.usc.edu/t/use-cmake-in-an-mpi-c-program/507/4
//https://stackoverflow.com/questions/11368215/loading-shared-library-in-open-mpi-mpi-run

//https://p4est.github.io/api/p4est-latest/installing_p4est.html
//https://education.molssi.org/parallel-programming/04-distributed-examples.html

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

#include <mpi.h>
#include <deal.II/base/mpi.h>

//import OpenIFEM libraries
//solid linear elastic solver
#include "mpi_shared_linear_elasticity.h"
//fluid slightly compressible navier stokes solver, used cause it seems more stable for simulations?
#include "mpi_scnsim.h"
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
extern template class Fluid::MPI::SCnsIM<2>;
extern template class Fluid::MPI::SCnsIM<3>;
// extern template class Fluid::MPI::InsIM<2>;
// extern template class Fluid::MPI::InsIM<3>;

//fluid-solid interface objects
extern template class MPI::FSI<2>;
extern template class MPI::FSI<3>;

using namespace dealii;

int main(int argc, char *argv[]){
  // input mesh names for fluid and solid meshes here, using arrays to automate mesh refinement studies or other meshes as long as parameters match
  const std::string simMeshSolid[] = {"M5_BC_Half"};
  const std::string simMeshFluid[] = {"Fluid_Half_Gradient"};
  const std::string meshPath = "meshes/";
  const std::string paramsPath = "parameters_2D_BC_Half.prm";
  //read parameters file to determine the dimensions present
  Parameters::AllParameters params(paramsPath);
  GridOut gridOut;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  //iterate through each fluid mesh that was given
  for(const std::string &meshFluid : simMeshFluid){
    //iterate through each solid mesh
    for(const std::string &meshSolid : simMeshSolid){
      //This section has to be hard coded, since the creation of the Sim object requires a constant variable input
      //the value of ‘dims’ is not usable in a constant expression
      if (params.dimension == 2){
        //solid is still regular triangulation object, not distributed
        Triangulation<2> triaSolid;
        parallel::distributed::Triangulation<2> triaFluid(MPI_COMM_WORLD);
        DoFHandler<2> dof_handler;
        GridIn<2> gridIn;

        if(params.simulation_type == "Solid" || params.simulation_type == "FSI"){
          // Dynamically define path for solid mesh location
          std::ifstream solidPath(meshPath + meshSolid + ".msh");
          // Check if given mesh path is valid
          if (!solidPath){
            std::cerr << "----------------------------------------------------" << "\n"
                      << "ERROR FINDING SOLID MESH FILE " << meshSolid  << "\n"
                      << "----------------------------------------------------" << "\n";
            //exit the program
            exit(0);
          }
          //define GridIn object to receive 2d mesh
          gridIn.attach_triangulation(triaSolid);
          //imports mesh from selected area
          gridIn.read_msh(solidPath);
        }

        if (params.simulation_type == "Fluid" || params.simulation_type == "FSI"){
          // Dynamically define path for fluid mesh location
          std::ifstream fluidPath(meshPath + meshFluid + ".msh");
          // Check if given mesh path is valid
          if (!fluidPath){
            std::cerr << "----------------------------------------------------" << "\n"
                      << "ERROR FINDING FLUID MESH FILE " << meshFluid << "\n"
                      << "----------------------------------------------------" << "\n";
            exit(0);
          }
          //define GridIn object to receive fluid mesh
          gridIn.attach_triangulation(triaFluid);
          //imports the fluid mesh from the valid file path
          gridIn.read_msh(fluidPath);
        }

        //sim.loadMesh(meshSolid, meshFluid);
          if (params.simulation_type == "Solid"){
            Solid::MPI::SharedLinearElasticity<2> solid(triaSolid, params);
            solid.run();
          }else if(params.simulation_type == "Fluid"){
            Fluid::MPI::SCnsIM<2> fluid(triaFluid, params);
            fluid.run();
          }else if(params.simulation_type == "FSI"){
            //combine solid and fluid meshes to make FSI simulation
            Solid::MPI::SharedLinearElasticity<2> solid(triaSolid, params);
            Fluid::MPI::SCnsIM<2> fluid(triaFluid, params);
            
            // Define penetration criterion, incompressible plane along mid plane
            auto penetration_criterion = [](const Point<2> &p, const Point<2> &disp) -> double {
              double midplane = 8.5;
              double gap = 0.1;
              // Check if point is between min and max collision zone
              if ((p[1] >  (midplane - gap/2)) && (p[1] < (midplane + gap/2))){
                if (disp[1] < 0){
                  return (midplane + gap/2 - p[1]);
                }else{
                  return (p[1] - (midplane-gap/2));
                }
              }
              return (0);
            };

            // keeping vertex point data for futureproofing mpi_fsi class, could also overload function instead
            auto penetration_direction = [](const Point<2> &p, const Point<2> &disp) -> Tensor<1, 2> {
              return (Tensor<1,2>({0,-disp[1]}));
            };

            double PMLlength = 50, SigmaMax = 340000;
            auto sigma_pml_field =
            [PMLlength, SigmaMax](const Point<2> &p, const unsigned int component) {
              (void)component;
              double SigmaPML = 0.0, L = 100;
              // For tube acoustics
              if (p[0] > L - PMLlength)
                {
                  // A quadratic increasing function from boundary-PMLlength to the
                  // boundary
                  SigmaPML =
                    SigmaMax * pow((p[0] + PMLlength - L) / PMLlength, 4);
                }
              return SigmaPML;
            };
            
            fluid.set_sigma_pml_field(sigma_pml_field);

            MPI::FSI<2> fsi(fluid, solid, params, true);
            
            //apply penetration criterion to simulation, can only set one penetration criterion for the whole model
            fsi.set_penetration_criterion(penetration_criterion);
            fsi.set_penetration_direction(penetration_direction);

            fsi.run();
          }else{
            //error occured
            exit(0);
          }

      } else if (params.dimension == 3){
        Triangulation<3> triaSolid;
        parallel::distributed::Triangulation<3> triaFluid(MPI_COMM_WORLD);
        DoFHandler<3> dof_handler;
        GridIn<3> gridIn;

        if(params.simulation_type == "Solid" || params.simulation_type == "FSI"){
          // Dynamically define path for solid mesh location
          std::ifstream solidPath(meshPath + meshSolid + ".msh");
          // Check if given mesh path is valid
          if (!solidPath){
            std::cerr << "----------------------------------------------------" << "\n"
                      << "ERROR FINDING SOLID MESH FILE " << meshSolid  << "\n"
                      << "----------------------------------------------------" << "\n";
            //exit the program
            exit(0);
          }
          //define GridIn object to receive 2d mesh
          gridIn.attach_triangulation(triaSolid);
          //imports mesh from selected area
          gridIn.read_msh(solidPath);
        }
        if (params.simulation_type == "Fluid" || params.simulation_type == "FSI"){
          // Dynamically define path for fluid mesh location
          std::ifstream fluidPath(meshPath + meshFluid + ".msh");
          // Check if given mesh path is valid
          if (!fluidPath){
            std::cerr << "----------------------------------------------------" << "\n"
                      << "ERROR FINDING FLUID MESH FILE " << meshFluid << "\n"
                      << "----------------------------------------------------" << "\n";
            exit(0);
          }
          //define GridIn object to receive fluid mesh
          gridIn.attach_triangulation(triaFluid);
          //imports the fluid mesh from the valid file path
          gridIn.read_msh(fluidPath);
        }

        //sim.loadMesh(meshSolid, meshFluid);
          if (params.simulation_type == "Solid"){
            Solid::MPI::SharedLinearElasticity<3> solid(triaSolid, params);
            solid.run();
          }else if(params.simulation_type == "Fluid"){
            Fluid::MPI::SCnsIM<3> fluid(triaFluid, params);
            fluid.run();
          }else if(params.simulation_type == "FSI"){
            //combine solid and fluid meshes to make FSI simulation
            Solid::MPI::SharedLinearElasticity<3> solid(triaSolid, params);
            Fluid::MPI::SCnsIM<3> fluid(triaFluid, params);
            MPI::FSI<3> fsi(fluid, solid, params, true);
            fsi.run();
          }else{
            //error occured
            exit(0);
          }

        
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
      //moving file system info is broken for mpi, maybe need to stop mpi connection first?
      //create folder with a title corresponding to the current solid/fluid mesh names
      std::filesystem::create_directory(p / outputFolder);

      //iterate through each file in the main directory
      for(const auto& dirEntry : std::filesystem::directory_iterator(p)){
        //checks if each file is relevant to simulation results/output
          //vtu -> info from separate segmented meshes, one for each processor being used
          //pvtu -> joins vtu files together for a single timestep, only needed for parallel processes
          //pvd -> joins pvtu/vtu files together through whole simulation
          
        //If files are not moved, then simulations will be overwritten with following simulations
        if (dirEntry.path().extension() == ".vtu" || dirEntry.path().extension() == ".pvd" || dirEntry.path().extension() == ".pvtu"){
          
          //moves the "selected" outputs to the new folder corresponding to the fluid mesh name
          std::filesystem::rename(p / dirEntry.path().filename(), p / outputFolder / dirEntry.path().filename());
        }
      }
    }
  }
}
