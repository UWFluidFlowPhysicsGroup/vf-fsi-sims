// Base hello world template taken from:
// https://education.molssi.org/parallel-programming/04-distributed-examples.html

#include <iostream>
#include <mpi.h>
#include <deal.II/base/mpi.h>

using namespace dealii;

const double L = 4, H = 1, a = 0.1, b = 0.4, h = 0.05, U = 1.5;

int main(int argc, char **argv) {
  // Using anything related to dealii requires libp4est.so.3, which causes the missing library error
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // std::string infile("parameters.prm");
  //     if (argc > 1)
  //       {
  //         infile = argv[1];
  //       }
  // Parameters::AllParameters params(infile);
  // std::cout << "test";

  //   // Initialize MPI
  //   // This must always be called before any other MPI functions
  //   MPI_Init(&argc, &argv);
  // // Finalize MPI
  // // This must always be called after all other MPI functions
  // MPI_Finalize();

  return 0;
}
