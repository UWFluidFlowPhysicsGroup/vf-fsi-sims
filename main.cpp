#include "mpi_fsi.h"
#include "mpi_scnsim.h"
#include "mpi_shared_linear_elasticity.h"
#include "parameters.h"
#include "utilities.h"
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

extern template class Fluid::MPI::SCnsIM<2>;
extern template class Fluid::MPI::SCnsIM<3>;
extern template class Solid::MPI::SharedLinearElasticity<2>;
extern template class Solid::MPI::SharedLinearElasticity<3>;
extern template class MPI::FSI<2>;
extern template class MPI::FSI<3>;

using namespace dealii;

// template <int dim>
// class SigmaPMLField : public Function<dim>
// {
// public:
//   SigmaPMLField(double sig, double l)
//     : Function<dim>(), SigmaPMLMax(sig), PMLLength(l)
//   {
//   }
//   virtual double value(const Point<dim> &p,
//                        const unsigned int component = 0) const;
//   virtual void value_list(const std::vector<Point<dim>> &points,
//                           std::vector<double> &values,
//                           const unsigned int component = 0) const;

// private:
//   double SigmaPMLMax;
//   double PMLLength;
// };

// template <int dim>
// class ArtificialBF : public TensorFunction<1, dim>
// {
// public:
// 	ArtificialBF() : TensorFunction<1, dim>() {}
// 	virtual Tensor<1, dim> value(const Point<dim> &p) const;
// 	virtual void value_list(const std::vector< Point< dim >> &points, 
// 							std::vector<Tensor<1, dim>> &values) const;
// };

// template <int dim>
// double SigmaPMLField<dim>::value(const Point<dim> &p,
//                                  const unsigned int component) const
// {
//   (void)component;
//   (void)p;
//   double SigmaPML = 0.0;
//   std::vector<double> boundary = {0.0, 51.5, -20.0};
//   std::vector<unsigned int> boundary_dir = {0, 0, 1};
//   // For tube acoustics
//   for (unsigned int i = 0; i < boundary.size(); ++i)
//      {
//        if (std::abs(p[boundary_dir[i]] - boundary[i]) < PMLLength)
//        // A quadratic increasing function from boundary-PMLlength to the boundary
//        SigmaPML = SigmaPMLMax * pow((PMLLength + boundary[i] - p[boundary_dir[i]]) / PMLLength, 4);
//      }
//   return SigmaPML;
// }

// template <int dim>
// void SigmaPMLField<dim>::value_list(const std::vector<Point<dim>> &points,
//                                     std::vector<double> &values,
//                                     const unsigned int component) const
// {
//   (void)component;
//   for (unsigned int i = 0; i < points.size(); ++i)
//     values[i] = this->value(points[i]);
// }

// template <int dim>
// Tensor<1, dim> ArtificialBF<dim>::value(const Point<dim> &p) const
// {
// 	Tensor<1, dim> value;
// 	double rho = 1.3e-3;
// 	double bf = 1.0e4/rho;
// 	if (p[0] > 4.0 - 5e-4 && p[0] < 5.0 + 5e-4)
// 		value[0] = bf;
// 	return value;
// }

// template <int dim>
// void ArtificialBF<dim>::value_list(const std::vector< Point< dim >> &points, 
// 							std::vector<Tensor<1, dim>> &values) const
// {
// 	for (unsigned int i = 0; i < points.size(); ++i)
//     values[i] = this->value(points[i]);
// }

int main(int argc, char *argv[])
{
  using namespace dealii;

  // double PMLlength = 3.0, SigmaMax = 340000;

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      std::string infile("parameters.prm");
      if (argc > 1)
        {
          infile = argv[1];
        }
      Parameters::AllParameters params(infile);

      if (params.dimension == 2)
        {
          // Read fluid mesh
          parallel::distributed::Triangulation<2> tria_fluid(MPI_COMM_WORLD);
          GridIn<2> grid_in;
          std::ifstream input_fluid("vocal_tract.inp");
          grid_in.attach_triangulation(tria_fluid);
          grid_in.read_abaqus(input_fluid);

	  // locally refine the glottal area
	  for (auto & cell : tria_fluid)
	  {
	    for (unsigned int i : {0, 1, 2, 3})
	    {
	      double x = cell.vertex(i)(0);
	      if (x > 12.7 && x < 15.0)
	      {
                cell.set_refine_flag();
	      }
            }
	  }
          tria_fluid.execute_coarsening_and_refinement();

          // Read solid mesh
          // Triangulation<2> tria_solid;
          // std::ifstream input_solid("vocal_fold.inp");
          // grid_in.attach_triangulation(tria_solid);
          // grid_in.read_abaqus(input_solid);

          Triangulation<2> tria_solid;
          std::ifstream input_solid("BLC_Ellipse_cm.msh");
          grid_in.attach_triangulation(tria_solid);
          grid_in.read_msh(input_solid);

          // Translate solid mesh
          Tensor<1, 2> offset({12.7, 0});
          GridTools::shift(offset, tria_solid);

          // // placeholder for hard code BC
          // std::shared_ptr<Functions::ZeroFunction<2>> ptr =
          //   std::make_shared<Functions::ZeroFunction<2>>(
          //     Functions::ZeroFunction<2>(3));
          // // initialize the pml field
          // auto pml = std::make_shared<SigmaPMLField<2>>(
          //   SigmaPMLField<2>(SigmaMax, PMLlength));
          // // artificial body force
          // auto bf_ptr = std::make_shared<ArtificialBF<2>>(
          // 	ArtificialBF<2>());
          
          // Fluid::MPI::SCnsIM<2> fluid(tria_fluid, params, ptr, pml, bf_ptr);
          Fluid::MPI::SCnsIM<2> fluid(tria_fluid, params);
          Solid::MPI::SharedLinearElasticity<2> solid(tria_solid, params);
          
          auto pml = [](const Point<2> &p,
              const unsigned int component) -> double {
            (void)component;
            // (void)p;
            double SigmaPML = 0.0;
            double SigmaPMLMax = 340000;
            double PMLLength = 3.0;
            
            std::vector<double> boundary = {0.0, 51.5, -20.0};
            std::vector<unsigned int> boundary_dir = {0, 0, 1};
            // For tube acoustics
            for (unsigned int i = 0; i < boundary.size(); ++i)
              {
                if (std::abs(p[boundary_dir[i]] - boundary[i]) < PMLLength)
                // A quadratic increasing function from boundary-PMLlength to the boundary
                SigmaPML = SigmaPMLMax * pow((PMLLength + boundary[i] - p[boundary_dir[i]]) / PMLLength, 4);
              }
            return SigmaPML;
          };

          fluid.set_sigma_pml_field(pml);

          auto bf = [](const Point<2> &p,
                               const unsigned int component) -> double {
            // Tensor<1, dim> value;
            double rho = 1.3e-3;
            double bf = 1.0e4/rho;

            if (p[0] > 4.0 - 5e-4 && p[0] < 5.0 + 5e-4 && component == 0){
              return bf;
            }
            return 0.0;
          };

          fluid.set_body_force(bf);

          MPI::FSI<2> fsi(fluid, solid, params);
          fsi.run();
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
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
