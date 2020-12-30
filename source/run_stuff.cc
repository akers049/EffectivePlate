#include <fstream>
#include "../include/effectivePlate.h"

using namespace dealii;
int main (int argc, char** argv)
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  effective_plate::PlateProblem  pp;

  char fileName[MAXLINE];
  std::cout << "Please enter an input file: " << std::endl;
  std::cin >> fileName;
  pp.read_input_file(fileName);

  pp.create_mesh();

  pp.setup_system();

  std::cout << "\n   Number of active cells:       "
            << pp.get_number_active_cells()
            << std::endl;


  std::cout << "   Number of degrees of freedom: "
            << pp.get_n_dofs()
            << std::endl << std::endl;

  // output inital mesh
  pp.output_results (0);



//  ep.rhs_numerical_deriv(1e-5);
//  ep.compute_and_compare_second_numer_deriv(1e-5);


//  pp.rhs_numerical_deriv(1e-9);
//  pp.compute_and_compare_second_numer_deriv(1e-9);

  pp.solve_forward_problem();



  return 0;

}
