/*
 * NeoHookean_Newton_CompressedStrip.h
 *
 *  Created on: Aug 10, 2017
 *      Author: andrew
 */

#ifndef EFFECTIVE_PLATE_H_
#define EFFECTIVE_PLATE_H_

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/arpack_solver.h>



#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/slepc_solver.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Constituitive.h"


#define MAXLINE 1024
#define DIM 2

namespace effective_plate
{
  using namespace dealii;
  /****************************************************************
                       Class Declarations
  ****************************************************************/


  class WX : public Function<DIM>
  {

  public:
      WX () : Function<DIM>() {}
    virtual ~WX (){}

    // const double PI = std::atan(1.0)*4;

    virtual double value (const Point<DIM> &p,
                          const unsigned int  component = 0) const;

    virtual void value_list(const std::vector< Point< DIM > > &  points,
                             std::vector< double > &   values,
                             const unsigned int  component = 0 )  const;

    void set_param_values(double lx_, double ly_, double delta_);


    double lx = 1.0;
    double ly = 1.0;
    double delta = 1.0;

  };

  class WY : public Function<DIM>
  {

  public:
      WY () : Function<DIM>() {}
    virtual ~WY (){}

    // const double PI = std::atan(1.0)*4;

    virtual double value (const Point<DIM> &p,
                          const unsigned int  component = 0) const;

    virtual void value_list(const std::vector< Point< DIM > > &  points,
                             std::vector< double > &   values,
                             const unsigned int  component = 0 )   const;

    void set_param_values(double lx_, double ly_, double delta_);

    double lx = 1.0;
    double ly = 1.0;
    double delta = 1.0;

  };

  /****  ElasticProblem  *****
   * This is the primary class used, with all the dealii stuff
   */
  class PlateProblem
  {
  public:
    PlateProblem();
    ~PlateProblem();

    void create_mesh();
    void deform_mesh();
    void setup_system ();
    void newton_iterate();

    void output_results(const unsigned int cycle) const;

    void read_input_file(char* filename);

    void save_current_state(unsigned int indx);
    void load_state(unsigned int indx);

    void rhs_numerical_deriv(double delta);

    void solve_forward_problem();

    void compute_and_compare_second_numer_deriv(double epsilon);

    // get methods for important constants

    unsigned int get_n_dofs(){return dof_handler.n_dofs();};
    unsigned int get_number_active_cells(){return triangulation.n_active_cells();};

    Vector<double>       present_solution;
    Vector<double>       evaluation_point;



  private:


    void right_hand_side (const std::vector<Point<DIM> > &points,
                          std::vector<Tensor<1, DIM> >   &values);

    void setup_system_constraints();

    void make_symmetry_constraints();

    void apply_boundaries_to_rhs(Vector<double> *rhs, std::vector<bool> *homogenous_dirichlet_dofs);

    double line_search_and_add_step_length(double last_residual, std::vector<bool> *homogenous_dirichlet_dofs);

    void assemble_system_energy();
    void assemble_system_rhs();
    void assemble_system_matrix();

    void assemble_wv_matrices();

    unsigned int get_system_eigenvalues(const int cycle);

    void apply_boundaries_and_constraints_system_matrix(std::vector<bool> *homogenous_dirichlet_dofs);
    void apply_boundaries_and_eigen_constraints_system_matrix(std::vector<bool> *homogenous_dirichlet_dofs);

    void apply_boundaries_to_matrix(SparseMatrix<double> &mat,
                         std::vector<bool> *homogenous_dirichlet_dofs,
                         bool zeroFlag = false);


    void solve();

    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                            int const maxSize, int* const endOfFileFlag);

    void add_small_perturbations(double amplitude, bool firstTime);

    void output_matrix_csv();

    void set_displacement(double value)
    {
      for(unsigned int i = 0; i < dof_handler.n_dofs(); i++)
        if(load_dofs[i] == true)
          present_solution[i] = value * load_factors[i];
    }




    Triangulation<DIM,DIM>   triangulation;
    DoFHandler<DIM>      dof_handler;

    DoFHandler<DIM>      dof_handler_wv;

    FESystem<DIM>        fe;

    FESystem<DIM>     fe_wv;

    Quadrature<DIM>      *quadrature_formula = NULL;

    ConstraintMatrix     constraints;
    ConstraintMatrix     constraints_wv;


    ConstraintMatrix     constraints_eigen;


    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    SparsityPattern      sparsity_pattern_wv;
    SparseMatrix<double> K_wlam;
    SparseMatrix<double> K_vlam;
    SparseMatrix<double> K_vv;


    PETScWrappers::SparseMatrix    system_matrix_petsc;


    plate_energy PE;


    Vector<double>       newton_update;
    Vector<double>       system_rhs;

    char output_directory[MAXLINE];

    Vector<double>       previous_solution;

    std::vector<bool> homo_dofs;
    std::vector<bool> homo_dofs_v;
    std::vector<bool> homo_dofs_w;

    std::vector<bool> load_dofs;

    std::vector<double> load_factors;

    std::vector<unsigned int>  grid_dimensions;
    std::vector<double> domain_dimensions;

    WX eval_wx;
    WY eval_wy;


    unsigned int iter = 0;

    unsigned int N = 1;
    unsigned int qx = 2;
    unsigned int qy = 2;

    double lx = 1.0;
    double ly = 1.0;
    double delta = 1.0;

    double mu = 1.0;
    double lambda = 1.0;
    double B = 1.0;
    double eta = 1.0;

    double load_val = 0.01;
    unsigned int load_steps = 10;

    double newton_tol = 1e-8;
    unsigned int newton_max_iter = 300;

    double system_energy = 0.0;
    unsigned int middle_dof = 0;
    double eps_w = 1.0e-5;
    double r_w = 300.0;


    bool firstFlag = true;

    std::vector<unsigned int> corner_dofs;

  };
}

#endif /* EFFECTIVE_PLATE_H_ */
