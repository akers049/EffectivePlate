/* ---------------------------------------------------------------------
 *
 *
 *
 *
 *
 * ---------------------------------------------------------------------
 *
 *
 * Author: Andrew Akerson
 */

#ifndef EFFECTIVE_PLATE_CC_
#define EFFECTIVE_PLATE_CC_
#include "../include/effectivePlate.h"

#include <fstream>
#include <iostream>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>


#define DIM 2



namespace effective_plate
{
  using namespace dealii;



  inline
  void get_F(Tensor<2,DIM> &grad_u, Tensor<2,DIM> &F)
  {
    F = grad_u;
    for(unsigned int i = 0; i < DIM; i ++)
      F[i][i] += 1.0;
  }


  double WX::value (const Point<DIM>  &p,
                                   const unsigned int  component) const
  {
    // this is a scalar function, so make sure component is zero...
    Assert (component == 0, ExcNotImplemented());
    Assert (p.dimension == 2, ExcNotImplemented())

    // Put your function for WX. p(0) is x1 value, p(1) is the x2 value

    double wx = (lx - delta)/2.0;

    return wx;
  }

  void WX::value_list(const std::vector< Point< DIM > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = WX::value(points[i], component);

  }

  void WX::set_param_values(double lx_, double ly_, double delta_)
  {
    lx = lx_;
    ly = ly_;
    delta = delta_;
  }

  double WY::value (const Point<DIM>  &p,
                    const unsigned int  component) const
  {
    // this is a scalar function, so make sure component is zero...
    Assert (component == 0, ExcNotImplemented());
    Assert (p.dimension == 2, ExcNotImplemented())

    // Put your function for WY. p(0) is x1 value, p(1) is the x2 value
    double x1 = p[0];
    double x2 = p[1];

    double wy = ((ly - delta)/2.0)*(1.0 - sin(x2*M_PI));

    return wy;
  }

  void WY::value_list(const std::vector< Point< DIM > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = WY::value(points[i], component);

  }

  void WY::set_param_values(double lx_, double ly_, double delta_)
  {
    lx = lx_;
    ly = ly_;
    delta = delta_;
  }



  PlateProblem::PlateProblem ()
    :
    dof_handler (triangulation),
    fe (FESystem<DIM>(FE_Q<DIM>(1), DIM), 1,
         FE_Q<DIM>(2), 1)
  {}




  PlateProblem::~PlateProblem ()
  {
    dof_handler.clear ();
    delete quadrature_formula;

  }


  void PlateProblem::create_mesh()
  {

    // creates our strip.
    Point<DIM> corner1, corner2;
    corner1(0) =  0;
    corner1(1) =  -domain_dimensions[1]/2.0;
    corner2(0) =  domain_dimensions[0];
    corner2(1) =  domain_dimensions[1]/2.0;

    GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, true);

    // Make sure to renumber the boundaries
    renumber_boundary_ids();

  }



  void PlateProblem::setup_system ()
  {
    // Sets up system. Makes the constraint matrix, and reinitializes the
    // vectors and matricies used throughout to be the proper size.

    N = triangulation.n_active_cells();

    dof_handler.distribute_dofs (fe);

    QGauss<1> quad_x(qx);
    QGauss<1> quad_y(qy);
    quadrature_formula = new QAnisotropic<DIM>(quad_x, quad_y);



    PE.set_moduli(mu, lambda, B, eta);

    present_solution.reinit (dof_handler.n_dofs());

    setup_system_constraints();

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);


    sparsity_pattern.copy_from (dsp);


//    GridTools::distort_random(0.4, triangulation, true);


    // MAKE SURE TO SET UP THE DOFS WE WILL APPLY HOMOGENOUS CONDITIONS TO
//    system_matrix.reinit (sparsity_pattern);
//
//    evaluation_point = present_solution;
//
//    // get the dofs that we will apply dirichlet condition to
//    homo_dofs.resize(dof_handler.n_dofs(), false);
//    load_dofs.resize(dof_handler.n_dofs(), false);
//
//    std::set< types::boundary_id > boundary_id_12;
//    boundary_id_12.insert(1);
//    boundary_id_12.insert(2);
//
//    std::set< types::boundary_id > boundary_id_2;
//    boundary_id_2.insert(2);
//
//    std::vector<bool> x1_comp = {true, false, false, false};
//    ComponentMask x1_mask(x1_comp);
//
//    std::vector<bool> all_components = {true, true, true, true};
//    ComponentMask all_mask(all_components);
//
//    DoFTools::extract_boundary_dofs(dof_handler,
//                                       all_mask,
//                                       homo_dofs,
//                                       boundary_id_12);
//
//    DoFTools::extract_boundary_dofs(dof_handler,
//                                       x1_mask,
//                                       load_dofs,
//                                       boundary_id_2);
//
//    a_lambda_dofs.resize(dof_handler.n_dofs(), false);
//    std::vector<bool> a_lambda_comp = {false, false, true, true};
//    ComponentMask a_lambda_mask(a_lambda_comp);
//    DoFTools::extract_dofs(dof_handler, a_lambda_mask, a_lambda_dofs);

  }


  void PlateProblem::setup_system_constraints ()
  {

    constraints.clear ();

    constraints.close ();

    // now do hanging nodes. Because some of the constraints might refer to the same dof
    // for both the symmetry constraint and the hanging node constraint, we will make them
    // separate, then merge them, giving precedence to the hanging node constraints;
    ConstraintMatrix hanging_node_constraints;
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    hanging_node_constraints.close();

    constraints.merge(hanging_node_constraints, ConstraintMatrix::MergeConflictBehavior::right_object_wins);
  }


  void PlateProblem::solve_forward_problem()
  {

    double du = load_val/(1.0*load_steps);
    double current_load_value = 0.0;

    for(unsigned int i = 0; i < load_steps; i ++)
    {
      current_load_value += du;
      set_displacement(current_load_value);

      std::cout << "  Iteration : " << i+1 << std::endl;
      newton_iterate();
      output_results(i);
    }



  }


  void PlateProblem::newton_iterate()
  {
    /* This function executes the newton iterations until it converges to a solution
     * or it exceeds the maximum number of iterations.
     */

    double current_residual;
    unsigned int iteration = 0;

    // Starts by getting the residual in the current configuration.
    evaluation_point = present_solution;

    assemble_system_rhs();
    apply_boundaries_to_rhs(&system_rhs, &homo_dofs);

//    current_residual = system_rhs.l2_norm();

    current_residual = system_rhs.l2_norm();

    // Loops until coverge or go over max iterations
    while((current_residual > newton_tol) &&
             (iteration < newton_max_iter))
    {
      evaluation_point = present_solution;

      apply_boundaries_to_rhs(&system_rhs, &homo_dofs);

//      std::cout << std::endl;

      // Assemble the stiffness matrix
      assemble_system_matrix();
      apply_boundaries_and_constraints_system_matrix(&homo_dofs);

      // solve for the newton step
      solve();

      // Find the step length and add it to the current solution.
      // This function also calls assemble_system_rhs() so we don't need to
      // do another rhs call.

      apply_boundaries_to_rhs(&newton_update, &homo_dofs);

      current_residual = line_search_and_add_step_length(current_residual, &homo_dofs);

      std::cout << "      Iteration "<< iteration << "     Current Residual : " << current_residual << std::endl;

      iteration ++;
    }
    present_solution = evaluation_point;

    // output iterations for convergance.
    std::cout << "    Converging Iterations : " << iteration << "\n";

  }



  double PlateProblem::line_search_and_add_step_length(double last_residual, std::vector<bool> *homogenous_dirichlet_dofs)
  {
   /* this function makes sure that the step sizes we are taking with
    * the newton iteration are making the residual get smaller.
    * Something very similar is used in the dealii step 57?. for now
    * it doesn't do much but we might need to give it more capabilites in the future.
    */

    double current_residual = 0.0;
    double lineSearchLength;
    for(lineSearchLength = 0.9; lineSearchLength > 1e-5; lineSearchLength *=0.5)
    {
      evaluation_point = present_solution;
      evaluation_point.add(lineSearchLength, newton_update);

      assemble_system_rhs();
      apply_boundaries_to_rhs(&system_rhs, homogenous_dirichlet_dofs);

      current_residual = 0.0;

      current_residual = system_rhs.l2_norm();

      if(current_residual < last_residual)
      {
        break;
      }

    }
    if(lineSearchLength < 1e-5)
    {
      std::cout << "THERE WAS A PROBLEM!!" << std::endl;
    }
    present_solution = evaluation_point;

    return current_residual;
  }

  void PlateProblem::assemble_system_energy()
  {
    system_energy = 0.0;

    FEValues<DIM> fe_values (fe, *quadrature_formula,
                             update_values   | update_gradients | update_hessians |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    // unsigned int   n_q_points    = quadrature_formula.size();

    unsigned int n_q_points = quadrature_formula->size();

    std::vector<Tensor<2,DIM>> displacement_gradients(n_q_points);
    std::vector<Tensor<1, DIM>> w_gradients(n_q_points);
    std::vector<Tensor<2, DIM>> w_hessians(n_q_points);


    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector displacements (0);
    const FEValuesExtractors::Scalar w (DIM);

    std::vector<double> wx(n_q_points);
    std::vector<double> wy(n_q_points);

    Tensor<2,DIM> F;

    typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      fe_values[displacements].get_function_gradients(evaluation_point, displacement_gradients);
      fe_values[w].get_function_gradients(evaluation_point, w_gradients);
      fe_values[w].get_function_hessians(evaluation_point, w_hessians);

      eval_wx.value_list (fe_values.get_quadrature_points(), wx);
      eval_wy.value_list (fe_values.get_quadrature_points(), wy);


      unsigned int cell_index = cell->active_cell_index();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        PE.set_params(wx[q_point], wy[q_point], lx, ly, delta);
        get_F(displacement_gradients[q_point], F);
        double lap_w = trace(w_hessians[q_point]);
        system_energy += PE.get_E(F, w_gradients[q_point], lap_w)*fe_values.JxW(q_point);
      }
    }

  }

  void PlateProblem::assemble_system_rhs()
  {
    // Assembling the system rhs. I choose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_rhs = 0.0;

    FEValues<DIM> fe_values (fe, *quadrature_formula,
                             update_values   | update_gradients | update_hessians |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    // unsigned int   n_q_points    = quadrature_formula.size();

    unsigned int n_q_points = quadrature_formula->size();

    std::vector<Tensor<2,DIM>> displacement_gradients(n_q_points);
    std::vector<Tensor<1, DIM>> w_gradients(n_q_points);
    std::vector<Tensor<2, DIM>> w_hessians(n_q_points);

    std::vector<double> wx(n_q_points);
    std::vector<double> wy(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector displacements (0);
    const FEValuesExtractors::Scalar w (DIM);

    Tensor<2,DIM> F;

    DE_data de_dat;

    Vector<double> cell_rhs(dofs_per_cell);

    typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0.0;

      fe_values.reinit (cell);

      fe_values[displacements].get_function_gradients(evaluation_point, displacement_gradients);
      fe_values[w].get_function_gradients(evaluation_point, w_gradients);
      fe_values[w].get_function_hessians(evaluation_point, w_hessians);

      eval_wx.value_list (fe_values.get_quadrature_points(), wx);
      eval_wy.value_list (fe_values.get_quadrature_points(), wy);

      unsigned int cell_index = cell->active_cell_index();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        PE.set_params(wx[q_point], wy[q_point], lx, ly, delta);
        get_F(displacement_gradients[q_point], F);
        double lap_w = trace(w_hessians[q_point]);

        double JxW = fe_values.JxW(q_point);


        PE.get_DE(F, w_gradients[q_point], lap_w, de_dat);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          cell_rhs(n) -= de_dat.dW_dgrad_w*fe_values[w].gradient(n, q_point)*JxW;

          cell_rhs(n) -= de_dat.dW_dlap_w*trace(fe_values[w].hessian(n, q_point))*JxW;

          for(unsigned int i = 0; i<DIM; ++i)
            for(unsigned int j = 0; j<DIM; ++j)
            {
              cell_rhs(n) -= de_dat.dW_dF[i][j]
                             *fe_values[displacements].gradient(n, q_point)[i][j]
                             *JxW;
            }

        }

      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        system_rhs(local_dof_indices[n]) += cell_rhs(n);

    }

    constraints.condense (system_rhs);

  }

  void PlateProblem::assemble_system_matrix()
  {

    time_t start, end;
    time(&start);

    // Assembling the system matrix. I chose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_matrix = 0.0;

    FEValues<DIM> fe_values (fe, *quadrature_formula,
                             update_values   | update_gradients | update_hessians |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    // unsigned int   n_q_points    = quadrature_formula.size();

    unsigned int n_q_points = quadrature_formula->size();

    std::vector<Tensor<2,DIM>> displacement_gradients(n_q_points);
    std::vector<Tensor<1, DIM>> w_gradients(n_q_points);
    std::vector<Tensor<2, DIM>> w_hessians(n_q_points);


    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector displacements (0);
    const FEValuesExtractors::Scalar w (DIM);

    Tensor<2,DIM> F;

    std::vector<double> wx(n_q_points);
    std::vector<double> wy(n_q_points);


    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

//    Tensor<2,DIM> F_inv;
    DDE_data dde_dat;

    typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {

      cell_matrix = 0.0;

      fe_values.reinit (cell);

      fe_values[displacements].get_function_gradients(evaluation_point, displacement_gradients);
      fe_values[w].get_function_gradients(evaluation_point, w_gradients);
      fe_values[w].get_function_hessians(evaluation_point, w_hessians);

      eval_wx.value_list (fe_values.get_quadrature_points(), wx);
      eval_wy.value_list (fe_values.get_quadrature_points(), wy);

      unsigned int cell_index = cell->active_cell_index();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        PE.set_params(wx[q_point], wy[q_point], lx, ly, delta);
        get_F(displacement_gradients[q_point], F);
        double lap_w = trace(w_hessians[q_point]);

        double JxW = fe_values.JxW(q_point);


        PE.get_DDE(F, w_gradients[q_point], lap_w, dde_dat);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
          for (unsigned int m = 0; m < dofs_per_cell; ++m)
          {
            cell_matrix(n,m) += dde_dat.d2W_d2lap_w*trace(fe_values[w].hessian(m, q_point))
                                *trace(fe_values[w].hessian(n, q_point))*JxW;

            for(unsigned int i = 0; i < DIM; i ++)
              for(unsigned int j = 0; j < DIM; j ++)
               {
                  cell_matrix(n,m) += dde_dat.d2W_d2grad_w[i][j]
                                      *fe_values[w].gradient(n, q_point)[i]
                                      *fe_values[w].gradient(m, q_point)[j]*JxW;
                  for(unsigned int k = 0; k < DIM; k ++)
                  {
                    cell_matrix(n,m) += dde_dat.d2W_dF_dgrad_w[i][j][k]
                                        *fe_values[displacements].gradient(m, q_point)[i][j]
                                        *fe_values[w].gradient(n, q_point)[k]*JxW;
                    cell_matrix(n,m) += dde_dat.d2W_dF_dgrad_w[i][j][k]
                                         *fe_values[displacements].gradient(n, q_point)[i][j]
                                         *fe_values[w].gradient(m, q_point)[k]*JxW;
                    for (unsigned int l = 0; l<DIM; l ++)
                    {
                      cell_matrix(n,m) += dde_dat.d2W_d2F[i][j][k][l]
                                          *fe_values[displacements].gradient(n, q_point)[i][j]
                                          *fe_values[displacements].gradient(m, q_point)[k][l]*JxW;


                    }
                  }
               }
          }
      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        for (unsigned int m=0; m<dofs_per_cell; ++m)
        {
          system_matrix.add (local_dof_indices[n],
                             local_dof_indices[m],
                             cell_matrix(n,m));
        }
    }

    time(&end);

    double time_taken = double(end - start);
//    std::cout << "    Time taken by asseble system : " << time_taken << " seconds." << std::endl;
  }


  void PlateProblem::apply_boundaries_to_rhs(Vector<double> *rhs, std::vector<bool> *homogenous_dirichlet_dofs)
  {
    for (unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {
      if ((*homogenous_dirichlet_dofs)[i] == true)
      {
        (*rhs)[i] = 0.0;
      }
    }
  }

  void PlateProblem::apply_boundaries_and_constraints_system_matrix(std::vector<bool> *homogenous_dofs)
  {
    constraints.condense (system_matrix);

//    std::map<types::global_dof_index,double> boundary_values;
//
//    std::vector<bool> side1_components = {true, true, true };
//
//    ComponentMask side1_mask(side1_components);
//
//    VectorTools::interpolate_boundary_values (dof_handler,
//                                              1,
//                                              ZeroFunction<DIM, double>(DIM + 1),
//                                              boundary_values,
//                                              side1_mask);
//    VectorTools::interpolate_boundary_values (dof_handler,
//                                              2,
//                                              ZeroFunction<DIM, double>(DIM + 1),
//                                              boundary_values,
//                                              side1_mask);
//
//    MatrixTools::apply_boundary_values (boundary_values,
//                                        system_matrix,
//                                        newton_update,
//                                        system_rhs);

    unsigned int m = system_matrix.m();
    // set values on the diagonal to the first diagonal element,
    // or 1 if it is nonexistent
    // This all follows the dealii built in apply_boundaries closely
    double first_nonzero_diagonal_entry = 1.0;
    for (unsigned int i=0; i<m; ++i)
    {
      if (system_matrix.diag_element(i) != 0.0)
      {
        first_nonzero_diagonal_entry = fabs(system_matrix.diag_element(i));
        break;
      }
    }
    // now march through matrix, zeroing out rows and columns.
    // If there is a current value on the diagonal of the constrained
    // boundary dof, don't touch it. If there is not one, then we can
    // just set it equal to the first nonzero entry we just found
    for (unsigned int row = 0; row < m; ++row)
    {

      const typename SparseMatrix<double>::iterator end_row = system_matrix.end(row);
      for (typename SparseMatrix<double>::iterator entry = system_matrix.begin(row);
                    entry != end_row; ++entry)
      {
        if(((*homogenous_dofs)[row] == true || (*homogenous_dofs)[entry->column()] == true)
            && (row != entry->column()))
        {
          entry->value() = 0.0;
        }
        else if((*homogenous_dofs)[row] == true
          && (row == entry->column()))
        {
          entry->value() = first_nonzero_diagonal_entry;
        }

      }

    }

  }

  void PlateProblem::solve ()
  {
    // direct solver for the system
    if(dof_handler.n_dofs() < 10000)
    {
      SparseDirectUMFPACK  A_direct;
      A_direct.initialize(system_matrix);
      A_direct.vmult (newton_update, system_rhs);
    }
    else
    {
      SolverControl solver_control(dof_handler.n_dofs(), 1e-11);
      SolverCG<> solver(solver_control);
      solver.solve(system_matrix, newton_update, system_rhs, PreconditionIdentity());
    }

    constraints.distribute (newton_update);

  }


  void PlateProblem::output_results (const unsigned int cycle) const
  {

    std::vector<std::string> solution_names;
    switch (DIM)
      {
      case 1:
        solution_names.push_back ("displacement");
        break;
      case 2:
        solution_names.push_back ("x1_displacement");
        solution_names.push_back ("x2_displacement");
        break;
      case 3:
        solution_names.push_back ("x1_displacement");
        solution_names.push_back ("x2_displacement");
        solution_names.push_back ("x3_displacement");
        break;
      default:
        Assert (false, ExcNotImplemented());
        break;
      }


    std::vector<std::string> solutionNames(DIM, "u");
    solution_names.emplace_back("w");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(DIM,
                     DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    // output the total displacements. this requires adding in the uniform solution on top of the displacements

    std::string filename0(output_directory);
    filename0 += "/lagrangian_solution";

    // see if the directory exists...
    struct stat st;
    if (stat(filename0.c_str(), &st) == -1)
      mkdir(filename0.c_str(), 0700);

    filename0 += "/lagrangian_solution-";
    filename0 += std::to_string(cycle);
    filename0 += ".vtk";
    std::ofstream output_lagrangian_solution (filename0.c_str());

    DataOut<DIM> data_out_lagrangian;

    data_out_lagrangian.add_data_vector(dof_handler,
                                       present_solution,
                                       solution_names,
                                       interpretation);

    data_out_lagrangian.build_patches ();
    data_out_lagrangian.write_vtk (output_lagrangian_solution);


    // Now output the w on deformed meshes

    std::string filename1(output_directory);
    filename1 += "/eulerian_solution_u";

    // see if the directory exists...
    if (stat(filename1.c_str(), &st) == -1)
      mkdir(filename1.c_str(), 0700);

    filename1 += "/eulerian_solution_u1-";
    filename1 += std::to_string(cycle);
    filename1 += ".vtk";
    std::ofstream output_eulerian_u1 (filename1.c_str());


    DataOut<DIM> data_out_output_eulerian_u1;

    data_out_output_eulerian_u1.add_data_vector(dof_handler,
                                       present_solution,
                                       solution_names,
                                       interpretation);


    MappingQEulerian<DIM> q_mapping(1,  dof_handler, present_solution);
    data_out_output_eulerian_u1.build_patches(q_mapping, 1);
    data_out_output_eulerian_u1.write_vtk (output_eulerian_u1);

  }



  void PlateProblem::read_input_file(char* filename)
  {
    FILE* fid;
    int endOfFileFlag;
    char nextLine[MAXLINE];

    int valuesWritten;
    bool fileReadErrorFlag = false;

    grid_dimensions.resize(DIM);
    domain_dimensions.resize(DIM);

    fid = std::fopen(filename, "r");
    if (fid == NULL)
    {
      std::cout << "Unable to open file \"" << filename  << "\"" <<  std::endl;
      fileReadErrorFlag = true;
    }
    else
    {

      // Read in the output name
      char directory_name[MAXLINE];
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%s", directory_name);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      sprintf(output_directory, "output/");
      strcat(output_directory, directory_name);


//      if(objective_type == 0) load_val = 1.0;


      // Read in the grid dimensions
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u %u", &grid_dimensions[0], &grid_dimensions[1]);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // Read in the domain dimensions
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg", &domain_dimensions[0], &domain_dimensions[1]);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in some good parameters for now
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg", &lx, &ly, &delta);
      if(valuesWritten != 3)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in some good  moduii parameters for now
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg", &mu, &lambda, &B, &eta);
      if(valuesWritten != 4)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the number of guass points in the x and y direction
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u  %u", &qx, &qy);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the newton raphson tolerance and the max number of steps
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg  %u", &newton_tol, &newton_max_iter);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }


      fileClose:
      {
        fclose(fid);
      }
    }

    if (fileReadErrorFlag)
    {
      // default parameter values
      std::cout << "Error reading input file, Exiting.\n" << std::endl;
      exit(1);
    }
    else
      std::cout << "Input file successfully read" << std::endl;

    // make the output directory
    struct stat st;
    if (stat("./output", &st) == -1)
       mkdir("./output", 0700);

    if (stat(output_directory, &st) == -1)
      mkdir(output_directory, 0700);

  }

  void PlateProblem::getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                          int const maxSize, int* const endOfFileFlag)
  {
    *endOfFileFlag = 0;
    do
    {
      if(fgets(nextLinePtr, maxSize, filePtr) == NULL)
      {
        *endOfFileFlag = 1;
        break;
      }
      while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
             (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
      {
        nextLinePtr = (nextLinePtr + 1);
      }
    }
    while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
  }

  void PlateProblem::save_current_state(unsigned int indx)
  {
    // create the output directory if it doesnt exist

    char saved_state_dir[MAXLINE];
    strcpy(saved_state_dir, output_directory);
    strcat(saved_state_dir, "/saved_state");

    // see if the directory exists
    struct stat st;
    if (stat(saved_state_dir, &st) == -1)
         mkdir(saved_state_dir, 0700);


    char index_char[32];
    sprintf(index_char, "%u", indx);

    // Triangulation
    char triag_file[MAXLINE];
    strcpy(triag_file, saved_state_dir);
    strcat(triag_file, "/triag_");
    strcat(triag_file, index_char);
    strcat(triag_file, ".dat");
    std::ofstream triag_out (triag_file);
    boost::archive::text_oarchive triag_ar(triag_out);
    triangulation.save(triag_ar, 1);

    // dof handler
    char dof_file[MAXLINE];
    strcpy(dof_file, saved_state_dir);
    strcat(dof_file, "/dof_");
    strcat(dof_file, index_char);
    strcat(dof_file, ".dat");
    std::ofstream dof_out (dof_file);
    boost::archive::text_oarchive dof_ar(dof_out);
    dof_handler.save(dof_ar, 1);

    // Present Solution
    char present_solution_file[MAXLINE];
    strcpy(present_solution_file, saved_state_dir);
    strcat(present_solution_file, "/present_solution_");
    strcat(present_solution_file, index_char);
    strcat(present_solution_file, ".dat");
    std::ofstream solution_out (present_solution_file);
    boost::archive::text_oarchive solution_ar(solution_out);
    present_solution.save(solution_ar, 1);

  }

  void PlateProblem::load_state(unsigned int indx)
  {
    // create the output directory

    char input_dir_path[MAXLINE];
    strcpy(input_dir_path, output_directory);
    strcat(input_dir_path, "/saved_state");
    struct stat st;
    if (stat(input_dir_path, &st) == -1)
    {
      std::cout << "Could not find the directory : " << input_dir_path << "\nExiting." <<std::endl;
      exit(-1);
    }

    char index_char[32];
    sprintf(index_char, "%u", indx);

    // Triangulation
    char triag_file[MAXLINE];
    strcpy(triag_file, input_dir_path);
    strcat(triag_file, "/triag_");
    strcat(triag_file, index_char);
    strcat(triag_file, ".dat");
    std::ifstream triag_in (triag_file);
    boost::archive::text_iarchive triag_ar(triag_in);
    triangulation.load(triag_ar, 1);

    // df_handler
    dof_handler.distribute_dofs(fe);
    char dof_file[MAXLINE];
    strcpy(dof_file, input_dir_path);
    strcat(dof_file, "/dof_");
    strcat(dof_file, index_char);
    strcat(dof_file, ".dat");
    std::ifstream dof_in (dof_file);
    boost::archive::text_iarchive dof_ar(dof_in);
    dof_handler.load(dof_ar, 1);

    // Present Solution
    char present_solution_file[MAXLINE];
    strcpy(present_solution_file, input_dir_path);
    strcat(present_solution_file, "/present_solution_");
    strcat(present_solution_file, index_char);
    strcat(present_solution_file, ".dat");
    std::ifstream solution_in (present_solution_file);
    boost::archive::text_iarchive solution_ar(solution_in);
    present_solution.load(solution_ar, 1);
    evaluation_point = present_solution;


  }


  void PlateProblem::renumber_boundary_ids()
  {

    // renumber boundary ids because they have problems being saved for nonuniform mesh.
    typename Triangulation<DIM>::active_cell_iterator cell =
     triangulation.begin_active(), endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (unsigned int f = 0; f < GeometryInfo<DIM>::faces_per_cell; ++f)
      {

        const Point<DIM> face_center = cell->face(f)->center();
        if (fabs(face_center[0]) < 1e-4)  // left
        {
          cell->face(f)->set_boundary_id (1);
        }
        else if (fabs(face_center[0] - domain_dimensions[0]) < 1e-4)
        {// right
          cell->face(f)->set_boundary_id (2);
        }
        else if (fabs(face_center[1] + domain_dimensions[1]/2.0) < 1e-6) //bottom
        {
          cell->face(f)->set_boundary_id (3);
        }
        else if (fabs(face_center[1] - domain_dimensions[1]/2.0) < 1e-6) //top
        {
         cell->face(f)->set_boundary_id (4);

        }

      }
  }

  void PlateProblem::add_small_perturbations(double amplitude, bool firstTime)
  {
    // This is just to make the solution vector non-zero so it doesn't start at the right answer.

    if(firstTime)
      std::srand(5); //some seed

    for(unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {

      double random_val =  0.5*(2.0*(std::rand()/double(RAND_MAX)) - 1.0)*amplitude;
      present_solution[i] += random_val;

    }

    // constraints.distribute(present_solution);

    evaluation_point = present_solution;

  }

  void PlateProblem::rhs_numerical_deriv(double delta)
  {
    add_small_perturbations(0.000174, true);
    assemble_system_energy();

//    stuffOn = true;
    assemble_system_rhs();
    output_results(100000);
    double e0 = system_energy;
    Vector<double> numer_rhs(dof_handler.n_dofs());
    Vector<double> numer_rhs_diff(dof_handler.n_dofs());

    std::cout << "Total Energy : " << e0 << std::endl;
    std::cout << "rhs difference : " << std::endl;
    for(unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {
      evaluation_point[i] += delta;
      assemble_system_energy();
      numer_rhs[i] = (system_energy - e0)/(delta);
      evaluation_point[i] -= delta;
    }

//    constraints.condense(system_rhs);
//    constraints.condense(numer_rhs);
    for(unsigned int i = 0 ; i < dof_handler.n_dofs(); i++)
    {
//      if(fabs(system_rhs[i]) >= 1e-6 )
//        numer_rhs_diff[i] = (numer_rhs[i] + system_rhs[i])/system_rhs[i];
//      else
        numer_rhs_diff[i] = (numer_rhs[i] + system_rhs[i]);


      std::cout << "   ";
      std::cout << std::setw(10) << numer_rhs[i] << "   ";
      std::cout << std::setw(10) << -system_rhs[i] << "   ";
      std::cout << std::setw(10) << numer_rhs_diff[i] << std::endl;

    }

    std::cout << std::endl;
//    constraints.condense(numer_rhs_diff);

    std::cout <<"\n\n\nL2 Norm of the difference : " << numer_rhs_diff.l2_norm() << std::endl;

  }

  void PlateProblem::compute_and_compare_second_numer_deriv(double epsilon)
  {
    unsigned int numberDofs = dof_handler.n_dofs();
//    present_solution = 0.0;

    add_small_perturbations(0.000174, true);

    evaluation_point = present_solution;
    assemble_system_rhs();
    assemble_system_matrix();
    system_rhs *= -1.0;
    Vector<double> residual = system_rhs;
    std::vector<std::vector<double>> numericalStiffness(numberDofs);
    for(unsigned int i = 0; i < numberDofs; i++)
      numericalStiffness[i].resize(numberDofs);

    for(unsigned int j = 0; j < numberDofs; j++)
    {
      evaluation_point[j] += epsilon;
      assemble_system_rhs();
      evaluation_point[j] -= epsilon;

      system_rhs *= -1.0;
      for(unsigned int i = 0; i < numberDofs; i++)
      {
        numericalStiffness[i][j] = (system_rhs[i] - residual[i])/epsilon;
      }
    }

    std::ofstream out("secondDeriv.dat");
    out << "" << std::endl;
    out << "Stiffness Matrix" << std::endl;


    for(unsigned int i = 0; i < numberDofs; i++)
    {
      for(unsigned int j = 0; j < numberDofs; j++)
      {
        out << system_matrix.el(i,j) << " ";
      }
      out << std::endl;
    }

    out << "\n\n\nNumerical" << std::endl;


    for(unsigned int i = 0; i < numberDofs; i++)
    {
      for(unsigned int j = 0; j < numberDofs; j++)
      {
        out << numericalStiffness[i][j]  << " ";
      }
      out << std::endl;
    }

    out << "\n\n\nDifference" << std::endl;
    eta = eta_;
    double diffNorm = 0;
    for(unsigned int i = 0; i < numberDofs; i++)
    {
      for(unsigned int j = 0; j < numberDofs; j++)
      {
        double diff;
        diff = (system_matrix.el(i,j) - numericalStiffness[i][j]);
        diffNorm += diff*diff;
        out << diff << " ";
      }
      out << std::endl;
    }
    diffNorm = sqrt(diffNorm);

    out << std::endl;
    out << "Norm of the difference: " << diffNorm << std::endl;

    out.close();


  }



}

#endif // EFFECTIVE_PLATE_CC_
