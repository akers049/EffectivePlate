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
#include <deal.II/lac/slepc_solver.h>



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
         FE_Q<DIM>(1), 1, FE_Q<DIM>(1), 1, FE_Q<DIM>(1), 1)
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
    corner1(1) =  0;
    corner2(0) =  domain_dimensions[0];
    corner2(1) =  domain_dimensions[1];

    GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, true);

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


    double number_dofs = dof_handler.n_dofs();

    PE.set_moduli(mu, lambda, B, eta);

    eval_wx.set_param_values(lx, ly, delta);
    eval_wy.set_param_values(lx, ly, delta);

    present_solution.reinit (dof_handler.n_dofs());


    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    setup_system_constraints();
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);


    sparsity_pattern.copy_from (dsp);


//    GridTools::distort_random(0.4, triangulation, true);


    system_matrix.reinit (sparsity_pattern);

    system_matrix_petsc.reinit (sparsity_pattern);


    evaluation_point = present_solution;

    homo_dofs.resize(dof_handler.n_dofs(), false);
    load_dofs.resize(dof_handler.n_dofs(), false);


    std::vector<Point<DIM>> support_points(dof_handler.n_dofs());
    MappingQ1<DIM> mapping;
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

    std::vector<bool> x1_components = {true, false, false, false, false};
    ComponentMask x1_mask(x1_components);
    std::vector<bool> is_x1_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);

    std::vector<bool> x2_components = {false, true, false, false, false};
    ComponentMask x2_mask(x2_components);
    std::vector<bool> is_x2_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);

    std::vector<bool> w_components = {false, false, true, false , false};
    ComponentMask w_mask(w_components);
    std::vector<bool> is_w_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, w_mask, is_w_comp);

    std::vector<bool> v_components = {false, false, false, true , false};
    ComponentMask v_mask(v_components);
    std::vector<bool> is_v_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, v_mask, is_v_comp);

    std::vector<bool> lam_components = {false, false, false, false , true};
    ComponentMask lam_mask(lam_components);
    std::vector<bool> is_lam_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, lam_mask, is_lam_comp);

    std::vector<bool> all_components = {true, true, true, true , true};
    ComponentMask all_mask(all_components);
    std::vector<bool> boundary_dof (number_dofs, false);


    DoFTools::extract_boundary_dofs(dof_handler,
        all_mask,
        boundary_dof);


     for(unsigned int  i = 0; i < number_dofs; i ++)
     {
       if(fabs( support_points[i](1) - domain_dimensions[1]/2.0) < 1.0e-6)
       {
         if(fabs(support_points[i](0)) < 1.0e-6)
         {
           if(is_x1_comp[i] || is_x2_comp[i] || is_w_comp[i])
             homo_dofs[i] = true;
         }

         if(fabs(support_points[i](0) - domain_dimensions[0]) < 1.0e-6)
         {
           if(is_x1_comp[i] == true)
             load_dofs[i] = true;

           if(is_x1_comp[i] || is_x2_comp[i] || is_w_comp[i])
             homo_dofs[i] = true;
         }

       }

       // do the extra points we are constraining
       if(fabs( support_points[i](0) - domain_dimensions[0]/2.0) < 1.0e-6 &&
           (fabs(support_points[i](1)) < 1.0e-6  || fabs(support_points[i](1) - domain_dimensions[1]) < 1.0e-6)
               && is_w_comp[i] == true)
       {
         homo_dofs[i] = true;
       }

       if(boundary_dof[i] && (is_v_comp[i] || is_lam_comp[i]))
         homo_dofs[i] = true;

     }

//    assemble_system_matrix();
//    setup_system_constraints();
  }


  void PlateProblem::setup_system_constraints ()
  {

    const unsigned int   number_dofs = dof_handler.n_dofs();

    // now do constraints that the average w displacement is zero
//    std::vector<bool> w_components = {false, false, true};
//    ComponentMask w_mask(w_components);
//
//    std::vector<bool> boundary_dof_w (number_dofs, false);
//
//
//    DoFTools::extract_boundary_dofs(dof_handler,
//        w_mask,
//        boundary_dof_w);
//
//    unsigned int first_boundary_dof = 0;
//    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
//    {
//      if ((boundary_dof_w[i] == true))
//      {
//        first_boundary_dof = i;
//        break;
//      }
//    }
//    constraints.add_line (first_boundary_dof);
//    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
//    {
//      if (i == first_boundary_dof)
//        continue;
//
//      if(boundary_dof_w[i] == true)
//        constraints.add_entry (first_boundary_dof, i, -1);
//
//    }

//    std::vector<Point<DIM>> support_points(dof_handler.n_dofs());
//    MappingQ1<DIM> mapping;
//    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);
//
//
//    std::vector<bool> x1_components = {true, false, false, false, false};
//    ComponentMask x1_mask(x1_components);
//    std::vector<bool> is_x1_comp(number_dofs);
//    DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);
//
//    std::vector<bool> x2_components = {false, true, false, false, false};
//    ComponentMask x2_mask(x2_components);
//    std::vector<bool> is_x2_comp(number_dofs);
//    DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);
//
//    std::vector<bool> w_components = {false, false, true, false , false};
//    ComponentMask w_mask(w_components);
//    std::vector<bool> is_w_comp(number_dofs);
//    DoFTools::extract_dofs(dof_handler, w_mask, is_w_comp);
//
//    std::vector<bool> v_components = {false, false, false, true, false};
//    ComponentMask v_mask(v_components);
//    std::vector<bool> is_v_comp(number_dofs);
//    DoFTools::extract_dofs(dof_handler, v_mask, is_v_comp);
//
//    std::vector<bool> lam_components = {false, false, false, false, true};
//    ComponentMask lam_mask(lam_components);
//    std::vector<bool> is_lam_comp(number_dofs);
//    DoFTools::extract_dofs(dof_handler, lam_mask, is_lam_comp);
//
//    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
//    {
//      if(is_v_comp[i] == true)
//      {
//        for (unsigned int k = 0; k < dof_handler.n_dofs(); ++k)
//        {
//          if( is_lam_comp[k] == true &&
//              (support_points[i](0) - support_points[k](0))*(support_points[i](0) - support_points[k](0)) < 1.0e-6 &&
//              (support_points[i](1) - support_points[k](1))*(support_points[i](1) - support_points[k](1)) < 1.0e-6)
//            {
//              constraints.add_line (i);
//              constraints.add_entry (i, k, B);
//              break;
//            }
//        }
//      }
//    }
//
//    constraints.close ();

    std::vector<bool> x1_components = {true, false, false, false, false};
    ComponentMask x1_mask(x1_components);
    std::vector<bool> is_x1_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);

    std::vector<bool> x2_components = {false, true, false, false, false};
    ComponentMask x2_mask(x2_components);
    std::vector<bool> is_x2_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);

    std::vector<bool> w_components = {false, false, true, false , false};
    ComponentMask w_mask(w_components);
    std::vector<bool> is_w_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, w_mask, is_w_comp);

    std::vector<bool> v_components = {false, false, false, true, false};
    ComponentMask v_mask(v_components);
    std::vector<bool> is_v_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, v_mask, is_v_comp);

    std::vector<bool> lam_components = {false, false, false, false, true};
    ComponentMask lam_mask(lam_components);
    std::vector<bool> is_lam_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, lam_mask, is_lam_comp);

    constraints_eigen.clear();

    std::vector<Point<DIM>> support_points(dof_handler.n_dofs());
    MappingQ1<DIM> mapping;
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
//      if(is_lam_comp[i] == true)
//      {
//        for (unsigned int k = 0; k < dof_handler.n_dofs(); ++k)
//        {
//          if( is_v_comp[k] == true &&
//              (support_points[i](0) - support_points[k](0))*(support_points[i](0) - support_points[k](0)) < 1.0e-6 &&
//              (support_points[i](1) - support_points[k](1))*(support_points[i](1) - support_points[k](1)) < 1.0e-6)
//          {
////            constraints_eigen.add_line (i);
////            constraints_eigen.add_entry (i, k, -B);
////            constraints.add_line (i);
////            constraints.add_entry (i, k, -B);
//            break;
//          }
//        }

//      if(is_w_comp[i] == true && support_points[i](1) > 1.0e-6 + domain_dimensions[1]/2.0)
//      {
//        for (unsigned int k = 0; k < dof_handler.n_dofs(); ++k)
//        {
//          if( is_w_comp[k] == true &&
//              (support_points[i](0) - support_points[k](0))*(support_points[i](0) - support_points[k](0)) < 1.0e-6 &&
//              ( ((support_points[i](1) - domain_dimensions[1]/2.0) - (domain_dimensions[1]/2.0 - support_points[k](1)))
//                 *((support_points[i](1) - domain_dimensions[1]/2.0) - (domain_dimensions[1]/2.0 - support_points[k](1)))
//                   < 1.0e-6))
//          {
//                        constraints_eigen.add_line (i);
//                        constraints_eigen.add_entry (i, k, 1.0);
//                        constraints.add_line (i);
//                        constraints.add_entry (i, k, 1.0);
//            break;
//          }
//        }
//      }
    }



//    for (unsigned int row = 0; row < system_matrix.m(); ++row)
//    {
//      if(is_lam_comp[row] == true)
//      {
//        double first_entry = 0.0;
//        unsigned int constraint_dof = 0;
////        bool firstFlag = false;
//
//        first_entry = system_matrix.el(row, row - 2);
//        constraint_dof = row - 2;
//        constraints_eigen.add_line(constraint_dof);
//        const typename SparseMatrix<double>::const_iterator end_row = system_matrix.end(row);
//        for (typename SparseMatrix<double>::const_iterator entry = system_matrix.begin(row);
//                               entry != end_row; ++entry)
//         {
//
////           if((is_w_comp[entry->column()] || is_v_comp[entry->column()]) && fabs(entry->value()) > 1.0e-12 )
////           {
////             if(firstFlag == false)
////             {
////               first_entry = entry->value();
////               constraint_dof = entry->column();
////               constraints_eigen.add_line(constraint_dof);
////               firstFlag = true;
////             }
//          if((is_w_comp[entry->column()] || is_v_comp[entry->column()])) // && fabs(entry->value()) > 1.0e-12 )
//          {
//             if(entry->column() == row - 2)
//               continue;
//             else
//             {
//               constraints_eigen.add_entry (constraint_dof, entry->column(), -entry->value()/first_entry);
//             }
//           }
//         }
//      }
//    }
    constraints_eigen.close ();

    constraints.close();

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

    double E0 = 0.0;
    assemble_system_energy();
    E0 = system_energy;
    previous_solution = present_solution;
    Vector<double> dpred = present_solution;

    bool updateFlag = true;

    for(unsigned int i = 0; i < load_steps; i ++)
    {

      current_load_value = du*i;
      set_displacement(current_load_value);

      previous_solution = present_solution;

//      present_solution += dpred;

      std::cout << "  Iteration : " << i+1 << std::endl;
      newton_iterate();

      dpred = present_solution;
      dpred -= previous_solution;

      unsigned int num_neg_eigs;

//      output_matrix_csv();
//      exit(-1);
//      if(i == 0)
      if(updateFlag)
      {
        num_neg_eigs = get_system_eigenvalues(-1);
        std::cout << "      Number of negative eigenvalues is : " << num_neg_eigs << std::endl;
      }
      if(num_neg_eigs > 0)
        updateFlag = false;

      output_results(i);
//      output_matrix_csv();
//      if(num_neg_eigs > 0)
//        break;


//      present_solution *= 2.0;
//      present_solution -= previous_solution;
//      exit(-1);
    }

//    evaluation_point = present_solution;
//    assemble_system_energy();
//    double E1 = system_energy;
//
//    std::cout << "Old Energy : " << E0 << std::endl;
//    std::cout << "New Energy : " << E1 << std::endl;
//
//    newton_iterate();
//    output_results(9999);



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
    for(lineSearchLength = 1.0; lineSearchLength > 1.0e-3; lineSearchLength *=0.5)
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
    if(lineSearchLength < 1e-3)
    {
      std::cout << "THERE WAS A PROBLEM!!" << std::endl;
    }
    present_solution = evaluation_point;

    return current_residual;
  }

  unsigned int PlateProblem::get_system_eigenvalues(const int cycle)
  {
    // get the system's eigenvalues. I really don't need to have it take in lambda_eval
    // and could just have it use the present_lambda, but its fine. the cycle is the
    // the output file will have appended on its name. -1 for no output. Outputs the number
    // of negative eigenvalues in the range -10 < lambda < 0.3

    evaluation_point = present_solution;
    assemble_system_matrix();
//    apply_boundaries_and_eigen_constraints_system_matrix(&homo_dofs);
    apply_boundaries_and_eigen_constraints_system_matrix(&homo_dofs);

    LAPACKFullMatrix<double> system_matrix_full;
    system_matrix_full.copy_from(system_matrix);

    Vector<double> eigenvalues;
    FullMatrix<double> eigenvectors;

    double tol = 1.0e-9;

    system_matrix_full.compute_eigenvalues_symmetric(-10, 0.3, tol, eigenvalues, eigenvectors);

    unsigned int num_neg_eigs = 0;
    for(unsigned int i = 0; i < eigenvalues.size(); i ++)
    {
//      std::cout << eigenvalues[i] << std::endl;

      if(eigenvalues[i] < -1.0e-6)
        num_neg_eigs ++;
    }

    double number_dofs = dof_handler.n_dofs();
    std::vector<bool> lam_components = {false, false, false, false, true};
    ComponentMask lam_mask(lam_components);
    std::vector<bool> is_lam_comp(number_dofs);
    DoFTools::extract_dofs(dof_handler, lam_mask, is_lam_comp);
    unsigned int total_found = 0;
    unsigned int need_to_find = 0;
    for(unsigned int i = 0; i < dof_handler.n_dofs(); i ++)
      if(is_lam_comp[i])
        need_to_find++;
    if(num_neg_eigs > 0)
    {
      Vector<double> dest(dof_handler.n_dofs());
      Vector<double> next_eig(dof_handler.n_dofs());

      for(unsigned int k = 0; k < num_neg_eigs; k ++)
      {
        double next_tot = 0.0;
        for(unsigned int i = 0; i < dof_handler.n_dofs(); i ++)
        {
          next_eig[i] = eigenvectors[i][k];
        }
        apply_boundaries_to_rhs(&next_eig, &homo_dofs);
        constraints_eigen.distribute(next_eig);

        system_matrix.vmult(dest, next_eig);

        for(unsigned int i = 0; i < number_dofs; i ++)
        {
          if(lam_components[i] == true)
          {
            next_tot += dest[i]*dest[i];

          }
        }
//        std::cout << "Eigenvalue : " << k << "  Value : " << eigenvalues[k] << "Total : " << sqrt(next_tot) << std::endl;
        if(sqrt(next_tot) < 1.0e-6)
        {
          total_found ++;
          present_solution.add(0.1, next_eig);
          std::cout << "  Negaitve Eigenvalue : " << eigenvalues[k] << std::endl;
          break;
//          std::cout << "     <<<<<<< Found one >>>>>>>" << std::endl;
        }
      }
    }

    num_neg_eigs = total_found;
//    std::cout << "\n\n TOTAL FOUND : " << total_found << "   NEEDE TO FIND : " << need_to_find << std::endl;
//    std::cout << "COWS : " << present_solution[538] << " " << present_solution[539] << std::endl;
//    std::cout << "COWS : " << present_solution[538] << " " << present_solution[539] << std::endl;


//    std::vector<PETScWrappers::MPI::Vector> eigenfunctions;
//    std::vector<double>                     eigenvalues;
//    IndexSet eigenfunction_index_set(dof_handler.n_dofs());
//    eigenfunction_index_set.add_range(0, dof_handler.n_dofs());
//    eigenfunctions.resize (219);
//    for (unsigned int i=0; i<eigenfunctions.size (); ++i)
//      eigenfunctions[i].reinit (eigenfunction_index_set, MPI_COMM_WORLD);
//
//    eigenvalues.resize (eigenfunctions.size());
//    system_matrix_petsc = 0.0;
//    for (unsigned int row = 0; row < system_matrix.m(); ++row)
//    {
//      const typename SparseMatrix<double>::const_iterator end_row = system_matrix.end(row);
//      for (typename SparseMatrix<double>::const_iterator entry = system_matrix.begin(row);
//                             entry != end_row; ++entry)
//       {
//         if(fabs(entry->value()) > 1e-10 )
//           system_matrix_petsc.set(row, entry->column(),entry->value());
//       }
//    }
//
//    system_matrix_petsc.compress(VectorOperation::insert);
//    SolverControl solver_control (dof_handler.n_dofs()*100, 1e-5);
////    SLEPcWrappers::SolverKrylovSchur eigensolver (solver_control);
//    SLEPcWrappers::SolverLAPACK eigensolver (solver_control);
//
//    // Before we actually solve for the eigenfunctions and -values, we have to
//    // also select which set of eigenvalues to solve for. Lets select those
//    // eigenvalues and corresponding eigenfunctions with the smallest real
//    // part (in fact, the problem we solve here is symmetric and so the
//    // eigenvalues are purely real). After that, we can actually let SLEPc do
//    // its work:
////    eigensolver.set_target_eigenvalue(-10.0);
//    eigensolver.set_which_eigenpairs (EPS_SMALLEST_REAL);
//
//    eigensolver.set_problem_type (EPS_HEP);
//
//    eigensolver.solve (system_matrix_petsc,
//                       eigenvalues, eigenfunctions,
//                       eigenfunctions.size());
//
//
////    for(unsigned int i = 0; i < dof_handler.n_dofs(); i ++)
////    {
////      present_solution[i] = 0.1*eigenfunctions[0][i];
////    }
//
//    unsigned int num_neg_eigs = 0;
//    if (cycle != -1)
//    {
//      std::string filename(output_directory);
//          filename += "/eigenvalues";
//
//      // see if the directory exists
//      struct stat st;
//      if (stat(filename.c_str(), &st) == -1)
//        mkdir(filename.c_str(), 0700);
//
//      filename += "/eigenvalues-";
//      filename += std::to_string(cycle);
//
//      std::ofstream outputFile;
//      outputFile.open(filename.c_str());
//
//      //outputFile << "# eigenvalues of the system matrix" << std::endl;
//
//      for (unsigned int i = 0 ; i < eigenvalues.size(); i ++)
//      {
//        double nextEigenVal = eigenvalues[i];
//
//        outputFile << std::setprecision(15) << nextEigenVal << std::endl;
//
//        if (nextEigenVal < 0.0)
//        {
//          num_neg_eigs ++;
//        }
//
//      }
//
//      // outputFile << "\nIs positive definite : " << num_neg_eigs << std::endl;
//      outputFile.close();
//    }
//    else
//    {
//      for (unsigned int i = 0; i < eigenvalues.size(); i++)
//      {
//        std::cout << eigenvalues[i] << std::endl;
//        if (eigenvalues[i] < 0.0)
//        {
//          num_neg_eigs ++;
////          break;
//        }
//      }
//    }

    return num_neg_eigs;
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
    std::vector<Tensor<1, DIM>> lam_gradients(n_q_points);
    std::vector<double> v_values(n_q_points);
    std::vector<double> lam_values(n_q_points);

//    std::vector<Tensor<2, DIM>> w_hessians(n_q_points);


    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector displacements (0);
    const FEValuesExtractors::Scalar w (DIM);
    const FEValuesExtractors::Scalar v (DIM+1);
    const FEValuesExtractors::Scalar lam (DIM+2);

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
      fe_values[lam].get_function_gradients(evaluation_point, lam_gradients);
      fe_values[v].get_function_values(evaluation_point, v_values);
      fe_values[lam].get_function_values(evaluation_point, lam_values);

//      fe_values[w].get_function_hessians(evaluation_point, w_hessians);

      eval_wx.value_list (fe_values.get_quadrature_points(), wx);
      eval_wy.value_list (fe_values.get_quadrature_points(), wy);


      unsigned int cell_index = cell->active_cell_index();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        PE.set_params(wx[q_point], wy[q_point], lx, ly, delta);
        get_F(displacement_gradients[q_point], F);
        double lap_w = v_values[q_point]; //trace(w_hessians[q_point]);
        system_energy += PE.get_E(F, w_gradients[q_point], lap_w)*fe_values.JxW(q_point);
        system_energy += (lam_gradients[q_point]*w_gradients[q_point] + lam_values[q_point]*v_values[q_point])*fe_values.JxW(q_point);

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
    std::vector<Tensor<1, DIM>> lam_gradients(n_q_points);
    std::vector<double> v_values(n_q_points);
    std::vector<double> lam_values(n_q_points);

//    std::vector<Tensor<2, DIM>> w_hessians(n_q_points);

    std::vector<double> wx(n_q_points);
    std::vector<double> wy(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector displacements (0);
    const FEValuesExtractors::Scalar w (DIM);
    const FEValuesExtractors::Scalar v (DIM+1);
    const FEValuesExtractors::Scalar lam (DIM+2);


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
      fe_values[lam].get_function_gradients(evaluation_point, lam_gradients);
      fe_values[v].get_function_values(evaluation_point, v_values);
      fe_values[lam].get_function_values(evaluation_point, lam_values);

//      fe_values[w].get_function_hessians(evaluation_point, w_hessians);

      eval_wx.value_list (fe_values.get_quadrature_points(), wx);
      eval_wy.value_list (fe_values.get_quadrature_points(), wy);

      unsigned int cell_index = cell->active_cell_index();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {


        PE.set_params(wx[q_point], wy[q_point], lx, ly, delta);
        get_F(displacement_gradients[q_point], F);
        double lap_w = 0.0; // trace(w_hessians[q_point]);

        double JxW = fe_values.JxW(q_point);

//        if(cell_index == 69)
//          PE.numerical_deriv_internal(F, w_gradients[q_point], lap_w);

        PE.get_DE(F, w_gradients[q_point], lap_w, de_dat);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          cell_rhs(n) -= de_dat.dW_dgrad_w*fe_values[w].gradient(n, q_point)*JxW;
          cell_rhs(n) -= lam_gradients[q_point]*fe_values[w].gradient(n, q_point)*JxW;

          cell_rhs(n) -= (B*v_values[q_point] + lam_values[q_point])*fe_values[v].value(n, q_point)*JxW;
          cell_rhs(n) -= w_gradients[q_point]*fe_values[lam].gradient(n, q_point)*JxW;
          cell_rhs(n) -= v_values[q_point]*fe_values[lam].value(n, q_point)*JxW;

//          cell_rhs(n) -= de_dat.dW_dlap_w*trace(fe_values[w].hessian(n, q_point))*JxW;

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
//    std::vector<Tensor<2, DIM>> w_hessians(n_q_points);


    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector displacements (0);
    const FEValuesExtractors::Scalar w (DIM);
    const FEValuesExtractors::Scalar v (DIM+1);
    const FEValuesExtractors::Scalar lam (DIM+2);

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
//      fe_values[w].get_function_hessians(evaluation_point, w_hessians);

      eval_wx.value_list (fe_values.get_quadrature_points(), wx);
      eval_wy.value_list (fe_values.get_quadrature_points(), wy);

      unsigned int cell_index = cell->active_cell_index();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        PE.set_params(wx[q_point], wy[q_point], lx, ly, delta);
        get_F(displacement_gradients[q_point], F);
        double lap_w = 0.0; //trace(w_hessians[q_point]);

        double JxW = fe_values.JxW(q_point);

        PE.get_DDE(F, w_gradients[q_point], lap_w, dde_dat);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
          for (unsigned int m = 0; m < dofs_per_cell; ++m)
          {
//            cell_matrix(n,m) += dde_dat.d2W_d2lap_w*trace(fe_values[w].hessian(m, q_point))
//                                *trace(fe_values[w].hessian(n, q_point))*JxW;

            cell_matrix(n,m) += B*fe_values[v].value(m, q_point)*fe_values[v].value(n, q_point)*JxW;

            cell_matrix(n,m) += fe_values[lam].value(m, q_point)*fe_values[v].value(n, q_point)*JxW;
            cell_matrix(n,m) += fe_values[lam].value(n, q_point)*fe_values[v].value(m, q_point)*JxW;

            cell_matrix(n,m) += fe_values[lam].gradient(n, q_point)*fe_values[w].gradient(m, q_point)*JxW;
            cell_matrix(n,m) += fe_values[lam].gradient(m, q_point)*fe_values[w].gradient(n, q_point)*JxW;


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

  void PlateProblem::apply_boundaries_and_eigen_constraints_system_matrix(std::vector<bool> *homogenous_dofs)
  {
    constraints_eigen.condense (system_matrix);


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
      SolverControl solver_control(dof_handler.n_dofs(), 1e-8);
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
    solution_names.push_back("w");
    solution_names.push_back("v");
    solution_names.push_back("lam");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(DIM,
                     DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);
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

      // ds and numSteps
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %u", &load_val, &load_steps);
      if(valuesWritten != 2)
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
//    apply_boundaries_to_rhs(&present_solution, &homo_dofs);
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

  void PlateProblem::output_matrix_csv()
  {
    assemble_system_matrix();
    apply_boundaries_and_constraints_system_matrix(&homo_dofs);

    std::ofstream out("system_matrix.csv");
    unsigned int numberDofs = dof_handler.n_dofs();


    for(unsigned int i = 0; i < numberDofs; i++)
    {
      for(unsigned int j = 0; j < numberDofs; j++)
      {
        out << system_matrix.el(i,j);
        if(j != numberDofs - 1)
          out << ", ";
      }
      out << std::endl;
    }

    out.close();


  }



}

#endif // EFFECTIVE_PLATE_CC_
