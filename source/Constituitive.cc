#ifndef CONSTITUITIVE_CC_
#define CONSTITUITIVE_CC_

#include "Constituitive.h"

#define DIM 2

using namespace dealii;



  void plate_energy::set_params(double a_, double b_,  double c_)
  {
    a = a_;
    b = b_;
    c = c_;
  }


  double plate_energy::get_E(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w)
  {

    return 0.0;
  }

  void plate_energy::get_DE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DE_data &de)
  {

  }

  void plate_energy::get_DDE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DDE_data &dde)
  {

  }






  #endif
