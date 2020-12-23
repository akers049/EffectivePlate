#ifndef CONSTITUITIVE_H_
#define CONSTITUITIVE_H_

#include <iostream>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/tensor.h>

#define DIM 2

using namespace dealii;

  class DE_data
  {
    public:
    DE_data(){};
    void clear()
    {
      dW_dF = 0.0;
      dW_dgrad_w = 0.0;
      dW_dlap_w = 0.0;
    };

    Tensor<2,DIM> dW_dF;
    Tensor<1,DIM> dW_dgrad_w;
    double dW_dlap_w = 0.0;
  };

  class DDE_data
  {
    public:
    DDE_data(){};
    void clear()
    {
      d2W_d2F = 0.0;
      d2W_dF_dgrad_w = 0.0;
      d2W_d2grad_w = 0.0;
      d2W_d2lap_w = 0.0;
    };

    Tensor<4,DIM> d2W_d2F;
    Tensor<3,DIM> d2W_dF_dgrad_w;
    Tensor<2,DIM> d2W_d2grad_w;
    double d2W_d2lap_w = 0.0;
  };

  class plate_energy
  {
    public:
    plate_energy(){};
    void set_params(double a_, double b_,  double c_);


    double get_E(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w);
    void get_DE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DE_data &de);
    void get_DDE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DDE_data &dde);

    double a = 1.0;
    double b = 1.0;
    double c = 1.0;
  };


  #endif
