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
    void set_params(double wx_, double wy_,  double lx_, double ly_, double delta_);
    void set_moduli(double mu_, double lambda_, double B_, double eta_)
    {
      mu = mu_;
      lambda = lambda_;
      B = B_;
      eta = eta_;
    };


    double get_E(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w);
    void get_DE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DE_data &de);
    void get_DDE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DDE_data &dde);

    void numerical_deriv_internal(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w);

    private:

    void get_dPhi_dC(Tensor<2,DIM> &C, Tensor<2,DIM> &dPhi_dC);

    double dirac(unsigned int i, unsigned int j)
    {
      double ret_val = 0.0;
      (i == j ? ret_val = 1.0: ret_val = 0.0);
      return ret_val;
    };

    double wx = 1.0;
    double wy = 1.0;
    double lx = 1.0;
    double ly = 1.0;
    double delta = 1.0;
    double dh = 1.0;
    double dv = 1.0;

    double gamma = 1.0;
    double B = 1.0;
    double mu = 1.0;
    double lambda = 1.0;
    double eta = 1.0;

  };


  #endif
