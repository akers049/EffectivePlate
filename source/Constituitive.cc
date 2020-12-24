#ifndef CONSTITUITIVE_CC_
#define CONSTITUITIVE_CC_

#include "Constituitive.h"

#define DIM 2

using namespace dealii;



  void plate_energy::set_params(double wx_, double wy_,  double lx_, double ly_, double delta_)
  {
    wx = wx_;
    wy = wy_;
    lx = lx_;
    ly = ly_;
    delta = delta_;

    dh = sqrt(lx*lx + (ly - 2.0*wy - delta)*(ly - 2.0*wy - delta));
    dv = sqrt(ly*ly + (lx - 2.0*wx - delta)*(lx - 2.0*wx - delta));
    gamma = M_PI/2.0 - atan((ly - 2.0*wy - delta)/lx) - atan((-lx + 2.0*wx + delta)/ly);
  }


  double plate_energy::get_E(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w)
  {

    Tensor<2,DIM> C;
    C = transpose(F)*F;
    for(unsigned int i = 0; i < DIM; i++)
      for(unsigned int j = 0; j < DIM; j++)
        C[i][j] += grad_w[i]*grad_w[j];

    double J = sqrt(determinant(C));
    double Ibar1 = trace(C)/J;

    double Phi_m = (1.0/(2.0*eta))*(
                pow( (lx*lx*C[0][0]/(dh*dh)) + (ly*ly*C[1][1]/(dv*dv)) - (2.0*sin(gamma)*lx*ly/(dv*dh))*J - cos(gamma)*cos(gamma), 2.0)
               + C[1][0]*C[1][0]) + 0.5*mu*(Ibar1 - 2.0) + 0.5*lambda*(J - 1.0)*(J - 1.0);

    double Phi_b = B*lap_w*lap_w/2.0;
    return Phi_m + Phi_b;
  }

  void plate_energy::get_DE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DE_data &de)
  {

    Tensor<2,DIM> C;
    C = transpose(F)*F;
    for(unsigned int i = 0; i < DIM; i++)
      for(unsigned int j = 0; j < DIM; j++)
        C[i][j] += grad_w[i]*grad_w[j];

    double detC = determinant(C);
    double J = sqrt(detC);
    double Ibar1 = trace(C)/J;

    Tensor<4,DIM> dC_dF;
    Tensor<3,DIM> dC_dgradw;

    for (unsigned int i = 0; i < DIM; i ++)
      for (unsigned int j = 0; j < DIM; j ++)
        for (unsigned int k = 0; k < DIM; k ++)
        {
          dC_dgradw[i][j][k] = grad_w[i]*(j == k ? 1.0 : 0.0) + grad_w[j]*(i == k? 1.0 : 0.0);
          for (unsigned int l = 0; l < DIM; l ++)
            dC_dF[i][j][k][l] = (l == i? 1.0: 0.0)*F[k][l] + (j == l ? 1.0 : 0.0)*F[k][i];
        }

    Tensor<2,DIM> C_inv = invert(C);
    Tensor<2,DIM> dPhi_dC;

    Tensor<2,DIM> LD;
    LD[0][0] = lx*lx/(dh*dh);
    LD[1][1] = ly*ly/(dv*dv);

    double parentisis_quoties = (lx*lx*C[0][0]/(dh*dh)) + (ly*ly*C[1][1]/(dv*dv)) - (2.0*sin(gamma)*lx*ly/(dv*dh))*J - cos(gamma)*cos(gamma);

    Tensor<2,DIM> dPhi_dC;

    dPhi_dC += (1.0/eta)*parentisis_quoties*(LD - sin(gamma)*(lx*ly/(dv*dh))*J*C_inv);
    dPhi_dC[0][1] += (1.0/(2.0*eta))*C[0][1];
    dPhi_dC[1][0] += (1.0/(2.0*eta))*C[0][1];

    dPhi_dC += -0.25*mu*trace(C)*(1.0/J)*C_inv;
    for (unsigned int i = 0; i < DIM; i ++)
      dPhi_dC[i][i] += 0.5*mu/J;

    dPhi_dC += lambda*(J - 1.0)*0.5*J*C_inv;

    de.dW_dF = double_contract<0,0,1,1>(dPhi_dC, dC_dF);

    de.dW_dgrad_w = double_contract<0,0,1,1>(dPhi_dC, dC_dgradw);



  }

  void plate_energy::get_DDE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DDE_data &dde)
  {

  }






  #endif
