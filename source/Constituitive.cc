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


    de.clear();

    Tensor<2,DIM> C;
    C = transpose(F)*F;
    for(unsigned int i = 0; i < DIM; i++)
      for(unsigned int j = 0; j < DIM; j++)
        C[i][j] += grad_w[i]*grad_w[j];

    Tensor<4,DIM> dC_dF;
    Tensor<3,DIM> dC_dgradw;

    for (unsigned int i = 0; i < DIM; i ++)
      for (unsigned int j = 0; j < DIM; j ++)
        for (unsigned int k = 0; k < DIM; k ++)
        {
          dC_dgradw[i][j][k] = grad_w[i]*(j == k ? 1.0 : 0.0) + grad_w[j]*(i == k? 1.0 : 0.0);
          for (unsigned int l = 0; l < DIM; l ++)
            dC_dF[i][j][k][l] = (l == i? 1.0: 0.0)*F[k][j] + (j == l ? 1.0 : 0.0)*F[k][i];
        }


    Tensor<2,DIM> dPhi_dC;
    get_dPhi_dC(C, dPhi_dC);


    de.dW_dF = double_contract<0,0,1,1>(dPhi_dC, dC_dF);

    de.dW_dgrad_w = double_contract<0,0,1,1>(dPhi_dC, dC_dgradw);

    de.dW_dlap_w = B*lap_w;



  }

  void plate_energy::get_DDE(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w, DDE_data &dde)
  {

    dde.clear();

    Tensor<2,DIM> C;
    C = transpose(F)*F;
    for(unsigned int i = 0; i < DIM; i++)
      for(unsigned int j = 0; j < DIM; j++)
        C[i][j] += grad_w[i]*grad_w[j];

    double detC = determinant(C);
    double J = sqrt(detC);
    double Ibar1 = trace(C)/J;

    Tensor<2,DIM> C_inv = invert(C);

    Tensor<2,DIM> LD;
    LD[0][0] = lx*lx/(dh*dh);
    LD[1][1] = ly*ly/(dv*dv);

    double by_eta = 1.0/eta;

    double alpha = J*sin(gamma)*lx*ly/(dh*dv);

    double LD_dd_C = double_contract<0,0,1,1>(LD, C);

    double by_J = 1.0/J;

    // penalty term
    Tensor<4,DIM> d2Phi_d2C;
    for (unsigned int p = 0; p < DIM; p ++)
      for (unsigned int q = 0; q < DIM; q ++)
        for (unsigned int l = 0; l < DIM; l ++)
          for (unsigned int m = 0; m < DIM; m ++)
           {
              d2Phi_d2C[p][q][l][m] += (LD[l][m] - alpha*C_inv[l][m])
                                          *(LD[p][q] - alpha*C_inv[p][q]);
              d2Phi_d2C[p][q][l][m] += (LD_dd_C - 2.0*alpha - cos(gamma)*cos(gamma))*alpha
                                       *(-0.5*C_inv[l][m]*C_inv[p][q] + 0.5*(C_inv[p][l]*C_inv[m][q] + C_inv[p][m]*C_inv[l][q]));

              d2Phi_d2C[p][q][l][m] += 0.25*(dirac(q, 0)*dirac(p, 1) + dirac(p, 0)*dirac(q, 1))
                                          *(dirac(m, 0)*dirac(l, 1) + dirac(l, 0)*dirac(m, 1));
           }

    d2Phi_d2C *= by_eta;

    // elastic energy
    for (unsigned int p = 0; p < DIM; p ++)
      for (unsigned int q = 0; q < DIM; q ++)
        for (unsigned int l = 0; l < DIM; l ++)
          for (unsigned int m = 0; m < DIM; m ++)
           {
             d2Phi_d2C[p][q][l][m] += 0.5*mu*( -0.5*by_J*C_inv[l][m]*dirac(p, q)
                                              - 0.5*( dirac(l, m)*by_J*C_inv[p][q] - 0.5*Ibar1*C_inv[l][m]*C_inv[p][q]
                                                      -0.5*Ibar1*(C_inv[p][l]*C_inv[m][q] + C_inv[p][m]*C_inv[l][q])));

             d2Phi_d2C[p][q][l][m] += 0.5*lambda*((detC*C_inv[l][m] - 0.5*J*C_inv[l][m])*C_inv[p][q] +
                                                   -0.5*(detC - J)*(C_inv[p][l]*C_inv[m][q] + C_inv[p][m]*C_inv[l][q]));
           }


    Tensor<2,DIM> dPhi_dC;
    get_dPhi_dC(C, dPhi_dC);

    Tensor<4,DIM> dC_dF;
    Tensor<3,DIM> dC_dgradw;

    for (unsigned int i = 0; i < DIM; i ++)
      for (unsigned int j = 0; j < DIM; j ++)
        for (unsigned int k = 0; k < DIM; k ++)
        {
          dC_dgradw[i][j][k] = grad_w[i]*(j == k ? 1.0 : 0.0) + grad_w[j]*(i == k? 1.0 : 0.0);
          for (unsigned int l = 0; l < DIM; l ++)
            dC_dF[i][j][k][l] = (l == i? 1.0: 0.0)*F[k][j] + (j == l ? 1.0 : 0.0)*F[k][i];
        }


    Tensor<6, DIM> d2C_d2F;
    Tensor<4, DIM> d2C_d2w;

    for (unsigned int i = 0; i < DIM; i ++)
      for (unsigned int j = 0; j < DIM; j ++)
        for (unsigned int r = 0; r < DIM; r ++)
          for (unsigned int s = 0; s < DIM; s ++)
          {
            d2C_d2w[i][j][r][s] += dirac(j, r)*dirac(i, s) + dirac(j, s)*dirac(i, r);
            for (unsigned int p = 0; p < DIM; p ++)
              for (unsigned int q = 0; q < DIM; q ++)
                d2C_d2F[i][j][r][s][p][q] += dirac(s, i)*dirac(r, p)*dirac(j, q) +
                                              dirac(p, r)*dirac(q, i)*dirac(j, s);

          }



    dde.d2W_d2F = double_contract<2, 0, 3, 1>(double_contract<0, 0, 1, 1>(dC_dF, d2Phi_d2C), dC_dF); // probably wrong
    dde.d2W_d2F += double_contract<0, 0, 1, 1>(dPhi_dC, d2C_d2F);

    dde.d2W_dF_dgrad_w = double_contract<2, 0, 3, 1>(double_contract<0, 0, 1, 1>(dC_dF, d2Phi_d2C), dC_dgradw);

    dde.d2W_d2grad_w =  double_contract<1, 0, 2, 1>(double_contract<0, 0, 1, 1>(dC_dgradw, d2Phi_d2C), dC_dgradw);

    dde.d2W_d2grad_w += double_contract<0, 0, 1, 1>(dPhi_dC, d2C_d2w);

    dde.d2W_d2lap_w = B;


  }

  void plate_energy::get_dPhi_dC(Tensor<2,DIM> &C, Tensor<2,DIM> &dPhi_dC)
  {

    dPhi_dC = 0.0;

    double detC = determinant(C);
    double J = sqrt(detC);
    double Ibar1 = trace(C)/J;
    Tensor<2,DIM> C_inv = invert(C);

    Tensor<2,DIM> LD;
    LD[0][0] = lx*lx/(dh*dh);
    LD[1][1] = ly*ly/(dv*dv);

    double parentisis_quoties = (lx*lx*C[0][0]/(dh*dh)) + (ly*ly*C[1][1]/(dv*dv)) - (2.0*sin(gamma)*lx*ly/(dv*dh))*J - cos(gamma)*cos(gamma);


    dPhi_dC += (1.0/eta)*parentisis_quoties*(LD - sin(gamma)*(lx*ly/(dv*dh))*J*C_inv);
    dPhi_dC[0][1] += (1.0/(2.0*eta))*C[0][1];
    dPhi_dC[1][0] += (1.0/(2.0*eta))*C[0][1];

    dPhi_dC += -0.25*mu*trace(C)*(1.0/J)*C_inv;
    for (unsigned int i = 0; i < DIM; i ++)
      dPhi_dC[i][i] += 0.5*mu/J;

    dPhi_dC += lambda*(J - 1.0)*0.5*J*C_inv;
  }

  void plate_energy::numerical_deriv_internal(Tensor<2,DIM> &F, Tensor<1,DIM> grad_w, double &lap_w)
   {

      double epsilon = 1.0e-7;

      double E_old = get_E(F, grad_w, lap_w);

      DE_data de_dat;
      get_DE(F, grad_w, lap_w, de_dat);

      Tensor<2,DIM> dE_dF_numer;

      std::cout << "     NUMERICAL DERIVATIVE SHIT : \n";
      for (unsigned int i = 0; i < DIM; i ++)
        for (unsigned int j = 0; j < DIM; j ++)
        {
          F[i][j] += epsilon;
          double E = get_E(F, grad_w, lap_w);

          dE_dF_numer[i][j] = (E - E_old)/epsilon;
          F[i][j] -= epsilon;

          std::cout << de_dat.dW_dF[i][j] << " " << dE_dF_numer[i][j] << " " << fabs(dE_dF_numer[i][j] - de_dat.dW_dF[i][j]) << std::endl;
        }

      Tensor<1,DIM> dE_dgradw_numer;
      for (unsigned int i = 0; i < DIM; i ++)
      {
        grad_w[i] += epsilon;
        double E = get_E(F, grad_w, lap_w);

        dE_dgradw_numer[i] = (E - E_old)/epsilon;
        grad_w[i] -= epsilon;
        std::cout << "     " << de_dat.dW_dgrad_w[i] << " " << dE_dgradw_numer[i] << " " << fabs(dE_dgradw_numer[i] - de_dat.dW_dgrad_w[i]) << std::endl;

      }


      std::cout << "\n\n";


   }






  #endif
