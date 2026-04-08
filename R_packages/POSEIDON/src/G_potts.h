#ifndef G_potts
#define G_potts

#include <RcppArmadillo.h>



double EQ66( arma::mat RHO, arma::field<arma::uvec> NN_inds);

arma::mat p_MF(double beta,
               arma::mat antiRHO,
               arma::rowvec log_e_pi);

arma::mat create_pMF(double beta,
                     arma::mat antiRHO,
                     arma::rowvec log_e_pi);

double bisection_LU2(arma::mat RHO,
                     arma::field<arma::uvec> NN_inds,
                     arma::rowvec log_e_pi,
                     double low, double upp, int N_sims= 30);
double bisection_LU3(arma::mat RHO,
                     arma::field<arma::uvec> NN_inds,
                     arma::rowvec log_e_pi,
                     double low, double upp, int N_sims= 30);

double est_itemp_PL2(arma::mat RHO,
                     arma::mat antiRHO,
                     arma::vec itemp_grid);
double est_itemp_PL2_withpi(arma::mat RHO,
                            arma::mat antiRHO,
                            arma::rowvec ElogPi,
                            arma::vec itemp_grid);
#endif
