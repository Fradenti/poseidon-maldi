#ifndef C_NORM_CONST
#define C_NORM_CONST

#include <RcppArmadillo.h>

double lbeta_normconst_cpp(double a, double b);
arma::rowvec lbeta_normconst_row(arma::rowvec a, arma::rowvec b);
arma::colvec ldirichlet_normconst_col(arma::mat beta_lk);
arma::mat lbeta_normconst_mat_cpp(arma::mat a, arma::mat b);
double       ldirichlet_normconst_row2double(arma::rowvec alpha_k);

#endif
