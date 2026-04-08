#ifndef C_NORM_CONST
#define C_NORM_CONST

#include <RcppArmadillo.h>


arma::rowvec lbeta_normconst_row(const arma::rowvec& a,
                                 const arma::rowvec& b);

arma::rowvec lbeta_normconst_row(const arma::rowvec& a,
                                 const arma::rowvec& b);

arma::colvec ldirichlet_normconst_col(const arma::mat& beta_lk);

double ldirichlet_normconst_row2double(const arma::rowvec& alpha_k);

#endif
