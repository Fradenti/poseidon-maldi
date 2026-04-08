#include "C_NORMCONST.h"


// -----------------------------------------------------------------------------

arma::rowvec lbeta_normconst_row(const arma::rowvec& a,
                                 const arma::rowvec& b){
  arma::rowvec C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - (C_ab));
}

// -----------------------------------------------------------------------------

arma::colvec ldirichlet_normconst_col(const arma::mat& beta_lk){

  int K = beta_lk.n_cols;
  arma::colvec C(K);

  for(int k=0; k<K; k++){
    C(k) = lgamma(arma::accu(beta_lk.col(k))) -
           arma::accu(lgamma(beta_lk.col(k)));
  }
  return(C);
}

// -----------------------------------------------------------------------------

double ldirichlet_normconst_row2double(const arma::rowvec& alpha_k){

  double  C = lgamma(arma::accu(alpha_k)) - arma::accu(lgamma(alpha_k));

  return(C);
}
