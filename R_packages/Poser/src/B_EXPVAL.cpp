#include "B_EXPVAL.h"


arma::colvec E_log_beta_col(arma::colvec a,
                        arma::colvec b){

  int n = a.n_elem;
  arma::colvec res(n);
  for(int i=0; i<n; i++){
    res[i] = R::digamma(a[i]) - R::digamma(a[i] + b[i]);
  }
  return(res);
}

// -----------------------------------------------------------------------------

arma::rowvec E_log_beta_row(arma::rowvec a,
                            arma::rowvec b){

  int n = a.n_elem;
  arma::rowvec res(n);
  for(int i=0; i<n; i++){
    res[i] = R::digamma(a[i]) - R::digamma(a[i] + b[i]);
  }
  return(res);
}

// -----------------------------------------------------------------------------

arma::colvec E_log_DIR_col(arma::colvec a){

  int n = a.n_elem;
  arma::colvec res(n);
  double Sum_a = arma::accu(a);

  for(int i = 0; i < n; i++){
    res[i] = R::digamma( a[i] );
  }

  return(res - R::digamma(Sum_a));
}

// -----------------------------------------------------------------------------

arma::rowvec E_log_DIR_row(arma::rowvec a){

  int n = a.n_elem;
  arma::rowvec res(n);
  double Sum_a = arma::accu(a);

  for(int i = 0; i < n; i++){
    res[i] = R::digamma( a[i] );
  }

  return(res - R::digamma(Sum_a));
}

// -----------------------------------------------------------------------------

arma::colvec E_log_IG(arma::colvec a,
                          arma::colvec b){
  int n = a.n_elem;
  arma::colvec res(n);

  for(int i = 0; i < n; i++){
    res[i] = R::digamma( a[i] );
  }

  arma::colvec ares = log(b) - res;
  return(ares);
}


// -----------------------------------------------------------------------------
// Some different versions of E[log_lik]
// -----------------------------------------------------------------------------
arma::mat E_log_f_mat(const arma::mat& Y,
                      double ml,
                      double kl,
                      double cl,
                      double dl){

  return(
    - .5 * (
        std::log(dl) - R::digamma(cl)   +
        1.0/kl + ( (cl/dl) * (Y-ml) % (Y-ml) )  )
  );
}

// -----------------------------------------------------------------------------

arma::rowvec E_log_f_row(const arma::rowvec& Y,
                          double ml,
                          double kl,
                          double cl,
                          double dl){
  return(
    - .5 * (
        std::log(dl) - R::digamma(cl)   +
        1.0/kl + ( (cl/dl) * (Y-ml) % (Y-ml) )  )
  );
}
