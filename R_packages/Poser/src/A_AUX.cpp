#include "A_AUX.h"
#include "B_EXPVAL.h"

double col_LogSumExp_cpp(arma::vec logX){
  double a = (logX).max();
  return(  a + log(accu( exp( logX-a ) )));
}

// -----------------------------------------------------------------------------

double row_LogSumExp_cpp(arma::rowvec logX){
  double a = (logX).max();
  return(  a + log(accu( exp( logX-a ) )));
}

// -----------------------------------------------------------------------------

arma::rowvec reverse_cumsum_row(arma::rowvec X){
  return( accu(X) - arma::cumsum(X));
}

// --------------------------------------- needed for faster computations ------
arma::mat utils_SUM_il_X_ElogLik(arma::mat Y,
                                 arma::cube XI,
                                 arma::mat mkcd_star,
                                 int J, int N, int K){


  arma::mat RES(K,J,arma::fill::zeros);
  int L = mkcd_star.n_rows;
  arma::mat Res(L, J,arma::fill::zeros);

  for(int i=0; i<N; i++){
    for(int l=0; l < L; l++){
      Res.row(l) = E_log_f_row(Y.row(i),
                               mkcd_star(l,0),
                               mkcd_star(l,1),
                               mkcd_star(l,2),
                               mkcd_star(l,3));
    }

    for(int k=0; k<K; k++){
      const arma::mat Temp = XI.col(k);
      RES.row(k) +=   Temp.row(i) * Res ; //NxJ
    }
  }
  return(arma::trans(RES));
}

// --------------------------------------- needed for faster computations ------

arma::cube utils_SUM_il_X_ElogLik_fusion(arma::field<arma::mat> Y, // T matrices N_t x J
                                         arma::field<arma::cube> X, // T cubes N_t x K x L
                                         arma::cube theta, // T matrices L x 4
                                         size_t J, size_t K, size_t T){

  arma::cube RES(J,K,T);

  for(size_t t=0; t<T; t++){

    arma::mat subRES(K,J,arma::fill::zeros);
    arma::mat subY = Y(t);
    size_t N = subY.n_rows;
    arma::vec ml  = (theta.slice(t)).col(0);
    arma::vec kl  = (theta.slice(t)).col(1);
    arma::vec al  = (theta.slice(t)).col(2);
    arma::vec bl  = (theta.slice(t)).col(3);

    for(size_t i=0; i<N; i++){
      arma::mat Res(ml.n_elem,J);
      for(size_t l=0; l<ml.n_elem; l++){
        Res.row(l) = E_log_f_row(subY.row(i),ml(l),kl(l),al(l),bl(l));
      }

      for(size_t k=0; k<K; k++){
        const arma::mat Temp = (X(t)).col(k);
        subRES.row(k) +=   Temp.row(i) * Res ; //KxJ
      }
    }
    RES.slice(t) = arma::trans(subRES);
  }

  return(RES);
}

// -----------------------------------------------------------------------------

int mod(int a, int n){
  return a - floor(a/n)*n;
}


arma::colvec reverse_cumsum_col(arma::colvec X){
  return( accu(X) - arma::cumsum(X));
}

// -----------------------------------------------------------------------------


// [[Rcpp::export]]
arma::rowvec log_EXP_pi(arma::rowvec a_k, arma::rowvec b_k){
  arma::rowvec E_tau = a_k/(a_k+b_k);
  arma::rowvec log_E_pi = log(E_tau);
  arma::rowvec sh = arma::shift( arma::cumsum(log(1-E_tau)) , +1);
  sh(0) = 0;
  log_E_pi = log_E_pi+sh;
  return(log_E_pi);
}

// antiRHO: matrix having the sum of the NN weights relative to each row
// [[Rcpp::export]]
arma::mat create_antiRHO(arma::mat RHO,
                         arma::field<arma::uvec> NN_inds){
  // to be called once at every iteration of the VB algo
  arma::mat antiRHO(RHO.n_rows, RHO.n_cols);
  for(size_t j=0; j<RHO.n_rows; j++){
    arma::uvec inds = NN_inds[j];
    antiRHO.row(j) = arma::sum( RHO.rows(inds), 0 );
  }
  return(antiRHO);
}



// [[Rcpp::export]]
double log_beta_normconts_PL(arma::mat antiRHO, double E_beta){
  arma::vec s = arma::sum(exp(antiRHO * E_beta),1);
  return( arma::accu(log(s)) );
}

