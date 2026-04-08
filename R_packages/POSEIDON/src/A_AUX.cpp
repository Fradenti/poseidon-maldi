#include "A_AUX.h"
#include "B_EXPVAL.h"

double col_LogSumExp_cpp(arma::vec logX){
  double a = max(logX);
  return(  a + log(arma::accu( exp( logX-a ) )));
}

// -----------------------------------------------------------------------------

double row_LogSumExp_cpp(arma::rowvec logX){
  double a = max(logX);
  return(  a + log(arma::accu( exp( logX - a ) )));
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


  arma::mat RES(K,J, arma::fill::zeros);
  int L = mkcd_star.n_rows;

  arma::mat Res(L, J, arma::fill::zeros);
  for(int i=0; i<N; i++){
  Res.zeros();
    for(int l=0; l < L; l++){
      Res.row(l) = E_log_f_row(Y.row(i),
                               mkcd_star(l,0),
                               mkcd_star(l,1),
                               mkcd_star(l,2),
                               mkcd_star(l,3));
    }

    for(int k=0; k<K; k++){
      arma::mat Temp = XI.col(k);
      RES.row(k) +=   Temp.row(i) * Res ; //NxJ
    }
  }
  return(arma::trans(RES));
}

// --------------------------------------- needed for faster computations ------

arma::cube utils_SUM_il_X_ElogLik_fusion(arma::field<arma::mat> Y, // D matrices N_t x J
                                         arma::field<arma::cube> XI, // D cubes N_t x K x L
                                         arma::cube mkcd_star, // D matrices L x 4
                                         int const J, int const K, int const D){

  arma::cube RES(J,K,D, arma::fill::zeros);
  for(int t=0; t<D; t++){
    int N = Y(t).n_rows;
    RES.slice(t) =utils_SUM_il_X_ElogLik(Y(t),
                                         XI(t),
                                         mkcd_star.slice(t),
                                         J, N, K);
  }
  return(RES);
}

// -----------------------------------------------------------------------------

int mod(int a, int n){
  return a - floor(a/n)*n;
}

// antiRHO: matrix having the sum of the NN weights relative to each row
// [[Rcpp::export]]
arma::mat create_antiRHO(arma::mat RHO,
                         arma::field<arma::uvec> NN_inds){
  // to be called once at every iteration of the VB algo
  arma::mat antiRHO = RHO;
  for(size_t j=0; j<RHO.n_rows; j++){
    arma::uvec inds = NN_inds[j];
    antiRHO.row(j) = arma::sum( RHO.rows(inds), 0 );
  }
  return(antiRHO);
}



// [[Rcpp::export]]
double log_beta_normconts_PL(arma::mat antiRHO, double E_beta){
  arma::vec s = arma::sum(exp( antiRHO * E_beta ),1);
  return( arma::accu(log(s)) );
 }

// [[Rcpp::export]]
double log_beta_normconts_PL2(arma::mat antiRHO, double E_beta){
  arma::vec temp = antiRHO * E_beta;
  double m = temp.max();
  return arma::accu( m + arma::log(arma::sum(arma::exp(temp - m),1)) );
  
}
