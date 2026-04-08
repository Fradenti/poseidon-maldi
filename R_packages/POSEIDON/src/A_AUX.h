#ifndef A_AUX
#define A_AUX

#include <RcppArmadillo.h>

double col_LogSumExp_cpp(arma::vec logX);
double row_LogSumExp_cpp(arma::rowvec logX);
arma::rowvec reverse_cumsum_row(arma::rowvec X);
arma::mat utils_SUM_il_X_ElogLik(arma::mat Y,
                                  arma::cube XI,
                                  arma::mat mkcd_star,
                                  int J, int N, int K);
arma::cube utils_SUM_il_X_ElogLik_fusion(arma::field<arma::mat> Y,
                                         arma::field<arma::cube> X,
                                         arma::cube theta,
                                         int J, int K, int T);
int mod(int a, int n);
arma::mat create_antiRHO(arma::mat RHO,
                         arma::field<arma::uvec> NN_inds);
double log_beta_normconts_PL(arma::mat antiRHO, double E_beta);

#endif
