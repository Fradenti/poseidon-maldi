#ifndef B_EXPVAL
#define B_EXPVAL

#include <RcppArmadillo.h>


arma::colvec E_log_beta_col(arma::colvec a,
                        arma::colvec b);
arma::rowvec E_log_beta_row(arma::rowvec a,
                            arma::rowvec b);
arma::colvec E_log_DIR_col(arma::colvec a);

arma::rowvec E_log_DIR_row(arma::rowvec a);

arma::colvec E_log_IG(arma::colvec a,
                      arma::colvec b);

arma::mat E_log_f_mat(const arma::mat& Y,
                      double ml,
                      double kl,
                      double cl,
                      double dl);

arma::rowvec E_log_f_row(const arma::rowvec& Y,
                         double ml,
                         double kl,
                         double cl,
                         double dl);


#endif
