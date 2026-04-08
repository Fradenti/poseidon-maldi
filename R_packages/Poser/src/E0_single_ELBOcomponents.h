#ifndef E0_single_ELBOcomponents
#define E0_single_ELBOcomponents

#include <RcppArmadillo.h>
#include "A_AUX.h"
#include "B_EXPVAL.h"
#include "C_NORMCONST.h"

double elbo_p_THETA(const arma::mat& mkcd_0, const arma::mat& mkcd_star);
double elbo_q_THETA(const arma::mat& mkcd_star);


double elbo_diff_omega(arma::mat B0, arma::mat B_star);

double elbo_diff_R(arma::cube XI, arma::mat B_star);
double elbo_diff_R2(arma::cube XI, arma::mat B_star);

double elbo_diff_pi(arma::rowvec A0,arma::rowvec A_star);

double elbo_p_Y(arma::mat Y, arma::cube XI, arma::mat RHO, arma::mat mkcd_star);

double elbo_diff_vk(arma::rowvec a_tilde_Vk, arma::rowvec b_tilde_Vk,
                    arma::rowvec conc_star);
double elbo_diff_U(arma::mat a_ulk_bar,
                      arma::mat b_ulk_bar,
                      arma::rowvec conc_obser_star);

double elbo_diff_alpha_DP(arma::rowvec conc_hyper,
                          arma::rowvec conc_star);
double elbo_diff_beta_DP(arma::rowvec conc_hyper,
                          arma::rowvec conc_star);

double elbo_p_C(arma::mat R, arma::field<arma::uvec> NN_inds,
                arma::rowvec E_log_pi_k, double beta_potts);
double elbo_p_C_onlyPi(arma::mat RHO,
                      arma::rowvec E_log_pi_k);
double elbo_p_C_onlyPotts(arma::mat R, arma::mat antiRHO, double beta_potts);

double elbo_q_C_v2(arma::mat R);

double elbo_diff_itemp(arma::mat RHO, arma::mat antiRHO, double itemp);
#endif
