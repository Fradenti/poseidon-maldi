#ifndef E1_fusion_ELBOcomponents
#define E1_fusion_ELBOcomponents

#include <RcppArmadillo.h>
#include "A_AUX.h"
#include "B_EXPVAL.h"
#include "C_NORMCONST.h"


// -----------------------------------------------------------------------------
double elbo_p_THETA(arma::mat mkcd_0,
                    arma::mat mkcd_star);
double elbo_q_THETA(arma::mat mkcd_star);
double elbo_diff_THETA_fusion(arma::cube mkcd0,
                              arma::cube mkcd_star_cube,
                              int const D);
// -----------------------------------------------------------------------------
double elbo_diff_omega(arma::mat B0,
                       arma::mat B_star);
double elbo_diff_omega_fusion(arma::cube B0,
                              arma::cube B_star,
                              int const D);
// -----------------------------------------------------------------------------
double elbo_diff_R2(arma::cube XI,
                    arma::mat B_star);
double elbo_diff_R_fusion(arma::field<arma::cube> XI,
                       arma::cube beta_star_cube,
                       int const D);
// -----------------------------------------------------------------------------
double elbo_p_Y(arma::mat Y,
                arma::cube XI,
                arma::mat RHO,
                arma::mat mkcd_star);
double elbo_p_Y_fusion(arma::field<arma::mat> Y,
                       arma::field<arma::cube> XI,
                       arma::mat RHO,
                       arma::cube mkcd_star,
                       int const D);
// -----------------------------------------------------------------------------
double elbo_diff_pi(arma::rowvec A0,
                    arma::rowvec A_star);
// -----------------------------------------------------------------------------
double elbo_p_C_fusion(arma::mat R,
                       arma::field<arma::uvec> NN_inds,
                       arma::rowvec E_log_pi_k, double beta_potts);

double elbo_q_C_fusion(arma::mat R);

double elbo_diff_C_fusion(arma::mat R,
                       arma::field<arma::uvec> NN_inds,
                       arma::rowvec E_log_pi_k, double beta_potts);


double elbo_p_v_CP(arma::colvec a_tilde_k,
                   arma::colvec b_tilde_k,
                   arma::colvec conc_star,
                   int const K);

double elbo_q_v_CP(arma::colvec a_tilde_k,
                   arma::colvec b_tilde_k,
                   int const K);

double elbo_diff_vk(arma::rowvec a_tilde_Vk,
                    arma::rowvec b_tilde_Vk,
                    arma::rowvec conc_star, int const K);

double elbo_diff_alpha_DP(arma::rowvec conc_hyper,
                          arma::rowvec conc_star);

double elbo_diff_C_fusion_onlyPotts(arma::mat RHO,
                                    arma::mat antiRHO,
                                    double itemp);

double elbo_diff_C_fusion_onlyPi(arma::mat RHO,
                              arma::rowvec E_log_pi_k);
#endif
