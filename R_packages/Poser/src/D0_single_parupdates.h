#ifndef D0_single_parupdates
#define D0_single_parupdates

#include <RcppArmadillo.h>
#include "B_EXPVAL.h"
#include "A_AUX.h"


arma::mat upd_THETA_cpp(const arma::mat& Y,
                        const arma::cube& XI,
                        const arma::mat& RHO,
                        const arma::mat& mkcd0);


arma::mat upd_Vk_cpp(const double& a_tilde0_DP,
                     const double& b_tilde0_DP, // if conc parameter is random this needs to be equal to s1/s2 and not to alpha_0
                     const arma::mat RHO);

arma::rowvec upd_alpha_DP(const arma::rowvec& a_tilde_Vk,
                          const arma::rowvec& b_tilde_Vk,
                          const arma::rowvec& conc_hyper);

arma::mat upd_Bstar(arma::cube XI, arma::mat B0);

arma::rowvec upd_astar_cpp(arma::mat RHO, arma::rowvec a0);

arma::cube upd_XI(const arma::mat& Y,
                  const arma::mat& RHO,
                  const arma::mat& E_log_omega,
                  const arma::mat& mkcd_star);

arma::mat upd_RHO(const arma::mat& Y,
                  const arma::mat& RHO,
                  const arma::cube& XI,
                  const arma::mat& mkcd_star,
                  arma::field<arma::uvec>& NN_inds,
                  double beta_potts,
                  const arma::rowvec& E_log_pi_k,
                  int reps);

arma::mat upd_RHO_onlyPi(const arma::mat& Y,
                         const arma::cube& XI,
                         const arma::mat& mkcd_star,
                         const arma::rowvec& E_log_pi_k);

arma::mat upd_RHO_onlyPotts(const arma::mat& Y,
                            const arma::mat& RHO,
                            const arma::cube& XI,
                            const arma::mat& mkcd_star,
                            const arma::field<arma::uvec>& NN_inds,
                            double beta_potts,
                            int reps);

arma::cube Update_Ulk_cpp(arma::cube XI, // n x k x l
                          double const& a_bar,
                          double const& b_bar);

arma::rowvec upd_beta_DP( arma::mat a_ulk_bar,
                         arma::mat b_ulk_bar,
                         arma::rowvec conc_hyper);

#endif
