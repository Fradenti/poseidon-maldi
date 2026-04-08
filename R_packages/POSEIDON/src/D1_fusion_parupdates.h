#ifndef D1_fusion_parupdates
#define D1_fusion_parupdates

#include <RcppArmadillo.h>
#include "B_EXPVAL.h"
#include "A_AUX.h"

arma::mat upd_single_THETA_cpp(arma::mat Y,
                        arma::cube XI,
                        arma::mat RHO,
                        arma::mat mkcd0);

arma::cube upd_THETA_cpp_fusion(arma::field<arma::mat> Y,
                                arma::field<arma::cube> X,
                                arma::mat R,
                                arma::cube mkcd0,
                                int const D);

// ------------------------------------------------------------------------------------------------

arma::mat upd_single_Bstar_cpp(arma::cube XI, arma::mat B0);
arma::cube upd_Bstar_cpp_fusion(arma::field<arma::cube> XI,
                                arma::cube B0,
                                int const L, int const K, int const D);

// ------------------------------------------------------------------------------------------------

arma::rowvec upd_astar_cpp_fusion(arma::mat R,
                                  arma::rowvec a0);


// ------------------------------------------------------------------------------------------------


arma::cube upd_single_XI_cpp(arma::mat Y,
                      arma::mat RHO,
                      arma::mat B_star,
                      arma::mat mkcd_star);

arma::field<arma::cube> upd_XI_cpp_fusion(arma::field<arma::mat> Y,
                                         arma::mat R,
                                         arma::cube bstarcube,
                                         arma::cube mkcd_star_cube, int const D);


// ------------------------------------------------------------------------------------------------

arma::mat upd_RHO_cpp_fusion_onlyPi(arma::field<arma::mat> Y,
                                    arma::field<arma::cube> X,
                                    arma::cube mkcd_starcube,
                                    arma::rowvec astar,
                                    int const D, int const J, int const K, int const L);

// ------------------------------------------------------------------------------------------------
arma::mat upd_RHO_cpp_fusion(arma::field<arma::mat> Y,
                             arma::mat RHO,
                             arma::field<arma::cube> X,
                             arma::cube mkcd_starcube,
                             arma::field<arma::uvec> NN_inds,
                             double itemp,
                             arma::rowvec astar,
                             int reps,
                             int const D, int const J, int const K, int const L);
//--------------------------------------------------------------------------------------------

arma::mat upd_RHO_cpp_fusion_onlyPotts(arma::field<arma::mat> Y,
                                       arma::mat RHO,
                                       arma::field<arma::cube> X,
                                       arma::cube mkcd_starcube,
                                       arma::field<arma::uvec> NN_inds,
                                       double itemp,
                                       int reps,
                                       int const D,
                                       int const J,
                                       int const K,
                                       int const L);

// ------------------------------------------------------------------------------------

arma::mat upd_Vk_cpp(double const a_tilde0_DP,
                     double const b_tilde0_DP, // if conc parameter is random this needs to be equal to s1/s2 and not to alpha_0
                     arma::mat RHO);

arma::rowvec upd_alpha_DP(arma::rowvec a_tilde_Vk,
                          arma::rowvec b_tilde_Vk,
                          arma::rowvec conc_hyper);


#endif
