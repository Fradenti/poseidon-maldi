#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "D1_fusion_parupdates.h"
#include "E1_fusion_ELBOcomponents.h"
#include "G_potts.h"


// [[Rcpp::export]]
Rcpp::List POSEIDON_cpp_fiSAN(const arma::field<arma::mat>& Y_list,
                              const arma::field<arma::uvec>& NN_inds,
                             double itemp,
                             arma::vec itemp_grid,
                             int const L, int const K, int const D,
                             arma::cube B0,
                             arma::cube mkcd_0,
                             arma::mat RHO_init,
                             arma::field<arma::cube> XI_init,
                             arma::rowvec conc_hyper, // (s1, s2)
                             int nsim = 100,
                             int compute_elbo_every = 1,
                             double epsilon = .001,
                             bool verbose = true,
                             int replicates_potts = 5,
                             bool estimate_itemp = false){
  int const J = Y_list(0).n_cols;
  arma::cube B(L,K,D);
  arma::field<arma::cube> XI = XI_init;
  arma::mat RHO = RHO_init;
  arma::mat antiRHO = RHO_init;
  arma::vec E(nsim);
  double old_E=0.0;
  arma::mat Elbo_components(nsim,7);
  arma::cube mkcd_star = mkcd_0;
  arma::rowvec conc_star(2);
  conc_star = conc_hyper;
  arma::rowvec a_vk_star(K);
  arma::rowvec b_vk_star(K);
  arma::rowvec E_log_pi_k(K);
  arma::mat Parameters_DP_alpha(3,K);
  arma::colvec estC(J);
  arma::vec iTEMP_collect(nsim);

  //-----------------------------------
  int Q = nsim;
  //Rcpp::Rcout << 1;
  for(int sim = 0; sim< nsim; sim++){
    Rcpp::checkUserInterrupt();

    // -------------------------------------------------------------------------
    mkcd_star = upd_THETA_cpp_fusion(Y_list, XI, RHO, mkcd_0, D);
    // -------------------------------------------------------------------------
    B    = upd_Bstar_cpp_fusion(XI, B0, L, K ,D);
    // -------------------------------------------------------------------------
    Parameters_DP_alpha = upd_Vk_cpp(1.0, conc_star[0]/conc_star[1], RHO);
    a_vk_star   = Parameters_DP_alpha.row(0);
    b_vk_star   = Parameters_DP_alpha.row(1);
    E_log_pi_k  = Parameters_DP_alpha.row(2);
    conc_star =  upd_alpha_DP(a_vk_star, b_vk_star, conc_hyper);
    // -------------------------------------------------------------------------
    XI   = upd_XI_cpp_fusion(Y_list, RHO, B,  mkcd_star, D);
    // -------------------------------------------------------------------------
    RHO  = upd_RHO_cpp_fusion(Y_list,
                              RHO,
                              XI,
                              mkcd_star,
                              NN_inds,
                              itemp,
                              E_log_pi_k,
                              replicates_potts,
                              D, J, K, L);
    antiRHO = create_antiRHO(RHO, NN_inds);
    // -------------------------------------------------------------------------
    if(estimate_itemp){


        itemp = est_itemp_PL2_withpi(RHO,antiRHO,E_log_pi_k,itemp_grid);

      if(verbose){Rcpp::Rcout << "Estimated inv. temp. = "<< itemp << "\n";}
      iTEMP_collect(sim) = itemp;
    }

    // -------------------------------------------------------------------------

    Elbo_components(sim,0) = elbo_p_Y_fusion(Y_list, XI, RHO, mkcd_star, D);
    Elbo_components(sim,1) = elbo_diff_omega_fusion(B0, B, D);
    Elbo_components(sim,2) = elbo_diff_THETA_fusion(mkcd_0, mkcd_star, D);
    Elbo_components(sim,3) = elbo_diff_R_fusion(XI, B, D);
    //Elbo_components(sim,4) = elbo_diff_C_fusion(RHO, NN_inds, E_log_pi_k, itemp);
    // this term simplifies with the item term E(log(q(beta))) in the elbo.
    // E(log(p(beta))) is fixed
    Elbo_components(sim,4) = elbo_q_C_fusion(RHO);
    Elbo_components(sim,5) = elbo_diff_vk(a_vk_star, b_vk_star, conc_star, K);
    Elbo_components(sim,6) = elbo_diff_alpha_DP(conc_hyper, conc_star);

    E[sim] =
      Elbo_components(sim,0) +
      Elbo_components(sim,1) +
      Elbo_components(sim,2) +
      Elbo_components(sim,3) -
      Elbo_components(sim,4) +
      Elbo_components(sim,5) +
      Elbo_components(sim,6);

    if(sim == 0){
      old_E = E[sim];
    }else if(sim >= 1) {
      //double diff = (E[sim]-old_E);
      double diff_ratio = (E[sim]-old_E)/(old_E);
      //if(verbose){Rcout << "Iteration #" << sim+1 << " - Elbo value: " << E[sim] << " - Delta: " << diff << "\n";}
      if(verbose){Rcpp::Rcout << "Iteration #" << sim+1 << " - Elbo value: " << E[sim] << " - Delta: " << diff_ratio*100 << "%\n";}
      old_E = E[sim];

      //if((diff >= 0) & (diff_ratio < epsilon) ) {
      if(std::fabs(diff_ratio) < epsilon)  {
        if(verbose){Rcpp::Rcout << "Elbo final value:" << E[sim] << "---- Last Diff:" << (diff_ratio)*100 <<"\n";}
        Q = sim;
        break;
      }
    }

  }
  if(Q == nsim){Q = Q-1;}
  arma::colvec E2 = E.rows(1,Q);
  arma::mat E3 = Elbo_components.rows(1,Q);
  arma::vec iTEMP_collect2 = iTEMP_collect.rows(1,Q);

  List results = List::create(_["XI"] = XI,
                              _["RHO"] = RHO,
                              _["B_star"] = B,
                              _["Stick_pars"] = Parameters_DP_alpha,
                              _["alphaDP"] = conc_star[0]/conc_star[1],
                              _["itemp"] = itemp,
                              _["MKCD"] = mkcd_star,
                              _["Elbo_val"] = E2,
                              _["Elbo_comp"] = E3,
                              _["Y_list"] = Y_list,
                              _["iTEMP_collect"] = iTEMP_collect2);

  return(results);
}
