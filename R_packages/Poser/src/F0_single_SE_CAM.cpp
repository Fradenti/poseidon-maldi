#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "D0_single_parupdates.h"
#include "E0_single_ELBOcomponents.h"


// [[Rcpp::export]]
Rcpp::List SE_CAM(arma::mat Y,
                    int L, int K,
                    arma::mat mkcd_0,
                    arma::mat RHO_init,
                    arma::cube XI_init,
                    arma::rowvec conc_hyper,// (s1, s2, s3, s4)
                    int nsim = 100,
                    int compute_elbo_every=1,
                    double epsilon = .001,
                    bool verbose = true){

  //Rcpp::Rcout << "0\n";
  arma::cube XI = XI_init;
  arma::mat RHO = RHO_init;

  arma::mat mkcd_star(L,4);

  arma::vec E(nsim); E.zeros();
  arma::mat Elbo_components(nsim,10);
  double old_E=0.0;

  //arma::cube THETA_collect(L,4,nsim);
  //arma::cube RHO_collect(Y.n_cols,K,nsim);
  arma::vec iTEMP_collect(nsim);
  //RHO_collect.zeros();
  //THETA_collect.zeros();
  iTEMP_collect.zeros();
  arma::colvec estC(Y.n_cols);

  int Q = nsim;

  arma::mat Parameters_DP_distr(3,K);
  arma::cube Parameters_DP_obser(L,K,3);

  arma::rowvec conc_distr_star(2);
  arma::rowvec conc_distr_0(2);
  conc_distr_0[0] = conc_hyper[0];
  conc_distr_star[0] = conc_hyper[0];
  conc_distr_0[1] = conc_hyper[1];
  conc_distr_star[1] = conc_hyper[1];

  arma::rowvec conc_obser_star(2);
  arma::rowvec conc_obser_0(2);
  conc_obser_star[0] = conc_hyper[2];
  conc_obser_star[1] = conc_hyper[3];
  conc_obser_0[0] = conc_hyper[2];
  conc_obser_0[1] = conc_hyper[3];


  arma::rowvec a_vk_star(K);
  arma::rowvec b_vk_star(K);
  arma::rowvec E_log_pi_k(K);
  mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);
  arma::mat a_ulk_bar(L,K);
  arma::mat b_ulk_bar(L,K);
  arma::mat E_log_omega(L,K);

  //Rcpp::Rcout << "1\n" << conc_obser_star << conc_distr_star << "\n";

  for(int sim = 0; sim< nsim; sim++){

    Rcpp::checkUserInterrupt();
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    Parameters_DP_obser   = Update_Ulk_cpp(XI, 1.0, conc_obser_star[0]/conc_obser_star[1] );
    a_ulk_bar = Parameters_DP_obser.slice(0);
    b_ulk_bar = Parameters_DP_obser.slice(1);
    E_log_omega = Parameters_DP_obser.slice(2);
    // -------------------------------------------------------------------------
    Parameters_DP_distr = upd_Vk_cpp(1.0, conc_distr_star[0]/conc_distr_star[1], RHO);
    a_vk_star   = Parameters_DP_distr.row(0);
    b_vk_star   = Parameters_DP_distr.row(1);
    E_log_pi_k  = Parameters_DP_distr.row(2);
    // -------------------------------------------------------------------------
    conc_distr_star =  upd_alpha_DP(a_vk_star, b_vk_star, conc_distr_0);
    conc_obser_star =  upd_beta_DP(a_ulk_bar, b_ulk_bar, conc_obser_0);
    // -------------------------------------------------------------------------
    //Rcpp::Rcout << "5\n";
    XI   = upd_XI(Y, RHO, E_log_omega,  mkcd_star);
    // -------------------------------------------------------------------------
    //Rcpp::Rcout << "6\n";
    RHO   = upd_RHO_onlyPi(Y,
                    XI,
                    mkcd_star,
                    E_log_pi_k);
    // -------------------------------------------------------------------------
    //Rcpp::Rcout << "7\n";
    mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);
    // -------------------------------------------------------------------------
    // monitor evolution of the parameters
    //RHO_collect.slice(sim) = RHO;
    //THETA_collect.slice(sim) = mkcd_star;
    // -------------------------------------------------------------------------


    if(mod(sim, compute_elbo_every) == 0){
      Elbo_components(sim,0) = elbo_p_Y(Y,XI,RHO,mkcd_star);
      Elbo_components(sim,1) = elbo_diff_U(a_ulk_bar, b_ulk_bar, conc_obser_star);
      Elbo_components(sim,2) = elbo_diff_R2(XI,E_log_omega);
      Elbo_components(sim,3) = elbo_p_THETA(mkcd_0,mkcd_star);
      Elbo_components(sim,4) = elbo_q_THETA(mkcd_star);
      Elbo_components(sim,5) = elbo_p_C_onlyPi(RHO, E_log_pi_k);
      Elbo_components(sim,6) = elbo_q_C_v2(RHO);
      //  Rcpp::Rcout << 4;
      Elbo_components(sim,7)  = elbo_diff_vk(a_vk_star, b_vk_star, conc_distr_star);
      // Rcpp::Rcout << 5;
      Elbo_components(sim,8)  = elbo_diff_alpha_DP(conc_distr_0, conc_distr_star);
      Elbo_components(sim,9)  = elbo_diff_beta_DP(conc_obser_0, conc_obser_star);
      // Rcpp::Rcout << 6;

      E[sim] =
        Elbo_components(sim,0) +
        Elbo_components(sim,1) + Elbo_components(sim,2) +
        Elbo_components(sim,3) - Elbo_components(sim,4) +
        Elbo_components(sim,5) - Elbo_components(sim,6) +
        Elbo_components(sim,7) + Elbo_components(sim,8) +
        Elbo_components(sim,9) ;

      if(sim == 0 ){
        old_E = E[sim];
      }else if(sim >= 1) {
        double diff = (E[sim]-old_E);
        if(verbose){Rcout << "Iteration #" << sim+1 << " - Elbo value: " << E[sim] << " - Delta: " << diff << "\n";}
        old_E = E[sim];

        if( (diff >= 0) & (diff < epsilon) ) {
          if(verbose){Rcout << "Elbo final value:" << E[sim] << "---- Last Diff:" << (diff) <<"\n";}
          Q = sim;
          break;
        }}
    }

  }


  if(Q == nsim){Q = Q-1;}
  arma::colvec E2 = E.rows(1,Q);
  arma::mat E3 = Elbo_components.rows(1,Q);
  //arma::cube RHO_collect2 = RHO_collect.slices(1,Q);
  //arma::cube THETA_collect2 = THETA_collect.slices(1,Q);

  List results = List::create(_["XI"] = XI ,
                              _["RHO"] = RHO,
                              _["MKCD"] = mkcd_star,
                              _["a_tilde_k"] = a_vk_star,
                              _["b_tilde_k"] = b_vk_star,
                              _["a_bar_lk"] = a_ulk_bar,
                              _["b_bar_lk"] = b_ulk_bar,
                              _["alphaDP"] = conc_distr_star[0]/conc_distr_star[1],
                              _["betaDP"]  = conc_obser_star[0]/conc_obser_star[1],
                              _["Elbo_val"] = E2,
                              _["Elbo_comp"] = E3,
                              _["Y"] = Y,
                              //_["THETA"] = THETA_collect2,
                              //_["RHO_collect"] = RHO_collect2,
                              _["estC"] = estC);

  return(results);
}

