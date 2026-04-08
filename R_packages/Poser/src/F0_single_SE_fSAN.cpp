#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "D0_single_parupdates.h"
#include "E0_single_ELBOcomponents.h"


// [[Rcpp::export]]
Rcpp::List SE_fSAN(arma::mat Y,
                      int L, int K,
                      arma::rowvec A0,
                      arma::mat B0,
                      arma::mat mkcd_0,
                      arma::mat RHO_init,
                      arma::cube XI_init,
                      int nsim = 100,
                      int compute_elbo_every = 1,
                      double epsilon = .001,
                      bool verbose = true){

  arma::mat B_star(L,K);
  arma::rowvec A_star(K);

  arma::cube XI = XI_init;
  arma::mat RHO = RHO_init;

  arma::mat mkcd_star(L,4);

  arma::vec E(nsim); E.zeros();
  arma::mat Elbo_components(nsim,8);
  double old_E=0.0;

  //arma::cube THETA_collect(L,4,nsim);
  //arma::cube RHO_collect(Y.n_cols,K,nsim);
  //RHO_collect.zeros();
  //THETA_collect.zeros();

  int Q = nsim;


  arma::rowvec E_log_pi_k(K);
  A_star   = upd_astar_cpp(RHO, A0);
  E_log_pi_k = E_log_DIR_row(A_star);

  mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);
  arma::mat E_log_omega(L,K);
  // start the timer
  Timer timer;

  for(int sim = 0; sim< nsim; sim++){
    //Rcout << 1;
    Rcpp::checkUserInterrupt();
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    B_star   = upd_Bstar(XI, B0);
    for(int k1=0; k1<K; k1++){
      E_log_omega.col(k1) = E_log_DIR_col(B_star.col(k1));
    }
    // -------------------------------------------------------------------------
    A_star   = upd_astar_cpp(RHO, A0);
    E_log_pi_k = E_log_DIR_row(A_star);
    // -------------------------------------------------------------------------
    XI   = upd_XI(Y, RHO, E_log_omega,  mkcd_star);
    // -------------------------------------------------------------------------
    //Rcout << 2;
    RHO   = upd_RHO_onlyPi(Y,
                            XI,
                            mkcd_star,
                            E_log_pi_k);

    // -------------------------------------------------------------------------
    mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);
    // -------------------------------------------------------------------------
    //Rcout << 3;

    // monitor evolution of the parameters
    //RHO_collect.slice(sim) = RHO;
    //THETA_collect.slice(sim) = mkcd_star;
    // -------------------------------------------------------------------------


    if(mod(sim, compute_elbo_every) == 0){
      Elbo_components(sim,0) = elbo_p_Y(Y,XI,RHO,mkcd_star);
      Elbo_components(sim,1) = elbo_diff_omega(B0,B_star);
      Elbo_components(sim,2) = elbo_diff_R2(XI,E_log_omega);
      Elbo_components(sim,3) = elbo_p_THETA(mkcd_0,mkcd_star);
      Elbo_components(sim,4) = elbo_q_THETA(mkcd_star);
      Elbo_components(sim,5) = elbo_p_C_onlyPi(RHO, E_log_pi_k);
      Elbo_components(sim,6) = elbo_q_C_v2(RHO);
      Elbo_components(sim,7) = elbo_diff_pi(A0, A_star);

      E[sim] =
        Elbo_components(sim,0) +
        Elbo_components(sim,1) + Elbo_components(sim,2) +
        Elbo_components(sim,3) - Elbo_components(sim,4) +
        Elbo_components(sim,5) - Elbo_components(sim,6) +
        Elbo_components(sim,7);

      if(sim == 0){
        old_E = E[sim];
      }else if(sim >= 1) {
        double diff = (E[sim]-old_E);
        double diff_ratio = (E[sim]-old_E)/(old_E);
        //if(verbose){Rcout << "Iteration #" << sim+1 << " - Elbo value: " << E[sim] << " - Delta: " << diff << "\n";}
        if(verbose){Rcout << "Iteration #" << sim+1 << " - Elbo value: " << E[sim] << " - Delta: " << diff_ratio*100 << "%\n";}
        old_E = E[sim];

        //if((diff >= 0) & (diff_ratio < epsilon) ) {
        if(std::fabs(diff_ratio) < epsilon)  {
          if(verbose){Rcout << "Elbo final value:" << E[sim] << "---- Last Diff:" << (diff_ratio) <<"\n";}
          Q = sim;
          break;
        }
      }

    }

    }
  timer.step("end");        // record the starting point


  if(Q == nsim){Q = Q-1;}
  arma::colvec E2 = E.rows(1,Q);
  arma::mat E3 = Elbo_components.rows(1,Q);
  //arma::cube RHO_collect2 = RHO_collect.slices(1,Q);
  //arma::cube THETA_collect2 = THETA_collect.slices(1,Q);

  List results = List::create(_["XI"] = XI ,
                              _["RHO"] = RHO,
                              _["A_star"] = A_star,
                              _["B_star"] = B_star,
                              _["MKCD"] = mkcd_star,
                              _["Elbo_val"] = E2,
                              _["Elbo_comp"] = E3,
                              _["Y"] = Y,
                              //_["THETA"] = THETA_collect2,
                              //_["RHO_collect"] = RHO_collect2,
                              _["time"] = timer);

  return(results);
}

