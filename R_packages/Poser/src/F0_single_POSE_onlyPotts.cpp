#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "D0_single_parupdates.h"
#include "E0_single_ELBOcomponents.h"
#include "G_potts.h"

// [[Rcpp::export]]
Rcpp::List POSE_onlyPotts( arma::mat Y,
                          arma::field<arma::uvec> NN_inds,
                          double itemp,
                          arma::vec itemp_grid,    // grid can be changed to runif set at each iteration
                          int L, int K,
                          arma::mat B0,
                          arma::mat mkcd_0,
                          arma::mat RHO_init,
                          arma::cube XI_init,
                          int nsim = 100,
                          int compute_elbo_every=1,
                          double epsilon = .001,
                          int replicates_potts = 5,
                          bool verbose = 1,
                          int algo_itemp = 0, // 0 for mcgregory, 1 for Lu+Arbel
                          bool estimate_itemp=0){

  arma::mat B_star(L,K);
  arma::cube XI = XI_init;
  arma::mat RHO = RHO_init;
  arma::mat antiRHO = RHO;
  arma::mat mkcd_star(L,4);
  arma::vec E(nsim); E.zeros();
  arma::mat Elbo_components(nsim,8);
  double old_E=0.0;
  //arma::cube THETA_collect(L,4,nsim);
  //arma::cube RHO_collect(Y.n_cols,K,nsim);
  arma::vec iTEMP_collect(nsim);
  //RHO_collect.zeros();
  //THETA_collect.zeros();
  iTEMP_collect.zeros();
  arma::colvec estC(Y.n_cols);
  int Q = nsim;

  arma::rowvec E_log_pi_k(K);
  E_log_pi_k.fill(0);//log(1./K));
  arma::rowvec logEXPpi(K);
  logEXPpi.zeros();
  // to initialize
  mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);
  // start the timer
  Timer timer;

  timer.step("start");        // record the starting point
  for(int sim = 0; sim< nsim; sim++){

    Rcpp::checkUserInterrupt();
    // -------------------------------------------------------------------------
    B_star   = upd_Bstar(XI, B0);
    // -------------------------------------------------------------------------
    XI   = upd_XI(Y, RHO, B_star,  mkcd_star);
    // -------------------------------------------------------------------------
    RHO   = upd_RHO(Y,
                    RHO,
                    XI,
                    mkcd_star,
                    NN_inds,
                    itemp,
                    E_log_pi_k,
                    replicates_potts);
    antiRHO = create_antiRHO(RHO, NN_inds);
    // -------------------------------------------------------------------------
    mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);


    if(estimate_itemp){
      if(algo_itemp== 0){

        itemp = est_itemp_PL2(RHO,antiRHO,itemp_grid);

      }else{

        itemp = bisection_LU2(RHO,
                      NN_inds,
                      logEXPpi,
                      min(itemp_grid),
                      max(itemp_grid));

      }
     if(verbose){Rcpp::Rcout << "Estimated inv. temp. = "<<itemp << "\n";}
    }
    // monitor evolution of the parameters
    //RHO_collect.slice(sim) = RHO;
    //THETA_collect.slice(sim) = mkcd_star;
    iTEMP_collect(sim) = itemp;
    // -------------------------------------------------------------------------

    if(mod(sim, compute_elbo_every) == 0){
      Elbo_components(sim,0) = elbo_p_Y(Y,XI,RHO,mkcd_star);
      Elbo_components(sim,1) = elbo_diff_omega(B0,B_star);
      Elbo_components(sim,2) = elbo_diff_R2(XI,B_star);
      Elbo_components(sim,3) = elbo_p_THETA(mkcd_0,mkcd_star);
      Elbo_components(sim,4) = elbo_q_THETA(mkcd_star);
      Elbo_components(sim,5) = 0; //elbo_p_C_onlyPotts(RHO, antiRHO, itemp);
      Elbo_components(sim,6) = elbo_q_C_v2(RHO);
      Elbo_components(sim,7) = 0; //elbo_diff_itemp(RHO, antiRHO, itemp);
      // Component in position 5 and 7 cancel out - I avoid to compute them

      E[sim] =
        Elbo_components(sim,0) +
        Elbo_components(sim,1) + Elbo_components(sim,2) +
        Elbo_components(sim,3) - Elbo_components(sim,4) +
        Elbo_components(sim,5) - Elbo_components(sim,6) + Elbo_components(sim,7);

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
  arma::vec iTEMP_collect2 = iTEMP_collect.rows(1,Q);

  List results = List::create(_["XI"] = XI ,
                              _["RHO"] = RHO,
                              _["estC"] = estC,
                              _["B_star"] = B_star,
                              _["MKCD"] = mkcd_star,
                              _["Elbo_val"] = E2,
                              _["Elbo_comp"] = E3,
                              _["Y"] = Y,
                              //_["THETA"] = THETA_collect2,
                              //_["RHO_collect"] = RHO_collect2,
                              _["itemp"] = itemp,
                              _["iTEMP_collect"] = iTEMP_collect2,
                              _["time"] = timer);
  return(results);
}
