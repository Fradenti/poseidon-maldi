#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "D0_single_parupdates.h"
#include "E0_single_ELBOcomponents.h"
#include "G_potts.h"


// [[Rcpp::export]]
Rcpp::List POSE_fSAN(arma::mat Y,
                          arma::field<arma::uvec> NN_inds,
                          double itemp,
                          arma::vec itemp_range,    // grid can be changed to runif set at each iteration
                          int L, int K,
                          arma::rowvec A0,
                          arma::mat B0,
                          arma::mat mkcd_0,
                          arma::mat RHO_init,
                          arma::cube XI_init,
                          int nsim = 100,
                          int compute_elbo_every=1,
                          double epsilon = .001,
                          int replicates_potts = 5,
                          bool verbose = true,
                          bool estimate_itemp=false){

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
  arma::vec iTEMP_collect(nsim);
  //RHO_collect.zeros();
  //THETA_collect.zeros();
  iTEMP_collect.zeros();
  arma::colvec estC(Y.n_cols);

  int Q = nsim;

  arma::rowvec E_log_pi_k(K);
  A_star   = upd_astar_cpp(RHO, A0);
  E_log_pi_k = E_log_DIR_row(A_star);

  mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);
  arma::mat E_log_omega(L,K);

  for(int sim = 0; sim< nsim; sim++){

    Rcpp::checkUserInterrupt();
    // -------------------------------------------------------------------------
    // mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);
    // -------------------------------------------------------------------------
    A_star   = upd_astar_cpp(RHO, A0);
    E_log_pi_k = E_log_DIR_row(A_star);
    B_star   = upd_Bstar(XI, B0);
    for(int k1=0; k1<K; k1++){
      E_log_omega.col(k1) = E_log_DIR_col(B_star.col(k1));
    }
    // -------------------------------------------------------------------------
    XI   = upd_XI(Y, RHO, E_log_omega,  mkcd_star);
    // -------------------------------------------------------------------------
    RHO   = upd_RHO(Y,
                    RHO,
                    XI,
                    mkcd_star,
                    NN_inds,
                    itemp,
                    E_log_pi_k,
                    replicates_potts);
    // -------------------------------------------------------------------------
    mkcd_star = upd_THETA_cpp(Y, XI, RHO, mkcd_0);
    if(estimate_itemp){
        arma::rowvec logEXPpi = A_star/arma::accu(A_star);
        itemp = bisection_LU2(RHO,
                              NN_inds,
                              logEXPpi,
                              min(itemp_range),
                              max(itemp_range));

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
      Elbo_components(sim,2) = elbo_diff_R2(XI,E_log_omega);
      Elbo_components(sim,3) = elbo_p_THETA(mkcd_0,mkcd_star);
      Elbo_components(sim,4) = elbo_q_THETA(mkcd_star);
      Elbo_components(sim,5) = elbo_p_C(RHO, NN_inds, E_log_pi_k, itemp);
      Elbo_components(sim,6) = elbo_q_C2(RHO);
      Elbo_components(sim,7) = elbo_diff_pi(A0, A_star);

      // MISSING iTEMP PART!

      E[sim] =
        Elbo_components(sim,0) +
        Elbo_components(sim,1) + Elbo_components(sim,2) +
        Elbo_components(sim,3) - Elbo_components(sim,4) +
        Elbo_components(sim,5) - Elbo_components(sim,6) + Elbo_components(sim,7);

      if(sim == 0 ){
        old_E = E[sim];
      }else if(sim >= 1) {

        double diff = (E[sim]-old_E);
        if(verbose){Rcout << "Iteration #" << sim+1 << " - Elbo value: " << E[sim] << " - Delta: " << diff << "\n";}
        old_E = E[sim];

        if((diff >= 0) & (diff < epsilon) ) {
          if(verbose){Rcout << "Elbo final value:" << E[sim] << "---- Last Diff:" << (diff) <<"\n";}
          Q = sim;
          break;
        }
        }
    }

  }

  if(Q == nsim){Q = Q-1;}
  arma::colvec E2 = E.rows(1,Q);
  arma::mat E3 = Elbo_components.rows(1,Q);
  //arma::cube RHO_collect2 = RHO_collect.slices(1,Q);
  //arma::cube THETA_collect2 = THETA_collect.slices(1,Q);
  arma::vec iTEMP_collect2 = iTEMP_collect.rows(1,Q);

  List results = List::create(_["XI"] = XI ,
                              _["RHO"] = RHO,
                              _["estC"] = estC,
                              _["A_star"] = A_star,
                              _["B_star"] = B_star,
                              _["MKCD"] = mkcd_star,
                              _["Elbo_val"] = E2,
                              _["Elbo_comp"] = E3,
                              _["Y"] = Y,
                              //_["THETA"] = THETA_collect2,
                              //_["RHO_collect"] = RHO_collect2,
                              _["itemp"] = itemp,
                              _["iTEMP_collect"] = iTEMP_collect2);
  return(results);
}

