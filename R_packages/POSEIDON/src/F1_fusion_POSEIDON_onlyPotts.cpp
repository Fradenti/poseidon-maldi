#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "D1_fusion_parupdates.h"
#include "E1_fusion_ELBOcomponents.h"
#include "G_potts.h"


// [[Rcpp::export]]
Rcpp::List POSEIDON_cpp_onlyPotts(const arma::field<arma::mat>& Y_list,
                                  const arma::field<arma::uvec>& NN_inds,
                             double itemp,
                             arma::vec itemp_grid,
                             int const L, int const K,
                             arma::cube B0,
                             arma::cube mkcd_0,
                             arma::mat RHO_init,
                             arma::field<arma::cube> XI_init,
                             int nsim = 100,
                             int compute_elbo_every = 1,
                             double epsilon = .001,
                             bool verbose = true,
                             int replicates_potts = 5,
                             bool estimate_itemp=false){
  int const J = Y_list(0).n_cols;
  int const D = Y_list.n_elem;

  Rcpp::Rcout << D;
  arma::cube B(L,K,D);
  arma::field<arma::cube> XI = XI_init;
  arma::mat RHO = RHO_init;
  arma::mat antiRHO = RHO;
  arma::vec E(nsim);
  double old_E=0.0;
  arma::mat Elbo_components(nsim,7);
  arma::cube mkcd_star = mkcd_0;
  int Q = nsim;
  arma::colvec estC(J);
  arma::rowvec E_log_pi_k(K);
  arma::vec iTEMP_collect(nsim);
  //Rcpp::Rcout <<1;
  // start the timer
  Timer timer;

  timer.step("start");        // record the starting point

  for(int sim = 0; sim< nsim; sim++){
    Rcpp::checkUserInterrupt();
    // -------------------------------------------------------------------------
    //Rcpp::Rcout <<1.1;
    mkcd_star = upd_THETA_cpp_fusion(Y_list, XI, RHO, mkcd_0, D);
      // -------------------------------------------------------------------------
      //Rcpp::Rcout <<1.2;
      B    = upd_Bstar_cpp_fusion(XI, B0, L, K ,D);
          // -------------------------------------------------------------------------
          //Rcpp::Rcout <<1.3;
          XI   = upd_XI_cpp_fusion(Y_list, RHO, B,  mkcd_star, D);
            // -------------------------------------------------------------------------
            //Rcpp::Rcout <<2;
            RHO = upd_RHO_cpp_fusion_onlyPotts(Y_list,
                                                 RHO,
                                                 XI,
                                                 mkcd_star,
                                                 NN_inds,
                                                 itemp,
                                                 replicates_potts,
                                                 D, J, K, L);
              // -------------------------------------------------------------------------
              antiRHO = create_antiRHO(RHO,NN_inds);

              if(estimate_itemp){

                  itemp = est_itemp_PL2(RHO,antiRHO,itemp_grid);

                if(verbose){Rcpp::Rcout << "Estimated inv. temp. = "<<itemp << "\n";}
              }
              iTEMP_collect(sim) = itemp;

              // -------------------------------------------------------------------------

              Elbo_components(sim,0) = elbo_p_Y_fusion(Y_list,XI,RHO,mkcd_star,D);
              Elbo_components(sim,1) = elbo_diff_omega_fusion(B0,B, D);
              Elbo_components(sim,2) = elbo_diff_THETA_fusion(mkcd_0,mkcd_star, D);
              Elbo_components(sim,3) = elbo_diff_R_fusion(XI,B,D);
              Elbo_components(sim,4) = 0;//elbo_diff_C_fusion_onlyPotts(RHO,
                                      //                              antiRHO,
                                      //                              itemp);
              Elbo_components(sim,5) = elbo_q_C_fusion(RHO);
              // Term number 4 simplifies with the itemp elbo term under the PL approximation.


              E[sim] =
                Elbo_components(sim,0) +
                Elbo_components(sim,1) +
                Elbo_components(sim,2) +
                Elbo_components(sim,3) +
                Elbo_components(sim,4) -
                Elbo_components(sim,5);

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
  timer.step("end");        // record the starting point

  if(Q == nsim){Q = Q-1;}
  arma::colvec E2 = E.rows(1,Q);
  arma::mat E3 = Elbo_components.rows(1,Q);
  arma::vec iTEMP_collect2 = iTEMP_collect.rows(1,Q);

  List results = List::create(_["XI"] = XI,
                              _["RHO"] = RHO,
                              _["B_star"] = B,
                              _["itemp"] = itemp,
                              _["MKCD"] = mkcd_star,
                              _["Elbo_val"] = E2,
                              _["Elbo_comp"] = E3,
                              _["Y_list"] = Y_list,
                              _["iTEMP_collect"] = iTEMP_collect2,
                              _["time"] = timer);


  return(results);
}
