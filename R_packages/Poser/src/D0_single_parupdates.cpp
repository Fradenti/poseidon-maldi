#include "D0_single_parupdates.h"


// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat upd_THETA_cpp(const arma::mat& Y,
                        const arma::cube& XI,
                        const arma::mat& RHO,
                        const arma::mat& mkcd0){ // (mu0,kappa0,tau0,gamma0)


  int L = XI.n_slices;
  int J = Y.n_cols;
  int N = Y.n_rows;
  arma::mat Y2 = Y%Y;
  arma::mat mkcd_star(L,4);
  arma::mat n_ijl_slice(N,J);
  arma::vec N_l(L);

  for(int l=0; l<L; l++){

    //arma::mat Xtemp = ; // N x K
    n_ijl_slice  = XI.slice(l) * arma::trans(RHO); // N x J
    N_l(l)          = arma::accu(n_ijl_slice);

    mkcd_star(l,1) = mkcd0(l,1) + N_l(l);      //k_l
    mkcd_star(l,2) = mkcd0(l,2) + N_l(l) * 0.5 ; //c_l

    double sum_Y_l  = arma::accu(n_ijl_slice % Y );
    double sum_Y2_l = arma::accu(n_ijl_slice % Y2 );

    double Ybar_l = 0;
    double S2_l   = 0;

    if(N_l[l] > 0){
      Ybar_l = sum_Y_l/N_l[l];
      S2_l   = sum_Y2_l - pow(Ybar_l,2) * N_l[l];
    }

    mkcd_star(l,0) = (mkcd0(l,0) * mkcd0(l,1) + sum_Y_l) / (mkcd_star(l,1)) ;

    mkcd_star(l,3) = mkcd0(l,3) +
                     0.5 * (
                      S2_l +
                      ((mkcd0(l,1) * N_l[l]) / (mkcd_star(l,1))) *
                         pow( (Ybar_l - mkcd0(l,0)), 2)
                      );
  }

  return(mkcd_star);
}

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat upd_Vk_cpp(const double& a_tilde0_DP,
                     const double& b_tilde0_DP, // if conc parameter is random this needs to be equal to s1/s2 and not to alpha_0
                     const arma::mat RHO){

  int K = RHO.n_cols;
  arma::rowvec M_k = arma::sum(RHO,0);

  arma::rowvec rev_cs_mk = reverse_cumsum_row(M_k);

  arma::rowvec a_tilde_vk      = M_k       + a_tilde0_DP;
  a_tilde_vk[K-1] = 1;
  arma::rowvec b_tilde_vk      = rev_cs_mk + b_tilde0_DP;
  b_tilde_vk[K-1] = 1e-8;

  arma::rowvec E_ln_Vk    = E_log_beta_row(a_tilde_vk, b_tilde_vk);
  arma::rowvec E_ln_1mVk  = E_log_beta_row(b_tilde_vk, a_tilde_vk);
  arma::rowvec sE_ln_1mVk = shift(E_ln_1mVk, +1);

  sE_ln_1mVk[0] = 0;
  arma::rowvec CS_E_ln_1mVk = arma::cumsum(sE_ln_1mVk);

  arma::mat results(3,K);

  results.row(0) = a_tilde_vk;
  results.row(1) = b_tilde_vk;
  results.row(2) = E_ln_Vk + CS_E_ln_1mVk ;

  return(results);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::rowvec upd_alpha_DP(const arma::rowvec& a_tilde_Vk,
                          const arma::rowvec& b_tilde_Vk,
                          const arma::rowvec& conc_hyper){
  int K = a_tilde_Vk.n_cols;
  arma::rowvec upd_par(2);
  arma::rowvec a_copy = a_tilde_Vk;
  a_copy.shed_col(K-1);
  arma::rowvec b_copy = b_tilde_Vk;
  b_copy.shed_col(K-1);

  upd_par[0] = conc_hyper[0] + (K - 1) ;
  upd_par[1] = conc_hyper[1] - arma::accu(E_log_beta_row(b_copy,
                                                         a_copy));

  return(upd_par);
}

// -----------------------------------------------------------------------------

arma::mat upd_Bstar(arma::cube XI, arma::mat B0){
  arma::mat N_kl = arma::sum(XI,0);
  return( arma::trans(N_kl) + B0);
}

// ----------------------------------------------------------------------------

arma::rowvec upd_astar_cpp(arma::mat RHO,
                           arma::rowvec a0){
  arma::rowvec M_k = (arma::sum(RHO,0));
  return( M_k + a0 );
}

// -----------------------------------------------------------------------------

arma::cube upd_XI(const arma::mat& Y,
                  const arma::mat& RHO,
                  const arma::mat& E_log_omega,
                  const arma::mat& mkcd_star){

  int K = RHO.n_cols;
  int L = mkcd_star.n_rows;
  int N = Y.n_rows;

  arma::cube log_XI(N,K,L);

  for(int l=0; l<L; l++){
    arma::mat Rik = E_log_f_mat(Y,
                                mkcd_star(l,0),mkcd_star(l,1),
                                mkcd_star(l,2),mkcd_star(l,3))  * RHO;
    arma::mat Gik = arma::repelem(E_log_omega.row(l),N,1);
    log_XI.slice(l) =  Rik + Gik;
  }

  for(int k = 0; k < K; k++){
    for(int i = 0; i < N; i++){
      log_XI.tube(i,k) = log_XI.tube(i,k) - col_LogSumExp_cpp(log_XI.tube(i,k));
    }}

  return(exp(log_XI));
}

// ----------------------------------------------------------------------------
arma::mat upd_RHO_onlyPi(const arma::mat& Y,
                         const arma::cube& XI,
                         const arma::mat& mkcd_star,
                         const arma::rowvec& E_log_pi_k){

  int J = Y.n_cols;
  int K = XI.n_cols;
  int N = XI.n_rows;

  arma::mat Temp_R(J,K, arma::fill::zeros);
  arma::mat SUM_il_X_ElogLik(J,K);


  SUM_il_X_ElogLik = utils_SUM_il_X_ElogLik(Y,
                                            XI,
                                            mkcd_star,
                                            J,N,K); // J x K


  arma::rowvec temp(K);
    for(int j = 0; j < J; j++){

      temp = SUM_il_X_ElogLik.row(j) +  E_log_pi_k;

      Temp_R.row(j) =   (exp( temp - row_LogSumExp_cpp(temp)) ) ;
    }

  return(Temp_R);
}

arma::mat upd_RHO_onlyPotts(const arma::mat& Y,
                  const arma::mat& RHO,
                  const arma::cube& XI,
                  const arma::mat& mkcd_star,
                  const arma::field<arma::uvec>& NN_inds,
                  double beta_potts,
                  int reps){

  int J = Y.n_cols;
  int K = XI.n_cols;
  int N = XI.n_rows;

  arma::mat Temp_R = RHO;
  arma::mat SUM_il_X_ElogLik(J,K);


  SUM_il_X_ElogLik = utils_SUM_il_X_ElogLik(Y,
                                            XI,
                                            mkcd_star,
                                            J,N,K); // J x K

  arma::rowvec hmrf=arma::zeros<arma::rowvec>(K);

  arma::rowvec temp(K);
  for(int rr = 0; rr < reps; rr++){
    for(int j = 0; j < J; j++){

      if(beta_potts > 0.0){
        arma::uvec inds = NN_inds[j] ;
        hmrf = arma::sum(Temp_R.rows(inds),0) * beta_potts; // 1 x K
      }

      temp = SUM_il_X_ElogLik.row(j)  + hmrf;

      Temp_R.row(j) =   (exp( temp - row_LogSumExp_cpp(temp)) ) ;
    }
  }

  return(Temp_R);
}

arma::mat upd_RHO(const arma::mat& Y,
                  const arma::mat& RHO,
                  const arma::cube& XI,
                  const arma::mat& mkcd_star,
                  arma::field<arma::uvec>& NN_inds,
                  double beta_potts,
                  const arma::rowvec& E_log_pi_k,
                  int reps){

  int J = Y.n_cols;
  int K = XI.n_cols;
  int N = XI.n_rows;

  arma::mat Temp_R = RHO;
  arma::mat SUM_il_X_ElogLik(J,K);


  SUM_il_X_ElogLik = utils_SUM_il_X_ElogLik(Y,
                                            XI,
                                            mkcd_star,
                                            J,N,K); // J x K

  arma::rowvec hmrf(K);
  hmrf.zeros();

  arma::rowvec temp(K);
  for(int rr = 0; rr < reps; rr++){
    for(int j = 0; j < J; j++){

      if(beta_potts > 0.0){
        arma::uvec inds = NN_inds[j] ;
        hmrf = arma::sum(Temp_R.rows(inds),0) * beta_potts; // 1 x K
      }

      temp = SUM_il_X_ElogLik.row(j)  + hmrf +  E_log_pi_k;

      Temp_R.row(j) =   (exp( temp - row_LogSumExp_cpp(temp)) ) ;
    }
  }

  return(Temp_R);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------


arma::cube Update_Ulk_cpp(arma::cube XI, // n x k x l
                          double const& a_bar,
                          double const& b_bar){

  int K = XI.n_cols;
  int L = XI.n_slices;
  //to improve
  arma::mat a_bar_Ulk(L,K);
  arma::mat b_bar_Ulk(L,K);

  arma::mat N_kl = (arma::sum(XI,0));
  arma::mat N_lk = arma::trans(N_kl);
  arma::mat E_lnOmega_lk(L,K);

  for(int k = 0; k<K; k++){

    arma::colvec rQl_k  = reverse_cumsum_col(N_lk.col(k));
    arma::colvec G1 = a_bar + N_lk.col(k);
    arma::colvec G2 = b_bar + rQl_k;

    G1[L-1] = 1;
    G2[L-1] = 1e-10;

    a_bar_Ulk.col(k) =  G1;
    b_bar_Ulk.col(k) =  G2;

    arma::colvec E_ln_Ul_k    = E_log_beta_col(G1, G2);
    arma::colvec E_ln_1mUl_k  = E_log_beta_col(G2, G1);
    arma::colvec sE_ln_1mUl_k = shift(E_ln_1mUl_k, +1);
    sE_ln_1mUl_k[0] = 0;

    arma::colvec CS_E_ln_1mUL_k = arma::cumsum(sE_ln_1mUl_k);
    E_lnOmega_lk.col(k) = CS_E_ln_1mUL_k + E_ln_Ul_k;
  }



  arma::cube results(L,K,3);
  results.slice(0) = a_bar_Ulk;
  results.slice(1) = b_bar_Ulk;
  results.slice(2) = E_lnOmega_lk;

  return(results);
}



arma::rowvec upd_beta_DP(arma::mat a_ulk_bar,
                        arma::mat b_ulk_bar,
                        arma::rowvec conc_hyper){

  int L = a_ulk_bar.n_rows;
  int K = a_ulk_bar.n_cols;

  arma::rowvec upd_par(2);
  a_ulk_bar.shed_row(L-1);
  b_ulk_bar.shed_row(L-1);

  arma::colvec Res(K);

  for(int k = 0; k < K; k ++){

    Res[k] =  arma::accu( E_log_beta_col(a_ulk_bar.col(k),b_ulk_bar.col(k)) );

  }

  upd_par[0] = conc_hyper[2] + K * (L - 1);
  upd_par[1] = conc_hyper[3] - arma::accu(Res);

 return(upd_par);
}
