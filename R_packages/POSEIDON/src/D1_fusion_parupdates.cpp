#include "D1_fusion_parupdates.h"

arma::mat upd_single_THETA_cpp(arma::mat Y,
                        arma::cube XI,
                        arma::mat RHO,
                        arma::mat mkcd0){ // (mu0,kappa0,tau0,gamma0)


  int L = XI.n_slices;
  int J = Y.n_cols;
  int N = Y.n_rows;
  arma::mat Y2 = Y % Y;
  arma::mat mkcd_star(L,4, arma::fill::zeros);
  arma::cube n_ijl(N,J,L, arma::fill::zeros);
  arma::vec N_l(L, arma::fill::zeros);

  for(int l=0; l<L; l++){

    //arma::mat Xtemp = ; // N x K
    n_ijl.slice(l)  = XI.slice(l) * arma::trans(RHO); // N x J
    N_l(l)          = arma::accu(n_ijl.slice(l));

    mkcd_star(l,1) = mkcd0(l,1) + N_l(l);      //k_l
    mkcd_star(l,2) = mkcd0(l,2) + N_l(l) * 0.5 ; //c_l

    double sum_Y_l  = arma::accu(n_ijl.slice(l) % Y );
    double sum_Y2_l = arma::accu(n_ijl.slice(l) % Y2 );

    double Ybar_l = 0.0;
    double S2_l   = 0.0;

    if(N_l[l] > 0){
      Ybar_l = sum_Y_l / N_l[l];
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

arma::cube upd_THETA_cpp_fusion(arma::field<arma::mat> Y,
                                arma::field<arma::cube> XI,
                                arma::mat RHO,
                                arma::cube mkcd0,
                                int const D){

  int L = XI(0).n_slices;
  arma::cube mkcd_star(L,4,D, arma::fill::zeros);

  for(int t=0; t<D; t++){

  mkcd_star.slice(t) =  upd_single_THETA_cpp(Y(t),
                                      XI(t),
                                      RHO,
                                      mkcd0.slice(t));

  }

  return(mkcd_star);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

arma::mat upd_single_Bstar_cpp(arma::cube XI, arma::mat B0){
  arma::mat N_kl = arma::sum(XI,0);
  return( arma::trans(N_kl) + B0);
}


arma::cube upd_Bstar_cpp_fusion(arma::field<arma::cube> XI,
                                arma::cube B0,
                                int const L, int const K, int const D){
  arma::cube RES(L,K,D, arma::fill::zeros);
  for( int t = 0; t<D; t++ ){
    RES.slice(t) = upd_single_Bstar_cpp(XI(t),B0.slice(t));
  }
  return( RES );

}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// identical to single case
arma::rowvec upd_astar_cpp_fusion(arma::mat R,
                           arma::rowvec a0){
  arma::rowvec Mk = arma::sum(R,0);
  return( Mk + a0);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

arma::cube upd_single_XI_cpp(arma::mat Y,
                  arma::mat RHO,
                  arma::mat B_star,
                  arma::mat mkcd_star){

  int K = RHO.n_cols;
  int L = mkcd_star.n_rows;
  int N = Y.n_rows;

  arma::cube log_XI(N,K,L, arma::fill::zeros);
  arma::mat Hl_k(L,K);

  for(int k1=0; k1<K; k1++){
    Hl_k.col(k1) = E_log_DIR_col(B_star.col(k1));
  }

  for(int l=0; l<L; l++){
    arma::mat Rik = E_log_f_mat(Y,
                                mkcd_star(l,0),mkcd_star(l,1),
                                mkcd_star(l,2),mkcd_star(l,3))  * RHO;
    arma::mat Gik = arma::repelem(Hl_k.row(l),N,1);
    log_XI.slice(l) =  Rik + Gik;
  }

  for(int k = 0; k < K; k++){
    for(int i = 0; i < N; i++){
      log_XI.tube(i,k) = log_XI.tube(i,k) - col_LogSumExp_cpp(log_XI.tube(i,k));
    }}

  return(exp(log_XI));
}


arma::field<arma::cube> upd_XI_cpp_fusion(arma::field<arma::mat> Y,
                                          arma::mat RHO,
                                          arma::cube bstarcube,
                                          arma::cube mkcd_star_cube,
                                          int T){

  arma::field<arma::cube> newX(T);
  for(int t=0; t<T; t++){

    newX(t) = upd_single_XI_cpp(Y(t),
                     RHO,
                     bstarcube.slice(t),
                     mkcd_star_cube.slice(t));

  }

  return(newX);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat upd_RHO_cpp_fusion(arma::field<arma::mat> Y,
                             arma::mat RHO,
                             arma::field<arma::cube> X,
                             arma::cube mkcd_starcube,
                             arma::field<arma::uvec> NN_inds,
                             double itemp,
                             arma::rowvec E_log_pi_k,
                             int reps,
                             int const D,
                             int const J,
                             int const K,
                             int const L){


  // arma::mat Temp_R(J,K); - ERRORE INIZIALIZZAVO CON GARBAGE!!!
  arma::mat Temp_R = RHO;
  arma::cube SUM_il_X_ElogLik(J,K,D);
  arma::mat E_log_omega_lk(L,K);

  //Rcpp::Rcout << "-7.0-";
  SUM_il_X_ElogLik = utils_SUM_il_X_ElogLik_fusion(Y,
                                                   X,
                                                   mkcd_starcube,
                                                   J,K,D);

  arma::mat SUM_il_X_ElogLik_collapsed = arma::sum(SUM_il_X_ElogLik, 2); // JxK
  // Rcpp::Rcout << SUM_il_X_ElogLik_collapsed << "\n----------\n";

  arma::rowvec temp(K);

  if(itemp == 0.0){
    for(int j = 0; j < J; j++){
      temp = SUM_il_X_ElogLik_collapsed.row(j) + E_log_pi_k;
      Temp_R.row(j) =   (exp(temp - row_LogSumExp_cpp(temp))) ;
      }
  }else{
    arma::rowvec Q(K);
    for(int rr = 0; rr < reps; rr++){
      for(int j = 0; j < J; j++){
      //Rcpp::Rcout << "-7.3-";
      arma::uvec inds = NN_inds[j] ;
      //Rcpp::Rcout << "-7.4-" << Temp_R.rows(inds) << "\n";
      Q = arma::sum(Temp_R.rows(inds),0) * itemp; // 1 x K
      //Rcpp::Rcout << "-7.5-";
      temp = SUM_il_X_ElogLik_collapsed.row(j)  +
        Q + E_log_pi_k;
      //Rcpp::Rcout << "-7.6-";
      Temp_R.row(j) =   (exp(temp - row_LogSumExp_cpp(temp))) ;
      //Rcpp::Rcout << "-7.7-";
      }
    }
  }
  return(Temp_R);
}



// [[Rcpp::export]]
arma::mat upd_RHO_cpp_fusion_onlyPi(arma::field<arma::mat> Y,
                                    arma::field<arma::cube> X,
                                    arma::cube mkcd_starcube,
                                    arma::rowvec E_log_pi_k,
                                    int const D,
                                    int const J,
                                    int const K,
                                    int const L){


  arma::mat Temp_R(J,K, arma::fill::zeros);
  arma::cube SUM_il_X_ElogLik(J,K,D);
  arma::mat E_log_omega_lk(L,K);

  //Rcpp::Rcout << "-7.0-";
  SUM_il_X_ElogLik = utils_SUM_il_X_ElogLik_fusion(Y,
                                                   X,
                                                   mkcd_starcube,
                                                   J,K,D);

  arma::mat SUM_il_X_ElogLik_collapsed = arma::sum(SUM_il_X_ElogLik, 2); // JxK

  arma::rowvec temp(K);

  for(int j = 0; j < J; j++){
      temp = SUM_il_X_ElogLik_collapsed.row(j) + E_log_pi_k;
      Temp_R.row(j) =   (exp(temp - row_LogSumExp_cpp(temp))) ;
  }

  return(Temp_R);
}

// [[Rcpp::export]]
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
                                      int const L){

  // arma::mat Temp_R(J,K); - ERRORE INIZIALIZZAVO CON GARBAGE!!!
  arma::mat Temp_R = RHO;
  arma::cube SUM_il_X_ElogLik(J,K,D);
  arma::mat E_log_omega_lk(L,K);

  //Rcpp::Rcout << "-7.0-";
  SUM_il_X_ElogLik = utils_SUM_il_X_ElogLik_fusion(Y,
                                                   X,
                                                   mkcd_starcube,
                                                   J,K,D);

  arma::mat SUM_il_X_ElogLik_collapsed = arma::sum(SUM_il_X_ElogLik, 2); // JxK

  arma::rowvec temp(K);

    arma::rowvec Q(K);
    for(int rr = 0; rr < reps; rr++){
      for(int j = 0; j < J; j++){
        //Rcpp::Rcout << "-7.3-";
        arma::uvec inds = NN_inds[j] ;
        //Rcpp::Rcout << "-7.4-" << Temp_R.rows(inds) << "\n";
        Q = arma::sum(Temp_R.rows(inds),0) * itemp; // 1 x K
        //Rcpp::Rcout << "-7.5-";
        temp = SUM_il_X_ElogLik_collapsed.row(j)  +
          Q;
        //Rcpp::Rcout << "-7.6-";
        Temp_R.row(j) =   (exp(temp - row_LogSumExp_cpp(temp))) ;
        //Rcpp::Rcout << "-7.7-";
      }
    }

  return(Temp_R);
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat upd_Vk_cpp(double const a_tilde0_DP,
                     double const b_tilde0_DP, // if conc parameter is random this needs to be equal to s1/s2 and not to alpha_0
                     arma::mat RHO){

  int K = RHO.n_cols;
  arma::rowvec M_k = arma::sum(RHO,0);

  arma::rowvec rev_cs_mk = reverse_cumsum_row(M_k);

  arma::rowvec a_tilde_vk      = M_k       + a_tilde0_DP;
  a_tilde_vk[K-1] = 1;
  arma::rowvec b_tilde_vk      = rev_cs_mk + b_tilde0_DP;
  b_tilde_vk[K-1] = 1e-10;

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
arma::rowvec upd_alpha_DP(arma::rowvec a_tilde_Vk,
                          arma::rowvec b_tilde_Vk,
                          arma::rowvec conc_hyper){
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
