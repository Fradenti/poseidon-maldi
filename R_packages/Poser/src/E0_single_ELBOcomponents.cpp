
#include "E0_single_ELBOcomponents.h"

double elbo_p_THETA(const arma::mat& mkcd_0,
                    const arma::mat& mkcd_star){

  double  Z = -
    arma::accu( (mkcd_0.col(2) + 1.5) % E_log_IG(  mkcd_star.col(2),  mkcd_star.col(3)) ) -
    arma::accu( mkcd_0.col(3) % (  mkcd_star.col(2) /   mkcd_star.col(3)) ) -
    arma::accu( (mkcd_0.col(1) / 2) % (
        1.0 /   mkcd_star.col(1) + (  mkcd_star.col(2) /   mkcd_star.col(3)) %
        (mkcd_star.col(0)-mkcd_0.col(0)) % (mkcd_star.col(0)-mkcd_0.col(0))));

  return(Z);
}

double elbo_q_THETA(const arma::mat& mkcd_star){

  double Z =
    arma::accu(
        mkcd_star.col(2) % log(mkcd_star.col(3)) -
        lgamma(mkcd_star.col(2)) +
        .5* log(mkcd_star.col(1)) -
        (mkcd_star.col(2) + 1.5) % E_log_IG(mkcd_star.col(2),mkcd_star.col(3)) -
        mkcd_star.col(2) );

  return(Z);
}

// ----------------------------------------------------------------------------

double elbo_diff_omega(arma::mat B0,
                       arma::mat B_star){
  int L = B_star.n_rows;
  int K = B_star.n_cols;
  arma::mat E_ln_omega(L,K);
  double Konsts = arma::accu(ldirichlet_normconst_col(B0) -
                             ldirichlet_normconst_col(B_star));
  for(int k=0; k<K; k++){
    E_ln_omega.col(k) = E_log_DIR_col(B_star.col(k)); // L x K
  }
  double G =  Konsts + arma::accu( E_ln_omega % (B0-B_star) );
  return(G);

}

// ----------------------------------------------------------------------------

double elbo_diff_R(arma::cube XI,
                arma::mat B_star){

  int L = B_star.n_rows;
  int K = B_star.n_cols;
  arma::mat E_ln_omega(L,K), N_kl(K,L);

  for(int k=0; k<K; k++){
    E_ln_omega.col(k) = E_log_DIR_col(B_star.col(k)); // L x K
    arma::mat XX = XI.col(k);
    N_kl.row(k) = arma::sum(XX,0); //K x L
  }


  return(arma::accu(N_kl % arma::trans(E_ln_omega))-
         arma::accu(XI % log(XI + 1e-12)));
}


double elbo_diff_R2(arma::cube XI,
                   arma::mat E_ln_omega){

  int L = E_ln_omega.n_rows;
  int K = E_ln_omega.n_cols;
  arma::mat N_kl(K,L);


  double top_g = 0.0;
  for(int k=0; k<K; k++){
  //  E_ln_omega.col(k) = E_log_DIR_col(B_star.col(k)); // L x K
    arma::mat XX = XI.col(k);
    N_kl.row(k) = arma::sum(XX,0); //K x L
    arma::umat ui = find(XX>0);
    top_g += arma::accu(XX.elem(ui) % log(XX.elem(ui)));
  }


  return(arma::accu(N_kl % arma::trans(E_ln_omega)) - top_g);
}

// ----------------------------------------------------------------------------

double elbo_diff_pi(arma::rowvec A0,
                    arma::rowvec A_star){

  double Konst1 = ldirichlet_normconst_row2double(A0);
  double Konst2 = ldirichlet_normconst_row2double(A_star);
  double G =  Konst1 -
              Konst2 +
              arma::accu( E_log_DIR_row(A_star) % (A0 - A_star) );
  return(G);

}

// ----------------------------------------------------------------------------

double elbo_p_Y(arma::mat Y,
                arma::cube XI,
                arma::mat RHO,
                arma::mat mkcd_star){

  double res = 0.0;
  int N = Y.n_rows;
  int L = mkcd_star.n_rows;
  int J = Y.n_cols;
  arma::mat M_jl_fixed_i(J,L);
  arma::mat ell_lj_fixed_i(L,J);

  for(int i=0; i<N; i++){

    for(int l=0; l<L; l++){
      ell_lj_fixed_i.row(l) = E_log_f_row(Y.row(i),
                         mkcd_star(l,0),mkcd_star(l,1),
                         mkcd_star(l,2),mkcd_star(l,3));
    }

    arma::mat XI_fixed_i = XI.row(i); // K x L
    M_jl_fixed_i = RHO * XI_fixed_i;
    res += arma::accu(arma::trans(M_jl_fixed_i) % ell_lj_fixed_i);

  }

  return(res);
}

// ----------------------------------------------------------------------------

double elbo_diff_alpha_DP(arma::rowvec conc_hyper,
                          arma::rowvec conc_star){
  double pa =

    ( conc_hyper[0] * log(conc_hyper[1]) - lgamma(conc_hyper[0]) ) -
    ( conc_star[0] * log(conc_star[1]) - lgamma(conc_star[0]) ) +
    ( conc_hyper[0] - conc_star[0] ) * (R::digamma(conc_star[0]) - log(conc_star[1])) -
    conc_star[0] * (conc_hyper[1] / conc_star[1] - 1);

  return( pa );
}

double elbo_diff_beta_DP(arma::rowvec conc_hyper,
                          arma::rowvec conc_star){
  return( elbo_diff_alpha_DP(conc_hyper,conc_star) );
}


// ----------------------------------------------------------------------------

double elbo_diff_vk(arma::rowvec a_tilde_Vk,
                    arma::rowvec b_tilde_Vk,
                    arma::rowvec conc_star){

  int K = a_tilde_Vk.n_cols;

  a_tilde_Vk.shed_col(K-1);
  b_tilde_Vk.shed_col(K-1);

  arma::rowvec Y1 = E_log_beta_row(a_tilde_Vk,b_tilde_Vk);
  arma::rowvec Y2 = E_log_beta_row(b_tilde_Vk,a_tilde_Vk);

  return(
         (K-1)  * ( R::digamma(conc_star[0]) - log(conc_star[1]) ) -
         arma::accu(lbeta_normconst_row(a_tilde_Vk,b_tilde_Vk)) -
         arma::accu(Y1 % (a_tilde_Vk - 1.0 )) +
         arma::accu(Y2 % (conc_star[0]/conc_star[1] - b_tilde_Vk ))
           );
}


// -----------------------------------------------------------------------------

double elbo_p_C(arma::mat RHO,
                 arma::field<arma::uvec> NN_inds,
                 arma::rowvec E_log_pi_k,
                 double beta_potts){

  int J = RHO.n_rows;
  int K = RHO.n_cols;
  arma::mat Q(J,K, arma::fill::zeros);
  double P_potts = 0.0;

  if(beta_potts > 0.0){
    for(int j=0; j<J; j++){
        arma::uvec inds = NN_inds[j] ;
        Q.row(j) = arma::sum(RHO.rows(inds),0);
      }
    P_potts = arma::accu( (Q * beta_potts) % RHO );
  }

  arma::rowvec M_k = arma::sum(RHO,0);
  // we are missing part coming from E(Z(beta, pi_k)), if weights and beta are random!
  return( P_potts + arma::accu(M_k % E_log_pi_k) );
}

double elbo_p_C_onlyPi(arma::mat RHO,
                      arma::rowvec E_log_pi_k){

  int J = RHO.n_rows;
  int K = RHO.n_cols;
  arma::mat Q(J,K);

  arma::rowvec M_k = arma::sum(RHO,0);
  return( arma::accu(M_k % E_log_pi_k) );
}

// [[Rcpp::export]]
double elbo_p_C_onlyPotts(arma::mat RHO,
                arma::mat antiRHO,
                double beta_potts){

  double P_potts = 0.0;
  P_potts = arma::accu( (antiRHO ) % RHO ) * beta_potts;
  double Z = log_beta_normconts_PL(antiRHO,beta_potts);
  return( P_potts - Z );
}


/*
double elbo_q_C(arma::mat RHO){
  double y = arma::accu(RHO % log(RHO + 1e-12));
  return(y);
}
*/

// [[Rcpp::export]]
double elbo_q_C_v2(arma::mat RHO){
  arma::umat ui = arma::find(RHO>0);
  double y = arma::accu(RHO.elem(ui) % log(RHO.elem(ui)));
  return(y);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double elbo_diff_U(arma::mat a_ulk_bar,
                   arma::mat b_ulk_bar,
                   arma::rowvec conc_obser_star){

  int L = a_ulk_bar.n_rows;
  int K = a_ulk_bar.n_cols;

  a_ulk_bar.shed_row(L-1);
  b_ulk_bar.shed_row(L-1);

  arma::colvec PP(K);
  arma::colvec QQ(K);
  double b_bar = conc_obser_star[0]/conc_obser_star[1];

  double C_N = (L-1) * K * ( R::digamma(conc_obser_star[0]) - log(conc_obser_star[1]) );
  double C_D  = arma::accu(lbeta_normconst_mat_cpp(a_ulk_bar,b_ulk_bar));

  for(int k = 0; k < K; k ++){

    PP[k] =  //arma::accu( (a_bar-1) * E_log_beta(a_ulk_bar.col(k),b_ulk_bar.col(k))) +
      arma::accu(  E_log_beta_col(b_ulk_bar.col(k),a_ulk_bar.col(k)));
    QQ[k] =  arma::accu( (a_ulk_bar.col(k)-1) % E_log_beta_col(a_ulk_bar.col(k),b_ulk_bar.col(k))) +
      arma::accu( (b_ulk_bar.col(k)-1) % E_log_beta_col(b_ulk_bar.col(k),a_ulk_bar.col(k)));

  }

  return(C_N + (b_bar - 1) * arma::accu(PP) - C_D - arma::accu(QQ));

}


// [[Rcpp::export]]
double elbo_diff_itemp(arma::mat RHO, arma::mat antiRHO, double itemp){

  double Z = + log_beta_normconts_PL(antiRHO,itemp) -
         arma::accu( (antiRHO ) % RHO )* itemp;

  return(Z);
}
