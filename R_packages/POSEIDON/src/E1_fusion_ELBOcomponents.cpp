#include "E1_fusion_ELBOcomponents.h"


double elbo_p_THETA(arma::mat mkcd_0,
                    arma::mat mkcd_star){
  double  Z = -
    arma::accu( (mkcd_0.col(2) + 1.5) % E_log_IG(  mkcd_star.col(2),  mkcd_star.col(3)) ) -
    arma::accu( mkcd_0.col(3) % (  mkcd_star.col(2) /   mkcd_star.col(3)) ) -
    arma::accu( (mkcd_0.col(1) / 2) % (
        1 /   mkcd_star.col(1) + (  mkcd_star.col(2) /   mkcd_star.col(3)) %
        (mkcd_star.col(0)-mkcd_0.col(0)) % (mkcd_star.col(0)-mkcd_0.col(0))));
 // Normalizing constant of prior is omitted as it is fixed
  return(Z);
}

double elbo_q_THETA(arma::mat mkcd_star){

  double Z =
    arma::accu(
      mkcd_star.col(2) % log(mkcd_star.col(3)) -
      lgamma(mkcd_star.col(2)) +
      .5* log(mkcd_star.col(1)) -
  (mkcd_star.col(2) + 1.5) % E_log_IG(mkcd_star.col(2),mkcd_star.col(3)) -
  mkcd_star.col(2) );
  // -LT/2, -.5 log(2pi) omitted as they are fixed
  return(Z);
}


double elbo_diff_THETA_fusion(arma::cube mkcd0,
                              arma::cube mkcd_star_cube,
                              int const D){

  double Z = 0.0;

  for(int t=0; t<D; t++){
    Z += elbo_p_THETA(mkcd0.slice(t), mkcd_star_cube.slice(t)) -
         elbo_q_THETA(mkcd_star_cube.slice(t));
  }

  return(Z);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------



double elbo_diff_omega(arma::mat B0,
                       arma::mat B_star){
  int const L = B_star.n_rows;
  int const K = B_star.n_cols;
  arma::mat E_ln_omega(L,K);

  double Konsts = arma::accu(ldirichlet_normconst_col(B0) -
                             ldirichlet_normconst_col(B_star));
  for(int k=0; k<K; k++){
    E_ln_omega.col(k) = E_log_DIR_col(B_star.col(k)); // L x K
  }
  double G =  Konsts + arma::accu( E_ln_omega % (B0-B_star) );
  return(G);

}

double elbo_diff_omega_fusion(arma::cube B0,
                              arma::cube B_star,
                              int const D){

  double G = 0.0;
  for(int t=0; t<D; t++){
    G += elbo_diff_omega(B0.slice(t), B_star.slice(t));
  }
  return(G);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_diff_pi(arma::rowvec A0,
                    arma::rowvec A_star){

  double Konst1 = ldirichlet_normconst_row2double(A0);
  double Konst2 = ldirichlet_normconst_row2double(A_star);
  double G =  Konst1 - Konst2 +
              arma::accu( E_log_DIR_row(A_star) % (A0 - A_star) );
  return(G);

}


// -----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_diff_R2(arma::cube XI,
                    arma::mat B_star){

  int L = B_star.n_rows;
  int K = B_star.n_cols;
  arma::mat E_ln_omega(L,K), N_kl(K,L);

  double top_g = 0.0;
  for(int k=0; k<K; k++){
    E_ln_omega.col(k) = E_log_DIR_col(B_star.col(k)); // L x K
    arma::mat XX = XI.col(k); //N x L
    N_kl.row(k) = arma::sum(XX,0); //1 x L
    arma::umat ui = find(XX>0);
    top_g += arma::accu( XX.elem(ui) % log(XX.elem(ui)) );
  }

  return( arma::accu(N_kl % arma::trans(E_ln_omega) ) - top_g );
}

// [[Rcpp::export]]
double elbo_diff_R_fusion(arma::field<arma::cube> XI,
                          arma::cube beta_star_cube,
                          int const D){

  double res = 0.0;

  for(int t=0; t<D; t++){


    res += elbo_diff_R2(XI(t),
                        beta_star_cube.slice(t)) ;

    }
  return(res);
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_p_Y(arma::mat Y,
                arma::cube XI,
                arma::mat RHO,
                arma::mat mkcd_star){

  double res =0.0;
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

// [[Rcpp::export]]
double elbo_p_Y_fusion(arma::field<arma::mat> Y,
                       arma::field<arma::cube> XI,
                       arma::mat RHO,
                       arma::cube mkcd_star,
                       int const D){
  double res = 0.0;

  for(int t=0; t<D; t++){

    res += elbo_p_Y(Y(t),
                    XI(t),
                    RHO,
                    mkcd_star.slice(t));

  }
  return(res);
}

// ----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_p_C_fusion(arma::mat RHO,
                arma::field<arma::uvec> NN_inds,
                arma::rowvec E_log_pi_k,
                double itemp){

  int J = RHO.n_rows;
  int K = RHO.n_cols;
  arma::mat Q(J,K);
  double P_potts = 0.0;

  if(itemp > 0){
    for(int j=0; j<J; j++){
      arma::uvec inds = NN_inds[j] ;
      Q.row(j) = arma::sum(RHO.rows(inds),0);
    }
    P_potts = arma::accu( (Q * itemp) % RHO );
  }

  arma::rowvec M_k = arma::sum(RHO,0);
  // we are missing part coming from E(Z(beta, pi_k)), if weights and beta are random!
  return( P_potts + arma::accu(M_k % E_log_pi_k) );
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_p_C_onlyPi(arma::mat RHO,
                       arma::rowvec E_log_pi_k){

  int J = RHO.n_rows;
  int K = RHO.n_cols;
  arma::mat Q(J,K);

  arma::rowvec M_k = arma::sum(RHO,0);
  // we are missing part coming from E(Z(beta, pi_k)), if weights and beta are random!
  return( arma::accu(M_k % E_log_pi_k) );
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_p_C_onlyPotts(arma::mat RHO,
                          arma::mat antiRHO,
                          double beta_potts){

  double P_potts = 0.0;
  P_potts = arma::accu( (antiRHO ) % RHO ) * beta_potts;
  double Z = log_beta_normconts_PL(antiRHO,beta_potts);
  return( P_potts - Z );
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_q_C_fusion(arma::mat RHO){
  arma::umat ui = arma::find(RHO>0.0);
  double y = arma::accu(RHO.elem(ui) % log(RHO.elem(ui)));
  return(y);
}


// ----------------------------------------------------------------------------


double elbo_diff_C_fusion(arma::mat RHO,
                       arma::field<arma::uvec> NN_inds,
                       arma::rowvec E_log_pi_k,
                       double itemp){

  double y = elbo_p_C_fusion(RHO,
                             NN_inds,
                             E_log_pi_k,
                             itemp) -
              elbo_q_C_fusion(RHO);
  return(y);
}


double elbo_diff_C_fusion_onlyPotts(arma::mat RHO,
                                    arma::mat antiRHO,
                                    double itemp){

  double y = elbo_p_C_onlyPotts(RHO,
                                       antiRHO,
                                       itemp) -
                               elbo_q_C_fusion(RHO);
  return(y);
}

double elbo_diff_C_fusion_onlyPi(arma::mat RHO,
                                    arma::rowvec E_log_pi_k){

  double y = elbo_p_C_onlyPi(RHO,E_log_pi_k) -
                                  elbo_q_C_fusion(RHO);
  return(y);
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_p_v_CP(arma::rowvec a_tilde_k,
                   arma::rowvec b_tilde_k,
                   arma::rowvec conc_star,
                   int const K){

  a_tilde_k.shed_col(K-1);
  b_tilde_k.shed_col(K-1);

  arma::rowvec Y =
    E_log_beta_row(b_tilde_k, a_tilde_k) * ( conc_star[0]/conc_star[1] - 1.0 );

  return( (K-1)* ( R::digamma(conc_star[0]) - log(conc_star[1]) ) + arma::accu(Y));
}

double elbo_q_v_CP(arma::rowvec a_tilde_k,
                   arma::rowvec b_tilde_k,
                   int const K){
  a_tilde_k.shed_col(K-1);
  b_tilde_k.shed_col(K-1);


  arma::rowvec Y = lbeta_normconst_row(a_tilde_k,b_tilde_k) +
    E_log_beta_row(a_tilde_k,b_tilde_k) % (a_tilde_k - 1.0 )+
    E_log_beta_row(b_tilde_k,a_tilde_k) % (b_tilde_k - 1.0 );

  return(arma::accu(Y));
}


double elbo_diff_vk(arma::rowvec a_tilde_Vk,
                    arma::rowvec b_tilde_Vk,
                    arma::rowvec conc_star,
                    int const K){

  double x = elbo_p_v_CP(a_tilde_Vk,
                         b_tilde_Vk,
                         conc_star,
                         K) -
             elbo_q_v_CP(a_tilde_Vk,
                         b_tilde_Vk,
                         K);
  return(x);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double elbo_diff_alpha_DP(arma::rowvec conc_hyper,
                          arma::rowvec conc_star){
  double pa =

    ( conc_hyper[0] * log(conc_hyper[1]) - lgamma(conc_hyper[0]) ) -
    ( conc_star[0] * log(conc_star[1]) - lgamma(conc_star[0]) ) +
    ( conc_hyper[0] - conc_star[0] ) * (R::digamma(conc_star[0]) - log(conc_star[1])) -
      conc_star[0] * (conc_hyper[1] / conc_star[1] - 1.0);

  return( pa );
}
