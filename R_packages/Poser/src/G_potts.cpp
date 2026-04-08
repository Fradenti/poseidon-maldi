#include "A_AUX.h"
#include "G_potts.h"



// -----------------------------------------------------------------------------
// beta as fixed parameter

// [[Rcpp::export]]
double EQ66( arma::mat RHO, arma::field<arma::uvec> NN_inds){
  // from EQ 66 in Lu
  arma::mat aRHO = create_antiRHO(RHO, NN_inds);
  return(arma::accu(RHO % aRHO));
}

// much faster version ---------------------------------------------------------

// [[Rcpp::export]]
arma::mat p_MF( double itemp, arma::mat antiRHO, arma::rowvec log_e_pi){ // eq. 68
      arma::mat replog = arma::repelem(log_e_pi, antiRHO.n_rows, 1); // J x K
      return( exp(replog + itemp * antiRHO ) );

}


// [[Rcpp::export]]
arma::mat create_pMF(double beta,
                         arma::mat antiRHO,
                         arma::rowvec log_e_pi){
  size_t J = antiRHO.n_rows;
  size_t K = antiRHO.n_cols;
  arma::mat RHO_tilde(J,K);

  RHO_tilde = p_MF(beta,  antiRHO, log_e_pi);

  arma::colvec norm = arma::sum(RHO_tilde,1);
  arma::mat NORM = arma::repelem(norm, 1, RHO_tilde.n_cols);

  return(RHO_tilde % (1.0/NORM));
}



// [[Rcpp::export]]
double bisection_LU2(arma::mat RHO,
                    arma::field<arma::uvec> NN_inds,
                    arma::rowvec log_e_pi,
                    double low, double upp, int N_sims) {

  arma::mat antiRHO = create_antiRHO(RHO, NN_inds);

  double mid, fnew, flow;//
  double valEQ66 = EQ66(RHO,NN_inds);
  int check;

  mid = (low+upp) * 0.5;
  arma::mat RT_mid = create_pMF(mid, // beta middle
                              antiRHO,
                              log_e_pi);
  arma::mat RT_low = create_pMF(low, //beta low
                                antiRHO,
                                log_e_pi);

  fnew = valEQ66 - EQ66(RT_mid,NN_inds);
  flow = valEQ66 - EQ66(RT_low,NN_inds);

  if( R::sign(flow) * R::sign(fnew) < 0){
    upp = mid;
    check = -1;
  }else{
    low = mid;
    flow = fnew;
    check = +1;
  }

  for(int i=1; i < N_sims; i++){

    mid = (low+upp) * 0.5;
    //Rcpp::Rcout << mid << "\n";

    if(check == -1){
      arma::mat RT_x = create_pMF(mid,
                                  antiRHO,
                                  log_e_pi);
      fnew = valEQ66 - EQ66(RT_x,NN_inds);
    }else if(check == +1){
      arma::mat RT_x = create_pMF(mid,
                                  antiRHO,
                                  log_e_pi);
      fnew = valEQ66 - EQ66(RT_x,NN_inds);
    }

    if( R::sign(flow) * R::sign(fnew) < 0){
      upp = mid;
      check = -1;
    }else{
      low = mid;
      flow = fnew;
      check = +1;
    }
    //Rcpp::Rcout<< fnew << "---" << flow << "---" << o <<"\n";
  }
  return(mid);
}

// -----------------------------------------------------------------------------
// beta as RV

// [[Rcpp::export]]
double est_itemp_PL2(arma::mat RHO,
                    arma::mat antiRHO,
                    arma::vec itemp_grid){

  int n_itemp = itemp_grid.n_elem;

  double  RaR = arma::accu(RHO % antiRHO);

  arma::rowvec X(n_itemp);


    for(int b = 0; b < n_itemp; b++){
      arma::vec br =  arma::sum( exp(itemp_grid[b] * antiRHO),1);
      X(b) =  (itemp_grid[b] * RaR) -
              arma::accu(log(br));
    }


  arma::mat res =  exp(X - row_LogSumExp_cpp(X)) * itemp_grid ;
  //Rcpp::Rcout<<  exp(X - row_LogSumExp_cpp(X));
  return( res[0] );

}

// [[Rcpp::export]]
double est_itemp_PL2_withpi(arma::mat RHO,
                            arma::mat antiRHO,
                            arma::rowvec ElogPi,
                            arma::vec itemp_grid){

  int n_itemp = itemp_grid.n_elem;

  double  RaR = arma::accu(RHO % antiRHO);

  arma::rowvec X(n_itemp);
  arma::mat ELOGPI = arma::repelem(ElogPi, RHO.n_rows,1 );


  for(int b = 0; b < n_itemp; b++){
    arma::vec br =  arma::sum( exp(itemp_grid[b] * antiRHO + ELOGPI),1);
    X(b) =  (itemp_grid[b] * RaR) -
      arma::accu(log(br));
  }


  arma::mat res =  exp(X - row_LogSumExp_cpp(X)) * itemp_grid ;
  //Rcpp::Rcout<<  exp(X - row_LogSumExp_cpp(X));
  return( res[0] );

}











