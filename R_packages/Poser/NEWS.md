# Poser 0.1.0

* In `A_AUX.cpp`
`arma::rowvec log_EXP_pi(arma::rowvec a_k, arma::rowvec b_k){ int K = a_k.n_rows;`
contained an harmless bug. Now in
`arma::rowvec log_EXP_pi(arma::rowvec a_k, arma::rowvec b_k){`
I remove allocating `K`
