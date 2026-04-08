hyperpar1 = list(
  a  = 1e-4,
  b  = 1e-4,
  m0 = 0,
  k0 = .01, # Changed
  c0 = 3,
  d0 = 2,
  s1 = 1,
  s2 = 1,
  s3 = 1,
  s4 = 1)

K <- 10
L <- 10
nrep <- 10 # starting points for Poseidon

K_upper <- 10
L_upper <- 20
Ls_choose <- c(10,20,30)
save_list <- c("i_sim", "i_met",
               "hyperpar1","K","L","nrep","K_upper","L_upper")
