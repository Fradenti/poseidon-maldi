sample_dahl <- function(sigma, lambda, alpha){
  n <- length(sigma)
  z <- rep(NA,n) # ordered with 1:n (not with sigma)
  q <- 0 # number of current clusters
  logpr <- 0
  for(t in 1:n){
    s <- sigma[t]
    if(t == 1){ # first element
      z[s] <- 1
      q <- q+1
    } else { # next
      # here I will save the probabilities to join one of the existing q clusters
      # or a new cluster.
      pr <- numeric(q+1) 
      pr[1:q] <- (t-1) / (alpha + t - 1)
      pr[q+1] <- alpha / (alpha + t - 1)
      # this is the normalizing constant for the lambda
      const <- sum(lambda[s,sigma[1:(t-1)]])
      # for each j I have some elements in that cluster and I compute the lambda
      for(j in 1:q){
        ind_j <- which(z == j)
        pr[j] <- pr[j] * sum(lambda[s,ind_j])/const
      }
      # now I need to sample
      z[s] <- sample(q+1,1, prob = pr)
      if(z[s] == q+1){
        q <- q+1
      }
      logpr <- logpr + log(pr[z[s]]) - log(sum(pr))
    }
  }  
  return(c(z,logpr))
}

sample_potts <- function(sigma, alpha, beta, w.sym){
  n <- length(sigma)
  z <- rep(NA,n) # ordered with 1:n (not with sigma)
  q <- 0 # number of current clusters
  logp <- 0
  for(t in 1:n){
    s <- sigma[t]
    if(t == 1){ # first element
      z[s] <- 1
      q <- q+1
    } else { # next
      pr <- numeric(q+1)
      pr[q+1] <- alpha / (alpha + t - 1)
      for(j in 1:q){
        ind_j <- which(z == j)
        mk <- sum(w.sym[s, ind_j])
        pr[j] <- length(ind_j)/(alpha + t - 1) * exp(beta * mk)
      }
      pr <- pr/sum(pr) # I don't normalize so I can compare with prob_Orbanz
      # now I need to sample
      z[s] <- sample(q+1,1, prob = pr)
      if(z[s] == q+1){
        q <- q+1
      }
      logp <- logp + log(pr[z[s]])
    }
  }
  return(c(z,logp))
}

withincl_z <- function(z, w){
  wc <- 0
  for(x in unique(z)){
    ind <- which(z == x)
    wc <- wc + sum(w[ind,ind])/2
  }
  wc
}
