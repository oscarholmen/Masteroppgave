

##########helper functions##########
f_gamma_l <- function(x, alpha, lambda){# calculates the density of the gamma distribution
  return(dgamma(x, shape = alpha, rate = lambda, log = TRUE))
}

S_gamma_l <- function(x, alpha, lambda){# calculates the survival function of the gamma distribution
  return(pgamma(x, shape = alpha, rate = lambda, lower.tail = FALSE, log.p = TRUE))
}

residual_time_pdf_dist_l <- function(x, alpha, lambda){ #calculates the residual time distribution of the gamma distribution
  return(pgamma(x, shape = alpha, rate = lambda, lower.tail = FALSE, log.p = TRUE) - log(alpha/lambda))
}



prev_time <- function(s, i ,j, y){ #returns the time of the previous event for process s at state (i,j)
  p1_prev <- s[i, j, 1]
  p2_prev <- s[i, j, 2]
  p1_prev_time <- ifelse(p1_prev == 0, 0, y[p1_prev])
  p2_prev_time <- ifelse(p2_prev == 0, 0, y[p2_prev])
  return(list(p1 = p1_prev_time, p2 = p2_prev_time))
}

make_new_arrays <- function(p_ij_arr, s_ij_arr, timestep){#creates new arrays for probabilities and states
  #dimensions of p_ij_arr is c(2, len_states)
  #dimensions of s_ij_arr is c(2, len_states, 2)
  len_states <- length(s_ij_arr[1,,1])
  
  
  #Step for updating possible states
  dim_new_s <- c(2, len_states + 1, 2) #Adding one more column for the new state
  
  new_s_ij_arr <- array(0, dim = dim_new_s)
  for (i in 1:len_states){
    new_s_ij_arr[1,i,1] <- s_ij_arr[1,i,1] + 1
    new_s_ij_arr[2,i,1] <- s_ij_arr[2,i,1]
  }
  new_s_ij_arr[1,len_states + 1,1] <- s_ij_arr[2,len_states,1] + 2
  new_s_ij_arr[1,len_states + 1,2] <- s_ij_arr[2,len_states,2]
  
  
  for (i in 1:len_states){
    new_s_ij_arr[2,i,2] <- s_ij_arr[2,i,2] + 1
    new_s_ij_arr[1,i,2] <- s_ij_arr[1,i,2]
  }
  new_s_ij_arr[2,len_states + 1,2] <- s_ij_arr[1,len_states,2] + 2
  new_s_ij_arr[2,len_states + 1,1] <- s_ij_arr[1,len_states,1]
  
  #Step for updating probabilities
  dim_new_p <- c(2, len_states + 1)
  new_p_ij_arr <- array(0, dim = dim_new_p)
  
  return(list(arr1 = new_p_ij_arr, arr2 = new_s_ij_arr))
}


#Prune states with very low probabilities
prune_states <- function(p_mat, s_arr, logthresh) { 
  cols <- apply(p_mat, 2, function(x) any(x > logthresh)) #get columns where any probability is above the threshold
  keep <- (cols)
  if (!any(keep)) { keep[length(keep)] <- TRUE } #ensure at least one state is kept
  new_p <- p_mat[, keep, drop = FALSE]
  new_s <- s_arr[, keep, , drop = FALSE]
  return(list(p = new_p, s = new_s)) 
  }

#log-sum-exp trick for numerical stability
logsumexp <- function(a) {
  m <- max(a)
  return(m + log(sum(exp(a - m))))
}


############################

# Kalman Filter Implementation
kalman_filter_method <- function(y, alpha1, lambda1, alpha2, lambda2, tmax, thresh = 1e-10){
  
  n <- length(y)
  logthresh <- log(thresh)
  
  #Initialize variables
  s_ij_arr <- array(c(0,0,0,0), dim = c(2,1,2))  #list of all possible last states
  p_ij_arr <- array(c(0,0), dim = c(2,1)) #initial probabilities in log scale
  loglik = 0 #initial likelihood
  
  for (t in 1:n) {
    y_t <- y[t]
    if(t == 1){#special case for the first event, since the initial s_ij_arr has state (0,0) in both rows
      prev1 <- prev_time(s_ij_arr, 1, 1, y)$p1
      prev2 <- prev_time(s_ij_arr, 1, 1, y)$p2

      delta1 <- y_t - prev1
      delta2 <- y_t - prev2
      log_p1 <- p_ij_arr[1,1] + residual_time_pdf_dist_l(delta1, alpha1, lambda1) + S_gamma_l(delta2, alpha2, lambda2)
      p_ij_arr[1,1] <- log_p1
      s_ij_arr[1,1,1] <- s_ij_arr[1,1,1] + 1

      prev1 <- prev_time(s_ij_arr, 2, 1, y)$p1
      prev2 <- prev_time(s_ij_arr, 2, 1, y)$p2

      delta1 <- y_t - prev1
      delta2 <- y_t - prev2
      log_p2 <- p_ij_arr[2,1] + residual_time_pdf_dist_l(delta2, alpha2, lambda2) + S_gamma_l(delta1, alpha1, lambda1)
      p_ij_arr[2,1] <- log_p2
      s_ij_arr[2,1,2] <- s_ij_arr[2,1,2] + 1
      loglik <- loglik + log_p1 + log_p2
      next #move to next iteration. From now the updates are the same for all t
    }
    
    #create new arrays for probabilities and states
    new_arrays <- make_new_arrays(p_ij_arr, s_ij_arr, timestep = t)
    new_p_ij_arr <- new_arrays$arr1
    new_s_ij_arr <- new_arrays$arr2
    
    #update probabilities
    #row by row 
    for(i in 1:2){
      #starting with all but the last of each row, since they have a different update rule
      for(j in 1:(length(new_p_ij_arr[i,])-1)){
        prev1 <- prev_time(s_ij_arr, i, j, y)$p1
        prev2 <- prev_time(s_ij_arr, i, j, y)$p2
          
        delta1 <- y_t - prev1
        delta2 <- y_t - prev2
          
        if(i == 1){ #if the event is from process 1
          log_p <- p_ij_arr[i,j] + f_gamma_l(delta1, alpha1, lambda1) + S_gamma_l(delta2, alpha2, lambda2)
          new_p_ij_arr[i,j] <- log_p
        } else { #if the event is from process 2
          log_p <- p_ij_arr[i,j] + f_gamma_l(delta2, alpha2, lambda2) + S_gamma_l(delta1, alpha1, lambda1)
          new_p_ij_arr[i,j] <- log_p
        }
      }
    
      #update the last of each row
      terms_vec <- numeric(length(p_ij_arr[i,]))
      for(j in (1:length(p_ij_arr[i,]))){
        prev1 <- prev_time(s_ij_arr, i, j, y)$p1
        prev2 <- prev_time(s_ij_arr, i, j, y)$p2
        if(i == 1){
          terms_vec[j] <- p_ij_arr[i,j] + f_gamma_l(y_t - prev2, alpha2, lambda2) + S_gamma_l(y_t - prev1, alpha1, lambda1)}
        else {
          terms_vec[j] <- p_ij_arr[i,j] + f_gamma_l(y_t - prev1, alpha1, lambda1) + S_gamma_l(y_t - prev2, alpha2, lambda2)
        }
      }
      if(i == 1){
        new_p_ij_arr[i+1,length(new_p_ij_arr[i,])] <- logsumexp(terms_vec)
      } else {
        new_p_ij_arr[i-1,length(new_p_ij_arr[i,])] <- logsumexp(terms_vec)
      }
    }

    normalizing_log <- logsumexp(as.vector(new_p_ij_arr))
    norm_p_ij_arr <- new_p_ij_arr - normalizing_log
    loglik <- loglik + normalizing_log
    #Prune states with very low probabilities
    pr <- prune_states(norm_p_ij_arr, new_s_ij_arr, logthresh)
    p_ij_arr <- pr$p
    s_ij_arr <- pr$s
    

  }  
  return(loglik)

}
#############

library(MCMCpack)
metropfunc <- function(params, timeseries) {
  if (any(params <= 0)) return(-Inf)
  
  likelihood <- kalman_filter_method(
    y = timeseries$sup_vec,
    alpha1 = params[1],
    lambda1 = params[2],
    alpha2 = params[3],
    lambda2 = params[4],
    tmax = timeseries$tmax,
    thresh = 1e-10
  )
  
  #specify priors
  log_prior_alpha1 <- dexp(params[1], rate = 1/4, log = TRUE)
  log_prior_lambda1 <- dgamma(params[2], shape = 0.01, rate = 0.01, log = TRUE)
  log_prior_alpha2 <- dexp(params[3], rate = 1/10, log = TRUE)
  log_prior_lambda2 <- dgamma(params[4], shape = 0.04, rate = 0.01, log = TRUE)
  
  logpost <- likelihood  + log_prior_lambda1 +  log_prior_lambda2 + log_prior_alpha1 + log_prior_alpha2
  return(logpost)
}

metropfunc_alt_prior <- function(params, timeseries) {
  if (any(params <= 0)) return(-Inf)
  
  likelihood <- kalman_filter_method(
    y = timeseries$sup_vec,
    alpha1 = params[1],
    lambda1 = params[2],
    alpha2 = params[3],
    lambda2 = params[4],
    tmax = timeseries$tmax,
    thresh = 1e-10
  )
  
  # params = c(alpha1, lambda, alpha2, w)
  lambda  <- sqrt(params[2]*params[4])  # lambda > 0
  w      <- sqrt(params[2]/params[4])   # w > 0
  
  # log-priors:
  log_prior_lambda <- -log(lambda)                 # improper Jeffreys 1/lambda
  log_prior_w      <- dlnorm(w, meanlog = log(4), sdlog = 0.5, log = TRUE)
  log_prior_alpha1 <- dexp(params[1], rate = 1/4, log = TRUE)
  log_prior_alpha2 <- dexp(params[3], rate = 1/10, log = TRUE)
  
  logpost <- likelihood + log_prior_lambda + log_prior_w + log_prior_alpha1 + log_prior_alpha2
  return(logpost)
}


