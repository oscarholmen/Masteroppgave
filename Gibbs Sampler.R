library(coda)

#Create Gibbs Sampler to separate data from superposition of two renewal processes
#Data is generated from a superposition of two renewal processes
#The two renewal processes are generated from Gamma distributions
#The parameters of the Gamma distributions are unknown
#The goal is to estimate the parameters of the two distributions

gibbs_superpos <- function(y, tmax, method = "one_dist",
                           mcmc_steps = 1e3, burnin = 1e2, d1 = 0.2, d2 = 0.2, 
                           init_guess_params = c(1, 1, 1, 1), p = 0.5, prior_lambda = c(1,1,1,1),
                           prior_alpha =c(1,1), block_size = 3, thinning = 5, CHECK_K_CORR = FALSE){
  
  total_steps = mcmc_steps + burnin
  chain_length = total_steps/thinning
  n = length(y)
  
  #initialize latent parameters k_i bearnoulli random variables for each observation
  k = rbinom(n, 1, 0.5) + 1
  
  
  #create storage for chain
  if(CHECK_K_CORR == TRUE){chain = matrix(NA, nrow = chain_length, ncol = n)}
  if(CHECK_K_CORR == FALSE){chain = matrix(NA, nrow = chain_length, ncol = (4 + 2 + 2))}
  
  
  alpha1 = init_guess_params[1];lambda1 = init_guess_params[2]
  alpha2 = init_guess_params[3];lambda2 = init_guess_params[4]
  
  a1 = prior_lambda[1]; b1 = prior_lambda[2]; a2 = prior_lambda[3]; b2 = prior_lambda[4]
  c1 = prior_alpha[1]; c2 = prior_alpha[2]
  
  #Create fictitious values at the beginning and end of y
  
  
  for (i in 1:total_steps) {
    if(i %% 1000 == 0) {print(i)}
    #####Create a vector for each process######
    y_p1 <- y[k == 1]
    y_p2 <- y[k == 2]
    ###############################


    ######Simulate latent variables after each process#######
    y_pluss_1 <- sim_after(y_p1, alpha1, lambda1, tmax)
    y_pluss_2 <- sim_after(y_p2, alpha2, lambda2, tmax)
    y_mod <- c(y, y_pluss_1, y_pluss_2)
    k_mod <- c(k, 1, 2)
    ###############################
    
    #######Simulate latent variables before each process#######
    y_minus_1 <- sim_before(y_p1, alpha1, lambda1, tmax)
    y_minus_2 <- sim_before(y_p2, alpha2, lambda2, tmax)
    y_mod <- c(y_minus_1, y_minus_2, y_mod)
    k_mod <- c(1, 2, k_mod)
    ###############################
    
    
    ######Gibbs step for k_i#######
    k_mod <- gibbs_step_k(y = y_mod, k = k_mod, alpha1 = alpha1, lambda1 = lambda1, alpha2 = alpha2, lambda2 = lambda2, block_size = block_size)
    ###############################
    
    
    #########Create y1 and y2 from k########
    y1 <- y_mod[k_mod == 1]
    y2 <- y_mod[k_mod == 2]
    ###############################
    
    
    ######Gibbs step for alpha, lambda#######
    switch(method, 
           "one_dist" = {
             process1 = gibbs_step_one_dist(y1, alpha = alpha1, lambda = lambda1, a = a1, b = b1, c = c1, d = d1)
             alpha1 = process1[1]
             lambda1 = process1[2]
             process2 = gibbs_step_one_dist(y2, alpha = alpha2, lambda = lambda2, a= a2, b = b2,c = c2, d = d2)
             alpha2 = process2[1]
             lambda2 = process2[2]
           },
           "block" = {
             process1 = gibbs_step_block(y1, alpha = alpha1, lambda = lambda1, a = a1, b = b1, c = c1, d = d1)
             alpha1 = process1[1]
             lambda1 = process1[2]
             process2 = gibbs_step_block(y2, alpha = alpha2, lambda = lambda2, a= a2, b = b2,c = c2, d = d2)
             alpha2 = process2[1]
             lambda2 = process2[2]
           },
           "reparam" = {
             process = gibbs_step_reparam(y1, y2, alpha1, lambda1, alpha2, lambda2, a1= a1 , b1 = b1, a2= a2 , b2 = b2, priors = prior_alpha, sd = c(d1, d2))
             alpha1 = process[1]
             lambda1 = process[2]
             alpha2 = process[3]
             lambda2 = process[4]
           },
           "block_reparam" = {
             process = gibbs_step_block_reparam(y1, y2, alpha1, lambda1, alpha2, lambda2, a1= a1 , b1 = b1, a2= a2 , b2 = b2, priors = prior_alpha, sd = c(d1, d2))
             alpha1 = process[1]
             lambda1 = process[2]
             alpha2 = process[3]
             lambda2 = process[4]
           })
    
    ###############################
    
    
    ########Switching place s.t. mean_1 < mean_2########
    if (alpha1/lambda1 > alpha2/lambda2) {
      temp_alpha = alpha1
      alpha1 = alpha2
      alpha2 = temp_alpha
      
      temp_lambda = lambda1
      lambda1 = lambda2
      lambda2 = temp_lambda
      
      temp_a = a1
      a1 = a2
      a2 = temp_a
      
      temp_b = b1
      b1 = b2
      b2 = temp_b
      #swith k_i
      index = which(k_mod == 1)
      k_mod[index] = 2
      k_mod[-index] = 1
    }
    ###############################
    
    ######Remove two first and last values of k_mod#######
    k_chain <- k_mod[3:(length(k_mod)-2)]
    k <- k_mod[3:(length(k_mod)-2)]
    ###############################
    
    ##Count number of k = 1 and k = 2###
    k1 <- sum(k == 1)
    k2 <- sum(k == 2)
    
    #########Store every 5th values in chain########
    if (i %% thinning == 0){
      if(CHECK_K_CORR == FALSE) {
        chain[i/thinning, ] <- c(alpha1, lambda1, alpha2, lambda2, k1, k2, k[100], k[200])}} #can also include: k_chain
    
    if(CHECK_K_CORR == TRUE){chain[i/thinning, ] <- c(k)}
    
  }
  chain  <- mcmc(chain[-(1:(burnin/thinning)),])
  #get number of k=1 and k=2 at last iteration
  k1 <- sum(k_chain == 1)
  k2 <- sum(k_chain == 2)
  print(c(k1, k2))
  return(chain)
}



#########Helper Functions############
#########Gibbs Step alpha, lambda########
gibbs_step_one_dist <- function(y, alpha = 1, lambda = 1, a=0 , b = 0, c = 1, d = 0.2){
  n <- length(y)
  P <- prod(diff(y))
  S <- (y[n] - y[1])
  
  #update lambda
  repeat {
    lambda <- rgamma(1, shape = alpha * (n - 1) + a, rate = b + S)
    if (lambda > 0.1) break
  }
  
  #update alpha using random walk
  u <- runif(1)
  alpha_new <- rnorm(1, mean = alpha, d)
  if (alpha_new > 0){
    #acc_crit <- min(1, exp(1/c*(alpha-alpha_new)) * lambda**((n-1)*(alpha_new-alpha))
    #                *((gamma(alpha)/gamma(alpha_new))**(n-1))
    #                  *prod(diff(y))**(alpha_new-alpha))
    
    acc_crit <- min(1, exp(1/c*(alpha-alpha_new) +((n-1)*(alpha_new-alpha))*log(lambda) +
                             (n-1) *lgamma(alpha) - (n-1) * lgamma(alpha_new) +
                             (alpha_new-alpha) *log(P)))
    
    if(u < acc_crit){
      alpha <- alpha_new
    }
  }
  
  return(c(alpha, lambda))
} 

#Gibbs step Block
gibbs_step_block <- function(y, alpha = 1, lambda = 1, a=0 , b = 0, c = 1, d = 0.2){
  n <- length(y)
  P <- prod(diff(y))
  S <- (y[n] - y[1])
  
  #update lambda
  
  #update alpha using random walk
  
  u <- runif(1)
  
  alpha_new <- rnorm(1, mean = alpha, d)
  if(alpha_new <= 0){
    acc_crit <- 0
  }
  else{
    lambda_new <- rgamma(1, shape = alpha*(n-1) + a , rate = b + S)
    acc_crit <- min(1, exp((alpha_new - alpha)*(-1/c + sum(log(diff(y)))-(n-1)*log(S)) 
                           - (n-1)*(lgamma(alpha_new) - lgamma(alpha)) + a*(log(lambda_new)-log(lambda)) 
                           - b*(lambda_new - lambda)
                           +lgamma(alpha_new*(n-1) + 1) - lgamma(alpha*(n-1) + 1)))
    if(u < acc_crit && lambda_new > 0.1){
      alpha <- alpha_new
      lambda <- lambda_new
    }
  }
  
  return(c(alpha, lambda))
} 

#Gibbs step Reparametrization
gibbs_step_reparam <- function(y1, y2, alpha1, lambda1, alpha2, lambda2, a1=0 , b1 = 0, a2=0 , b2 = 0, priors = c(2,1), sd = c(0.2, 0.3)){
  n1 <- length(y1)
  P1 <- prod(diff(y1))
  S1 <- (y1[n1] - y1[1])
  n2 <- length(y2)
  P2 <- prod(diff(y2))
  S2 <- (y2[n2] - y2[1])
  
  #update alpha using random walk
  u1 <- runif(1)
  u2 <- runif(1)
  
  alpha1_new <- rnorm(1, mean = alpha1, sd[1])
  alpha2_new <- rnorm(1, mean = alpha2, sd[2])
  if (alpha1_new > 0){
    acc_crit1 <- min(1, exp(1/priors[1]*(alpha1-alpha1_new) +((n1-1)*(alpha1_new-alpha1))*log(lambda1) +
                             (n1-1) *lgamma(alpha1) - (n1-1) * lgamma(alpha1_new) +
                             (alpha1_new-alpha1) *log(P1)))

    if(u1 < acc_crit1){
      alpha1 <- alpha1_new
    }
  }
  if (alpha2_new > 0){
    acc_crit2 <- min(1, exp(1/priors[2]*(alpha2-alpha2_new) +((n2-1)*(alpha2_new-alpha2))*log(lambda2) +
                             (n2-1) *lgamma(alpha2) - (n2-1) * lgamma(alpha2_new) +
                             (alpha2_new-alpha2) *log(P2)))

    if(u2 < acc_crit2){
      alpha2 <- alpha2_new
    }
  }
  #update lambdas
  lambda_old <- sqrt(lambda1 * lambda2)
  w_old <- sqrt(lambda1 / lambda2)
  
  lambda_prop <- rgamma(1, shape = (alpha1 * (n1 - 1) + alpha2 * (n2 - 1)-1), rate = w_old*S1 + (1/w_old)*S2)
  lambda_old <- lambda_prop
  w_prop <- rnorm(1, mean = w_old, sd = 2)
  acc_crit_w <- min(1, (w_prop/w_old)*(alpha1*n1 - alpha2*n2) * exp((lambda_prop*(S1*(w_old-w_prop) + S2*(1/w_old - 1/w_prop)) )))
  if(runif(1) < acc_crit_w){
    w_old <- w_prop
    
  }
  lambda1_new <- lambda_prop * w_old
  lambda2_new <- lambda_prop / w_old
  if(lambda1_new > 0.1 && lambda2_new > 0.1){
    lambda1 <- lambda1_new
    lambda2 <- lambda2_new
  }
  
  return(c(alpha1, lambda1, alpha2, lambda2))
}


gibbs_step_block_reparam <- function(y1, y2, alpha1, lambda1, alpha2, lambda2, a1=0 , b1 = 0, a2=0 , b2 = 0, priors = c(2,1), sd = c(0.2, 0.3)){
  n1 <- length(y1)
  P1 <- prod(diff(y1))
  S1 <- (y1[n1] - y1[1])
  c1 <- priors[1]
  n2 <- length(y2)
  P2 <- prod(diff(y2))
  S2 <- (y2[n2] - y2[1])
  c2 <- priors[2]
  
  #update alpha using random walk
  u1 <- runif(1)
  u2 <- runif(1)
  u3 <- runif(1)
  u4 <- runif(1)
  
  alpha1_new <- rnorm(1, mean = alpha1, sd[1])
  alpha2_new <- rnorm(1, mean = alpha2, sd[2])
  #update lambdas
  lambda_old <- sqrt(lambda1 * lambda2)
  w_old <- sqrt(lambda1 / lambda2)
  
  lambda_prop <- rgamma(1, shape = (alpha1 * (n1 - 1) + alpha2 * (n2 - 1)-1), rate = w_old*S1 + (1/w_old)*S2)
  lambda_old <- lambda_prop
  w_prop <- rnorm(1, mean = w_old, sd = 2)
  acc_crit_w <- min(1, (w_prop/w_old)*(alpha1*n1 - alpha2*n2) * exp((lambda_prop*(S1*(w_old-w_prop) + S2*(1/w_old - 1/w_prop)) )))
  if(runif(1) < acc_crit_w){
    w_old <- w_prop
    
  }
  lambda1_new <- lambda_prop * w_old
  lambda2_new <- lambda_prop / w_old
  if(lambda1_new > 0.1 && lambda2_new > 0.1){
    acc_crit1 <- min(1, exp((alpha1_new - alpha1)*(-1/c1 + sum(log(diff(y1)))-(n1-1)*log(S1)) 
                            - (n1-1)*(lgamma(alpha1_new) - lgamma(alpha1)) + a1*(log(lambda1_new)-log(lambda1)) 
                            - b1*(lambda1_new - lambda1)
                            +lgamma(alpha1_new*(n1-1) + 1) - lgamma(alpha1*(n1-1) + 1)))
    
    if(u3 < acc_crit1 && lambda1_new > 0.1 && alpha1_new > 0){
      alpha1 <- alpha1_new
      lambda1 <- lambda1_new
    }
    acc_crit2 <- min(1, exp((alpha2_new - alpha2)*(-1/c2 + sum(log(diff(y2)))-(n2-1)*log(S2)) 
                            - (n2-1)*(lgamma(alpha2_new) - lgamma(alpha2)) + a2*(log(lambda2_new)-log(lambda2)) 
                            - b2*(lambda2_new - lambda2)
                            +lgamma(alpha2_new*(n2-1) + 1) - lgamma(alpha2*(n2-1) + 1)))
    if(u4 < acc_crit2 && lambda2_new > 0.1 && alpha2_new > 0){
      alpha2 <- alpha2_new
      lambda2 <- lambda2_new
    }
  }

  
  return(c(alpha1, lambda1, alpha2, lambda2))
}


#########Gibbs Step k_i########
############functions for updating k in blocks###############
block_k <- function(block_start, prev_p1, prev_p2, next_p1, next_p2, alpha1, alpha2, lambda1, lambda2, block_size, y){
  
  #prev_p1: last event index for process 1 before the block
  #prev_p2: last event index for process 2 before the block
  #next_p1: first event index for process 1 after the block
  #next_p2: first event index for process 2 after the block
  #alpha1: shape parameter for process 1
  #alpha2: shape parameter for process 2
  #lambda1: rate parameter for process 1
  #lambda2: rate parameter for process 2
  #block_size: number of events in the block
  #block_start: start index of the block
  #y: vector of all events
  
  B <- binary_list(block_size) #Make a list of all binary vectors of length block_size
  P <- numeric(0) # Empty vector for storing probabilities
  k <- numeric(block_size) #Initialize the allocation vector
  
  #Loop through each binary vector
  for(b in B){
    temp_prev_p1 <- prev_p1
    temp_prev_p2 <- prev_p2
    p_B <- 1 #Initialize the probability for this binary vector
    index <- block_start #Index for the event times in y
    for(digit in b){
      if(digit == 0){#The digit being 0 means it belongs to process 1
        #Calculate the probability of the next event belonging to process 1
        prob <- dgamma(y[index] - y[temp_prev_p1], shape = alpha1, rate = lambda1)
        p_B <- p_B * prob
        temp_prev_p1 <- index #Update the last event time for process 1
      } 
      else {#The digit being 1 means it belongs to process 2
        #Calculate the probability of the next event belonging to process 2
        prob <- dgamma(y[index] - y[temp_prev_p2], shape = alpha2, rate = lambda2)
        p_B <- p_B * prob
        temp_prev_p2 <- index #Update the last event time for process 2
      }
      index <- index + 1
    }
    #After processing all events in the block, account for the time until the next known event
    p_B <- p_B * dgamma(y[next_p1] - y[temp_prev_p1], shape = alpha1, rate = lambda1)
    p_B <- p_B * dgamma(y[next_p2] - y[temp_prev_p2], shape = alpha2, rate = lambda2)
    
    P <- c(P, p_B) #Store the probability for this binary vector
  }
  #Sample the right configuration
  number <- sample(1:length(B), size = 1, prob = P)
  right_config <- B[[number]]
  #change the allocation of events in y according to the sampled binary vector
  for(i in 1:block_size){
    if(right_config[i] == 0){
      k[i] <- 1
    } else {
      k[i] <- 2
    }
  }
  return(k)
}

gibbs_step_k <- function(y,k, alpha1, lambda1, alpha2, lambda2, block_size){
  n <- length(y)
  n_blocks <- floor((n-2)/block_size)
  #select n_blocks starting points at random
  block_starts <- sample(3:(n-block_size-2), n_blocks, replace = TRUE)
  
  for(i in block_starts){
    i_minus_1 <- find_previous(y, k, i, 1)
    i_minus_2 <- find_previous(y, k, i, 2)
    i_plus_1 <- find_next(y, k, i + (block_size-1), 1)
    i_plus_2 <- find_next(y, k, i + (block_size-1), 2)

    new_k_block <- block_k(block_start = i, prev_p1 = i_minus_1, prev_p2 = i_minus_2, next_p1 = i_plus_1,
                           next_p2 = i_plus_2, alpha1 = alpha1, alpha2 = alpha2, lambda1 = lambda1, lambda2 = lambda2, block_size = block_size, y = y)
    k[i:(i + (block_size - 1))] <- new_k_block
  }
  return(k)
}
########################

#########Helper Functions to find previous and next values of k########
find_previous <- function(y, k, i, val){
  for (j in (i-1):1) {
    if (k[j] == val) {
      return(j)
    }
  }
}
find_next <- function(y, k, i, val){
  for (j in (i+1):length(y)) {
    if (k[j] == val) {
      return(j)
    }
  }
}
############################################

#######Helper functions simulate latent event before and after#########
sim_before <- function(y, alpha, lambda, tmax){
  repeat{
    u <- runif(1)
    sim_value <- qgamma(u, shape = alpha, rate = lambda)
    y_minus <- y[1] - sim_value
    if(y_minus <= 0) break
  }
  return(y_minus)
}
sim_after <- function(y, alpha, lambda, tmax){
  repeat{
    n <- length(y)
    u <- runif(1)
    sim_value <- qgamma(u, shape = alpha, rate = lambda)
    y_pluss <- y[n] + sim_value
    if(y_pluss >= tmax) break
  }
  return(y_pluss)
}
################

###########Make binary list for block_k function############
#Make a list of all binary vectors of length m
binary_list <- function(m) {
  n <- 2^m
  lapply(0:(n-1), function(x) {
    rev(as.integer(intToBits(x))[1:m])  # take first m bits, reversed to usual order
  })
}
##############################################

###########################
log_lik_reparam <- function(y, alpha, lambda){
  n <- length(y)
  ints <- diff(y)
  log_lik <- sum(dgamma(ints, shape = alpha, rate = lambda, log = TRUE))
  return(log_lik)
}


##########display plots functions##########
display_traceplot <- function(chain) {
  par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
  traceplot(chain[, 1], main = paste("Traceplot of Alpha1"))
  traceplot(chain[, 2], main = paste("Traceplot of Lambda1"))
  traceplot(chain[, 3], main = paste("Traceplot of Alpha2"))
  traceplot(chain[, 4], main = paste("Traceplot of Lambda2"))
}

display_densplot <- function(chain) {
  par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
  densplot(chain[, 1], main = paste("Density Plot of Alpha1"))
  densplot(chain[, 2], main = paste("Density Plot of Lambda1"))
  densplot(chain[, 3], main = paste("Density Plot of Alpha2"))
  densplot(chain[, 4], main = paste("Density Plot of Lambda2"))
}

display_acf <- function(chain) {
  par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
  acf(chain[, 1], main = paste("ACF of Alpha1"))
  acf(chain[, 2], main = paste("ACF of Lambda1"))
  acf(chain[, 3], main = paste("ACF of Alpha2"))
  acf(chain[, 4], main = paste("ACF of Lambda2"))
}
######################################



