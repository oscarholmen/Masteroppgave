source("~/Library/CloudStorage/OneDrive-NTNU/Masteroppgave/UpdatedGeneratingDataset.R")
source("~/Library/CloudStorage/OneDrive-NTNU/Masteroppgave/Gibbs Sampler.R")
source("~/Library/CloudStorage/OneDrive-NTNU/Masteroppgave/KalmanFilterMethod.R")


#Generating the data
set.seed(123) #For reproducibility
#set.seed(420) #testing other dataset
timeseries <- generate_vec(4, 2, 25, 5, tmax = 400)

#Testing the function via optim
neg_loglik <- function(params, y, tmax, thresh = 1e-10){
  alpha1 <- params[1]
  lambda1 <- params[2]
  alpha2 <- params[3]
  lambda2 <- params[4]
  if(any(params <= 0)){
    return(Inf)
  }
  -kalman_filter_method(y, alpha1, lambda1, alpha2, lambda2, tmax, thresh)
}

microbenchmark::microbenchmark(
  neg_loglik(c(5,1,10,3), y = timeseries$sup_vec, tmax = timeseries$tmax),times = 10
)

res <- optim(par = c(2,1,10,3), fn = neg_loglik, y = timeseries$sup_vec, tmax = timeseries$tmax,
             method = "L-BFGS-B", lower = c(1e-5, 1e-5, 1e-5, 1e-5),
             control = list(fnscale = 1, trace = 1, maxit = 1000),
             hessian = TRUE)
res$par
res$hessian
res$value
eigen(res$hessian)$values

#Making a 2D grid to visualize the likelihood surface along two parameters
lambda1_seq <- seq(1.5,3,length.out = 20)
lambda2_seq <- seq(3,6.5,length.out = 20)
loglik_matrix <- matrix(0, nrow = length(lambda1_seq), ncol = length(lambda2_seq))
for(i in 1:length(lambda1_seq)){
  for(j in 1:length(lambda2_seq)){
    loglik_matrix[i,j] <- metropfunc_alt_prior(c(alpha1 = 4,lambda1 = lambda1_seq[i],alpha2 = 25,lambda2 = lambda2_seq[j])
                                     ,timeseries = timeseries)
  }
}
library(reshape2)
library(ggplot2)
loglik_df <- melt(loglik_matrix)
colnames(loglik_df) <- c("lambda1_index", "lambda2_index", "loglik")
loglik_df$lambda1 <- lambda1_seq[loglik_df$lambda1_index]
loglik_df$lambda2 <- lambda2_seq[loglik_df$lambda2_index]
#Contour plot
ggplot(loglik_df, aes(x = lambda1, y = lambda2, z = loglik)) +
  geom_contour_filled(bins = 20) +
  labs(title = "Log-Likelihood Contour Plot", x = expression(lambda[1]), y = expression(lambda[2])) +
  theme_minimal()

#Checks to gradients and Hessians
library(numDeriv)
grad <- grad(func = neg_loglik, x = res$par, y = timeseries$sup_vec, tmax = timeseries$tmax)
norm_grad <- sqrt(sum(grad^2))
print(grad)
print(norm_grad)
hess <- hessian(func = neg_loglik, x = res$par, y = timeseries$sup_vec, tmax = timeseries$tmax)
print(hess)
eigen_hess <- eigen(hess)
print(eigen_hess$values)
kappa(hess)

#making confidence sets using the loglikelihood ratio test and asymptotic chi-squared distribution
alpha_level <- 0.05
chi_sq_cutoff <- qchisq(1 - alpha_level, df = 2) 
loglikelihood_cutoff <- -res$value - chi_sq_cutoff / 2
loglik_df$in_conf_set <- loglik_df$loglik >= loglikelihood_cutoff
#Contour plot with confidence set
ggplot(loglik_df, aes(x = alpha2, y = lambda2)) +
  geom_contour_filled(aes(z = loglik), bins = 20, alpha = 0.7) +
  geom_point(data = subset(loglik_df, in_conf_set), aes(x = alpha2, y = lambda2), color = "red", size = 0.5) +
  labs(title = "Log-Likelihood Contour Plot with Confidence Set", x = expression(alpha[1]), y = expression(lambda[1])) +
  theme_minimal()


#Using the Kalman filter method within MCMC as likelihood
microbenchmark::microbenchmark(
  metropfunc(c(4, 2, 25, 5), timeseries),times = 10
)

chain_filter <- MCMCmetrop1R(
  fun = metropfunc,
  theta.init = c(2, 1, 10, 3),
  timeseries = timeseries,
  mcmc = 5000,
  burnin = 500,
  tune = 1,
  verbose = 500,
  optim.method = "L-BFGS-B",
  logfun = TRUE,
  optim.lower = c(1e-5, 1e-5, 1e-5, 1e-5),
  optim.upper = c(1e3, 1e3, 1e3, 1e3),
)


ess_filter <- effectiveSize(chain_filter[, 1:4])
print(ess_filter)
display_traceplot(chain_filter)
display_densplot(chain_filter)
display_acf(chain_filter)
summary(chain_filter)

chain_filter_mat <- as.matrix(chain_filter[,1:4])
set.seed(1)
idx <- sample(1:nrow(chain_alt_prior_mat), 2000)
pairs(chain_filter_mat[idx,], main = "Scatterplot Filter Independent Priors",
      labels = c(expression(alpha[1]), expression(lambda[1]), expression(alpha[2]), expression(lambda[2])),
      pch = 15, cex = 0.5)


      

metropfunc_alt_prior(c(4, 2, 25, 5), timeseries)


chain_alt_prior <- MCMCmetrop1R(
  fun = metropfunc_alt_prior,
  theta.init = c(4, 2, 25, 5),
  timeseries = timeseries,
  mcmc = 5000,
  burnin = 500,
  tune = 1,
  verbose = 500,
  optim.method = "L-BFGS-B",
  logfun = TRUE,
  optim.lower = c(1e-5, 1e-5, 1e-5, 1e-5),
  optim.upper = c(1e3, 1e3, 1e3, 1e3),
)

ess_alt_prior <- effectiveSize(chain_alt_prior[, 1:4])
print(ess_alt_prior)
display_traceplot(chain_alt_prior)
display_densplot(chain_alt_prior)
display_acf(chain_alt_prior)
summary(chain_alt_prior)
#scatterplot 
chain_alt_prior_mat <- as.matrix(chain_alt_prior[,1:4])
set.seed(1)
idx <- sample(1:nrow(chain_alt_prior_mat), 2000)
pairs(chain_alt_prior_mat[idx,], main = "Scatterplot Filter Reparametrized Prior", 
      labels = c(expression(alpha[1]), expression(lambda[1]), expression(alpha[2]), expression(lambda[2])),
      pch = 15, cex = 0.5)


######### Visualising the expected length of the filtered state ################
kalman_filter_method(
  y = timeseries$sup_vec,
  alpha1 = 4,
  lambda1 = 2,
  alpha2 = 25,
  lambda2 = 5,
  tmax = timeseries$tmax,
  thresh = 1e-10
)
