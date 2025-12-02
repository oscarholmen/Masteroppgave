source("~/Library/CloudStorage/OneDrive-NTNU/Masteroppgave/Gibbs Sampler.R")
source("~/Library/CloudStorage/OneDrive-NTNU/Masteroppgave/UpdatedGeneratingDataset.R")

#Generating the data
set.seed(123) #For reproducibility
timeseries = generate_vec(4, 2, 25, 5, tmax = 400)

#Setting parameters for Gibbs sampler, no extra for better mixing####
block_size = 1
gibbs_steps = 100000
burnin_steps = 5000
thinning = 5
prior_lambda = c(0.01,0.01,0.05,0.01)
prior_alpha = c(2, 5)
initial_guess = c(4,2, 25,5)
initial_guess1 = c(1,1, 10,3)
initial_guess2 = c(5,2, 20,5)
initial_guess3 = c(3,1, 10,1)
########
########Finding best parameters for normal Gibbs sampler by tuning sd1 and sd2#########
sd_process1 = c(1.0, 1.25, 1.5, 1.75)
sd_process2 = c(2.0, 2.5, 3.0)
ess_results <- matrix(NA, nrow = length(sd_process1), ncol = length(sd_process2))
for (i in 1:length(sd_process1)) {
  for (j in 1:length(sd_process2)) {
    set.seed(123)
    chain <- gibbs_superpos(
      y = timeseries$sup_vec,
      tmax = timeseries$tmax,
      method = "one_dist",
      mcmc_steps = gibbs_steps,
      burnin = burnin_steps,
      init_guess_params = initial_guess,
      d1 = sd_process1[i],
      d2 = sd_process2[j],
      prior_lambda = prior_lambda,
      prior_alpha = prior_alpha,
      thinning = thinning
    )
    ess <- effectiveSize(chain[, 1:4])
    print(c(sd_process1[i], sd_process2[j], ess))
    ess_results[i, j] <- min(ess)# Store the minimum ESS among parameters
  }
}
#######Running normal Gibbs sampler#############
set.seed(123)
normal_chain <- gibbs_superpos(
  y = timeseries$sup_vec,
  tmax = timeseries$tmax,
  method = "one_dist",
  mcmc_steps = 100000,
  burnin = burnin_steps,
  init_guess_params = initial_guess,
  d1 = 1.25,
  d2 = 1.75,
  prior_lambda = prior_lambda,
  prior_alpha = prior_alpha,
  thinning = thinning
)


ess_normal <- effectiveSize(normal_chain[, 1:4])
print(ess_normal)
summary(normal_chain)
display_traceplot(normal_chain)
display_densplot(normal_chain)
display_acf(normal_chain)
idx <- sample(1:nrow(normal_chain), 2000)
pairs(normal_chain[idx,1:4], main = "Independent Priors",
      labels = c(expression(alpha[1]), expression(lambda[1]), expression(alpha[2]), expression(lambda[2])),
      pch = 15, cex = 0.5)
###############################################
###############Gelman Rubin Diagnostic for normal Gibbs sampler###############
for (i in 1:3) {
  chain_name <- paste0("chain", i)
  init_guess <- get(paste0("initial_guess", i))
  set.seed(123)
  assign(chain_name, gibbs_superpos(
    y = timeseries$sup_vec,
    tmax = timeseries$tmax,
    method = "one_dist",
    mcmc_steps = 10000,
    burnin = 2000,
    init_guess_params = init_guess,
    d1 = 1.25,
    d2 = 1.75,
    prior_lambda = prior_lambda,
    prior_alpha = prior_alpha,
    block_size = 1,
    thinning = thinning
  ))
}
gelman.diag(mcmc.list(chain1, chain2, chain3), autoburnin = FALSE, multivariate = FALSE)
gelman.plot(mcmc.list(chain1[,1:4], chain2[,1:4], chain3[,1:4]), autoburnin = FALSE)
##########################################




###########################Reparametrized Gibbs Sampler###################
sd_process1 = c(1.25)
sd_process2 = c(1.50)
ess_results_reparam <- matrix(NA, nrow = length(sd_process1), ncol = length(sd_process2))
for (i in 1:length(sd_process1)) {
  for (j in 1:length(sd_process2)) {
    set.seed(123)
    chain <- gibbs_superpos(
      y = timeseries$sup_vec,
      tmax = timeseries$tmax,
      method = "reparam",
      mcmc_steps = mcmc_steps,
      burnin = burnin_steps,
      init_guess_params = initial_guess,
      d1 = sd_process1[i],
      d2 = sd_process2[j],
      prior_lambda = prior_lambda,
      prior_alpha = prior_alpha,
      thinning = thinning
    )
    ess <- effectiveSize(chain[, 1:4])
    print(c(sd_process1[i], sd_process2[j], ess))
    ess_results_reparam[i, j] <- min(ess)# Store the minimum ESS among parameters
  }
}
ess_results_reparam

###########################
chain_reparam <- gibbs_superpos(
  y = timeseries$sup_vec,
  tmax = timeseries$tmax,
  method = "reparam",
  mcmc_steps = 100000,
  burnin = burnin_steps,
  init_guess_params = initial_guess,
  d1 = 1.25,
  d2 = 1.5,
  prior_lambda = prior_lambda,
  prior_alpha = prior_alpha,
  thinning = thinning
)
ess_reparam <- effectiveSize(chain_reparam[, 1:4])
print(ess_reparam)
summary(chain_reparam)
display_traceplot(chain_reparam)
display_densplot(chain_reparam)
display_acf(chain_reparam)

#Gelman rubin for reparametrized Gibbs sampler#
for (i in 1:3) {
  chain_name <- paste0("chain_reparam", i)
  init_guess <- get(paste0("initial_guess", i))
  set.seed(123)
  assign(chain_name, gibbs_superpos(
    y = timeseries$sup_vec,
    tmax = timeseries$tmax,
    method = "reparam",
    mcmc_steps = 10000,
    burnin = 15000,
    init_guess_params = init_guess,
    d1 = 1.25,
    d2 = 1.5,
    prior_lambda = prior_lambda,
    prior_alpha = prior_alpha,
    thinning = thinning
  ))
}
gelman.diag(mcmc.list(chain_reparam1, chain_reparam2, chain_reparam3),autoburnin = FALSE, multivariate = FALSE)
gelman.plot(mcmc.list(chain_reparam1[,1:4], chain_reparam2[,1:4], chain_reparam3[,1:4]), autoburnin = FALSE)
display_traceplot(chain_reparam3)
###########################


#Scatterplot to visualise correlation between parameters




################Using Block Gibbs Sampler#########################
#testing sd values for block Gibbs sampler
sd_process1 = c(1.5, 1.75)
sd_process2 = c(2.25, 2.5, 2.75)
ess_results_block <- matrix(NA, nrow = length(sd_process1), ncol = length(sd_process2))
for (i in 1:length(sd_process1)) {
  for (j in 1:length(sd_process2)) {
    set.seed(123)
    chain <- gibbs_superpos(
      y = timeseries$sup_vec,
      tmax = timeseries$tmax,
      method = "block",
      mcmc_steps = 25000,
      burnin = burnin_steps,
      init_guess_params = initial_guess,
      d1 = sd_process1[i],
      d2 = sd_process2[j],
      prior_lambda = prior_lambda,
      prior_alpha = prior_alpha,
      block_size = 1,
      thinning = thinning
    )
    ess <- effectiveSize(chain[, 1:4])
    print(c(sd_process1[i], sd_process2[j], ess))
    ess_results_block[i, j] <- min(ess)# Store the minimum ESS among parameters
  }
}
print(ess_results_block)

chain_block <- gibbs_superpos(
  y = timeseries$sup_vec,
  tmax = timeseries$tmax,
  method = "block",
  mcmc_steps = gibbs_steps,
  burnin = burnin_steps,
  init_guess_params = initial_guess,
  d1 = 0.5,
  d2 = 1.75,
  prior_lambda = prior_lambda,
  prior_alpha = prior_alpha,
  block_size = 2,
  thinning = thinning
)
ess_block <- effectiveSize(chain_block[, 1:4])
print(ess_block)
summary(chain_block)
display_traceplot(chain_block)
display_densplot(chain_block)
display_acf(chain_block)
pairs(chain_block[idx,1:4], main = "Block Gibbs Sampler",
      labels = c(expression(alpha[1]), expression(lambda[1]), expression(alpha[2]), expression(lambda[2])),
      pch = 15, cex = 0.5)

#Getting Gelman-Rubin diagnostic for block Gibbs sampler#
for (i in 1:3) {
  chain_name <- paste0("chain_block", i)
  init_guess <- get(paste0("initial_guess", i))
  set.seed(123)
  assign(chain_name, gibbs_superpos(
    y = timeseries$sup_vec,
    tmax = timeseries$tmax,
    method = "block",
    mcmc_steps = 10000,
    burnin = 5000,
    init_guess_params = init_guess,
    d1 = 1.5,
    d2 = 1.75,
    prior_lambda = prior_lambda,
    prior_alpha = prior_alpha,
    block_size = 1,
    thinning = thinning
  ))
}
gelman.diag(mcmc.list(chain_block1, chain_block2, chain_block3), autoburnin = FALSE, multivariate = FALSE)
gelman.plot(mcmc.list(chain_block1[,1:4], chain_block2[,1:4], chain_block3[,1:4]), autoburnin = FALSE)

###########################



#######Finding optimal block size for block Gibbs sampler#########
#Setting parameters to gridsearch to find optimal block size####
block_sizes = c(1,2,3,4)
ess_mat <- matrix(NA, nrow = 4, ncol = length(block_sizes))
time_elapsed = numeric(length(block_sizes))
efficency = numeric(length(block_sizes))
######

########calculate ESS and efficiency for each block size########
set.seed(123)
for (j in 1:length(block_sizes)) {
  bs <- block_sizes[j]
  
  time_elapsed[j] <- system.time(
  chain <- gibbs_superpos(
    y = timeseries$sup_vec,
    tmax = timeseries$tmax,
    method = "reparam",
    mcmc_steps = 50000,
    burnin = burnin_steps,
    init_guess_params = initial_guess,
    d1 = 1.25,
    d2 = 1.5,
    prior_lambda = prior_lambda,
    prior_alpha = prior_alpha,
    block_size = bs,
    thinning = thinning
  ))[["elapsed"]]
  
  # ensure chain is in matrix form before passing to as.mcmc
  chain_mcmc <- as.mcmc(as.matrix(chain[, 1:4]))
  
  ess_mat[, j] <- effectiveSize(chain_mcmc)
}

print(ess_mat[,4])
print(time_elapsed)
##########

####Results for block updates and indicator block size 3 ########
set.seed(123)
block_chain_reparam <- gibbs_superpos(
  y = timeseries$sup_vec,
  tmax = timeseries$tmax,
  method = "reparam",
  mcmc_steps = gibbs_steps,
  burnin = burnin_steps,
  init_guess_params = initial_guess,
  d1 = 1.25,
  d2 = 1.75,
  prior_lambda = prior_lambda,
  prior_alpha = prior_alpha,
  block_size = 2,
  thinning = thinning
)

ess_block_reparam <- effectiveSize(block_chain_reparam[, 1:4])
print(ess_block_reparam)
summary(block_chain_reparam)

display_traceplot(block_chain_reparam[, 1:4])
display_densplot(block_chain_reparam[, 1:4])
display_acf(block_chain_reparam[, 1:4])
pairs(block_chain_reparam[idx,1:4], main = "Reparametrized Prior",
      labels = c(expression(alpha[1]), expression(lambda[1]), expression(alpha[2]), expression(lambda[2])),
      pch = 15, cex = 0.5)


########


##### set up scheme to check convergence with Gelman-Rubin diagnostic #####
for (i in 1:3) {
  chain_name <- paste0("chain_block_norm", i)
  init_guess <- get(paste0("initial_guess", i))
  set.seed(123)
  assign(chain_name, gibbs_superpos(
    y = timeseries$sup_vec,
    tmax = timeseries$tmax,
    method = "reparam",
    mcmc_steps = 10000,
    burnin = 5000,
    init_guess_params = init_guess,
    d1 = sd_process1,
    d2 = sd_process2,
    prior_lambda = prior_lambda,
    prior_alpha = prior_alpha,
    block_size = 2,
    thinning = thinning
  ))
}

gelman.diag(mcmc.list(chain_block_norm1, chain_block_norm2, chain_block_norm3), autoburnin = FALSE, multivariate = FALSE)
gelman.plot(mcmc.list(chain_block_norm1[,1:4], chain_block_norm2[,1:4], chain_block_norm3[,1:4]), autoburnin = FALSE)
######




#####################Block k and reparametrized Gibbs sampler combined#####################
# Finding best sd1 and sd2 for block reparametrized Gibbs sampler
sd_process1 = c(0.5, 1.0 , 1.25)
sd_process2 = c(1.5, 2.0)
block_size = c(2)
set.seed(123)
for (i in 1:length(sd_process1)) {
  for (j in 1:length(sd_process2)) {
    for (k in 1:length(block_size)){
      set.seed(123)
      chain <- gibbs_superpos(
        y = timeseries$sup_vec,
        tmax = timeseries$tmax,
        method = "block_reparam",
        mcmc_steps = 20000,
        burnin = 0,
        init_guess_params = initial_guess,
        d1 = sd_process1[i],
        d2 = sd_process2[j],
        prior_lambda = prior_lambda,
        prior_alpha = prior_alpha,
        block_size = block_size[k],
        thinning = thinning
      )
      ess <- effectiveSize(chain[, 1:4])
      print(c(sd_process1[i], sd_process2[j], block_size[k], ess))
    }
  }
}


#
chain_block_reparam <- gibbs_superpos(
  y = timeseries$sup_vec,
  tmax = timeseries$tmax,
  method = "block_reparam",
  mcmc_steps = 100000,
  burnin = burnin_steps,
  init_guess_params = initial_guess,
  d1 = 1.25,
  d2 = 2.0,
  prior_lambda = prior_lambda,
  prior_alpha = prior_alpha,
  block_size = 2,
  thinning = thinning
)
ess_block_reparam <- effectiveSize(chain_block_reparam[, 1:4])
print(ess_block_reparam)
summary(chain_block_reparam)
display_traceplot(chain_block_reparam[, 1:4])
display_densplot(chain_block_reparam[, 1:4])
display_acf(chain_block_reparam[, 1:4])
#########




#Posterior correlation between parameters
cor_matrix <- cor(chain_mcmc[, 1:4])
print(cor_matrix)


########### Making densities so i can plot posterior densities together ##########
density_normal1 <- density(normal_chain[,1], bw = 0.1, from = 0, to = 8)
density_normal2 <- density(normal_chain[,2], bw = 0.05, from = 0, to = 4)
density_normal3 <- density(normal_chain[,3], bw = 1.2, from = 0, to = 50)
density_normal4 <- density(normal_chain[,4], bw = 0.25, from = 0, to = 10)
density_block1 <- density(chain_block[,1], bw = 0.1, from = 0, to = 8)
density_block2 <- density(chain_block[,2], bw = 0.05, from = 0, to = 4)
density_block3 <- density(chain_block[,3], bw = 1.2, from = 0, to = 50)
density_block4 <- density(chain_block[,4], bw = 0.25, from = 0, to = 10)
density_reparam1 <- density(block_chain_reparam[,1], bw = 0.1, from = 0, to = 8)
density_reparam2 <- density(block_chain_reparam[,2], bw = 0.05, from = 0, to = 4)
density_reparam3 <- density(block_chain_reparam[,3], bw = 1.2, from = 0, to = 50)
density_reparam4 <- density(block_chain_reparam[,4], bw = 0.25, from = 0, to = 10)
density_filter_indep1 <- density(chain_filter[,1], bw = 0.1, from = 0, to = 8)
density_filter_indep2 <- density(chain_filter[,2], bw = 0.05, from = 0, to = 4)
density_filter_indep3 <- density(chain_filter[,3], bw = 1.2, from = 0, to = 50)
density_filter_indep4 <- density(chain_filter[,4], bw = 0.25, from = 0, to = 10)
density_filter_reparam1 <- density(chain_alt_prior[,1], bw = 0.1, from = 0, to = 8)
density_filter_reparam2 <- density(chain_alt_prior[,2], bw = 0.05, from = 0, to = 4)
density_filter_reparam3 <- density(chain_alt_prior[,3], bw = 1.2, from = 0, to = 50)
density_filter_reparam4 <- density(chain_alt_prior[,4], bw = 0.25, from = 0, to = 10)
#making x axis for plotting
x_1 <- seq(0, 8, length.out = 512)
x_2 <- seq(0, 4, length.out = 512)
x_3 <- seq(0, 50, length.out = 512)
x_4 <- seq(0, 10, length.out = 512)

##########################
df1 <- data.frame(x = x_1,
                     Normal = density_normal1$y,
                     Reparam = density_reparam1$y,
                     Block = density_block1$y,
                     Filter_Indep = density_filter_indep1$y,
                     Filter_Reparam = density_filter_reparam1$y)
df2 <- data.frame(x = x_2,
                     Normal = density_normal2$y,
                     Reparam = density_reparam2$y,
                     Block = density_block2$y,
                     Filter_Indep = density_filter_indep2$y,
                     Filter_Reparam = density_filter_reparam2$y)
df3 <- data.frame(x = x_3,
                     Normal = density_normal3$y,
                     Reparam = density_reparam3$y,
                     Block = density_block3$y,
                     Filter_Indep = density_filter_indep3$y,
                     Filter_Reparam = density_filter_reparam3$y)
df4 <- data.frame(x = x_4,
                     Normal = density_normal4$y,
                     Reparam = density_reparam4$y,
                     Block = density_block4$y,
                     Filter_Indep = density_filter_indep4$y,
                     Filter_Reparam = density_filter_reparam4$y)

library(ggplot2)
ggplot(df1, aes(x = x)) +
  geom_line(aes(y = Normal, color = "Gibbs:Independent Priors"), size = 1) +
  geom_line(aes(y = Reparam, color = "Gibbs:Reparametrized Priors"), size = 1) +
  geom_line(aes(y = Block, color = "Gibbs:Block"), size = 1) +
  geom_line(aes(y = Filter_Indep, color = "Filter:Independent Priors"), size = 1) +
  geom_line(aes(y = Filter_Reparam, color = "Filter:Reparametrized Priors"), size = 1) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "black") +
  labs(y = "Posterior Density") +
  theme_bw() +
  theme(legend.position = c(0.17,0.85))+
  scale_y_continuous(labels = NULL)
ggplot(df2, aes(x = x)) +
  geom_line(aes(y = Normal, color = "Gibbs:Independent Priors"), size = 1) +
  geom_line(aes(y = Reparam, color = "Gibbs:Reparametrized Priors"), size = 1) +
  geom_line(aes(y = Block, color = "Gibbs:Block"), size = 1) +
  geom_line(aes(y = Filter_Indep, color = "Filter:Independent Priors"), size = 1) +
  geom_line(aes(y = Filter_Reparam, color = "Filter:Reparametrized Priors"), size = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  labs(y = "Posterior Density") +
  theme_bw() +
  theme(legend.position = c(0.17,0.85))+
  scale_y_continuous(labels = NULL)
ggplot(df3, aes(x = x)) +
  geom_line(aes(y = Normal, color = "Gibbs:Independent Priors"), size = 1) +
  geom_line(aes(y = Reparam, color = "Gibbs:Reparametrized Priors"), size = 1) +
  geom_line(aes(y = Block, color = "Gibbs:Block"), size = 1) +
  geom_line(aes(y = Filter_Indep, color = "Filter:Independent Priors"), size = 1) +
  geom_line(aes(y = Filter_Reparam, color = "Filter:Reparametrized Priors"), size = 1) +
  geom_vline(xintercept = 25, linetype = "dashed", color = "black") +
  labs(y = "Posterior Density") +
  theme_bw() +
  theme(legend.position = c(0.17,0.85))+
  scale_y_continuous(labels = NULL)
ggplot(df4, aes(x = x)) +
  geom_line(aes(y = Normal, color = "Gibbs:Independent Priors"), size = 1) +
  geom_line(aes(y = Reparam, color = "Gibbs:Reparametrized Priors"), size = 1) +
  geom_line(aes(y = Block, color = "Gibbs:Block"), size = 1) +
  geom_line(aes(y = Filter_Indep, color = "Filter:Independent Priors"), size = 1) +
  geom_line(aes(y = Filter_Reparam, color = "Filter:Reparametrized Priors"), size = 1) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "black") +
  labs(y = "Posterior Density") +
  theme_bw() +
  theme(legend.position = c(0.17,0.85))+
  scale_y_continuous(labels = NULL)




#These densityplots are not correct, need to fix them

display_densplot(normal_chain)
display_densplot(chain_block_reparam)
display_densplot(chain_reparam)
