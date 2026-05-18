#### Upload packages ####

# upload all relevant packages
source('code/upload_packages.R')
source('code/helpers.R')
source('code/model_code.R')


# load model input data
load('data/model_data_input.RData')

# data_input: list of data input for the models
# init: list of initial values for model parameters
# const: list of constant values measuring the range of the indices


# Number of chains and cores
nchains <- 4
ncores <- nchains

# Define initial values for each chain
init_list <- list(
  chain1 = init,  
  chain2 = init,
  chain3 = init,
  chain4 = init)  

# Create the parallel cluster
cl <- makeCluster(ncores)


# Export necessary objects to the workers
clusterExport(cl, c("model_nimble_inverse_gamma", "data_input", "const", "init_list", "run_single_chain"))
clusterEvalQ(cl, {
  library(nimble)
  library(nimbleHMC)
})

# Run each chain in parallel
chain_results <- parLapply(
  cl = cl,
  X = 1:nchains,
  fun = function(chain_id) {
    run_single_chain(
      chain_id = chain_id,
      model_code = model_nimble_new,
      data = dd_cov,
      consts = const_new,
      init = init_list,
      monitors = names(init_list$chain1),
      niter = 6000,
      nburnin = 3000,
      thin = 10
    )
  }
)

stopCluster(cl)


# Save the results

# all parameters

mcmc.out_all = MCMCsummary(
  chain_results,
  params = "all",
  excl = NULL,
  ISB = TRUE,
  exact = TRUE,
  probs = c(0.025, 0.5, 0.975),
  hpd_prob = 0.95,
  HPD = FALSE,
  pg0 = FALSE,
  digits = NULL,
  round = NULL,
  Rhat = TRUE,
  n.eff = TRUE,
  func = NULL,
  func_name = NULL
) 



# check convergence using Rhat
mcmc.out_all$Rhat[mcmc.out_all$Rhat>1.1]

ifelse(sum(mcmc.out_all$Rhat[mcmc.out_all$Rhat>=1.1])>0,
       'Rhat>1.1 for some parameters',
       'Rhat<1.1 for all parameters')


# Extract TFR posterior estimates

mcmc.out_TFR_final = MCMCsummary(
  chain_results,
  params = "TFR",
  excl = NULL,
  ISB = TRUE,
  exact = TRUE,
  probs = c(0.025, 0.5, 0.975),
  hpd_prob = 0.95,
  HPD = FALSE,
  pg0 = FALSE,
  digits = NULL,
  round = NULL,
  Rhat = TRUE,
  n.eff = TRUE,
  func = NULL,
  func_name = NULL
) 

# add relevant variables

mcmc.out_TFR_final$years = rep(1751:1910,8)

names(mcmc.out_TFR_final)[3:5] = c('TFR_lower','TFR_median','TFR_upper')

mcmc.out_TFR_final$country = rep(c("SWE","ENG","FRA",
                                   "NOR","DEN","NLD","USA","FIN"),each=160)




# Extract theta posterior estimates (log-scale)




mcmc.out_theta_final = MCMCsummary(
  chain_results,
  params = "theta",
  excl = NULL,
  ISB = TRUE,
  exact = TRUE,
  probs = c(0.025, 0.5, 0.975),
  hpd_prob = 0.95,
  HPD = FALSE,
  pg0 = FALSE,
  digits = NULL,
  round = NULL,
  Rhat = TRUE,
  n.eff = TRUE,
  func = NULL,
  func_name = NULL
) 

mcmc.out_theta_final$years = rep(1751:1910,8)

names(mcmc.out_theta_final)[3:5] = c('theta_lower','theta_median','theta_upper')

mcmc.out_theta_final$country = rep(c("SWE","ENG","FRA",
                                     "NOR","DEN","NLD","USA","FIN"),
                                   each=160)



# save results in a .RData file

save(mcmc.out_all,mcmc.out_theta_final,mcmc.out_TFR_final,
     file='results/bayesian_model_results.RData')





