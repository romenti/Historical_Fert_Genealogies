#### Upload packages ####

# upload all relevant packages
source('code/upload_packages.R')
source('code/helpers.R')
source('code/model_code.R')


# load model input data
load('data/model_data_input.RData')

# upload all relevant packages

source('code/upload_packages.R')
source('code/model_code.R')

# load model input data
load('data/model_data_input.RData')


#### Model without Sweden ####

model_nimble_inverse_gamma_no_sweden = nimble::nimbleCode({
  
  # Hyperparameters
  
  #sigma_beta_c ~ dunif(0,40)
  sigma_beta_c ~ dinvgamma(1,1)
  beta_dummy ~ dnorm(mean=0,sd=1)
  sd_sigma_beta_sp ~ dunif(0,40)
  mu_sigma_beta_sp ~ dnorm(mean=0,sd=10)
  for(j in 2:A){
    sigma_TFR[j] ~ dinvgamma(1,1)
    beta_c[j] ~ dnorm(mean=0,sd=sigma_beta_c)
    log_sigma_theta[j] ~  dnorm(0,sd=1)
    sigma_theta[j] ~ dinvgamma(1,1)
    log_sigma_beta_sp[j] ~ dnorm(mean=mu_sigma_beta_sp,sd=sd_sigma_beta_sp)
    sigma_beta_sp[j] <- exp(log_sigma_beta_sp[j])
  }
  
  
  for(j in 2:A){
    for(p in 1:18){
      beta_sp[j,p] ~ dnorm(mean=0,sd=sigma_beta_sp[j])
    }
  }
  
  for(i in 1:Ti){
    for(j in 2:A){
      r_splines[i,j] <- inprod(Z[i,1:18,j], beta_sp[j,1:18])
      c_splines[i,j] <- BG[i,j]*beta_c[j]
    }
  }
  
  for(i in 1:Ti){
    for(j in 2:A){
      theta_mean[i,j] <- r_splines[i,j]+c_splines[i,j]
    }
  }
  
  # Likelihood for theta
  for(i in 1:Ti){
    for(j in 2:A){
      theta[i,j] ~ dnorm(mean=theta_mean[i,j],sd=sigma_theta[j])
    }
  }
  
  for(j in 2:A){
    for(i in Index_TFR[j]:Ti){
      TFR[i,j] ~ dnorm(mean=TFR_nat[i,j],sd=sigma_TFR[j])
    }
  }
  
  for(j in 2:A){
    for(i in 1:(Index_TFR[j]-1)){
      TFR[Index_TFR[j]-i,j] ~ dnorm(mean=TFR[Index_TFR[j]-i+1,j],sd=sigma_TFR[j])
    }
  }
  
  for(i in 1:Ti){
    for(j in 2:A){
      
      # Priors for beta
      beta[i,1,j] ~ dnorm(mean=0,sd=1)
      beta[i,2,j] ~ dnorm(mean=0,sd=1)
      
      # Compute h from q5
      h[i,j] <- log(q5[i,j])
      
      # Likelihood for q5
      q5[i,j] ~ dbeta(q5_a[i,j], q5_b[i,j])
      
      # Priors for k
      k[i,j] ~ dnorm(mean=0,sd=1)
      
      # Compute mortality rates (mx)
      mx[i,1,j] <- exp(af[1] + bf[1] * h[i,j] + cf[1] * (h[i,j]^2) + vf[1] * k[i,j])
      mx[i,2,j] <- -0.25 * (mx[i,1,j] + log(1 - q5[i,j]))
      mx[i,3:11,j] <- exp(af[3:11] + bf[3:11] * h[i,j] + cf[3:11] * (h[i,j]^2) + vf[3:11] * k[i,j])
      
      
      # Compute lx
      
      lx[i,1,j] <- 1
      lx[i,2,j] <- lx[i,1,j] * exp(-mx[i,1,j])
      lx[i,3,j] <- lx[i,2,j] * exp(-4*mx[i,2,j])
      for(z in 4:12){
        lx[i,z,j] <- lx[i,z-1,j] * exp(-5*mx[i,z-2,j])
      }
      
      # Compute Lx
      Lx[i,1,j] <- 0.5 * (lx[i,1,j] + lx[i,2,j]) + 2 * (lx[i,2,j] + lx[i,3,j])
      Lx[i,2,j] <- 5 * (lx[i,3,j] + lx[i,4,j]) / 2
      for(z in 3:10){
        Lx[i,z,j] <- 5 * (lx[i,z+1,j] + lx[i,z+2,j]) / 2
      }
      
      # Compute gamma
      for(z in 1:7){
        gamma[i,z,j] <- m[z] + X[1,z] * beta[i,1,j] + X[2,z] * beta[i,2,j]
      }
      
      # Compute phi
      phi[i,1:7,j] <- exp(gamma[i,1:7,j]) / sum(exp(gamma[i,1:7,j]))
      
      # Compute Fx
      Fx[i,1,j] <- 0
      Fx[i,2:8,j] <- phi[i,1:7,j] * TFR[i,j] / 5
      
      # Compute Kx
      Kx[i,1:7,j] <- (Lx[i,3:9,j] / Lx[i,4:10,j] * Fx[i,1:7,j] + Fx[i,2:8,j]) * Lx[i,1,j] / 2
      Kx_cov[i,1:7,j] <- Kx[i,1:7,j]*exp((-1)*beta_dummy*dummy_usa[i,j])
      # Adjusted Kx with theta
      Kx_star[i,1:7,j] <- Kx_cov[i,1:7,j] * exp(theta[i,j])
    }
  }
  
  
  # Likelihood for C_true
  for(j in 2:A){
    for(i in Index[j]:Ti){
      C_expected_true[i,j] <- inprod(Kx_cov[i,1:7,j], W_true[i,1:7,j])
      C_true[i,j] ~ dpois(C_expected_true[i,j])
    }
  }
  
  # Likelihood for C
  for(j in 2:A){
    for(i in 1:Ti){
      C_expected[i,j] <- inprod(Kx_star[i,1:7,j], W_gen[i,1:7,j])
      C[i,j] ~ dpois(C_expected[i,j])
    }
  }
  
})


# Create the parallel cluster
cl <- makeCluster(ncores)

# Export necessary objects to the workers
clusterExport(cl, c("model_nimble_inverse_gamma_no_sweden", "data_input", "const", "init_list", "run_single_chain"))
clusterEvalQ(cl, {
  library(nimble)
  library(nimbleHMC)
})

# Run each chain in parallel

chain_results_no_sweden <- parallel::parLapply(
  cl = cl,
  X = 1:nchains,
  fun = function(chain_id) {
    run_single_chain(
      chain_id = chain_id,
      model_code = model_nimble_inverse_gamma_no_sweden,
      data = data_input,
      consts = const,
      init = init_list,
      monitors = c("TFR","theta"),
      niter = 4000,
      nburnin = 2000,
      thin = 10
    )
  }
)

stopCluster(cl)


#### Model without Sweden and England & Wales ####


# Create the parallel cluster
cl <- makeCluster(ncores)

# Export necessary objects to the workers
clusterExport(cl, c("model_nimble_inverse_gamma_no_sweden_eng", "data_input", "const", "init_list", "run_single_chain"))
clusterEvalQ(cl, {
  library(nimble)
  library(nimbleHMC)
})

# Run each chain in parallel

chain_results_no_sweden_eng <- parallel::parLapply(
  cl = cl,
  X = 1:nchains,
  fun = function(chain_id) {
    run_single_chain(
      chain_id = chain_id,
      model_code = model_nimble_inverse_gamma_no_sweden_eng,
      data = data_input,
      consts = const,
      init = init_list,
      monitors = c("TFR","theta"),
      niter = 4000,
      nburnin = 2000,
      thin = 10
    )
  }
)

stopCluster(cl)



#### Model without FamiLinx ####


#Remove parameters from the initilization list #

init_no_gen <- init[!(names(init) %in% c("mu_sigma_beta_sp", "sd_sigma_beta_sp", "beta_sp", 
                                         "log_sigma_theta","sigma_beta_c","log_sigma_beta_sp",
                                         "sigma_TFR",
                                         "beta_c","theta","sigma_theta"))]



# Define initial values for each chain
init_list_no_gen <- list(
  chain1 = init_no_gen,  
  chain2 = init_no_gen,
  chain3 = init_no_gen,
  chain4 = init_no_gen)  

const_no_gen = const
const_no_gen$Index_TFR = NULL

# Create the parallel cluster
cl <- makeCluster(6)


# Export necessary objects to the workers
clusterExport(cl, c("model_nimble_no_gen", "data_input", "const_no_gen", "init_list_no_gen", "run_single_chain"))
clusterEvalQ(cl, {
  library(nimble)
  library(nimbleHMC)
})

# Run each chain in parallel
chain_results_simple_model_no_gen <- parLapply(
  cl = cl,
  X = 1:nchains,
  fun = function(chain_id) {
    run_single_chain(
      chain_id = chain_id,
      model_code = model_nimble_no_gen,
      data = data_input,
      consts = const_no_gen,
      init = init_list_no_gen,
      monitors = c(names(init_list_no_gen$chain1)),
      niter = 4000,
      nburnin = 2000,
      thin = 5
    )
  }
)

stopCluster(cl)


# Store model results in a matrix

post_mat_no_gen <- do.call(rbind, lapply(chain_results_simple_model_no_gen, function(s) as.matrix(s)))


col_TFR =  function(i,j) sprintf("TFR[%d, %d]", i, j)
Years = seq(1751,1910,1)
country = c("SWE","ENG","FRA","NOR","DEN","NLD","USA","FIN")
A = 8
Ti = 160
data_TFR_no_gen = numeric()
for(j in 1:A){
  for (i in const$Index[j]:Ti) {
    temp = data.frame(
      year    = Years[i],
      country = country[j],
      TFR   = post_mat_no_gen[, col_TFR(i, j)]) 
    data_TFR_no_gen = rbind(data_TFR_no_gen,temp) 
  }
}

data_TFR_no_gen_median = data_TFR_no_gen %>%
  group_by(year,country) %>%
  summarise(TFR_median=median(TFR,na.rm=T))

save(data_TFR_no_gen_median,file="results/bayesian_model_result_no_gen.RData")











