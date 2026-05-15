#### Run Model in Parallel ####
run_single_chain <- function(chain_id, model_code, data, consts, 
                             init, monitors, niter, nburnin, thin) {
  library(nimble)
  library(nimbleHMC)
  
  # Initialize the model
  dynModel <- nimbleModel(
    code = model_code,
    data = data, 
    constants = consts, 
    buildDerivs = TRUE,
    inits = init[[chain_id]]  # Use chain-specific initial values
  )
  
  # Run HMC
  model_output <- nimbleHMC(
    dynModel, 
    data = data, 
    inits = init[[chain_id]],
    monitors = monitors, 
    thin = thin,
    niter = niter, 
    nburnin = nburnin,
    summary = TRUE,
    WAIC = FALSE,
    progressBar = getNimbleOption("MCMCprogressBar")
  )
  
  return(model_output$samples)  # Return MCMC samples
}


