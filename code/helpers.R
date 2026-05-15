#### Calculate Population Pyramids from FamiLinx ####

pop_pyramid_calculation = function(data,country_name,year_min=NA,year_max=NA){
  
  if(is.na(country_name) | is.null(data)){
    return(print("Enter the genealogical data set and the country name"))
  } else{
    pop_final = NULL
    pop_final = data %>%
      filter(!(is.na(gender)),
             gender!="unknown",
             country_birth_final==country_name) %>%
      mutate(years = purrr::map2(birth_year,death_year,~seq(.x,.y,by=1))) %>%
      unnest(years) %>%
      filter(years>=year_min,
             years<=year_max) %>%
      mutate(age = years-birth_year,
             age_class = case_when(age==0 ~ "0",
                                   age>=1 & age<=4 ~ "1-4",
                                   
                                   age>=5 & age<=9 ~ "5-9",
                                   age>=10 & age<=14 ~ "10-14",
                                   age>=15 & age<=19 ~ "15-19",
                                   age>=20 & age<=24 ~ "20-24",
                                   age>=25 & age<=29 ~ "25-29",
                                   age>=30 & age<=34 ~ "30-34",
                                   age>=35 & age<=39 ~ "35-39",
                                   age>=40 & age<=44 ~ "40-44",
                                   age>=45 & age<=49 ~ "45-49",
                                   age>=50 & age<=54 ~ "50-54",
                                   age>=55 & age<=59 ~ "55-59",
                                   age>=60 & age<=64 ~ "60-64",
                                   age>=65 & age<=69 ~ "65-69",
                                   age>=70 & age<=74  ~ "70-74",
                                   age>=75 & age<=79 ~ "75-79",
                                   age>=80  ~ "80+"
             )) %>%
      mutate(age_class = factor(age_class, levels = c("0","1-4","5-9","10-14","15-19",
                                                      "20-24","25-29","30-34","35-39",
                                                      "40-44","45-49","50-54","55-59",
                                                      "60-64","65-69","70-74","75-79",
                                                      "80+")),
             gender = factor(gender,levels = c("male","female"))) %>%
      group_by(years,age_class,gender,.drop=FALSE) %>%
      dplyr::tally() %>%
      #complete(years,age_class,gender, fill = list(n = 0)) %>%
      mutate(country=country_name) %>%
      pivot_wider(names_from = "gender", values_from = "n") %>%
      mutate(total = male+female) %>%
      pivot_longer(!c("years","country","age_class"),names_to = "gender", values_to = "counts")
    return(pop_final)
    
  }
  
}

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


