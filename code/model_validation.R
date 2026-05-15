#### Upload packages ####

# upload all relevant packages
source('code/upload_packages.R')
source('code/helpers.R')
source('code/model_code.R')


#### 30-Year Validation ####

# load model input data
load('data/model_validation_input.RData')

data_input_30year = data_input

data_input_30year$C_true[131:160,] = NA
data_input_30year$W_true[131:160,,] = NA

# Define initial values for each chain
init_list_half <- list(
  chain1 = init,  
  chain2 = init,
  chain3 = init,
  chain4 = init)  

# Create the parallel cluster
ncores = 4
cl <- makeCluster(ncores)


# Export necessary objects to the workers
clusterExport(cl, c("model_nimble_inverse_gamma", "data_input_30year", "const", "init_list", "run_single_chain"))
clusterEvalQ(cl, {
  library(nimble)
  library(nimbleHMC)
})

# Run each chain in parallel
chain_results_30year <- parallel::parLapply(
  cl = cl,
  X = 1:nchains,
  fun = function(chain_id) {
    run_single_chain(
      chain_id = chain_id,
      model_code = model_nimble_inverse_gamma,
      data = data_input_30year,
      consts = const,
      init = init_list,
      monitors = c("Kx_cov"),
      niter = 6000,
      nburnin = 3000,
      thin = 10
    )
  }
)

stopCluster(cl)

post_mat_30years <- do.call(rbind, lapply(chain_results_30year, function(s) as.matrix(s)))


col_Kx_cov   <- function(i,j,z) sprintf("Kx_cov[%d, %d, %d]", i, j, z)
Age = seq(15,45,5)
Years = seq(1751,1910,1)

country = c("SWE","ENG","FRA","NOR","DEN","NLD","USA","FIN")
data_pred_30years = numeric()


for(i in 131:160){
  for(j in 1:A){
    for(z in 1:7){
      temp = data.frame(year=Years[i],
                        country=country[j],
                        age=Age[z],
                        Kx=post_mat_30years[,col_Kx_cov(i,z,j)]*data_input$W_true)
      data_pred_30years = rbind(data_pred_30years,temp)
    }
  }
}

pred_tot_30year <- data_pred_30year %>%
  group_by(year, country, age) %>%
  mutate(draw = row_number()) %>%     # draw index within each (year,country,age)
  ungroup() %>%
  group_by(year, country, draw) %>%
  summarise(lambda = sum(Kx), .groups = "drop")





set.seed(1)
pred_tot_30year <- pred_tot_30year %>% mutate(C_rep = rpois(n(), lambda))


score_30year = pred_tot_30year %>%
  filter(year>1880) %>%
  group_by(year, country) %>%
  left_join(data_CWR_true) %>%
  mutate(CW_rep=C_rep/W_true) %>%
  summarise(
    CW_true = first(CW_true),
    r_med  = median(CW_rep),
    r_l80  = quantile(CW_rep, alpha80/2),        # 0.30
    r_u80  = quantile(CW_rep, 1 - alpha80/2),    # 0.90
    r_l90  = quantile(CW_rep, alpha90/2),        # 0.30
    r_u90  = quantile(CW_rep, 1 - alpha90/2),    # 0.90
    r_l95  = quantile(CW_rep, alpha95/2),        # 0.025
    r_u95  = quantile(CW_rep, 1 - alpha95/2),    # 0.975
    .groups="drop"
  ) %>%
  mutate(
    cover80 = as.integer(CW_true >= r_l80 & CW_true <= r_u80),
    cover90 = as.integer(CW_true >= r_l90 & CW_true <= r_u90),
    cover95 = as.integer(CW_true >= r_l95 & CW_true <= r_u95),
    rel_err = abs(CW_true - r_med) / r_med*100
  )

cover80_30year = mean(score_30year$cover80)*100
cover95_30year = mean(score_30year$cover95)*100
cover90_30year = mean(score_30year$cover90)*100
rel_err_30year = mean(score_30year$rel_err)


#### 20-Year Validation ####

# load model input data
load('data/model_data_input.RData')


data_input_20year = data_input

data_input_20year$C_true[141:160,] = NA
data_input_20year$W_true[141:160,,] = NA

# Define initial values for each chain
init_list_half <- list(
  chain1 = init,  
  chain2 = init,
  chain3 = init,
  chain4 = init)  

# Create the parallel cluster
ncores = 4
cl <- makeCluster(ncores)


# Export necessary objects to the workers
clusterExport(cl, c("model_nimble_inverse_gamma", "data_input_20year", "const", "init_list", "run_single_chain"))
clusterEvalQ(cl, {
  library(nimble)
  library(nimbleHMC)
})

# Run each chain in parallel
chain_results_20year <- parallel::parLapply(
  cl = cl,
  X = 1:nchains,
  fun = function(chain_id) {
    run_single_chain(
      chain_id = chain_id,
      model_code = model_nimble_inverse_gamma,
      data = data_input_20year,
      consts = const,
      init = init_list,
      monitors = c("Kx_cov"),
      niter = 6000,
      nburnin = 2000,
      thin = 10
    )
  }
)

stopCluster(cl)

post_mat_20years <- do.call(rbind, lapply(chain_results_20year, function(s) as.matrix(s)))



col_Kx   <- function(i,j,z) sprintf("Kx_cov[%d, %d, %d]", i, j, z)
Age = seq(15,45,5)
Years = seq(1751,1910,1)
country = c("SWE","ENG","FRA","NOR","DEN","NLD","USA","FIN")
data_pred_20year = numeric()


for(i in 141:160){
  for(j in 1:A){
    for(z in 1:7){
      temp = data.frame(year=Years[i],
                        country=country[j],
                        age=Age[z],
                        Kx=post_mat_20years[,col_Kx_cov(i,z,j)]*data_input$W_true)
      data_pred_20year = rbind(data_pred_20year,temp)
    }
  }
}



pred_tot_20year <- data_pred_20years %>%
  group_by(year, country, age) %>%
  mutate(draw = row_number()) %>%     # draw index within each (year,country,age)
  ungroup() %>%
  group_by(year, country, draw) %>%
  summarise(lambda = sum(Kx), .groups = "drop")





set.seed(1)
pred_tot_20year <- pred_tot_20year %>% mutate(C_rep = rpois(n(), lambda))



score_20year = pred_tot_20year %>%
  filter(year>1890) %>%
  group_by(year, country) %>%
  left_join(data_CWR_true) %>%
  mutate(CW_rep=C_rep/W_true) %>%
  summarise(
    CW_true = first(CW_true),
    r_med  = median(CW_rep),
    r_l80  = quantile(CW_rep, alpha80/2),        # 0.30
    r_u80  = quantile(CW_rep, 1 - alpha80/2),    # 0.90
    r_l90  = quantile(CW_rep, alpha90/2),        # 0.30
    r_u90  = quantile(CW_rep, 1 - alpha90/2),    # 0.90
    r_l95  = quantile(CW_rep, alpha95/2),        # 0.025
    r_u95  = quantile(CW_rep, 1 - alpha95/2),    # 0.975
    .groups="drop"
  ) %>%
  mutate(
    cover80 = as.integer(CW_true >= r_l80 & CW_true <= r_u80),
    cover90 = as.integer(CW_true >= r_l90 & CW_true <= r_u90),
    cover95 = as.integer(CW_true >= r_l95 & CW_true <= r_u95),
    rel_err = abs(CW_true - r_med) / r_med*100
  )

cover80_20year = mean(score_20year$cover80)*100
cover95_20year = mean(score_20year$cover95)*100
cover90_20year = mean(score_20year$cover90)*100
rel_err_20year = mean(score_20year$rel_err)


validation_table = data.frame(validation_set=c("30-year",
                            "20-year"),
           mare = c(rel_err_30year,
                    rel_err_20year),
           cov80 = c(cover80_30year,
                     cover80_20year),
           cov90 = c(cover90_30year,
                     cover90_20year),
           cov95 = c(cover95_30year,
                     cover95_20year)
           ) 
  



  
table1_a_latex = kable(validation_table, format = "latex", booktabs = F,
                         align = c("l", rep("r", 4)),
                         linesep='',
                         col.names = c("Validation Set",c("MARE", "80% Coverage", 
                                                          "90% Coverage", "95% Coverage"))) %>%
  kable_styling(latex_options = "hold_position")


writeLines(table1_a_latex, "tables/table1_a.tex")






