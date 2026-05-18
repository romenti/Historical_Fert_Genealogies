# Clean the R working environment
source('code/upload_packages.R')
source('code/helpers.R')

##### Read Child Mortality data ####
load("data/mortality_bounds.RData")
load("data/genealogical_data_cw.RData")
load("data/wilmoth_constants.RData")
load("data/data_accurate.RData")
load("data/svd_input.RData")
load("data/historical_tfr.RData")


#### Spline Construction ####

years <- 1751:1910

# Define knots at every 10 years from 1760 to 1900
knots <- seq(1760, 1900, by = 10)

# Generate the B-spline basis using cubic splines (degree = 3)
B <- splines2::bSpline(years,
                       knots = knots,
                       degree = 3,
                       intercept = TRUE,
                       Boundary.knots = c(1751, 1910))

# Convert B-spline matrix to data frame for plotting
B_df <- as.data.frame(B)
B_df$Year <- years

# Convert to long format for ggplot
B_long <- pivot_longer(B_df, cols = -Year, names_to = "Basis", values_to = "Value")


ktemp <- ncol(B)
D1 <- diff(diag(ktemp), differences = 1)

# Compute the transformation matrix Dcomb
Dcomb <- t(D1) %*% solve(D1 %*% t(D1))

# Compute the 'random' spline component by transforming the basis
Z <- B %*% Dcomb

# Compute a fixed spline component as a simple weighted sum (here weights = 1)
G <- rep(1, ktemp)
BG <- B %*% G


Z = array(rep(Z,8),dim=c(160,18,8))

BG = matrix(rep(BG,8),160,8)


data_input_new$X_sp


dummy_usa = matrix(0,nrow=160,ncol=8)
dummy_usa[112:116,7] = 1


#### Read historical TFR values ####

tfr_hist_sweden = hist_data %>%
  filter(country=="SWE") %>%
  arrange(years) %>%
  tidyr::complete(years = full_seq(years, 1)) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(., x = years, na.rm = FALSE))) %>%
  select(TFR_true) %>%
  as.matrix()

tfr_hist_eng = hist_data %>%
  filter(country=="ENG") %>%
  arrange(years) %>%
  tidyr::complete(years = full_seq(years, 1)) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(., x = years, na.rm = FALSE))) %>%
  select(TFR_true) %>%
  as.matrix()


tfr_hist_france = hist_data %>%
  filter(country=="FRA") %>%
  arrange(years) %>%
  tidyr::complete(years = full_seq(years, 1)) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(., x = years, na.rm = FALSE))) %>%
  select(TFR_true) %>%
  as.matrix()


tfr_hist_norway = hist_data %>%
  filter(country=="NOR") %>%
  arrange(years) %>%
  #tidyr::complete(years = full_seq(years, 1)) %>%
  #mutate(across(where(is.numeric), ~ zoo::na.approx(., x = years, na.rm = FALSE))) %>%
  select(TFR_true) %>%
  as.matrix()



tfr_hist_denmark = hist_data %>%
  filter(country=="DEN") %>%
  arrange(years) %>%
  tidyr::complete(years = full_seq(years, 1)) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(., x = years, na.rm = FALSE))) %>%
  select(TFR_true) %>%
  as.matrix()




tfr_hist_netherland = hist_data %>%
  filter(country=="NLD") %>%
  arrange(years) %>%
  tidyr::complete(years = full_seq(years, 1)) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(., x = years, na.rm = FALSE))) %>%
  select(TFR_true) %>%
  as.matrix()



tfr_hist_usa = hist_data %>%
  filter(country=="USA") %>%
  arrange(years) %>%
  tidyr::complete(years = full_seq(years, 1)) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(., x = years, na.rm = FALSE))) %>%
  select(TFR_true) %>%
  as.matrix()



tfr_hist_finland = hist_data %>%
  filter(country=="FIN") %>%
  arrange(years) %>%
  tidyr::complete(years = full_seq(years, 1)) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(., x = years, na.rm = FALSE))) %>%
  select(TFR_true) %>%
  as.matrix()




#### Data Input ####


data_input = list(C = matrix(cbind(as.integer(CW_sweden$C),
                              as.integer(CW_eng$C),
                              as.integer(CW_france$C),
                              as.integer(CW_norway$C),
                              as.integer(CW_denmark$C),
                              as.integer(CW_netherlands$C),
                              as.integer(CW_usa$C),
                              as.integer(CW_finland$C)),nrow=160,ncol=8),
                  W_gen = array(
                    rbind(
                      as.matrix(CW_sweden[, paste0('W', seq(15, 45, 5))]),
                      as.matrix(CW_eng[, paste0('W', seq(15, 45, 5))]),
                      as.matrix(CW_france[, paste0('W', seq(15, 45, 5))]),
                      as.matrix(CW_norway[, paste0('W', seq(15, 45, 5))]),
                      as.matrix(CW_denmark[, paste0('W', seq(15, 45, 5))]),
                      as.matrix(CW_netherlands[, paste0('W', seq(15, 45, 5))]),
                      as.matrix(CW_usa[, paste0('W', seq(15, 45, 5))]),
                      as.matrix(CW_finland[, paste0('W', seq(15, 45, 5))])
                    ),
                    dim = c(160, 7, 8)
                  ),
                  C_true=matrix(cbind(as.integer(C_sweden_true),
                                 as.integer(C_eng_true),
                                 c(rep(0,65),as.integer(C_france_true)),
                                 c(rep(0,95),as.integer(C_norway_true)),
                                 c(rep(0,84),as.integer(C_denmark_true)),
                                 c(rep(0,99),as.integer(C_netherlands_true)),
                                 c(rep(0,79),as.integer(C_usa_true)),
                                 c(rep(0,114),as.integer(C_finland_true))),
                                nrow=160,ncol=8),
                  W_true = array(
                    rbind(
                      W_sweden_true,
                      W_eng_true,
                      rbind(matrix(0,nrow=65,ncol=7),
                            W_france_true),
                      rbind(matrix(0,nrow=95,ncol=7),
                            W_norway_true),
                      rbind(matrix(0,nrow=84,ncol=7),
                            W_denmark_true),
                      rbind(matrix(0,nrow=99,ncol=7),
                            W_netherlands_true),
                      rbind(matrix(0,nrow=79,ncol=7),
                            W_usa_true),
                      rbind(matrix(0,nrow=114,ncol=7),
                            W_finland_true)
                    ),
                    dim = c(160, 7, 8)
                  ),
                  q5_a = matrix(cbind(q5_1_sweden,
                               q5_1_uk,
                               q5_1_france,
                               q5_1_norway,
                               q5_1_denmark,
                               q5_1_netherlands,
                               q5_1_usa,
                               q5_1_finland
                               ),nrow=160,ncol=8),
                  q5_b = matrix(cbind(q5_2_sweden,
                                      q5_2_uk,
                                      q5_2_france,
                                      q5_2_norway,
                                      q5_2_denmark,
                                      q5_2_netherlands,
                                      q5_2_usa,
                                      q5_2_finland
                  ),nrow=160,ncol=8),
                  Z=Z,
                  BG=BG,
                  dummy_usa=dummy_usa,
                  m=m,
                  X=X,
                  af=wilmoth$af[1:11],
                  bf=wilmoth$bf[1:11],
                  cf=wilmoth$cf[1:11],
                  vf=wilmoth$vf[1:11],
                  TFR_nat = matrix(cbind(tfr_hist_sweden,
                                         c(rep(0,49),tfr_hist_eng),
                                         c(rep(0,49),tfr_hist_france),
                                         c(rep(0,84),tfr_hist_norway),
                                         c(rep(0,99),tfr_hist_denmark),
                                         c(rep(0,89),tfr_hist_netherland),
                                         c(rep(0,49),tfr_hist_usa),
                                         c(rep(0,49),tfr_hist_finland)),
                                     nrow=160,ncol=8)
                  )


#### Constants for the Model ####

# A = number of countries

# Ti = number of years

# Index = first year with accurate population data by country


# Index = first year with historical TFR value by country


const = list(A = 8,
             Ti = 160,
             Index = c(1,1,66,96,85,100,80,115),
             Index_TFR = c(1,50,50,85,100,90,50,50))



#### Initialize Model Parameters #####
A = 8
Ti = 160

init = list(q5 = matrix(rbeta(Ti*A,shape1 = c(q5_1_sweden,
                                      q5_1_uk,
                                      q5_1_france,
                                      q5_1_norway,
                                      q5_1_denmark,
                                      q5_1_netherlands,
                                      q5_1_usa,
                                      q5_1_finland),
                              shape2 =c(q5_2_sweden,
                                        q5_2_uk,
                                        q5_2_france,
                                        q5_2_norway,
                                        q5_2_denmark,
                                        q5_2_netherlands,
                                        q5_2_usa,
                                        q5_2_finland)),nrow=Ti,ncol=A),
            TFR = matrix(pmax(.10, rnorm(Ti*A, 5, sd=.50)),nrow=Ti,ncol=A),
            sigma_TFR=rep(2,8),
            k = matrix(rnorm(Ti*A,mean=0,sd=0.1),nrow=Ti,ncol=A),
            beta = array(rep(matrix(rnorm(Ti*2,mean=0,sd=0.1),nrow=Ti,ncol=2),A),
                         dim=c(Ti,2,A)),
            theta = matrix(0,nrow=Ti,ncol=A),
            beta_c = rep(0,A),
            beta_dummy = -0.5,
            mu_sigma_beta_sp = 1,
            sd_sigma_beta_sp = 2,
            log_sigma_beta_sp = rep(2,A),
            beta_sp = matrix(0,nrow=A,ncol=18),
            sigma_theta = rep(2,A),
            sigma_beta_c = 2)


save(data_input,
     init,
     const,file="data/model_input_data.RData")
