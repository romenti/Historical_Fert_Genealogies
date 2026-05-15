#### upload all the relevant packages ####
source('code/upload_packages.R')

#### upload relevant data sets ####

# Bayesian model results
load('results/bayesian_model_result.RData')

# Bayesian model results without FamiLinx
load('results/bayesian_model_result_no_gen.RData')

# Indirect estimation results
load('results/indirect_TFR_results.RData')

# Historical TFR values
load('data/historical_TFR.RData')


table_tfr_accuracy = hist_data %>%
  left_join(data_TFR_no_gen_median %>%
              select(years=year,country,TFR_median)) %>%
  mutate(e=abs(TFR_median-TFR_true)/TFR_true*100) %>%
  group_by(country) %>%
  summarise(mae_bTFR_plus=mean(e,na.rm=T)) %>%
  left_join(hist_data %>%
              left_join(mcmc.out_TFR_inverse_gamma %>%
                          select(years,country,TFR_median)) %>%
              mutate(e=abs(TFR_median-TFR_true)/TFR_true*100) %>%
              group_by(country) %>%
              summarise(mae_bTFR_star=mean(e,na.rm=T))) %>%
  left_join(hist_data %>%
              left_join(data_iTFR) %>%
              mutate(e=abs(iTFR_star-TFR_true)/TFR_true*100) %>%
              group_by(country) %>%
              summarise(mae_iTFR_star=mean(e,na.rm=T))) %>%
  left_join(hist_data %>%
              left_join(data_xTFR) %>%
              mutate(e=abs(xTFR_star-TFR_true)/TFR_true*100) %>%
              group_by(country) %>%
              summarise(mae_xTFR_star=mean(e,na.rm=T))) %>%
  mutate(mae_bTFR_plus=round(mae_bTFR_plus,3),
         mae_bTFR_star=round(mae_bTFR_star,3),
         mae_iTFR_star=round(mae_iTFR_star,3),
         mae_xTFR_star=round(mae_xTFR_star,3)) %>%
  select(country,mae_bTFR_star,mae_bTFR_plus,
         mae_iTFR_star,mae_xTFR_star) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    col.names = c("Country", "bTFR+", "bTFR$^{*}$", "iTFR$^{*}$", "xTFR$^{*}$"),
    caption = "Mean absolute relative error (\\%) of fertility measures relative to historical TFR"
  ) %>%
  kable_styling(latex_options = c("hold_position"))



writeLines(table_tfr_accuracy, "tables/table1_b.tex")


