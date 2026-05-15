#### Paper Figures ####


#### Figure 1 ####

load('data/gen_population_pyramids.RData')

plot_pyramid = data_pyramid %>%
  filter(years %in% c(1751,1800,1850,1900),
         country %in% c('USA','SWE','ENG')) %>%
  mutate(count_gen=ifelse(gender=="female",count_gen*(-1),count_gen)) %>%
  ggplot()+  # default x-axis is age in years;
  # case data graph
  geom_col(mapping = aes(
    x = Age,
    y = count_gen,
    fill = gender),         
    colour = "white")+       # white around each bar
  geom_step(data =  data_pyramid   %>% filter(gender == "female",country %in% c('ENG','SWE','USA')), 
            aes(x = Age, y= expected, color = label),size=1) +
  geom_step(data =  data_pyramid   %>% filter(gender == "male",country %in% c('ENG','SWE','USA')), 
            aes(x = Age, y = (-1)*expected,color=label),size=1) +
  
  # flip the X and Y axes to make pyramid vertical
  coord_flip()+
  scale_fill_manual(name = "Sex",
                    values = c("female" = "#990099",
                               "male" = "#009900"),
                    labels = c("Female", "Male"))+
  scale_color_manual(values = c('true population' = "black"),name="",labels = c("Expected Population Structure"),na.translate = F)+
  ggthemes::theme_fivethirtyeight()+
  scale_y_continuous(labels = abs)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = 'bottom',
        rect = element_rect(fill = 'white'),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        strip.text.x = element_text(face="bold"), 
        strip.text.y = element_text(face="bold"),
        aspect.ratio = 0.66)+
  facet_wrap(country~years,scales = "free",nrow=3)+
  ylab("Genealogy Population Size")+
  xlab("Age")

ggsave("figures/Figure1.pdf",plot_pyramid,height = 15, width = 25, units = "cm",dpi=700)

#### Figure 2 ####


load('data/motivation_datasets.RData')

plot_missing_age = missing_age %>%
  ggplot(aes(x = years, y = perc*100)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = c(1751, seq(1770, 1910, 20)),
                     labels = c(1751, seq(1770, 1910, 20))) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +  # <- ADD THIS LINE
  facet_wrap(~country, nrow = 2) +
  ggthemes::theme_fivethirtyeight(base_size = 18) +
  ylab("% of Children under 5\nwith Missing Mother Age") +
  xlab("Year") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "bottom",
    legend.text = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 16),
    rect = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text.x = element_text(face = "bold", size = 16),
    strip.text.y = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("figures/Figure2.pdf",plot_missing_age,height = 15, width = 25, units = "cm",dpi=700)

#### Figure 4 ####

load("results/model_results.RData")

plot_TFR = mcmc.out_TFR_inverse_gamma %>%
  mutate(country = factor(country,levels = c("SWE","ENG","FRA","DEN",
                                             "NOR","NLD","USA","FIN"))) %>%
  ggplot()+
  geom_line(aes(x=years,y=TFR_median,shape="bTFR*",color="bTFR*",linetype = "bTFR*"),size=2)+
  geom_line(data= data_iTFR,aes(x=years,y=iTFR_star,color="iTFR*",linetype="iTFR*"),size=2)+
  #geom_line(data= data_xTFR,aes(x=years,y=xTFR_star,color="xTFR*",linetype="xTFR*"),size=2)+
  geom_point(data = hist_data %>%
               filter(years <= 1910,
                      years %% 5 == 0),aes(x=years,y=TFR_true,
                                           shape="Historical Estimates",color="Historical Estimates"),size=3)+
  geom_ribbon(data=mcmc.out_TFR_uniform ,aes(x=years,ymin = TFR_lower, ymax = TFR_upper),  linetype=2, alpha = 0.3,fill="grey")+
  
  scale_x_continuous(breaks=c(1751,seq(1770,1910,20)),labels=c(1751,seq(1770,1910,20)))+
  scale_colour_manual(name = "",
                      values = c("bTFR*"="#0072B2", "Historical Estimates"="#FF0000","iTFR*"="#009E73")) +  
  scale_linetype_manual(name="",
                        values=c("bTFR*"="solid","iTFR*"="twodash","xTFR*"='dotted'))+
  scale_shape_manual(name = "",
                     values = c("Historical Estimates"=4))+
  facet_wrap(~country,nrow=2)+
  ggthemes::theme_fivethirtyeight(base_size = 18)+ 
  #coord_cartesian(ylim=c(2,9))+
  ylab("TFR Estimates")+
  xlab("Year")+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=14),
        axis.text.y = element_text(size=14),
        legend.position = "none",
        legend.text = element_text(face="bold",size=16),
        legend.title  = element_text(face="bold",size=16),
        rect = element_rect(fill = 'white'),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        strip.background = element_rect(fill = "white", colour = "white"),
        aspect.ratio = 0.66,
        strip.text.x = element_text(face="bold",size=16),
        strip.text.y = element_text(face="bold",size=16),
        axis.title.y = element_text(face="bold",size=16),
        axis.title.x = element_text(face="bold",size=16),
        plot.title = element_text(hjust=0.5,face="bold")
  ) 



legend_data <- data.frame(
  x = c(1, 2,1,2, 1, 2, 1, 1),  # Duplicate x values for line groups
  y = c(1, 1,1,1, 1, 1, 1, 1),
  type = c("iTFR*", "iTFR*","xTFR*","xTFR*" ,"bTFR*", "bTFR*", "Historical Estimates", "Historical Estimates")
)


manual_legend <- ggplot(legend_data, aes(x, y, group = type, linetype = type, color = type, shape = type)) +
  geom_line(size = 1, data = subset(legend_data, type != "Historical Estimates")) +  # Lines only for relevant types
  geom_point(size = 3, data = subset(legend_data, type == "Historical Estimates")) +  # Points for historical estimates
  scale_color_manual(name = "",
                     values = c("bTFR*"="#0072B2", "Historical Estimates"="#FF0000","iTFR*"="#009E73")) +
  scale_linetype_manual(name = "",
                        values = c("iTFR*" = "twodash",
                                   "bTFR*" = "solid",
                                   "Historical Estimates" = "blank")) +
  scale_shape_manual(name = "",
                     values = c("iTFR*" = NA,
                                "bTFR*" = NA,
                                "Historical Estimates" = 4)) +  
  theme_void() +
  theme(legend.text = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=16),
        legend.key.width = unit(2, "cm"),
        legend.spacing.y = unit(2, "cm")
  ) +
  guides(
    color = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  )

legend_grob <- cowplot::get_legend(manual_legend)
legend_grob <- get_legend(manual_legend)


plot_TFR_final = plot_TFR/legend_grob+
  plot_layout(heights = c(10, 2))  


ggsave('figure/Figure4.pdf',plot_TFR_final,height = 30, width = 50, units = "cm")


#### Figure 5 ####

mcmc.out_theta_inverse_gamma = mcmc.out_theta_final

bias_plot = mcmc.out_theta_final %>%
  ggplot()+
  geom_line(aes(x=years,y=theta_median,
                color='Bayesian Bias-Adjustement Multiplier',
                linetype='Bayesian Bias-Adjustement Multiplier'),
            size=2)+
  geom_line(data=data_multiplier,aes(x=years,y=log(multiplier),
                                     color='Indirect Bias-Adjustment Multiplier',
                                     linetype='Indirect Bias-Adjustment Multiplier'),size=2)+
  geom_hline(aes(yintercept=0),linetype='dotted',size=2)+
  #geom_point(data=canada_true_fertility,aes(x=years,y=TFR_true,shape="Historical Estimates"))+
  geom_ribbon(aes(x=years,ymin = theta_lower, ymax = theta_upper),  linetype=2, alpha = 0.5,fill="#377EB8")+
  theme_classic()+
  scale_color_manual(name='',values=c('Bayesian Bias-Adjustement Multiplier'='#377EB8',
                                      'Indirect Bias-Adjustment Multiplier'="#FFB000"))+
  scale_linetype_manual(name='',values=c('Bayesian Bias-Adjustement Multiplier'='solid',
                                         'Indirect Bias-Adjustment Multiplier'="twodash"))+
  scale_x_continuous(breaks=c(1751,seq(1770,1910,20)),labels=c(1751,seq(1770,1910,20)))+
  coord_cartesian(ylim=c(-1,0.35))+
  xlab('Year')+
  ylab('Bias-adjustment Multiplier')+
  ggthemes::theme_fivethirtyeight(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=14),
        axis.text.y = element_text(size=14),
        legend.position = "bottom",
        aspect.ratio = 0.66,
        legend.text = element_text(face="bold",size=16),
        legend.title  = element_text(face="bold",size=16),
        #legend.background = element_rect(fill="lightgray", size=0.5, linetype="solid"),
        strip.text.x = element_text(face="bold",size=16),
        strip.text.y = element_text(face="bold",size=16),
        axis.title.y = element_text(face="bold",size=16),
        axis.title.x = element_text(face="bold",size=16),
        plot.title = element_text(hjust=0.5,face="bold"),
        rect = element_rect(fill = 'white'),
        plot.background = element_rect(fill = "white", colour = "white")
  ) +
  facet_wrap(~country,nrow=2)


ggsave('figure/Figure5.pdf',bias_plot,height = 30, width = 50, units = "cm")


