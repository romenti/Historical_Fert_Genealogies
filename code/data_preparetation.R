# Clean the R working environment
source('code/upload_packages.R')
source('code/helpers.R')

##### Read FamiLinx data ####
load("data/data_focal.RData")


#### Final sample selection ####

# Death Year>=1751
# Birth Year<=1910
# Age at Death>=0 & Age at Death<=110
# Cleaned data set taken from the publication of 
# Colasurdo and Omenti (2024)


data_final = indicator_demo %>%
  filter(country_birth_final==country_death_final) %>%
  filter(!(is.na(age_death)),
         !(is.na(gender)),
         gender!="unknown",
         country_birth_final %in% c("UK","FRANCE","SWEDEN",
                                    "USA","FINLAND","CANADA",
                                    "DENMARK","NORWAY","NETHERLANDS"),
         death_year>=1751,
         birth_year<=1910,
         (death_year-birth_year)>=0,
         (death_year-birth_year)<=110)



#### Create population pyramids #####

#### Population Pyramids: USA ####
usa_population = pop_pyramid_calculation(data_final,"USA",year_min=1751,year_max=1910)


#### Population Pyramids: UK ####
uk_population = pop_pyramid_calculation(data_final,"UK",year_min=1751,year_max=1910)


#### Population Pyramids: SWEDEN ####
sweden_population = pop_pyramid_calculation(data_final,"SWEDEN",year_min=1751,year_max=1910)


#### Population Pyramids: Denmark ####
denmark_population = pop_pyramid_calculation(data_final,"DENMARK",year_min=1751,year_max=1910)


#### Population Pyramids: Norway ####
norway_population = pop_pyramid_calculation(data_final,"NORWAY",year_min=1751,year_max=1910)



#### Population Pyramids: France ####
france_population = pop_pyramid_calculation(data_final,"FRANCE",year_min=1751,year_max=1910)


#### Population Pyramids: the Netherlands ####
netherlands_population = pop_pyramid_calculation(data_final,"NETHERLANDS",year_min=1751,year_max=1910)

#### Population Pyramids: Finland ####
finland_population = pop_pyramid_calculation(data_final,"FINLAND",year_min=1751,year_max=1910)




#### Child-Woman Ratios ####

#### Child-Woman Ratios: USA ####

CW_usa = CW_ratio_calculation(usa_population_exposure,year_min=1741,year_max=1910)
CW_usa = filter(CW_usa,years>=1751) %>% ungroup()

#### Child-Woman Ratios: Canada ####

CW_canada = CW_ratio_calculation(canada_population_exposure,year_min=1741,year_max=1910)
CW_canada = filter(CW_canada,years>=1751) %>% ungroup()
plot(CW_canada$C)



#### Child-Woman Ratios: Sweden ####
CW_sweden = CW_ratio_calculation(sweden_population,year_min=1751,year_max=1910)


#### Child-Woman Ratios: England ####
CW_eng = CW_ratio_calculation(eng_population,year_min=1741,year_max=1910)

#### Child-Woman Ratios: France ####
CW_france = CW_ratio_calculation(france_population,year_min=1751,year_max=1910)


#### Child-Woman Ratios: Denmark ####
CW_denmark = CW_ratio_calculation(denmark_population,year_min=1751,year_max=1910)

#### Child-Woman Ratios: Norway ####

CW_norway = CW_ratio_calculation(norway_population,year_min=1751,year_max=1910)

#### Child-Woman Ratios: Netherlands ####
CW_netherlands = CW_ratio_calculation(netherlands_population,year_min=1751,year_max=1910)

#### Child-Woman Ratios: Finland ####
CW_finland = CW_ratio_calculation(finland_population,year_min=1751,year_max=1910)

#### Child-Woman Ratios: USA ####
CW_usa = CW_ratio_calculation(usa_population,year_min=1751,year_max=1910)
data_input$C

save(CW_sweden,
     CW_eng,
     CW_usa,
     CW_norway,
     CW_denmark,
     CW_finland,
     CW_france,
     CW_netherlands,
     file="data/genealogical_data_cw.RData")





