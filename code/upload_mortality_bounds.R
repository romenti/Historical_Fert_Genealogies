# Clean the R working environment
source('code/upload_packages.R')
source('code/helpers.R')

##### Read FamiLinx data ####
load("data/child_mortality_data.RData")


france_mort = sapply(france_mortality_final$qx, function(this.q) {
  LearnBayes::beta.select( list(x= this.q*0.5, p=.05), list(x=this.q*2, p=.95))
})

sweden_mort = sapply(sweden_true_mort$qx, function(this.q) {
  LearnBayes::beta.select( list(x= this.q*0.5, p=.05), list(x=this.q*2, p=.95))
})


usa_mort = sapply(usa_mortality_final$qx, function(this.q) {
  LearnBayes::beta.select( list(x= this.q*0.9, p=.05), list(x=this.q*1.1, p=.95))
})

netherlands_mort = sapply(netherlands_mortality_final$qx, function(this.q) {
  LearnBayes::beta.select( list(x= this.q*0.9, p=.05), list(x=this.q*1.1, p=.95))
})


denmark_mort = sapply(denmark_mortality_final$qx, function(this.q) {
  LearnBayes::beta.select( list(x= this.q*0.9, p=.05), list(x=this.q*1.1, p=.95))
})

norway_mort = sapply(norway_mortality_final$qx, function(this.q) {
  LearnBayes::beta.select( list(x= this.q*0.9, p=.05), list(x=this.q*1.1, p=.95))
})

uk_mort = sapply(uk_mortality_final$qx, function(this.q) {
  LearnBayes::beta.select( list(x= this.q*0.9, p=.05), list(x=this.q*1.1, p=.95))
})

finland_mort = sapply(finland_mortality_final$qx, function(this.q) {
  LearnBayes::beta.select( list(x= this.q*0.9, p=.05), list(x=this.q*1.1, p=.95))
})


q5_1_sweden = sweden_mort[1,]
q5_2_sweden = sweden_mort[2,]

q5_1_uk = uk_mort[1,]
q5_2_uk = uk_mort[2,]

q5_1_usa = usa_mort[1,]
q5_2_usa = usa_mort[2,]



q5_1_finland = finland_mort[1,]
q5_2_finland = finland_mort[2,]


q5_1_netherlands = netherlands_mort[1,]
q5_2_netherlands = netherlands_mort[2,]

q5_1_denmark = denmark_mort[1,]
q5_2_denmark = denmark_mort[2,]

q5_1_netherlands = netherlands_mort[1,]
q5_2_netherlands = netherlands_mort[2,]

q5_1_norway = norway_mort[1,]
q5_2_norway = norway_mort[2,]


save(q5_1_sweden,q5_2_sweden,
     q5_1_uk,q5_2_uk,
     q5_1_france,q5_2_france,
     q5_1_norway,q5_2_norway,
     q5_1_denmark,q5_2_denmark,
     q5_1_usa,q5_2_usa,
     q5_1_finland,q5_2_finland,
     q5_1_netherlands,q5_2_netherlands,
     file='data/mortality_bounds.RData')



