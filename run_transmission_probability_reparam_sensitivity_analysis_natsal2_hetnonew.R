rm(list=ls())

# setwd('/Users/Joanna/OneDrive - Imperial College London/backup/papers/transmission_probs/ct_transmission_prob')

library(foreign)
library(rstan)
library(survey)
library(deSolve)
library(boot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

options(survey.lonely.psu="adjust")

# Read in Natsal-2 data

# Edit this line to point to your Natsal-2 dataset
n2 <- read.dta("/Users/Joanna/OneDrive - Imperial College London/backup/M_Genitalium/Natsal-2 data/stata11/natsal_2000_for_archive.dta")
n2 <- n2[(n2$c_result %in% c(0,1)),] # keep only respondents with a chlamydia test result

n2design <- svydesign(
  ids = ~ totalpsu,
  weights = ~ urine_wt,
  strata = ~ totalstr,
  data = n2
)

####################################
# tabulate number of new partners in the last year
####################################

##############
# men
##############

npart_m <- n2$hetnonew[
  (n2$rsex == 'male') & 
    (n2$het == 'het. sex. intercourse aged 13+') & 
    (n2$het1yr %in% 1:901) &
    (n2$hetnonew != 999)
  ]

ntab_m <- table(npart_m)
ntab_m <- data.frame(ntab_m)
names(ntab_m) <- c('n', 'freq')
ntab_m$prop <- ntab_m$freq/sum(ntab_m$freq)
ntab_m$low.95 <- ntab_m$high.95 <- ntab_m$prop_wt <- NA
for(i in 1:(nrow(ntab_m))){
  ntab_m$low.95[i] <- binom.test(ntab_m$freq[i], sum(ntab_m$freq))$conf.int[[1]]
  ntab_m$high.95[i] <- binom.test(ntab_m$freq[i], sum(ntab_m$freq))$conf.int[[2]]
  ntab_m$prop_wt[i] <- svymean(
    ~ (hetnonew) == as.numeric(as.character(ntab_m$n[i])),
    design = subset(n2design, 
                    (n2$rsex == 'male') & 
                      (n2$het == 'het. sex. intercourse aged 13+') & 
                      (n2$het1yr %in% 1:901)
    )
  )[2]
}

##############
# women
##############

npart_f <- n2$hetnonew[
  (n2$rsex == 'female') & 
    (n2$het == 'het. sex. intercourse aged 13+') & 
    (n2$het1yr %in% 1:901) &
    (n2$hetnonew != 999)
  ]
ntab_f <- table(npart_f)
ntab_f <- data.frame(ntab_f)
names(ntab_f) <- c('n', 'freq')
ntab_f$prop <- ntab_f$freq/sum(ntab_f$freq)
ntab_f$low.95 <- ntab_f$high.95 <- ntab_f$prop_wt <- NA
for(i in 1:(nrow(ntab_f))){
  ntab_f$low.95[i] <- binom.test(ntab_f$freq[i], sum(ntab_f$freq))$conf.int[[1]]
  ntab_f$high.95[i] <- binom.test(ntab_f$freq[i], sum(ntab_f$freq))$conf.int[[2]]
  ntab_f$prop_wt[i] <- svymean(
    ~ (hetnonew) == as.numeric(as.character(ntab_f$n[i])),
    design = subset(n2design, 
                    (n2$rsex == 'female') & 
                      (n2$het == 'het. sex. intercourse aged 13+') & 
                      (n2$het1yr %in% 1:901)
    )
  )[2]
}

####################################
# tabulate partner numbers and infection status
####################################

for(i in sort(unique(npart_m))){
  
  print(i)
  print(table(n2$c_result[(n2$rsex == 'male') & 
                          (n2$het == 'het. sex. intercourse aged 13+') & 
                          (n2$het1yr %in% 1:901) &
                          (n2$hetnonew == i)]))
  print(100*svyciprop(~(c_result == 1), 
                       design = subset(n2design, 
                                       (n2$rsex == 'male') & 
                                         (n2$het == 'het. sex. intercourse aged 13+') & 
                                         (n2$het1yr %in% 1:901) &
                                         (n2$hetnonew == i) & 
                                         (n2$c_result %in% c(0,1))
                       ),
                       method = "beta"
  )
  )
  print(100*attr(svyciprop(~(c_result == 1), 
                  design = subset(n2design, 
                                  (n2$rsex == 'male') & 
                                  (n2$het == 'het. sex. intercourse aged 13+') & 
                                  (n2$het1yr %in% 1:901) &
                                  (n2$hetnonew == i) & 
                                  (n2$c_result %in% c(0,1))
                                  ),
                  method = "beta"
                  ), 'ci'
        ))
  
}

for(i in sort(unique(npart_f))){
  
  print(i)
  print(table(n2$c_result[(n2$rsex == 'female') & 
                            (n2$het == 'het. sex. intercourse aged 13+') & 
                            (n2$het1yr %in% 1:901) &
                            (n2$hetnonew == i)]))
  print(100*svyciprop(~(c_result == 1), 
                      design = subset(n2design, 
                                      (n2$rsex == 'female') & 
                                        (n2$het == 'het. sex. intercourse aged 13+') & 
                                        (n2$het1yr %in% 1:901) &
                                        (n2$hetnonew == i) & 
                                        (n2$c_result %in% c(0,1))
                      ),
                      method = "beta"
  )
  )
  print(100*attr(svyciprop(~(c_result == 1), 
                           design = subset(n2design, 
                                           (n2$rsex == 'female') & 
                                             (n2$het == 'het. sex. intercourse aged 13+') & 
                                             (n2$het1yr %in% 1:901) &
                                             (n2$hetnonew == i) & 
                                             (n2$c_result %in% c(0,1))
                           ),
                           method = "beta"
  ), 'ci'
  ))
  
}


####################################
# data for Stan
####################################

dt <- list(
  
  # numbers of partners reported
  N_n_m = length(npart_m),
  n_m = npart_m,
  N_n_f = length(npart_f),
  n_f = npart_f,
  
  # survey weights
  wt_m = n2$urine_wt[(n2$rsex == 'male') & (n2$het == 'het. sex. intercourse aged 13+') & (n2$het1yr %in% 1:901) & (n2$hetnonew != 999)],
  wt_f = n2$urine_wt[(n2$rsex == 'female') & (n2$het == 'het. sex. intercourse aged 13+') & (n2$het1yr %in% 1:901) & (n2$hetnonew != 999)],
  
  # infection status 
  infect_m = n2$c_result[(n2$rsex == 'male') & (n2$het == 'het. sex. intercourse aged 13+') & (n2$het1yr %in% 1:901) & (n2$hetnonew != 999)],
  infect_f = n2$c_result[(n2$rsex == 'female') & (n2$het == 'het. sex. intercourse aged 13+') & (n2$het1yr %in% 1:901) & (n2$hetnonew != 999)],
  
  # clearance rate data - women
  studnum_f = 9, # number of studies 
  studnum_bytype_f = c(4, 3, 2), # number of studies 
  studobs_f = c(4, 5, 1, 1, 3, 1, 1, 5, 4), # number of observations (time periods) in each study
  cumobs_f = c(0,  4,  9, 10, 11, 14, 15, 16, 21), # cumulative number of observations at the start of each study
  Nobs_f = 25, # total number of observations (=sum(studobs))
  r_f = c(10, 7, 6, 6, 2, 7, 1, 0, 3, 23, 3, 17, 0, 0, 8, 3, 2, 2, 4, 0, 2, 44, 23, 7, 2), # number who cleared infection at each time point 
  n_test_f = c(23, 14, 14, 8, 12, 28, 4, 8, 6, 129, 15, 93, 1, 1, 13, 7, 20, 5, 15, 1, 13, 82, 37, 14, 6), # number tested at each time point 
  t_f = c(0.038, 0.058, 0.077, 0.125, 0.012, 0.030, 0.049, 0.088, 0.274, 0.045, 0.083, 0.250, 0.500, 0.750, 1.000, 1.375, 0.083, 0.500, 0.417, 0.917, 0.500, 1.000, 1.000, 1.000, 1.000), # estimated mean follow-up 
  seind_f = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), # did the study use culture as opposed to NAATs?
  T_f = c(999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 999.000, 0.000, 0.000, 0.083, 0.083, 0.500, 0.000, 1.000, 2.000, 3.000), # already followed for...
  
  # clearance rate data - men
  studnum_m = 8, # number of studies 
  studnum_bytype_m = c(6, 2, 0), # number of studies 
  studobs_m = c(1, 1, 4, 1, 5, 1, 4, 1), # number of observations (time periods) in each study
  cumobs_m = c(0, 1,  2,  6,  7, 12, 13, 17), # cumulative number of observations at the start of each study
  Nobs_m = 18, # total number of observations (=sum(studobs))
  r_m = c(0,  4,  3, 13,  3,  2,  7,  3,  2,  1,  0,  1,  5,  1,  0,  0,  0,  1), # number who cleared infection at each time point 
  n_test_m = c(10, 13, 17, 27,  6,  2, 21, 15,  9,  4,  4,  4, 14,  5,  2,  2,  1,  9), # number tested at each time point 
  t_m = c(0.019, 0.023, 0.019, 0.038, 0.058, 0.077, 0.077, 0.012, 0.030, 0.049, 0.088, 0.190, 0.045, 0.019, 0.038, 0.058, 0.077, 0.500), # estimated mean follow-up 
  seind_m = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0), # did the study use culture as opposed to NAATs?
  T_m = c(999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999), # already followed for...
  
  # population (men, women)
  population = c(sum(c(1849600, 1756500, 4292600, 4298200)),  sum(c(1782900, 1750200, 4353400, 4380000))) 
)

####################################
# run the MCMC sampling
####################################

##############
# set initial states for chains
##############

# function to help check the initial conditions are allowed
pi_int <- function(t, y, parms){ # parms = c(alpha, beta, A, lambda){
  dydt = exp(-parms[2]*t) * t^parms[1] / (t + (parms[4]/parms[3] ))
  return(list(dydt))
}

inits <- list(init1 = list(), init2 = list(), init3 = list(), init4 = list())
inits[[1]] <- list(mu_nb = 1, beta = c(1,1), p1 = c(0.23, 0.23), lambda_slow = c(0.74, 0.74), psi = 7/77, A = matrix(0.0001,ncol=11,nrow=2))

for(i in 2:4){ 
  
  inits[[i]]$mu_nb <- exp(log(inits$init1$mu_nb) + rnorm(1, 0, 1)) 
  inits[[i]]$beta <- exp(log(inits$init1$beta) + rnorm(2, 0, 1))
  inits[[i]]$p1 <- inv.logit(logit(inits$init1$p1) + rnorm(2, 0, 1))
  inits[[i]]$lambda_slow <- exp(log(inits$init1$lambda_slow) + rnorm(2, 0, 1))
  inits[[i]]$psi <- inv.logit(logit(inits$init1$psi) + rnorm(1, 0, 1))
  inits[[i]]$A <- inv.logit(logit(inits$init1$A) + matrix(rnorm(22, 0, 1), nrow=2))
  
}

fit_all <- stan(
  file = "transmission_probability_reparam_sens_analysis.stan",
  data = dt,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = inits,
  seed = 67890, refresh=10,
  control=list(adapt_delta=0.8)
)


op_1 <- extract(fit_all)

# to compare A = ppp*rho between the different numbers of partners
qs_m <- apply(op_1$A[,1,], 2, quantile, p=c(0.5, 0.025, 0.975))
plot(0:10, qs_m[1,], ylim=c(0,1), xlab="Number of new partners", ylab ="A (men)", main = "Natsal-2, men")
arrows(0:10, y0=qs_m[2,], y1=qs_m[3,], code=3, angle=90, length=0.05)
for(i in 0:10){
  h <- hist(op_1$A[,1,i+1], n=20, plot=FALSE)
  polygon(c(i, rep(i + 0.05*h$density, each=2), i), rep(h$breaks, each=2), border=NA, col=rgb(0,1,0,0.5))
  polygon(c(i, rep(i - 0.05*h$density, each=2), i), rep(h$breaks, each=2), border=NA, col=rgb(0,1,0,0.5))
}

qs_f <- apply(op_1$A[,2,], 2, quantile, p=c(0.5, 0.025, 0.975))
plot(0:10, qs_f[1,], ylim=c(0,1), xlab = "Number of new partners", ylab="A (women)", main = "Natsal-2, women")
arrows(0:10, y0=qs_f[2,], y1=qs_f[3,], code=3, angle=90, length=0.05)
for(i in 0:10){
  h <- hist(op_1$A[,2,i+1], n=20, plot=FALSE)
  polygon(c(i, rep(i + 0.015*h$density, each=2), i), rep(h$breaks, each=2), border=NA, col=rgb(0,1,0,0.5))
  polygon(c(i, rep(i - 0.015*h$density, each=2), i), rep(h$breaks, each=2), border=NA, col=rgb(0,1,0,0.5))
}


