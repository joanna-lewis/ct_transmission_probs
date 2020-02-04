rm(list=ls())

# setwd('/Users/Joanna/OneDrive - Imperial College London/backup/papers/transmission_probs/ct_transmission_prob')

library(foreign)
library(rstan)
library(survey)
library(deSolve)
library(boot)

set.seed(12345)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

options(survey.lonely.psu="adjust")

# Read in Natsal-2 data

# Edit this line to point to your Natsal-2 dataset
#n2 <- read.dta13("C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\M_Genitalium\\Natsal-2 data\\stata11\\natsal_2000_for_archive.dta")
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

mean(n2$hetnonew[(n2$rsex == 'male') & 
                   (n2$het == 'het. sex. intercourse aged 13+') & 
                   (n2$het1yr %in% 1:901) &
                   (n2$c_result %in% c(0,1))]
     )
svymean(~hetnonew, design = subset(n2design, 
                                   (n2$rsex == 'male') & 
                                     (n2$het == 'het. sex. intercourse aged 13+') & 
                                     (n2$het1yr %in% 1:901) &
                                     (n2$c_result %in% c(0,1)))
)

tapply(n2$hetnonew[(n2$rsex == 'male') &
                     (n2$het == 'het. sex. intercourse aged 13+') & 
                     (n2$het1yr %in% 1:901) &
                     (n2$c_result %in% c(0,1))],
       n2$c_result[(n2$rsex == 'male') &
                 (n2$het == 'het. sex. intercourse aged 13+') & 
                 (n2$het1yr %in% 1:901) &
                 (n2$c_result %in% c(0,1))],
       mean
)

svyby(~hetnonew, by = ~c_result, design = subset(n2design, 
                                   (n2$rsex == 'male') & 
                                     (n2$het == 'het. sex. intercourse aged 13+') & 
                                     (n2$het1yr %in% 1:901) & 
                                     (n2$c_result %in% c(0,1))
                                     ),
      svymean)

print(100*svyciprop(~(c_result == 1), 
                    design = subset(n2design, 
                                    (n2$rsex == 'male') & 
                                      (n2$het == 'het. sex. intercourse aged 13+') & 
                                      (n2$het1yr %in% 1:901) &
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
                                           (n2$c_result %in% c(0,1))
                         ),
                         method = "beta"
), 'ci'
))



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

mean(n2$hetnonew[(n2$rsex == 'female') & 
                   (n2$het == 'het. sex. intercourse aged 13+') & 
                   (n2$het1yr %in% 1:901) &
                   (n2$c_result %in% c(0,1))]
)
svymean(~hetnonew, design = subset(n2design, 
                                   (n2$rsex == 'female') & 
                                     (n2$het == 'het. sex. intercourse aged 13+') & 
                                     (n2$het1yr %in% 1:901) &
                                     (n2$c_result %in% c(0,1)))
)

tapply(n2$hetnonew[(n2$rsex == 'female') &
                     (n2$het == 'het. sex. intercourse aged 13+') & 
                     (n2$het1yr %in% 1:901) &
                     (n2$c_result %in% c(0,1))],
       n2$c_result[(n2$rsex == 'female') &
                     (n2$het == 'het. sex. intercourse aged 13+') & 
                     (n2$het1yr %in% 1:901) &
                     (n2$c_result %in% c(0,1))],
       mean
)
svyby(~hetnonew, by = ~c_result, design = subset(n2design, 
                                                 (n2$rsex == 'female') & 
                                                   (n2$het == 'het. sex. intercourse aged 13+') & 
                                                   (n2$het1yr %in% 1:901) & 
                                                   (n2$c_result %in% c(0,1))
),
svymean)

####################################
# overall prevalence
####################################

print(100*svyciprop(~(c_result == 1), 
                      design = subset(n2design,
                                      (n2$rsex == 'male') & 
                                      (n2$het == 'het. sex. intercourse aged 13+') & 
                                      (n2$het1yr %in% 1:901) &
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
                                      (n2$c_result %in% c(0,1))
                    ),
                    method = "beta"
), 'ci')
)

print(100*svyciprop(~(c_result == 1), 
                    design = subset(n2design,
                                    (n2$rsex == 'female') & 
                                      (n2$het == 'het. sex. intercourse aged 13+') & 
                                      (n2$het1yr %in% 1:901) &
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
                                           (n2$c_result %in% c(0,1))
                         ),
                         method = "beta"
), 'ci')
)
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

inits <- list(init1 = list(), init2 = list(), init3 = list(), init4 = list(), init5 = list())
inits[[1]] <- list(mu_nb = c(1,1), beta = c(1,1), p1 = c(0.23, 0.23), lambda_slow = c(0.74, 0.74), psi = 7/77, A = c(0.0001,0.0001))

for(i in 2:4){ 
  
    inits[[i]]$mu_nb <- exp(log(inits$init1$mu_nb) + rnorm(2, 0, 1)) 
    inits[[i]]$beta <- exp(log(inits$init1$beta) + rnorm(2, 0, 1))
    inits[[i]]$p1 <- inv.logit(logit(inits$init1$p1) + rnorm(2, 0, 1))
    inits[[i]]$lambda_slow <- exp(log(inits$init1$lambda_slow) + rnorm(2, 0, 1))
    inits[[i]]$psi <- inv.logit(logit(inits$init1$psi) + rnorm(1, 0, 1))
    inits[[i]]$A <- inv.logit(logit(inits$init1$A) + rnorm(2, 0, 1))
  
}

fit_all <- stan(
  file = "transmission_probability_reparam_different_means.stan",
  data = dt,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = inits,
  seed = 67890, refresh=10,
  control=list(adapt_delta=0.8)
)

#save(fit_all, file= 'transmission_probability_reparam_different_means_natsal2_hetnonew_1000.RData')

check_divergences(fit_all)
check_treedepth(fit_all) 
check_energy(fit_all)

sum_all <- summary(fit_all)
hist(sum_all$summary[,'n_eff']/80000)
table(sum_all$summary[,'n_eff']/80000 < 0.001)

hist(sum_all$summary[,'Rhat'])

op_0 <- extract(fit_all)
mli <- which(op_0$lp__ == max(op_0$lp__))

###################
# plots 
###################

plot(op_0$lp__, type='l')

#quartz()
par(mfrow=c(1,1))
plot(rep(h1$breaks, each=2), c(0, rep(h1$density, each=2), 0), type='l', col='#F98400', lwd = 3, xlim=c(0,1), ylim=c(0,20), xlab="Transmission probability per partnership", ylab='Density', main="Natsal-2")
lines(rep(h2$breaks, each=2), c(0, rep(h2$density, each=2), 0), col='#00A08A', lwd=3 )
legend('topright', lwd=3, col=c('#F98400', '#00A08A'), legend=c('Male-to-female', 'Female-to-male'), bty='n')
text(0, 20, 'A', cex=5, adj = c(0,1))

quantile(op_0$rho_mf, p=c(0.5, 0.025, 0.975))
quantile(op_0$rho_fm, p=c(0.5, 0.025, 0.975))

############
# number of partners 
############

# tabulate numbers of partners simulated for men and women
simtab_m <- table((matrix(rep( 1:nrow(op_0$n_m_sim), times=ncol(op_0$n_m_sim) ), nrow=nrow(op_0$n_m_sim))), op_0$n_m_sim)
simtab_f <- table((matrix(rep( 1:nrow(op_0$n_f_sim), times=ncol(op_0$n_f_sim) ), nrow=nrow(op_0$n_f_sim))), op_0$n_f_sim)

quartz(width=8, height=4.5)
par(mfrow=c(1,2))

plot(0,0,xlim=c(0,15), ylim=c(0.0001,1), pch='', xlab = "Number of partners", ylab = "Proportion of respondents", main='Natsal-2, Men', log='y', yaxt="n")
for(i in 1:length(simtab_m[,1]) )
	points(as.numeric(dimnames(simtab_m)[[2]]), simtab_m[i,]/sum(simtab_m[i,]), 
		pch=16, col=rgb(0,0,0,0.01), cex=0.5
		) 
lines(as.numeric(dimnames(simtab_m)[[2]]), apply(simtab_m/dt$N_n_m, 2, quantile, p=0.5))
lines(as.numeric(dimnames(simtab_m)[[2]]), apply(simtab_m/dt$N_n_m, 2, quantile, p=0.025), lty=2)
lines(as.numeric(dimnames(simtab_m)[[2]]), apply(simtab_m/dt$N_n_m, 2, quantile, p=0.975), lty=2)
axis(2, at=c(0.0001, 0.001, 0.01, 0.1, 1), label=c("0.0001", "0.001", "0.01", "0.1", "1"))

points(as.numeric(as.character(ntab_m$n)), ntab_m$prop_wt, pch='+', cex=1, col='red')

plot(0,0,xlim=c(0,15), ylim=c(0.0001,1), pch='', xlab = "Number of partners", ylab = "Proportion of respondents", main='Natsal-2, Women', log='y', yaxt="n")
for(i in 1:length(simtab_f[,1]) )
	points(as.numeric(dimnames(simtab_f)[[2]]), simtab_f[i,]/sum(simtab_f[i,]), 
		pch=16, col=rgb(0,0,0,0.01), cex=0.5
		) 
lines(as.numeric(dimnames(simtab_f)[[2]]), apply(simtab_f/dt$N_n_f, 2, quantile, p=0.5))
lines(as.numeric(dimnames(simtab_f)[[2]]), apply(simtab_f/dt$N_n_f, 2, quantile, p=0.025), lty=2)
lines(as.numeric(dimnames(simtab_f)[[2]]), apply(simtab_f/dt$N_n_f, 2, quantile, p=0.975), lty=2)
axis(2, at=c(0.0001, 0.001, 0.01, 0.1, 1), label=c("0.0001", "0.001", "0.01", "0.1", "1"))

points(as.numeric(as.character(ntab_f$n)), ntab_f$prop_wt, pch='+', cex=1, col='red')

#########################
# infection status with bars
#########################

# men
layout(matrix(1:4, nrow=2, byrow=TRUE), widths=c(1,1), heights=c(6,3))

par(mar=c(1,4,4,2))
plot(1,1, pch='', xlim=c(0,25), ylim=c(0.001,1), ylab="Prevalence", main="Natsal-2, Men", xaxt='n')
infect_table_m <- matrix(nrow=dim(op_0$sim_infect_m)[1], ncol=length(unique(dt$n_m)))

for(i in 1:nrow(infect_table_m)){
	infect_table_m[i,] <- table(dt$n_m, op_0$sim_infect_m[i,])[,2] / table(dt$n_m)
	points( as.numeric(names(table(dt$n_m))), infect_table_m[i,], pch=16, col=rgb(0,0,0,0.1) )
	}

lines(as.numeric(names(table(dt$n_m))), apply(infect_table_m, 2, quantile, p=0.5) )
lines(as.numeric(names(table(dt$n_m))), apply(infect_table_m, 2, quantile, p=0.025), lty=2 )
lines(as.numeric(names(table(dt$n_m))), apply(infect_table_m, 2, quantile, p=0.975), lty=2 )

points(as.numeric(names(table(dt$n_m))), tapply(dt$wt_m*dt$infect_m, dt$n_m, sum) / tapply(dt$wt_m, dt$n_m, sum), pch='+', cex=2, col='red')

plot(1,1, pch='', xlim=c(0,25), ylim=c(0.001,1), ylab="Prevalence", main="Natsal-2, Women", xaxt='n')
infect_table_f <- matrix(nrow=dim(op_0$sim_infect_f)[1], ncol=length(unique(dt$n_f)))

for(i in 1:nrow(infect_table_f)){
	infect_table_f[i,] <- table(dt$n_f, op_0$sim_infect_f[i,])[,2] / table(dt$n_f)
	points( as.numeric(names(table(dt$n_f))), infect_table_f[i,], pch=16, col=rgb(0,0,0,0.1) )
	}

lines(as.numeric(names(table(dt$n_f))), apply(infect_table_f, 2, quantile, p=0.5) )
lines(as.numeric(names(table(dt$n_f))), apply(infect_table_f, 2, quantile, p=0.025), lty=2 )
lines(as.numeric(names(table(dt$n_f))), apply(infect_table_f, 2, quantile, p=0.975), lty=2 )

points(as.numeric(names(table(dt$n_f))), tapply(dt$wt_f*dt$infect_f, dt$n_f, sum) / tapply(dt$wt_f, dt$n_f, sum), pch='+', cex=2, col='red')

par(mar=c(5,4,0,2))
plot(0, 0, xlim=c(0, 25), ylim=c(0, 1400), pch='', xlab='Number of partners', ylab='Frequency', yaxt='n')
for(i in 1:nrow(ntab_m)){
	rect(as.numeric(as.character(ntab_m$n[i]))-0.45, 0, as.numeric(as.character(ntab_m$n[i]))+0.45, ntab_m$freq[i])
	text(as.numeric(as.character(ntab_m$n[i])), 1300, as.character(ntab_m$freq[i]), adj=1, srt=90, cex=0.8)
	}
axis(2, at=c(0, 500, 1000))

plot(0, 0, xlim=c(0, 25), ylim=c(0, 1400), pch='', xlab='Number of partners', ylab='Frequency', yaxt='n')
for(i in 1:nrow(ntab_f)){
	rect(as.numeric(as.character(ntab_f$n[i]))-0.45, 0, as.numeric(as.character(ntab_f$n[i]))+0.45, ntab_f$freq[i])
	text(as.numeric(as.character(ntab_f$n[i])), 1300, as.character(ntab_f$freq[i]), adj=1, srt=90, cex=0.8)
	}
axis(2, at=c(0, 500, 1000), labels=c(0, 500, 1000))

#########################
# some posterior summaries
#########################

mean(op_0$rho_mf)
quantile(op_0$rho_mf, p=c(0.5, 0.025, 0.975))

mean(op_0$rho_fm)
quantile(op_0$rho_fm, p=c(0.5, 0.025, 0.975))

