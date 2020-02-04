rm(list=ls())

# setwd('/Users/Joanna/OneDrive - Imperial College London/backup/papers/transmission_probs/ct_transmission_prob')

library(foreign)
library(rstan)
library(survey)
library(deSolve)
library(boot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

alldata <- NULL # set up data frame

# NHANES 09-10 data - Windows address C:\Users\ucbpjle\OneDrive - Imperial College London\backup\data\NHANES

nhanes_09_10_dem <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2009-2010/DEMO_F.XPT')
#nhanes_09_10_dem <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2009-2010\\DEMO_F.XPT')
nhanes_09_10_ctgc <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2009-2010/CHLMDA_F.XPT')
#nhanes_09_10_ctgc <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2009-2010\\CHLMDA_F.XPT')
nhanes_09_10_beh <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2009-2010/SXQ_F.XPT')
#nhanes_09_10_beh <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2009-2010\\SXQ_F.XPT')

nhanes_09_10 <- merge(nhanes_09_10_dem, nhanes_09_10_ctgc, all=TRUE)
nhanes_09_10 <- merge(nhanes_09_10, nhanes_09_10_beh, all=TRUE)

data_09_10 <- data.frame(
  survey = "2009-2010",
  id = nhanes_09_10$SEQN,
  sex = nhanes_09_10$RIAGENDR,
  partners = ifelse(nhanes_09_10$RIAGENDR == 1, nhanes_09_10$SXQ827, nhanes_09_10$SXQ727),
  shift_partners = ifelse(nhanes_09_10$RIAGENDR == 1, nhanes_09_10$SXQ827, nhanes_09_10$SXQ727),
  ct = nhanes_09_10$URXUCL,
  weight = nhanes_09_10$WTMEC2YR,
  psu = nhanes_09_10$SDMVPSU,
  str = nhanes_09_10$SDMVSTRA
)

data_09_10$shift_partners[which(nhanes_09_10$SXQ648 == 2)] <- 0 # if state no new partners

# >= 1 new partner, and one partner total => one new partner
data_09_10$shift_partners[which((nhanes_09_10$RIAGENDR == 1) & (nhanes_09_10$SXQ648 == 1) & (nhanes_09_10$SXQ827 == 1))] <- 1
data_09_10$shift_partners[which((nhanes_09_10$RIAGENDR == 2) & (nhanes_09_10$SXQ648 == 1) & (nhanes_09_10$SXQ727 == 1))] <- 1

# >= 1 new partner, and >1 partner total => assume new partners = partners-1
data_09_10$shift_partners[which((nhanes_09_10$RIAGENDR == 1) & (nhanes_09_10$SXQ648 == 1) & (nhanes_09_10$SXQ827 > 1))] <- 
  nhanes_09_10$SXQ827[which((nhanes_09_10$RIAGENDR == 1) & (nhanes_09_10$SXQ648 == 1) & (nhanes_09_10$SXQ827 > 1))] - 1
data_09_10$shift_partners[which((nhanes_09_10$RIAGENDR == 2) & (nhanes_09_10$SXQ648 == 1) & (nhanes_09_10$SXQ727 > 1))] <- 
  nhanes_09_10$SXQ727[which((nhanes_09_10$RIAGENDR == 2) & (nhanes_09_10$SXQ648 == 1) & (nhanes_09_10$SXQ727 > 1))] - 1


data_09_10 <- data_09_10[which(complete.cases(data_09_10)),]
data_09_10 <- data_09_10[data_09_10$partners %in% (1:999),]
data_09_10$weight <- nrow(data_09_10) * data_09_10$weight / sum(data_09_10$weight)
alldata <- rbind(alldata, data_09_10)

# NHANES 11-12 data

nhanes_11_12_dem <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2011-2012/DEMO_G.XPT')
nhanes_11_12_ctgc <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2011-2012/CHLMDA_G.XPT')
nhanes_11_12_beh <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2011-2012/SXQ_G.XPT')

#nhanes_11_12_dem <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2011-2012\\DEMO_G.XPT')
#nhanes_11_12_ctgc <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2011-2012\\CHLMDA_G.XPT')
#nhanes_11_12_beh <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2011-2012\\SXQ_G.XPT')

nhanes_11_12 <- merge(nhanes_11_12_dem, nhanes_11_12_ctgc, all=TRUE)
nhanes_11_12 <- merge(nhanes_11_12, nhanes_11_12_beh, all=TRUE)

data_11_12 <- data.frame(
  survey = "2011-2012",
  id = nhanes_11_12$SEQN,
  sex = nhanes_11_12$RIAGENDR,
  partners = ifelse(nhanes_11_12$RIAGENDR == 1, nhanes_11_12$SXQ827, nhanes_11_12$SXQ727),
  shift_partners = NA, # ifelse(nhanes_11_12$RIAGENDR == 1, nhanes_11_12$SXQ827 - 1, nhanes_11_12$SXQ727 - 1),
  ct = nhanes_11_12$URXUCL,
  weight = nhanes_11_12$WTMEC2YR,
  psu = nhanes_11_12$SDMVPSU,
  str = nhanes_11_12$SDMVSTRA
)

data_11_12$shift_partners[which(nhanes_11_12$SXQ648 == 2)] <- 0 # if state no new partners

# >= 1 new partner, and one partner total => one new partner
data_11_12$shift_partners[which((nhanes_11_12$RIAGENDR == 1) & (nhanes_11_12$SXQ648 == 1) & (nhanes_11_12$SXQ827 == 1))] <- 1
data_11_12$shift_partners[which((nhanes_11_12$RIAGENDR == 2) & (nhanes_11_12$SXQ648 == 1) & (nhanes_11_12$SXQ727 == 1))] <- 1

# >= 1 new partner, and >1 partner total => assume new partners = partners-1
data_11_12$shift_partners[which((nhanes_11_12$RIAGENDR == 1) & (nhanes_11_12$SXQ648 == 1) & (nhanes_11_12$SXQ827 > 1))] <- 
  nhanes_11_12$SXQ827[which((nhanes_11_12$RIAGENDR == 1) & (nhanes_11_12$SXQ648 == 1) & (nhanes_11_12$SXQ827 > 1))] - 1
data_11_12$shift_partners[which((nhanes_11_12$RIAGENDR == 2) & (nhanes_11_12$SXQ648 == 1) & (nhanes_11_12$SXQ727 > 1))] <- 
  nhanes_11_12$SXQ727[which((nhanes_11_12$RIAGENDR == 2) & (nhanes_11_12$SXQ648 == 1) & (nhanes_11_12$SXQ727 > 1))] - 1

data_11_12 <- data_11_12[which(complete.cases(data_11_12)),]
data_11_12 <- data_11_12[data_11_12$partners %in% (1:999),]
data_11_12$weight <- nrow(data_11_12) * data_11_12$weight / sum(data_11_12$weight)
alldata <- rbind(alldata, data_11_12)

# NHANES 13-14 data

nhanes_13_14_dem <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2013-2014/DEMO_H.XPT')
nhanes_13_14_ctgc <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2013-2014/CHLMDA_H.XPT')
nhanes_13_14_beh <- read.xport('/Users/Joanna/OneDrive - Imperial College London/backup/data/NHANES/data/2013-2014/SXQ_H.XPT')

#nhanes_13_14_dem <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2013-2014\\DEMO_H.XPT')
#nhanes_13_14_ctgc <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2013-2014\\CHLMDA_H.XPT')
#nhanes_13_14_beh <- read.xport('C:\\Users\\ucbpjle\\OneDrive - Imperial College London\\backup\\data\\NHANES\\data\\2013-2014\\SXQ_H.XPT')

nhanes_13_14 <- merge(nhanes_13_14_dem, nhanes_13_14_ctgc, all=TRUE)
nhanes_13_14 <- merge(nhanes_13_14, nhanes_13_14_beh, all=TRUE)

data_13_14 <- data.frame(
  survey = "2013-2014",
  id = nhanes_13_14$SEQN,
  sex = nhanes_13_14$RIAGENDR,
  partners = ifelse(nhanes_13_14$RIAGENDR == 1, nhanes_13_14$SXQ827, nhanes_13_14$SXQ727),
  shift_partners = NA, #ifelse(nhanes_13_14$RIAGENDR == 1, nhanes_13_14$SXQ827 - 1, nhanes_13_14$SXQ727 - 1),
  ct = nhanes_13_14$URXUCL,
  weight = nhanes_13_14$WTMEC2YR,
  psu = nhanes_13_14$SDMVPSU,
  str = nhanes_13_14$SDMVSTRA
)

data_13_14$shift_partners[which(nhanes_13_14$SXQ648 == 2)] <- 0 # if state no new partners

# >= 1 new partner, and one partner total => one new partner
data_13_14$shift_partners[which((nhanes_13_14$RIAGENDR == 1) & (nhanes_13_14$SXQ648 == 1) & (nhanes_13_14$SXQ827 == 1))] <- 1
data_13_14$shift_partners[which((nhanes_13_14$RIAGENDR == 2) & (nhanes_13_14$SXQ648 == 1) & (nhanes_13_14$SXQ727 == 1))] <- 1

# >= 1 new partner, and >1 partner total => assume new partners = partners-1
data_13_14$shift_partners[which((nhanes_13_14$RIAGENDR == 1) & (nhanes_13_14$SXQ648 == 1) & (nhanes_13_14$SXQ827 > 1))] <- 
  nhanes_13_14$SXQ827[which((nhanes_13_14$RIAGENDR == 1) & (nhanes_13_14$SXQ648 == 1) & (nhanes_13_14$SXQ827 > 1))] - 1
data_13_14$shift_partners[which((nhanes_13_14$RIAGENDR == 2) & (nhanes_13_14$SXQ648 == 1) & (nhanes_13_14$SXQ727 > 1))] <- 
  nhanes_13_14$SXQ727[which((nhanes_13_14$RIAGENDR == 2) & (nhanes_13_14$SXQ648 == 1) & (nhanes_13_14$SXQ727 > 1))] - 1


data_13_14 <- data_13_14[which(complete.cases(data_13_14)),]
data_13_14 <- data_13_14[data_13_14$partners %in% (1:999),]
data_13_14$weight <- nrow(data_13_14) * data_13_14$weight / sum(data_13_14$weight)
alldata <- rbind(alldata, data_13_14)

head(alldata)

####################################
# number of new partners in the last year
####################################

######
# men
######

npart_m <- alldata$shift_partners[(alldata$sex == 1)]
ntab_m <- table(npart_m)
ntab_m <- data.frame(ntab_m)
names(ntab_m) <- c('n', 'freq')
ntab_m$prop <- ntab_m$freq/sum(ntab_m$freq)
ntab_m$low.95 <- ntab_m$high.95 <- ntab_m$prop_wt <- NA

sub <- subset(alldata, alldata$sex == 1)
for(i in 1:(nrow(ntab_m))){
  ntab_m$low.95[i] <- binom.test(ntab_m$freq[i], sum(ntab_m$freq))$conf.int[[1]]
  ntab_m$high.95[i] <- binom.test(ntab_m$freq[i], sum(ntab_m$freq))$conf.int[[2]]
  ntab_m$prop_wt[i] <- sum(sub$weight[sub$shift_partners == as.numeric(as.character(ntab_m$n[i]))]) / sum(sub$weight)
}

######
# women
######

npart_f <- alldata$shift_partners[(alldata$sex == 2)]
ntab_f <- table(npart_f)
ntab_f <- data.frame(ntab_f)
names(ntab_f) <- c('n', 'freq')
ntab_f$prop <- ntab_f$freq/sum(ntab_f$freq)
ntab_f$low.95 <- ntab_f$high.95 <- ntab_f$prop_wt <- NA

sub <- subset(alldata, alldata$sex == 2)
for(i in 1:(nrow(ntab_f))){
  ntab_f$low.95[i] <- binom.test(ntab_f$freq[i], sum(ntab_f$freq))$conf.int[[1]]
  ntab_f$high.95[i] <- binom.test(ntab_f$freq[i], sum(ntab_f$freq))$conf.int[[2]]
  ntab_f$prop_wt[i] <- sum(sub$weight[sub$shift_partners == as.numeric(as.character(ntab_f$n[i]))]) / sum(sub$weight)
}

nhanesdesign <- svydesign(
  ids = ~ psu,
  weights = ~ weight,
  strata = ~ str,
  data = alldata,
  nest=TRUE
)

####################################
# tabulate partner numbers and infection status
####################################

for(i in sort(unique(npart_m))){
  
  print(i)
  print(table(alldata$ct[(alldata$sex == 1) & 
                         (alldata$shift_partners == i)]))
  print(100*svyciprop(~(ct == 1), 
                      design = subset(nhanesdesign, 
                                      (alldata$sex == 1) & 
                                      (alldata$shift_partners == i) & 
                                      (alldata$ct %in% c(1,2))
                      ),
                      method = "beta"
  ))
  print(100*attr(svyciprop(~(ct == 1), 
                           design = subset(nhanesdesign, 
                                          (alldata$sex == 1) & 
                                          (alldata$shift_partners == i) & 
                                          (alldata$ct %in% c(1,2))
                           ),
                           method = "beta"
  ), 'ci'
  ))
  
}

mean(alldata$shift_partners[(alldata$sex == 1) & 
                   (alldata$ct %in% c(1,2))]
)
svymean(~shift_partners, design = subset(nhanesdesign, 
                                   (alldata$sex == 1) & 
                                   (alldata$ct %in% c(1,2)))
)

tapply(alldata$shift_partners[(alldata$sex == 1) &
                     (alldata$ct %in% c(1,2))],
       alldata$ct[(alldata$sex == 1) &
                     (alldata$ct %in% c(1,2))],
       mean
)

svyby(~shift_partners, by = ~ct, design = subset(nhanesdesign, 
                                                 (alldata$sex == 1) & 
                                                 (alldata$ct %in% c(1,2))
),
svymean)

for(i in sort(unique(npart_f))){
  
  print(i)
  print(table(alldata$ct[(alldata$sex == 2) & 
                           (alldata$shift_partners == i)]))
  print(100*svyciprop(~(ct == 1), 
                      design = subset(nhanesdesign, 
                                      (alldata$sex == 2) & 
                                        (alldata$shift_partners == i) & 
                                        (alldata$ct %in% c(1,2))
                      ),
                      method = "beta"
  ))
  print(100*attr(svyciprop(~(ct == 1), 
                           design = subset(nhanesdesign, 
                                           (alldata$sex == 2) & 
                                             (alldata$shift_partners == i) & 
                                             (alldata$ct %in% c(1,2))
                           ),
                           method = "beta"
  ), 'ci'
  ))
  
}

mean(alldata$shift_partners[(alldata$sex == 2) & 
                              (alldata$ct %in% c(1,2))]
)
svymean(~shift_partners, design = subset(nhanesdesign, 
                                         (alldata$sex == 2) & 
                                           (alldata$ct %in% c(1,2)))
)

tapply(alldata$shift_partners[(alldata$sex == 2) &
                                (alldata$ct %in% c(1,2))],
       alldata$ct[(alldata$sex == 2) &
                    (alldata$ct %in% c(1,2))],
       mean
)

svyby(~shift_partners, by = ~ct, design = subset(nhanesdesign, 
                                                 (alldata$sex == 2) & 
                                                   (alldata$ct %in% c(1,2))
),
svymean)

####################################
# overall prevalence
####################################

print(100*svyciprop(~(ct == 1), 
                    design = subset(nhanesdesign,
                                    (alldata$sex == 1) & 
                                      (alldata$ct %in% c(1,2))
                    ),
                    method = "beta"
)
)

print(100*attr(svyciprop(~(ct == 1), 
                         design = subset(nhanesdesign,
                                         (alldata$sex == 1) & 
                                           (alldata$ct %in% c(1,2))
                         ),
                         method = "beta"
), 'ci')
)

print(100*svyciprop(~(ct == 1), 
                    design = subset(nhanesdesign,
                                    (alldata$sex == 2) & 
                                      (alldata$ct %in% c(1,2))
                    ),
                    method = "beta"
)
)

print(100*attr(svyciprop(~(ct == 1), 
                         design = subset(nhanesdesign,
                                         (alldata$sex == 2) & 
                                           (alldata$ct %in% c(1,2))
                         ),
                         method = "beta"
), 'ci')
)

################
# data for STAN
################

dt <- list(
  
  # numbers of partners reported
  N_n_m = length(npart_m),
  n_m = npart_m,
  N_n_f = length(npart_f),
  n_f = npart_f,
  
  # survey weights
  wt_m = alldata$weight[alldata$sex == 1],
  wt_f = alldata$weight[alldata$sex == 2],
  
  # infection status 
  infect_m = 2 - alldata$ct[alldata$sex == 1], # 1 = infected; 0 = susceptible
  infect_f = 2 - alldata$ct[alldata$sex == 2], # 1 = infected; 0 = susceptible
  
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
  
  diag = c(307340, 878708),
  population = c(52984487, 51823583) # population aged 15-39, 2009
)

####################################
# run the MCMC sampling
####################################

##############
# set initial states for chains
##############

inits <- list(init1 = list(), init2 = list(), init3 = list(), init4 = list(), init5 = list())
inits[[1]] <- list(mu_nb = 1, beta = c(1,1), p1 = c(0.23, 0.23), lambda_slow = c(0.74, 0.74), psi = 7/77, A = c(0.0001,0.0001))

for(i in 2:4){ 
  
  inits[[i]]$mu_nb <- exp(log(inits$init1$mu_nb) + rnorm(1, 0, 1)) 
  inits[[i]]$beta <- exp(log(inits$init1$beta) + rnorm(2, 0, 1))
  inits[[i]]$p1 <- inv.logit(logit(inits$init1$p1) + rnorm(2, 0, 1))
  inits[[i]]$lambda_slow <- exp(log(inits$init1$lambda_slow) + rnorm(2, 0, 1))
  inits[[i]]$psi <- inv.logit(logit(inits$init1$psi) + rnorm(1, 0, 1))
  inits[[i]]$A <- inv.logit(logit(inits$init1$A) + rnorm(2, 0, 1))
  
}

fit_all <- stan(
  file = "transmission_probability_reparam.stan",
  data = dt,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = inits,
  seed = 67890, refresh=10,
  control=list(adapt_delta=0.8)
)

#save(fit_all, file='transmission_probability_reparam_nhanes_hetnonew_1000.RData')

check_divergences(fit_all)
check_treedepth(fit_all)
check_energy(fit_all)

sum_all <- summary(fit_all)
hist(sum_all$summary[,'n_eff']/25000)
table(sum_all$summary[,'n_eff']/25000 < 0.001)

hist(sum_all$summary[,'Rhat'])

op_0 <- extract(fit_all)
mli <- which(op_0$lp__ == max(op_0$lp__))

###################
# plots 
###################

plot(op_0$lp__, type='l')

h1 <- hist(op_0$rho_mf, breaks = c(seq(0,1,0.05), max(op_0$rho_mf)), plot=FALSE)
h2 <- hist(op_0$rho_fm, breaks = c(seq(0,1,0.05), max(op_0$rho_fm)), plot=FALSE)
#quartz()
par(mfrow=c(1,1))
plot(rep(h1$breaks, each=2), c(0, rep(h1$density, each=2), 0), 
     type='l', col='#F98400', lwd = 3, 
     xlim=c(0,1), ylim=c(0,12), xaxt = "n",
     xlab="Transmission probability per partnership (%)", ylab='Density', main="NHANES opposite-sex partnerships"
)
lines(rep(h2$breaks, each=2), c(0, rep(h2$density, each=2), 0), col='#00A08A', lwd=3 )
axis(1, at = seq(0,1,0.2), labels=seq(0,100,20))
legend('topright', inset=c(0,0.2), lwd=3, col=c('#F98400', '#00A08A'), legend=c('Male-to-female', 'Female-to-male'), bty='n')
text(1, 12, 'B', cex=5, adj = c(1,1))

############
# number of partners 
############

# tabulate numbers of partners simulated for men and women
simtab_m <- table((matrix(rep( 1:nrow(op_0$n_m_sim), times=ncol(op_0$n_m_sim) ), nrow=nrow(op_0$n_m_sim))), op_0$n_m_sim)
simtab_f <- table((matrix(rep( 1:nrow(op_0$n_f_sim), times=ncol(op_0$n_f_sim) ), nrow=nrow(op_0$n_f_sim))), op_0$n_f_sim)

par(mfrow=c(1,2))

plot(0,0,xlim=c(0,15), ylim=c(0.0001,1), pch='', xlab = "Number of partners", ylab = "Proportion of respondents", main='NHANES, Men', log='y')
for(i in 1:length(simtab_m[,1]) )
	points(as.numeric(dimnames(simtab_m)[[2]]), simtab_m[i,]/sum(simtab_m[i,]), 
		pch=16, col=rgb(0,0,0,0.01), cex=0.5
		) 
lines(as.numeric(dimnames(simtab_m)[[2]]), apply(simtab_m/dt$N_n_m, 2, quantile, p=0.5))
lines(as.numeric(dimnames(simtab_m)[[2]]), apply(simtab_m/dt$N_n_m, 2, quantile, p=0.025), lty=2)
lines(as.numeric(dimnames(simtab_m)[[2]]), apply(simtab_m/dt$N_n_m, 2, quantile, p=0.975), lty=2)

points(as.numeric(as.character(ntab_m$n)), ntab_m$prop_wt, pch='+', cex=1, col='red')

plot(0,0,xlim=c(0,15), ylim=c(0.0001,1), pch='', xlab = "Number of partners", ylab = "Proportion of respondents", main='NHANES, Women', log='y')
for(i in 1:length(simtab_f[,1]) )
	points(as.numeric(dimnames(simtab_f)[[2]]), simtab_f[i,]/sum(simtab_f[i,]), 
		pch=16, col=rgb(0,0,0,0.01), cex=0.5
		) 
lines(as.numeric(dimnames(simtab_f)[[2]]), apply(simtab_f/dt$N_n_f, 2, quantile, p=0.5))
lines(as.numeric(dimnames(simtab_f)[[2]]), apply(simtab_f/dt$N_n_f, 2, quantile, p=0.025), lty=2)
lines(as.numeric(dimnames(simtab_f)[[2]]), apply(simtab_f/dt$N_n_f, 2, quantile, p=0.975), lty=2)

points(as.numeric(as.character(ntab_f$n)), ntab_f$prop_wt, pch='+', cex=1, col='red')

#########################
# infection status with bars
#########################

# quartz(width=8.27, height=5.83)

# men
layout(matrix(1:4, nrow=2, byrow=TRUE), widths=c(1,1), heights=c(6,3))

par(mar=c(1,4,4,2))
plot(1,1, pch='', xlim=c(0,25), ylim=c(0.001,1), ylab="Prevalence", main="NHANES, Men", xaxt='n')
infect_table_m <- matrix(nrow=dim(op_0$sim_infect_m)[1], ncol=length(unique(dt$n_m)))

for(i in 1:nrow(infect_table_m)){
  infect_table_m[i,] <- table(dt$n_m, op_0$sim_infect_m[i,])[,2] / table(dt$n_m)
  points( as.numeric(names(table(dt$n_m))), infect_table_m[i,], pch=16, col=rgb(0,0,0,0.1) )
}

lines(as.numeric(names(table(dt$n_m))), apply(infect_table_m, 2, quantile, p=0.5) )
lines(as.numeric(names(table(dt$n_m))), apply(infect_table_m, 2, quantile, p=0.025), lty=2 )
lines(as.numeric(names(table(dt$n_m))), apply(infect_table_m, 2, quantile, p=0.975), lty=2 )

points(as.numeric(names(table(dt$n_m))), tapply(dt$wt_m*dt$infect_m, dt$n_m, sum) / tapply(dt$wt_m, dt$n_m, sum), pch='+', cex=2, col='red')

plot(1,1, pch='', xlim=c(0,25), ylim=c(0.001,1), ylab="Prevalence", main="NHANES, Women", xaxt='n')
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
plot(0, 0, xlim=c(0, 25), ylim=c(0, 2000), pch='', xlab='Number of partners', ylab='Frequency', yaxt='n')
for(i in 1:nrow(ntab_m)){
  rect(as.numeric(as.character(ntab_m$n[i]))-0.45, 0, as.numeric(as.character(ntab_m$n[i]))+0.45, ntab_m$freq[i])
  text(as.numeric(as.character(ntab_m$n[i])), 1500, as.character(ntab_m$freq[i]), adj=1, srt=90, cex=0.8)
}
axis(2, at=c(0, 500, 1000, 1500, 2000), labels=c(0, "", 1000, "", ""))
mtext("1000", side=2, line=1, at=1000, cex=0.8)
mtext("2000", side=2, line=1, at=2000, cex=0.8)

plot(0, 0, xlim=c(0, 25), ylim=c(0, 2000), pch='', xlab='Number of partners', ylab='Frequency', yaxt='n')
for(i in 1:nrow(ntab_f)){
  rect(as.numeric(as.character(ntab_f$n[i]))-0.45, 0, as.numeric(as.character(ntab_f$n[i]))+0.45, ntab_f$freq[i])
  text(as.numeric(as.character(ntab_f$n[i])), 1500, as.character(ntab_f$freq[i]), adj=1, srt=90, cex=0.8)
}
axis(2, at=c(0, 500, 1000, 1500, 2000), labels=c(0, "", 1000, "", ""))
mtext("1000", side=2, line=1, at=1000, cex=0.8)
mtext("2000", side=2, line=1, at=2000, cex=0.8)

#########################
# some posterior summaries
#########################

mean(op_0$rho_mf)
quantile(op_0$rho_mf, p=c(0.5, 0.025, 0.975))

mean(op_0$rho_fm)
quantile(op_0$rho_fm, p=c(0.5, 0.025, 0.975))

