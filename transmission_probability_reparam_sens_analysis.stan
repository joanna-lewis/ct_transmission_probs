/////////////////
// Model for inferring male-to-female and female-to-male per-partnership chlamydia transmission probabilities
// Sensitivity analysis
// Joanna Lewis
// Imperial College London 
// April 2019
/////////////////

functions {

  // integrator to help with calculating per-person prevalence
  real[] prev_int(real t,
                 real[] y,
                 real[] theta, // (alpha, beta, A, lambda)
                 real[] x_r,
                 int[] x_i) {
        real dydt[1];
        dydt[1] = theta[3]*exp(-theta[2]*t) * t^theta[1] / (theta[3]*t + theta[4]);
        return dydt;
  }

} 

data {
  
  // number of partners reported and infection status
  int<lower=0> N_n_m; // number of male survey participants 
  int<lower=0> n_m[N_n_m]; // number of partners reported by each man
  
  int<lower=0> N_n_f; // number of female survey participants 
  int<lower=0> n_f[N_n_f]; // number of partners reported by each woman

  // survey weights
  real<lower=0> wt_m[N_n_m]; // survey weights (men)
  real<lower=0> wt_f[N_n_f]; // survey weights (women)
  
  // infection status (chlamydia test result)
  int<lower=0> infect_m[N_n_m]; // infection status (men)
  int<lower=0> infect_f[N_n_f]; // infection status (women)
  
  // clearance rate data from literature - women (see Price et al. Stat. Med. 2013)
  int<lower=0> studnum_f; // number of studies 
  int<lower=0> studnum_bytype_f[3]; // number of studies 
  int<lower=0> studobs_f[studnum_f]; // number of observations (time periods) in each study
  int<lower=0> cumobs_f[studnum_f]; // cumulative number of observations at the start of each study
  int<lower=0> Nobs_f; // total number of observations (=sum(studobs))
  int<lower=0> r_f[Nobs_f]; // number who cleared infection at each time point 
  int<lower=0> n_test_f[Nobs_f]; // number tested at each time point 
  real<lower=0> t_f[Nobs_f]; // estimated mean follow-up 
  int<lower=0> seind_f[Nobs_f]; // did the study use culture as opposed to NAATs?
  real<lower=0> T_f[Nobs_f]; // already followed for...

  // clearance rate data from literature - women (see Lewis et al. JID 2017)
  int<lower=0> studnum_m; // number of studies 
  int<lower=0> studnum_bytype_m[3]; // number of studies 
  int<lower=0> studobs_m[studnum_m]; // number of observations (time periods) in each study
  int<lower=0> cumobs_m[studnum_m]; // cumulative number of observations at the start of each study
  int<lower=0> Nobs_m; // total number of observations (=sum(studobs))
  int<lower=0> r_m[Nobs_m]; // number who cleared infection at each time point 
  int<lower=0> n_test_m[Nobs_m]; // number tested at each time point 
  real<lower=0> t_m[Nobs_m]; // estimated mean follow-up 
  int<lower=0> seind_m[Nobs_m]; // did the study use culture as opposed to NAATs?
  real<lower=0> T_m[Nobs_m]; // already followed for...

  // population size, for simulating number of diagnoses
  int<lower=0> population[2]; // (male, female)

}

transformed data {

  real x_r[0]; // dummy for ode solver
  int x_i[0]; // dummy for ode solver
  real s_max[1]; // upper limit of integration for ode solver
  
  s_max[1] = 5000.0;
  
}


parameters {

  // hyperparameters
  real<lower=0> mu_nb; // parameterises gamma distribution
  real<lower=0> beta[2]; // parameterises gamma distribution (male, female)
 
  // clearance rate
  real<lower=0,upper=1> p1[2]; // probability of clearing fast (male, female)
  real<lower=0> lambda_slow[2]; // slow clearance rates (male, female)
  real<lower=0,upper=1> psi; // sensitivity of culture given initial positive culture (assume no sex difference)

  // transmission probability
  real<lower=0,upper=1> A[2,11]; // force of infection per partner (men, women): per-partnership prevalence in opposite sex, multiplied by transmission probability
  
  // observed prevalence
  real<lower=0,upper=1> pi_m_obs[11]; // prevalence in men reporting each number
  real<lower=0,upper=1> pi_f_obs[11]; // prevalence in women reporting each number

  
}

transformed parameters {

  real<lower=0> alpha[2]; // parameterises gamma distribution (male, female)
  
  real<lower=0,upper=1> w1[2]; //  expected proportion of fast clearers in a prevalent population (men,women)
  real theta_m[Nobs_m]; // weighted average of clearance probabilities for two classes (men)
  real theta_f[Nobs_f]; // weighted average of clearance probabilities for two classes (women)
  real<lower=0> lambda[2,2]; // (men, women)
  real<lower=0,upper=1> pk[2]; // temporary variable, containing probability of this k (men, women)
  
  real<lower=0,upper=1> pi_m[11]; // prevalence in men reporting each number
  real<lower=0,upper=1> pi_f[11]; // prevalence in women reporting each number
  
  real theta[2,4]; // parameters for ODE solver  
  real init[1]; // initial condition for ODE solver
  
//  real<lower=0,upper=1> ppp[2,11]; // per-partnership prevalence in partnerships offered by (men, women)
  
  alpha[1] = mu_nb * beta[1];
  alpha[2] = mu_nb * beta[2];
  
  /////////////////
  // for ODE solver
  /////////////////

  theta[1,1] = alpha[1];
  theta[1,2] = beta[1];
  theta[1,3] = A[1,1];
  theta[1,4] = lambda_slow[1];
  theta[2,1] = alpha[2];
  theta[2,2] = beta[2];
  theta[2,3] = A[2,1];
  theta[2,4] = lambda_slow[2];
  init[1] = 0.0;

  /////////////////
  // clearance probabilities - see Price et al. Stat. Med. 2013
  /////////////////

  lambda[1,1] = 49; // men
  lambda[1,2] = lambda_slow[1]; // men
  lambda[2,1] = 120; // women
  lambda[2,2] = lambda_slow[2]; // women

  //////////
  // men
  //////////
  
  // proportion of participants in each study time point expected to have recovered.
  // clinic studies
  for (i in 1:studnum_bytype_m[1]) { 
    for (j in 1:studobs_m[i]) {
      theta_m[cumobs_m[i]+j] =  0;
      for(k in 0:n_test_m[cumobs_m[i]+j]){
        theta_m[cumobs_m[i]+j] = theta_m[cumobs_m[i]+j] + exp( lchoose(n_test_m[cumobs_m[i]+j],k) ) * p1[1]^k * (1-p1[1])^(n_test_m[cumobs_m[i]+j]-k) *  ((k / (1.0*n_test_m[cumobs_m[i]+j])) * (1 - exp(-lambda[1,1] * t_m[cumobs_m[i]+j])) + (1 - (k / (1.0*n_test_m[cumobs_m[i]+j]))) * (1 - exp(-lambda[1,2] * t_m[cumobs_m[i]+j])));
        }
      if(seind_m[cumobs_m[i]+j])
        theta_m[cumobs_m[i]+j] = 1 + psi*(theta_m[cumobs_m[i]+j] - 1);
      }
    }
  // left-truncated; single observation
  w1[1] = (p1[1] / lambda[1,1]) / (p1[1] / lambda[1,1] + (1 - p1[1]) / lambda[1,2]);
  for (i in (studnum_bytype_m[1]+1):(studnum_bytype_m[1]+studnum_bytype_m[2])) { // should find a way of using input data to state study categories
    for (j in 1:studobs_m[i]) {
      theta_m[cumobs_m[i]+j] = 0;
      for(k in 0:n_test_m[cumobs_m[i]+j]){
        theta_m[cumobs_m[i]+j] = theta_m[cumobs_m[i]+j] + exp( lchoose(n_test_m[cumobs_m[i]+j],k) ) * w1[1]^k * (1-w1[1])^(n_test_m[cumobs_m[i]+j]-k) *  ((k / (1.0*n_test_m[cumobs_m[i]+j])) * (1 - exp(-lambda[1,1] * t_m[cumobs_m[i]+j])) + (1 - (k / (1.0*n_test_m[cumobs_m[i]+j]))) * (1 - exp(-lambda[1,2] * t_m[cumobs_m[i]+j])));
        }
      if(seind_m[cumobs_m[i]+j])
        theta_m[cumobs_m[i]+j] = 1 + psi*(theta_m[cumobs_m[i]+j] - 1);
      }
    }
  // left-truncated; repeat observations
  pk[1] = 0; // set pk, to avoid error messages if no repeat-observation studies
  if(studnum_bytype_m[3] != 0){
    for (i in (studnum_bytype_m[1]+studnum_bytype_m[2]+1):studnum_m) {
      for (j in 1:studobs_m[i]) {
          theta_m[cumobs_m[i]+j] = 0;
          for (k in 0:n_test_m[cumobs_m[i]+j]){
          
            pk[1] = exp( lchoose(n_test_m[cumobs_m[i]+j],k) ) * w1[1]^k * (1-w1[1])^(n_test_m[cumobs_m[i]+j]-k); // temporary variable, containing probability of this k

            theta_m[cumobs_m[i]+j] = theta_m[cumobs_m[i]+j] + pk[1] * (k / (1.0*n_test_m[cumobs_m[i]+1])) * exp(-lambda[1,1]*T_m[cumobs_m[i]+j]) * (1 - exp(-lambda[1,1]*t_m[cumobs_m[i]+j])) / ( (k / (1.0*n_test_m[cumobs_m[i]+1])) * exp(-lambda[1,1]*T_m[cumobs_m[i]+j]) + (1 - k / (1.0*n_test_m[cumobs_m[i]+1])) * exp(-lambda[1,2]*T_m[cumobs_m[i]+j]) );   // fast category

            theta_m[cumobs_m[i]+j] = theta_m[cumobs_m[i]+j] + pk[1] * (1 - k / (1.0*n_test_m[cumobs_m[i]+1])) * exp(-lambda[1,2]*T_m[cumobs_m[i]+j]) * (1 - exp(-lambda[1,2]*t_m[cumobs_m[i]+j])) / ( (k / (1.0*n_test_m[cumobs_m[i]+1])) * exp(-lambda[1,1]*T_m[cumobs_m[i]+j]) + (1 - k / (1.0*n_test_m[cumobs_m[i]+1])) * exp(-lambda[1,2]*T_m[cumobs_m[i]+j]) ) ;  // slow category
          
            }
        if(seind_m[cumobs_m[i]+j])
          theta_m[cumobs_m[i]+j] = 1 + psi*(theta_m[cumobs_m[i]+j] - 1);    
        }
      }
    }

  //////////
  // women
  //////////
  
  // proportion of participants in each study time point expected to have recovered.
  // clinic studies
  for (i in 1:studnum_bytype_f[1]) { 
    for (j in 1:studobs_f[i]) {
      theta_f[cumobs_f[i]+j] =  0;
      for(k in 0:n_test_f[cumobs_f[i]+j]){
        theta_f[cumobs_f[i]+j] = theta_f[cumobs_f[i]+j] + exp( lchoose(n_test_f[cumobs_f[i]+j],k) ) * p1[2]^k * (1-p1[2])^(n_test_f[cumobs_f[i]+j]-k) *  ((k / (1.0*n_test_f[cumobs_f[i]+j])) * (1 - exp(-lambda[2,1] * t_f[cumobs_f[i]+j])) + (1 - (k / (1.0*n_test_f[cumobs_f[i]+j]))) * (1 - exp(-lambda[2,2] * t_f[cumobs_f[i]+j])));
        }
      if(seind_f[cumobs_f[i]+j])
        theta_f[cumobs_f[i]+j] = 1 + psi*(theta_f[cumobs_f[i]+j] - 1);
      }
    }
  // left-truncated; single observation
  w1[2] = (p1[2] / lambda[2,1]) / (p1[2] / lambda[2,1] + (1 - p1[2]) / lambda[2,2]);
  for (i in (studnum_bytype_f[1]+1):(studnum_bytype_f[1]+studnum_bytype_f[2])) { // should find a way of using input data to state study categories
    for (j in 1:studobs_f[i]) {
      theta_f[cumobs_f[i]+j] = 0;
      for(k in 0:n_test_f[cumobs_f[i]+j]){
        theta_f[cumobs_f[i]+j] = theta_f[cumobs_f[i]+j] + exp( lchoose(n_test_f[cumobs_f[i]+j],k) ) * w1[2]^k * (1-w1[2])^(n_test_f[cumobs_f[i]+j]-k) *  ((k / (1.0*n_test_f[cumobs_f[i]+j])) * (1 - exp(-lambda[2,1] * t_f[cumobs_f[i]+j])) + (1 - (k / (1.0*n_test_f[cumobs_f[i]+j]))) * (1 - exp(-lambda[2,2] * t_f[cumobs_f[i]+j])));
        }
      if(seind_f[cumobs_f[i]+j])
        theta_f[cumobs_f[i]+j] = 1 + psi*(theta_f[cumobs_f[i]+j] - 1);
      }
    }
  // left-truncated; repeat observations
  pk[2] = 0; // set pk, to avoid error messages if no repeat-observation studies
  if(studnum_bytype_f[3] != 0){
    for (i in (studnum_bytype_f[1]+studnum_bytype_f[2]+1):studnum_f) {
      for (j in 1:studobs_f[i]) {
          theta_f[cumobs_f[i]+j] = 0;
          for (k in 0:n_test_f[cumobs_f[i]+j]){
          
            pk[2] = exp( lchoose(n_test_f[cumobs_f[i]+j],k) ) * w1[2]^k * (1-w1[2])^(n_test_f[cumobs_f[i]+j]-k); // temporary variable, containing probability of this k

            theta_f[cumobs_f[i]+j] = theta_f[cumobs_f[i]+j] + pk[2] * (k / (1.0*n_test_f[cumobs_f[i]+1])) * exp(-lambda[2,1]*T_f[cumobs_f[i]+j]) * (1 - exp(-lambda[2,1]*t_f[cumobs_f[i]+j])) / ( (k / (1.0*n_test_f[cumobs_f[i]+1])) * exp(-lambda[2,1]*T_f[cumobs_f[i]+j]) + (1 - k / (1.0*n_test_f[cumobs_f[i]+1])) * exp(-lambda[2,2]*T_f[cumobs_f[i]+j]) );   // fast category
            theta_f[cumobs_f[i]+j] = theta_f[cumobs_f[i]+j] + pk[2] * (1 - k / (1.0*n_test_f[cumobs_f[i]+1])) * exp(-lambda[2,2]*T_f[cumobs_f[i]+j]) * (1 - exp(-lambda[2,2]*t_f[cumobs_f[i]+j])) / ( (k / (1.0*n_test_f[cumobs_f[i]+1])) * exp(-lambda[2,1]*T_f[cumobs_f[i]+j]) + (1 - k / (1.0*n_test_f[cumobs_f[i]+1])) * exp(-lambda[2,2]*T_f[cumobs_f[i]+j]) ) ;  // slow category
          
            }
        if(seind_f[cumobs_f[i]+j])
          theta_f[cumobs_f[i]+j] = 1 + psi*(theta_f[cumobs_f[i]+j] - 1);    
        }
      }
    }

  /////////////////
  // per-partnership prevalence
  /////////////////

//  for(i in 1:11){
//    theta[1,1] = alpha[1]+i-1;
//    theta[1,2] = beta[1] + 1;
//    theta[1,3] = A[1,i];
//    ppp[1,i] = (beta[1]+1)^(alpha[1]+i) * integrate_ode_rk45(ppp_int, init, 0.0, s_max, theta[1], x_r, x_i)[1,1] / ((alpha[1]+i-1)*tgamma(alpha[1]+i-1) );
//    if(is_nan(ppp[1,i]))
//      ppp[1,i] = -999;

//    theta[2,1] = alpha[2]+i-1;
//    theta[2,2] = beta[2] + 1;
//    theta[2,3] = A[2,i];
//    ppp[2,i] = (beta[2]+1)^(alpha[2]+i) * integrate_ode_rk45(ppp_int, init, 0.0, s_max, theta[2], x_r, x_i)[1,1] / ((alpha[2]+i-1)*tgamma(alpha[2]+i-1));
//    if(is_nan(ppp[2,i]))
//      ppp[2,i] = -999;
//  }
    
//  print(ppp)
    
  /////////////////
  // prevalence
  /////////////////

  // prevalence in men reporting different numbers of partners
  for(i in 1:11){
    theta[1,1] = alpha[1]+i-1;
    theta[1,2] = beta[1] + 1;
    theta[1,3] = A[1,i];
    pi_m[i] = (beta[1]+1)^(alpha[1]+i-1) * integrate_ode_rk45(prev_int, init, 0.0, s_max, theta[1], x_r, x_i)[1,1] / tgamma(alpha[1]+i-1);
  }
	  
  // prevalence in women reporting different numbers of partners
  for(i in 1:11){
    theta[2,1] = alpha[2]+i-1;
    theta[2,2] = beta[2] + 1;
    theta[2,3] = A[2,i];
    pi_f[i] = (beta[2]+1)^(alpha[2]+i-1) * integrate_ode_rk45(prev_int, init, 0.0, s_max, theta[2], x_r, x_i)[1,1] / tgamma(alpha[2]+i-1);
  }

}

model {

  ///////////
  // priors
  ///////////

  mu_nb ~ exponential(0.1); // parameterises negative binomial
  beta ~ exponential(0.1); // parameterises negative binomial 
  p1 ~ beta(1,1);
  lambda_slow ~ exponential(0.001);
  psi ~ beta(78, 8);

//  print(A);
//  print(pi_f);
//  print(pi_m);

  for(i in 1:11){
    A[1,i] ~ exponential(0.001); //uniform(0,1); //ppp[2,i]); // force of infection / sigma
    A[2,i] ~ exponential(0.001); //uniform(0,1); //ppp[1,i]); 
  }

  pi_m_obs ~ beta(1,1);
  pi_f_obs ~ beta(1,1);

  ///////////
  // partner number
  ///////////
    
  // men   
  for (i in 1:N_n_m) 
    target += wt_m[i] * neg_binomial_lpmf(n_m[i] | alpha[1], beta[1]);
        
  // women
  for (i in 1:N_n_f)
    target += wt_f[i] * neg_binomial_lpmf(n_f[i] | alpha[2], beta[2]);
    
  ///////////
  // clearance rate
  ///////////

  for (i in 1:studnum_f) {
    for (j in 1:studobs_f[i]) {
      r_f[cumobs_f[i]+j] ~ binomial(n_test_f[cumobs_f[i]+j], theta_f[cumobs_f[i]+j]);
      }
    }
    
  for (i in 1:studnum_m) {
    for (j in 1:studobs_m[i]) {
      r_m[cumobs_m[i]+j] ~ binomial(n_test_m[cumobs_m[i]+j], theta_m[cumobs_m[i]+j]);
      }
    }

  ///////////
  // infection status
  ///////////
  
  // men
  for(i in 1:N_n_m){
    if(n_m[i] <=10){
      target += wt_m[i] * bernoulli_lpmf(infect_m[i] | pi_m[n_m[i]+1]);  
      target += wt_m[i] * bernoulli_lpmf(infect_m[i] | pi_m_obs[n_m[i]+1]);  
      }
    }

  // women
  for (i in 1:N_n_f){ 
    if(n_f[i] <= 10){
      target += wt_f[i] * bernoulli_lpmf(infect_f[i] | pi_f[n_f[i]+1]); 
      target += wt_f[i] * bernoulli_lpmf(infect_f[i] | pi_f_obs[n_f[i]+1]); 
      }
    }
    
}
  
generated quantities {

//  int<lower=0> n_m_sim[N_n_m]; // simulate numbers reported by men
//  int<lower=0> n_f_sim[N_n_f]; // simulate numbers reported by women

//  real<lower=0> sim_infect_m[N_n_m]; // simulated infection status (men)
//  real<lower=0> sim_infect_f[N_n_f]; // simulated infection status (women)

  real rho_mf[11]; // male-to-female transmission probability
  real rho_fm[11]; // female-to-male transmission probability

//  real<lower=0,upper=1> phi_sim[2]; // proportion of infections prompting fast treatment
//  int<lower=0> diag_sim[2]; // simulated number of diagnoses
  
  // transmission probability
  for(i in 1:11){
    rho_mf[i] = A[2,i]/(pi_m_obs[i]);
    rho_fm[i] = A[1,i]/(pi_f_obs[i]);
  }
  
  //////////
  // simulate partnership numbers
  //////////
  
  // men
//  for(i in 1:N_n_m)
//    n_m_sim[i] = neg_binomial_rng(alpha[1], beta[1]);
  
  // women  
//  for(i in 1:N_n_f)
//    n_f_sim[i] = neg_binomial_rng(alpha[2], beta[2]);
    
  //////////
  // simulate infection status
  //////////
  
//  for(i in 1:3)
//    sim_infect_m[i] = bernoulli_rng( pi_m[n_m[i]+1] );
  
//  for(i in 1:3)
//    sim_infect_f[i] = bernoulli_rng( pi_f[n_f[i]+1] );
    
//  //////////
//  // simulate number of diagnoses
//  //////////

//  phi_sim[1] = beta_rng(11, 5); // based on Geisler 2008 (all men positive at baseline)
//  phi_sim[2] = beta_rng(27, 90); // based on Geisler 2008 (only women re-testing as positive)
  
//  diag_sim[1] = poisson_rng( population[1] * A[1] * alpha[1] * phi_sim[1] / beta[1] ); 
//  diag_sim[2] = poisson_rng( population[2] * A[2] * alpha[2] * phi_sim[2] / beta[2] ); 
    
}

