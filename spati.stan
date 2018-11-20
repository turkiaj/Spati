data { 
  int<lower=0> N;   // number of observations
  int<lower=1> p;   // number of predictors
  int<lower=1> J;   // number of groups in data 
  int<lower=1> k;   // number of group-level predictors
  int<lower=1,upper=J> group[N]; //group indicator
  matrix[N,p] X;    // fixed-effect design matrix
  matrix[N,k] Z;    // random-effect design matrix
  vector[N] Y;      // response
} 

transformed data { 
  int Pc; 
  matrix[N, p - 1] Xc;    // Centered X    
  matrix[N, p - 1] Xp;    // Original X without intercept
  vector[p - 1] means_X;  // column means of X before centering 
  
  Pc = p - 1;  // the intercept is removed from the design matrix 
  for (i in 2:p) { 
    means_X[i - 1] = mean(X[, i]); 
    Xc[, i - 1] = X[, i] - means_X[i - 1]; 
    Xp[, i - 1] = X[, i]; 
  } 
}

parameters { 

  vector[Pc] beta;                // fixed effect coefficients
  real temp_Intercept;           
  cholesky_factor_corr[k] L;      // Cholesky factor of group ranef corr matrix
  vector<lower=0>[k] sigma_b;     // group-level random-effect standard deviations
  real<lower=0> sigma_e;          // residual standard deviations 
  vector[k] z[J];                 // unscaled group-level effects
  real<lower=-1,upper=1> ar1;     // autoregressive effect, AR(1)
}

transformed parameters {
  
  vector[k] b[J];                 // random effects
  matrix[k, k] D;		              // variance-covariance matrix of random-effects
  
  // Premultiply diagonal matrix [sigma_b] with the Cholesky decomposition L of
  // the correlation matrix D to get variance-covariance matrix of group-level effects

  // diag(sigma_b) * L
  D = diag_pre_multiply(sigma_b, L); 
  
  // Group-level effects are generated by multipying D with z 
  // that has standard normal distribution
    
  for(j in 1:J) 
    b[j] = D * z[j];
}

model { 
  vector[N] mu;
  
  // Residuals
  vector[N] e;
  
  // Group variables for AR
  int group_size;
  int current_group;
 
  // Priors
  sigma_e ~ student_t(3, 0, 10);
  sigma_b ~ student_t(3, 0, 10);
  L ~ lkj_corr_cholesky(1); 
  
  // Standard normal prior for random effects
  for (j in 1:J)
    z[j] ~ normal(0,1);

  // Likelihood 
  mu = temp_Intercept + Xc * beta;
  
  // - add group effects
  for (n in 1:N) 
  {
     mu[n] += Z[n] * b[group[n]]; 
  }

  Y ~ normal(mu, sigma_e);
}

generated quantities {
  real beta_Intercept;             
  vector[N] Y_rep;               // Repeated response
  real mu;

  // Correlation matrix of random-effects, C = L'L
  corr_matrix[k] C; 
  C = multiply_lower_tri_self_transpose(L); 
  
  beta_Intercept = temp_Intercept - dot_product(means_X, beta);
  
  // Posterior predictive distribution for model checking

  for (n in 1:N) 
  {
    // In Bayesian statistics, group variation is considered as part of mu, not variance
    mu = beta_Intercept + Xp[n] * beta + Z[n] * b[group[n]];
    Y_rep[n] = normal_rng(mu, sigma_e);
  }
  
} 