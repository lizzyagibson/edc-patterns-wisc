// Noninformative priors on everything

// The input data.
data {
  int<lower=0> N;   // number of observations
  int<lower=0> C;   // number of covariates
  matrix[N, C] x;   // covariate matrix (minus sex)
  vector[N] y;      // outcome vector
  vector[N] sex;    // sex vector (separate for interaction)
  
  vector<lower=0>[N] ewa;    // mu for pattern
  vector<lower=0>[N] sd_ewa; // std dev for pattern
  
  matrix[N, C] x_const; // covariates for prediction -- this lets me get post predicted values holding other things const
}

// The parameters accepted by the model.
parameters {
  real<lower=0> alpha;     // intercept
  real<lower=0> sigma;     // error scale
  
  vector[C] beta_c;    // coefficients for covariates
  real beta_int;       // coefficient for interaction
  real beta_sex;       // coefficient for sex
  real beta_p;         // coefficients for patterns
  
  vector<lower=0>[N] WA;  // pattern scores with uncertainty
}

// interaction term here
// Unknown but are known given the values of the objects in the parameters block
// Saved in the output and hence should be of interest to the researcher
// Are usually the arguments to the log-likelihood function that is evaluated in the model block, 
// although in hierarchical models the line between the prior and the likelihood can be drawn in multiple ways
transformed parameters {
    vector<lower=0>[N] inter; // pattern score * sex
    inter = WA .* sex;
}

// The model to be estimated.
//  With no prior in the model block, the effect is an improper prior on all real numbers. 
model {

  for (n in 1:N) { // draw patterns scores from this distribution
      WA[n] ~ normal(ewa[n], sd_ewa[n]);
    }
  
  y ~ normal(((WA * beta_p) + (inter * beta_int) + 
              (sex * beta_sex) + (x * beta_c) + alpha), sigma);  // likelihood
}

// to get predicted values
generated quantities {
  real y_tilde[N] = normal_rng(((WA * beta_p) + (inter * beta_int) + 
                                (sex * beta_sex) + (x * beta_c) + alpha), sigma);
                                
  real y_f[N] = normal_rng(((WA * beta_p) + (WA*beta_int) + 
                                (beta_sex) + (x_const * beta_c) + alpha), sigma);
  real y_m[N] = normal_rng(((WA * beta_p) + 
                                (x_const * beta_c) + alpha), sigma);
}




