// The input data.
data {
  int<lower=0> N;   // number of observations
  int<lower=0> C;   // number of covariates
  int<lower=0> K;   // number of patterns
  matrix[N, C] x;   // covariate matrix (minus sex)
  vector[N] y;      // outcome vector
  vector[N] sex;    // sex vector (separate for interaction)
  
  matrix<lower=0>[N,K] ewa;    // mu for patterns
  matrix<lower=0>[N,K] sd_ewa; // std dev for patterns
}

// The parameters accepted by the model.
parameters {
  real alpha;              // intercept
  real<lower=0> sigma;     // error scale
  
  vector[C] beta_c;        // coefficients for covariates
  real beta_int;      // coefficient for interaction
  real beta_sex;      // coefficient for sex
  
  vector[K] beta_p;         // coefficients for patterns
  matrix<lower=0>[N,K] WA;  // pattern score with uncertainty
}

// interaction term here
// Unknown but are known given the values of the objects in the parameters block
// Saved in the output and hence should be of interest to the researcher
// Are usually the arguments to the log-likelihood function that is evaluated in the model block, 
// although in hierarchical models the line between the prior and the likelihood can be drawn in multiple ways
transformed parameters {
    vector<lower=0>[N] P1;
    vector<lower=0>[N] inter; // pattern score * sex
    P1 = WA[,1]; // extract for interaction
    inter = P1 .* sex;
}

// The model to be estimated.
//  With no prior in the model block, the effect is an improper prior on all real numbers. 
model {
  // sigma ~ inv_gamma(0.001, 0.001); // prior on error
  // Cite: stat.columbia.edu/~gelman/research/published/taumain.pdf
  
  // alpha ~ normal(100, 15); // prior on alpha = IQ

  // beta coefficient priors
  // student t parameters: degrees of freedom nu, location mu, and scale sigma
  // smaller nu, fatter tails
  // beta_c ~ normal(0,5); // student_t(1, 0, 2.5);
  //  beta_p ~ normal(0,5); // student_t(1, 0, 2.5);
  
  for (n in 1:N) { // prior on data
    for (k in 1:K) {
      WA[n,k] ~ normal(ewa[n,k], sd_ewa[n,k]);
    }}
  
  y ~ normal(((WA * beta_p) + (inter * beta_int) + 
              (sex * beta_sex) + (x * beta_c) + alpha), sigma);  // likelihood

}

// to get predicted values
generated quantities {
  real y_tilde[N] = normal_rng(((WA * beta_p) + (inter * beta_int) + 
                                (sex * beta_sex) + (x * beta_c) + alpha), sigma);
}




