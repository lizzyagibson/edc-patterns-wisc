// Noninformative priors on everything
// Include both patterns in single model

// The input data.
data {
  int<lower=0> N;   // number of observations
  int<lower=0> C;   // number of covariates
  matrix[N, C] x;   // covariate matrix (minus sex)
  vector[N] y;      // outcome vector
  vector[N] sex;    // sex vector (separate for interaction)
  
  // pattern 1
  vector<lower=0>[N] alpha_w_1;    // gamma likelihood
  vector<lower=0>[N] beta_w_1;     // gamma likelihood
  
  real<lower=0> alpha_a_1;    // gamma likelihood
  real<lower=0> beta_a_1;     // gamma likelihood
  
  // pattern 2
  vector<lower=0>[N] alpha_w_2;    // gamma likelihood
  vector<lower=0>[N] beta_w_2;     // gamma likelihood
  
  real<lower=0> alpha_a_2;    // gamma likelihood
  real<lower=0> beta_a_2;     // gamma likelihood
}

// The parameters accepted by the model.
parameters {
  real<lower=0> alpha;     // intercept // Robbie: Won't affect analysis but called beta in the manuscript
  real<lower=0> sigma;     // error scale
  
  vector[C] beta_c;    // coefficients for covariates
  
  // interaction term for each pattern with sex
  real beta_int_1;       // coefficient for interaction with pattern 1
  real beta_int_2;       // coefficient for interaction with pattern 2
  
  real beta_sex;       // coefficient for sex
  
  // coefficient for each pattern
  real beta_p_1;         // coefficients for pattern 1
  real beta_p_2;         // coefficients for pattern 2

  // Pattern scores for each pattern  
  vector<lower=0>[N] W_1;  // pattern scores (unscaled) with uncertainty
  real<lower=0> a_1;       // individual value of sparse vector (with uncertainty) to scale scores
  
  vector<lower=0>[N] W_2;  // pattern scores (unscaled) with uncertainty
  real<lower=0> a_2;       // individual value of sparse vector (with uncertainty) to scale scores
  
  // Robbie: Is this the sparse vector a in the manuscript? If so perhaps make a bit more clear here?
  // Lizzy: Yes!
}

// interaction term here
// Unknown but are known given the values of the objects in the parameters block
// Are usually the arguments to the log-likelihood function that is evaluated in the model block, 
// although in hierarchical models the line between the prior and the likelihood can be drawn in multiple ways
transformed parameters {
  
    // create transformed parameters for each pattern and interaction
    vector<lower=0>[N] WA_1; // Wa drawn from distribution
    vector<lower=0>[N] inter_1; // pattern score * sex
    vector<lower=0>[N] WA_2; // Wa drawn from distribution
    vector<lower=0>[N] inter_2; // pattern score * sex
    
    // pattern 1
    WA_1    = W_1 * a_1;
    inter_1 = WA_1 .* sex; 
    
    // pattern 2
    WA_2    = W_2 * a_2;
    inter_2 = WA_2 .* sex; 
    
    // Robbie: What does the '.' in front of the '*' mean sorry? I understand (I think) that this is the interaction term between 
    // sex and individual pattern score.
    // Lizzy: .* makes it element-wise multiplication. For W * a, a is a scalar value (single number), so W * a multiplies every value
    // in the W vector by a.
    // For WA .* sex, WA and sex are both vectors, so WA * sex would give vector multiplication, which we don't want.
    // WA .* sex gives each individual's WA value times sex, which is what we want for the interaction term.
}

// The model to be estimated.
//  With no prior in the model block, the effect is an improper prior on all real numbers. 
model {
    
  // Draw pattern scores for each pattern
  for (n in 1:N) { // draw patterns scores from this distribution
      W_1[n] ~ gamma(alpha_w_1[n], beta_w_1[n]);
  }
  
  for (n in 1:N) { // draw patterns scores from this distribution
      W_2[n] ~ gamma(alpha_w_2[n], beta_w_2[n]);
  }

  a_1 ~ gamma(alpha_a_1, beta_a_1); // single a distribution, all W are multiplied by this 
  a_2 ~ gamma(alpha_a_2, beta_a_2); // single a distribution, all W are multiplied by this 
  // distribution for individual value of sparse vector (with uncertainty) to scale scores  
  // Robbie: See above. This is the sparse matrix to penalise larger number of patterns correct?
  // Lizzy: Yes!
  
  // Include both patterns and interactions with sex
  y ~ normal(((WA_1 * beta_p_1) + (inter_1 * beta_int_1) + (sex * beta_sex) + 
              (WA_2 * beta_p_2) + (inter_2 * beta_int_2) + 
              (x * beta_c) + alpha), sigma);  // likelihood
}

// to get predicted values
generated quantities {
  real y_tilde[N] = normal_rng(((WA_1 * beta_p_1) + (inter_1 * beta_int_1) + (sex * beta_sex) + 
                                (WA_2 * beta_p_2) + (inter_2 * beta_int_2) +
                                (x * beta_c) + alpha), sigma);
}




