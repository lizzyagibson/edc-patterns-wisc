// Noninformative priors on everything

// The input data.
data {
  int<lower=0> N;   // number of observations
  int<lower=0> C;   // number of covariates
  matrix[N, C] x;   // covariate matrix (minus sex)
  vector[N] y;      // outcome vector
  vector[N] sex;    // sex vector (separate for interaction)
  
  vector<lower=0>[N] alpha_w;    // gamma likelihood
  vector<lower=0>[N] beta_w;     // gamma likelihood
  
  real<lower=0> alpha_a;    // gamma likelihood
  real<lower=0> beta_a;     // gamma likelihood
}

// The parameters accepted by the model.
parameters {
  real<lower=0> alpha;     // intercept // Robbie: Won't affect analysis but called beta in the manuscript
  real<lower=0> sigma;     // error scale
  
  vector[C] beta_c;    // coefficients for covariates
  real beta_int;       // coefficient for interaction
  real beta_sex;       // coefficient for sex
  real beta_p;         // coefficients for patterns
  
  vector<lower=0>[N] W;  // pattern scores (unscaled) with uncertainty
  real<lower=0> a;       // individual value of sparse vector (with uncertainty) to scale scores
  // Robbie: Is this the sparse vector a in the manuscript? If so perhaps make a bit more clear here?
  // Lizzy: Yes!
}

// interaction term here
// Unknown but are known given the values of the objects in the parameters block
// Are usually the arguments to the log-likelihood function that is evaluated in the model block, 
// although in hierarchical models the line between the prior and the likelihood can be drawn in multiple ways
transformed parameters {
    vector<lower=0>[N] WA; // Wa drawn from distribution
    vector<lower=0>[N] inter; // pattern score * sex
    WA    = W * a;
    inter = WA .* sex; 
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
    
  for (n in 1:N) { // draw patterns scores from this distribution
      W[n] ~ gamma(alpha_w[n], beta_w[n]);
  }
  
  a ~ gamma(alpha_a, beta_a); // single a distribution, all W are multiplied by this 
  // distribution for individual value of sparse vector (with uncertainty) to scale scores  
  // Robbie: See above. This is the sparse matrix to penalise larger number of patterns correct?
  // Lizzy: Yes!
  y ~ normal(((WA * beta_p) + (inter * beta_int) + 
              (sex * beta_sex) + (x * beta_c) + alpha), sigma);  // likelihood
}

// to get predicted values
generated quantities {
  real y_tilde[N] = normal_rng(((WA * beta_p) + (inter * beta_int) + 
                                (sex * beta_sex) + (x * beta_c) + alpha), sigma);
}




