# edc-patterns-wisc
Patterns of phthalates, phenols, and parabens measured in pregnant women in Mothers &amp; Newborns cohort

## in this repo:

[codebook.md](codebook.md) explains variables.

### sourced from other scripts

1. `packages.R` has all the packages used and a few functions I wrote.  
2. `read_data.R` reads in the data. :)

### analysis 

1. `summary_ppp.Rmd` has descriptive stats, table 1, table 2, OLS health models.
  * You cannot run this file.
  * HTML version included.
2. `bayes_inf_and_non.Rmd` has Bayesian health models.
  * this calls Stan files
  * You can run this file with the Stan model fits.
  * HTML version included.
3. `sensitivity.Rmd` makes the supplemental figure for the sensitivity analysis.
  * You cannot run this file.
  * HTML version included.

### Stan

1. `bn2mf_score_reg_int_gamma.stan` has the completely non-informative prior model (main model).
2. `bn2mf_score_reg_int_infp1_nonc.stan` and `bn2mf_score_reg_int_infp2_nonc.stan` have informative priors on everything but covariates.
3. `bn2mf_score_reg_int_infp1.stan` and `bn2mf_score_reg_int_infp2.stan` have informative priors on everything.

### ignore

1. `misc` folder.
2. `matlab` folder.
* Nothing in these affects current analysis.
