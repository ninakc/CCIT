# load libraries
library(CCIT)


# Step 1: Generate a synthetic data set.
## Step 1a: Generate covariates and treatment
syn_data_covs <- generate_syn_data_covs(sample_size = 10000, gps_spec = 2)

## Step 1b: Determine covariate balance
org_cor <- CausalGPS::absolute_corr_fun(
  data.table::as.data.table(syn_data_covs$treat),
  data.table::as.data.table(syn_data_covs[c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6')])
)

print(org_cor)

## Step 1c: Add outcome to generated data
syn_data <- generate_syn_data_outcome(cf = syn_data_covs[, c("cf1", "cf2",
                                                             "cf3", "cf4",
                                                             "cf5", "cf6")],
                                      em = syn_data_covs[, c("em1", "em2",
                                                             "em3", "em4")],
                                      treat = syn_data_covs[, c("treat")],
                                      outcome_sd = 1,
                                      em_spec = 1,
                                      heterogenous_intercept = FALSE,
                                      beta = 1)


# Step 2: Split Data into Exploration, Validation, and Inference
split_syn_data <- split_dataset(syn_data, num_exposure_cats = 10)


# Step 3,4: Match the data on covariates for each stratum and combine them.
matched_syn_data <- stratified_GPS_matching(data = split_syn_data,
                                            delta = 0.1,
                                            exposure_name = 'treat',
                                            confounders_names = c('cf1', 'cf2',
                                                                  'cf3', 'cf4',
                                                                  'cf5', 'cf6'),
                                            names_of_strata_vars = c('em1',
                                                                     'em2',
                                                                     'em3',
                                                                     'em4',
                                                                     'subsample'),
                                            outcome_name = 'Y')


# Step 5: Conduct conditional average treatment effect.
CCIT_result <- CCIT(
  filter(matched_syn_data, subsample == 'exploration'),
  matched.validation.sample.outcomes = filter(matched_syn_data, subsample == 'validation'),
  matched.inference.sample.outcomes = filter(matched_syn_data, subsample == 'inference'),
  lambdas=c(1),
  stopping.rule = TRUE
)
