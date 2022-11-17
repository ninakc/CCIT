

# Step 1: Generate a synthetic data set.
## Step 1a: Generate covariates and treatment
synth_data_covs <- generate_syn_data_covs(sample_size = 10000, gps_spec = 2)

## Step 1b: Determine covariate balance
absolute_corr_fun(as.data.table(synth_data_covs$treat), as.data.table(synth_data_covs[c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6')]))


## Step 1c: Add outcome to generated data
synth_data <-
  do.call(
    generate_syn_data_outcome,
    c(as.list(synth_data_covs), beta = 1, em_spec = 1, outcome_sd = 1, heterogenous_intercept = FALSE)
  )

# Step 2: Split Data into Exploration, Validation, and Inference
split_synth_data <- split_dataset(synth_data, num_exposure_cats = 10)

# Step 3: Match the data on covariates.
matched_synth_data <- stratified_GPS_matching(
  split_synth_data, 0.1, 'treat', c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'), c('em1', 'em2', 'em3', 'em4', 'subsample'), 'Y'
)

# Step 4: Combine generated matched data.

# Step 5: Conduct conditional average treatment effect.
CCIT(
  filter(matched_synth_data, subsample == 'exploration'),
  matched.validation.sample.outcomes = filter(matched_synth_data, subsample == 'validation'),
  matched.inference.sample.outcomes = filter(matched_synth_data, subsample == 'inference'),
  lambdas=c(1),
  stopping.rule = TRUE
  )
