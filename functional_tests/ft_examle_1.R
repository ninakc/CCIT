# set this to the directory with source code
source_dir <- '..'
source(paste0(source_dir, '/R/generate_synthetic_data_covs.R'))
source(paste0(source_dir, '/R/generate_synthetic_data_outcome.R'))
source(paste0(source_dir, '/_R/split_dataset.R'))
source(paste0(source_dir, '/_R/stratified_GPS_matching.R'))
source(paste0(source_dir, '/_R/CCIT.R'))
source(paste0(source_dir, '/_R/rpart_funcs.R'))
source(paste0(source_dir, '/_R/create.sequence.R'))
source(paste0(source_dir, '/_R/evaluate.sequence.R'))

library(truncnorm)
library(CausalGPS)
library(dplyr)
library(data.table)
library(caret)
library(rpart)

# Step 1: Generate a synthetic data set.
## Step 1a: Generate covariates and treatment
synth_data_covs <- generate_syn_data_covs(sample_size = 10000, gps_spec = 2)

## Step 1b: Determine covariate balance
org_cor <- CausalGPS::absolute_corr_fun(
  data.table::as.data.table(synth_data_covs$treat),
  data.table::as.data.table(synth_data_covs[c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6')])
)

print(org_cor)

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
CCIT_result <- CCIT(
  filter(matched_synth_data, subsample == 'exploration'),
  matched.validation.sample.outcomes = filter(matched_synth_data, subsample == 'validation'),
  matched.inference.sample.outcomes = filter(matched_synth_data, subsample == 'inference'),
  lambdas=c(1),
  stopping.rule = TRUE
  )
