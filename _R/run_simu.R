
# wrapper function to simulate data, fit tree, and evaluate outcomes
run_simu <- function(gps_spec = 1, num_exposure_cats, sample_size = 20000,
                     em_spec = 1, heterogenous_intercept = FALSE, beta = NULL, outcome_sd = 1, correct_splits = NULL, true_trt_effect_func = NULL, noise.var,
                     n_trials, lambdas, stopping.rule,
                     exploration.sample_covs = NULL, inference.sample_covs = NULL, val.sample_covs = NULL,
                     matched.exploration.sample = NULL, matched.validation.sample = NULL, matched.inference.sample = NULL) {

  # if the matched covariate dataset (everything but the outcome) is not passed in, generate one
  if (is.null(matched.exploration.sample) | is.null(matched.validation.sample) | is.null(matched.inference.sample) |
    is.null(exploration.sample_covs) | is.null(val.sample_covs) | is.null(inference.sample_covs)) {
    # generate covariate data
    synth_data_covs <- generate_syn_data_covs(sample_size = sample_size, gps_spec = gps_spec)

    # discretize treatment values for gps matching
    a.vals <- seq(min(synth_data_covs$treat), max(synth_data_covs$treat), length.out = num_exposure_cats)
    delta_n <- a.vals[2] - a.vals[1]

    synth_data_covs <-
      synth_data_covs %>%
      mutate(treat_level = cut(treat, breaks = a.vals)) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)

    # split data into subsamples, stratifying on exposure level bins
    synth_data_covs <-
      split_dataset(data = synth_data_covs, exposure_bins = synth_data_covs$treat_level)

    exploration.sample_covs <-
      synth_data_covs %>%
      filter(subsample == "exploration") %>%
      tibble::rowid_to_column("orig_id")

    val.sample_covs <-
      synth_data_covs %>%
      filter(subsample == "validation") %>%
      tibble::rowid_to_column("orig_id")

    inference.sample_covs <-
      synth_data_covs %>%
      filter(subsample == "inference") %>%
      tibble::rowid_to_column("orig_id")

    # GPS matching within subsamples and effect modifier groups
    matched.exploration.sample <-
      stratified_GPS_matching(exploration.sample_covs, delta_n,
        bin_seq = a.vals, exposure_name = "treat",
        confounders_names = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
        em_names = c("em1", "em2", "em3", "em4"),
        outcome_name = NA
      ) %>%
      mutate(id = row_number()) %>%
      # rename(w = treat) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    # print('here')
    # print(names(matched.exploration.sample))

    matched.validation.sample <-
      stratified_GPS_matching(val.sample_covs, delta_n,
        bin_seq = a.vals, exposure_name = "treat",
        confounders_names = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
        em_names = c("em1", "em2", "em3", "em4"),
        outcome_name = NA
      ) %>%
      mutate(id = row_number()) %>%
      # rename(w = treat) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)

    matched.inference.sample <-
      stratified_GPS_matching(inference.sample_covs, delta_n,
        bin_seq = a.vals, exposure_name = "treat",
        confounders_names = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
        em_names = c("em1", "em2", "em3", "em4"),
        outcome_name = NA
      ) %>%
      mutate(id = row_number()) %>%
      # rename(w = treat) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
  }

  results <- data.frame()
  if (n_trials > 0) {
    for (i in 1:n_trials) {
      # generate outcome data
      exploration.sample_outcome <-
        do.call(
          generate_syn_data_outcome,
          c(as.list(exploration.sample_covs %>% select(-treat_level, -subsample, -orig_id)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd, heterogenous_intercept = heterogenous_intercept)
        ) %>%
        tibble::rowid_to_column()

      validation.sample_outcome <-
        do.call(
          generate_syn_data_outcome,
          c(as.list(val.sample_covs %>% select(-treat_level, -subsample, -orig_id)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd, heterogenous_intercept = heterogenous_intercept)
        ) %>%
        tibble::rowid_to_column()

      inference.sample_outcome <-
        do.call(
          generate_syn_data_outcome,
          c(as.list(inference.sample_covs %>% select(-treat_level, -subsample, -orig_id)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd, heterogenous_intercept = heterogenous_intercept)
        ) %>%
        tibble::rowid_to_column()

      # add outcome data into matched samples
      matched.exploration.sample.outcomes <-
        matched.exploration.sample %>%
        select(-Y, -row_index) %>%
        left_join(exploration.sample_outcome, by = c(
          "orig_id" = "rowid", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "em1", "em2", "em3", "em4",
          "treat"
        ))

      matched.validation.sample.outcomes <-
        matched.validation.sample %>%
        select(-Y, -row_index) %>%
        left_join(validation.sample_outcome, by = c(
          "orig_id" = "rowid", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "em1", "em2", "em3", "em4",
          "treat"
        ))
      matched.inference.sample.outcomes <-
        matched.inference.sample %>%
        select(-Y, -row_index) %>%
        left_join(inference.sample_outcome, by = c(
          "orig_id" = "rowid", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "em1", "em2", "em3", "em4",
          "treat"
        ))

      CCIT_results <- CCIT(
        matched.exploration.sample.outcomes, matched.validation.sample.outcomes, matched.inference.sample.outcomes, lambdas, stopping.rule)

      est.treatment.effects <- CCIT_results$est.treatment.effects
      
      selected.trees <- CCIT_results$selected.trees
      
      tree.list <- CCIT_results$tree.list
      
      selected.tree.size <- CCIT_results$selected.tree.size
    
      true_trt_effects <- true_trt_effect_func(est.treatment.effects[[1]])

      mse <- sapply(est.treatment.effects, function(est.treatment.effect) {
        mean((est.treatment.effect$est - true_trt_effects$eff)^2)
      })

      bias <- sapply(est.treatment.effects, function(est.treatment.effect) {
        mean((est.treatment.effect$est - true_trt_effects$eff))
      })

      # TRUE if tree is exactly correct - CIT paper uses the correct number of splits, but shouldn't order effect the correct number??
      selected.correct.splits <- sapply(selected.trees, function(t) {
        list(row.names(t$splits)) %in% correct_splits
      })

      # TRUE if one of the trees in the generated sequence is correct
      correct.tree.in.sequence <- any(sapply(tree.list, function(t) {
        list(row.names(t$splits)) %in% correct_splits
      }))

      # Number of noise variables selected
      numb.noise <- sapply(selected.trees, function(t) {
        sum(t$frame$var %in% noise.var)
      })

      iter_results <-
        selected.tree.size %>%
        mutate(
          selected.correct.splits = selected.correct.splits,
          correct.tree.in.sequence = correct.tree.in.sequence,
          mse = mse,
          bias = bias,
          numb.noise = numb.noise,
          iter = i
        )
      results <- rbind(results, iter_results)
    }
  }
  return(list(
    results = results, matched.exploration.sample = matched.exploration.sample,
    matched.inference.sample = matched.inference.sample,
    matched.validation.sample = matched.validation.sample,
    inference.sample_covs = inference.sample_covs,
    exploration.sample_covs = exploration.sample_covs,
    val.sample_covs = val.sample_covs
  ))
}

# pps calculation doesn't complete
# pps = 1
# for(i in 1:(nrow(inference.sample)-1)){
#   for(j in (i+1):nrow(inference.sample)){
#     a=b=0
#     if(true_trt_effects$eff[i] == true_trt_effects$eff[j]){a = 1}
#
#     # Treat observations in the 1-observation terminal node be in the same node
#     if(is.na(est.treatment.effects$pred[i]) | is.na(est.treatment.effects$pred[j])) {
#       b = 1
#     } else if (est.treatment.effects$pred[i] == est.treatment.effects$pred[j]) {
#       b = 1
#     }
#     pps = pps - abs(a-b)/choose(nrow(inference.sample), 2)
#   }
# }
