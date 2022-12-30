
# wrapper function to simulate data, fit tree, and evaluate outcomes
run_simu <- function(gps_spec = 1,
                     num_exposure_cats,
                     sample_size = 20000,
                     em_spec = 1,
                     heterogenous_intercept = FALSE,
                     beta = NULL,
                     outcome_sd = 1,
                     correct_splits = NULL,
                     true_trt_effect_func = NULL,
                     noise_var,
                     n_trials,
                     lambdas,
                     stopping_rule,
                     exploration_sample_covs = NULL,
                     inference_sample_covs = NULL,
                     val_sample_covs = NULL,
                     matched_exploration_sample = NULL,
                     matched_validation_sample = NULL,
                     matched_inference_sample = NULL) {

  # if the matched covariate dataset (everything but the outcome) is not passed in, generate one
  if (is.null(matched_exploration_sample) | is.null(matched_validation_sample) | is.null(matched_inference_sample) |
    is.null(exploration_sample_covs) | is.null(val_sample_covs) | is.null(inference_sample_covs)) {
    # generate covariate data
    synth_data_covs <- generate_syn_data_covs(sample_size = sample_size, gps_spec = gps_spec)

    # discretize treatment values for gps matching
    a_vals <- seq(min(synth_data_covs$treat), max(synth_data_covs$treat), length.out = num_exposure_cats)
    delta_n <- a_vals[2] - a_vals[1]

    synth_data_covs <-
      synth_data_covs %>%
      mutate(treat_level = cut(treat, breaks = a_vals)) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)

    # split data into subsamples, stratifying on exposure level bins
    synth_data_covs <-
      split_dataset(data = synth_data_covs, exposure_bins = synth_data_covs$treat_level)

    exploration_sample_covs <-
      synth_data_covs %>%
      filter(subsample == "exploration") %>%
      tibble::rowid_to_column("orig_id")

    val_sample_covs <-
      synth_data_covs %>%
      filter(subsample == "validation") %>%
      tibble::rowid_to_column("orig_id")

    inference_sample_covs <-
      synth_data_covs %>%
      filter(subsample == "inference") %>%
      tibble::rowid_to_column("orig_id")

    # GPS matching within subsamples and effect modifier groups
    matched_exploration_sample <-
      stratified_GPS_matching(exploration_sample_covs, delta_n,
        bin_seq = a_vals, exposure_name = "treat",
        confounders_names = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
        em_names = c("em1", "em2", "em3", "em4"),
        outcome_name = NA
      ) %>%
      mutate(id = row_number()) %>%
      # rename(w = treat) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    # print('here')
    # print(names(matched_exploration_sample))

    matched_validation_sample <-
      stratified_GPS_matching(val_sample_covs, delta_n,
        bin_seq = a_vals, exposure_name = "treat",
        confounders_names = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
        em_names = c("em1", "em2", "em3", "em4"),
        outcome_name = NA
      ) %>%
      mutate(id = row_number()) %>%
      # rename(w = treat) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)

    matched_inference_sample <-
      stratified_GPS_matching(inference_sample_covs, delta_n,
        bin_seq = a_vals, exposure_name = "treat",
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
      exploration_sample_outcome <-
        do.call(
          generate_syn_data_outcome,
          c(as.list(exploration_sample_covs %>% select(-treat_level, -subsample, -orig_id)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd, heterogenous_intercept = heterogenous_intercept)
        ) %>%
        tibble::rowid_to_column()

      validation_sample_outcome <-
        do.call(
          generate_syn_data_outcome,
          c(as.list(val_sample_covs %>% select(-treat_level, -subsample, -orig_id)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd, heterogenous_intercept = heterogenous_intercept)
        ) %>%
        tibble::rowid_to_column()

      inference_sample_outcome <-
        do.call(
          generate_syn_data_outcome,
          c(as.list(inference_sample_covs %>% select(-treat_level, -subsample, -orig_id)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd, heterogenous_intercept = heterogenous_intercept)
        ) %>%
        tibble::rowid_to_column()

      # add outcome data into matched samples
      matched_exploration_sample_outcomes <-
        matched_exploration_sample %>%
        select(-Y, -row_index) %>%
        left_join(exploration_sample_outcome, by = c(
          "orig_id" = "rowid", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "em1", "em2", "em3", "em4",
          "treat"
        ))

      matched_validation_sample_outcomes <-
        matched_validation_sample %>%
        select(-Y, -row_index) %>%
        left_join(validation_sample_outcome, by = c(
          "orig_id" = "rowid", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "em1", "em2", "em3", "em4",
          "treat"
        ))
      matched_inference_sample_outcomes <-
        matched_inference_sample %>%
        select(-Y, -row_index) %>%
        left_join(inference_sample_outcome, by = c(
          "orig_id" = "rowid", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "em1", "em2", "em3", "em4",
          "treat"
        ))

      CCIT_results <- CCIT(
        matched_exploration_sample_outcomes, matched_validation_sample_outcomes, matched_inference_sample_outcomes, lambdas, stopping_rule)

      est_treatment_effects <- CCIT_results$est_treatment_effects

      selected_trees <- CCIT_results$selected_trees

      tree_list <- CCIT_results$tree_list

      selected_tree_size <- CCIT_results$selected_tree_size

      true_trt_effects <- true_trt_effect_func(est_treatment_effects[[1]])

      mse <- sapply(est_treatment_effects, function(est_treatment_effects) {
        mean((est_treatment_effects$est - true_trt_effects$eff)^2)
      })

      bias <- sapply(est_treatment_effects, function(est_treatment_effects) {
        mean((est_treatment_effects$est - true_trt_effects$eff))
      })

      # TRUE if tree is exactly correct - CIT paper uses the correct number of splits, but shouldn't order effect the correct number??
      selected_correct_splits <- sapply(selected_trees, function(t) {
        list(row.names(t$splits)) %in% correct_splits
      })

      # TRUE if one of the trees in the generated sequence is correct
      correct_tree_in_sequence <- any(sapply(tree_list, function(t) {
        list(row.names(t$splits)) %in% correct_splits
      }))

      # Number of noise variables selected
      numb_noise <- sapply(selected_trees, function(t) {
        sum(t$frame$var %in% noise_var)
      })

      iter_results <-
        selected_tree_size %>%
        mutate(
          selected_correct_splits = selected_correct_splits,
          correct_tree_in_sequence = correct_tree_in_sequence,
          mse = mse,
          bias = bias,
          numb_noise = numb_noise,
          iter = i
        )
      results <- rbind(results, iter_results)
    }
  }
  return(list(
    results = results, matched_exploration_sample = matched_exploration_sample,
    matched_inference_sample = matched_inference_sample,
    matched_validation_sample = matched_validation_sample,
    inference_sample_covs = inference_sample_covs,
    exploration_sample_covs = exploration_sample_covs,
    val_sample_covs = val_sample_covs
  ))
}
