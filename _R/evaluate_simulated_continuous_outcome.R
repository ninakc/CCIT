library(dplyr)
library(truncnorm)
library(CausalGPS)
library(caret)
library(rpart)

dir_out = '/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/'

# wrapper function to simulate data, fit tree, and evaluate outcomes
evaluate_simulated_continuous_outcome <- function(lambdas, correct_splits, noise_var, num_exposure_cats, gps_spec = 1, em_spec = 1,
                                                  sample_size = 20000, heterogenous_intercept = FALSE) {
  synth_data <- generate_syn_data_het(sample_size = 20000, outcome_type = 'continuous',
                                       gps_spec = 1, em_spec = em_spec, cova_spec = 1, heterogenous_intercept = heterogenous_intercept,
                                       em_as_confounder = FALSE, outcome_sd = 1, beta = 0.3)

  a_vals <- seq(min(synth_data$treat), max(synth_data$treat), length.out = num_exposure_cats)

  delta_n <- a_vals[2] - a_vals[1]

  synth_data <-
    synth_data %>%
    mutate(treat_level = cut(treat, breaks = a_vals))

  synth_data <-
    split_dataset(data = synth_data, exposure_bins = synth_data$treat_level)

  exploration_sample <-
    synth_data %>% filter(subsample == 'exploration')

  val_sample <-
    synth_data %>% filter(subsample == 'validation')

  inference_sample <-
    synth_data %>% filter(subsample == 'inference')

  matched_exploration_sample <-
    stratified_GPS_matching(exploration_sample, delta_n) %>%
    mutate(id = row_number()) %>%
    #rename(w = treat) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)

  matched_validation_sample <-
    stratified_GPS_matching(val_sample, delta_n) %>%#
    mutate(id = row_number()) %>%
    #rename(w = treat) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)

  overall_effect1 <- lm(Y ~ w + 1, data = matched_exploration_sample)$coefficients['w']

  parms_used1 <- list(w  = matched_exploration_sample$w,
                     # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                     Y   = matched_exploration_sample$Y,
                     overall_effect = overall_effect1)

  lists1 <- create_sequence(matched_exploration_sample, ulist_used, parms_used1)

  tree_list1 <- lists1[[1]]
  g_h_list1 <- lists1[[2]]
  tree_sizes1 <- sapply(tree_list1, function(t) {ifelse(is.null(t$splits), 0,nrow(t$splits))})

  # tree size selected via the validation sample (no cross-fitting for now)
  selected_tree_sizes <- data.frame()
  for (b in 1:1) {
    boot_sample_index <- sample(1:nrow(matched_validation_sample), replace = TRUE)
    boot_matched_validation_sample <- matched_validation_sample[boot_sample_index, ]

    complex_vals1 <- evaluate_sequence(tree_list1, boot_matched_validation_sample, matched_exploration_sample, lambdas)

    selected_tree_size1 <- complex_vals1 %>% group_by(lambda) %>%
      filter(complex_val == max(complex_val)) %>%
      select(tree_size) %>%
      ungroup()
    selected_tree_sizes <- rbind(selected_tree_sizes, selected_tree_size1)

  }

  selected_tree_size <- selected_tree_sizes %>%
    group_by(lambda) %>%
    summarise_at(vars(tree_size), function(s) {round(mean(s))}) %>%
    ungroup()

  print(selected_tree_size$tree_size)
  print(tree_sizes1)
  selected_trees <- lapply(selected_tree_size$tree_size, function(s) {tree_list1[[which(tree_sizes1 == s)]]})
  names(selected_trees) <- lambdas
  est_treatment_effects <-  lapply(selected_trees,
                                   function(t) {mutate(inference_sample, pred = predict(t, newdata = inference_sample))})


  true_trt_effects <- true_trt_effect_func(inference_sample)

  mse <- sapply(est_treatment_effects, function(est_treatment_effects) { mean((est_treatment_effects$pred - true_trt_effects$eff)^2)})

  bias <- sapply(est_treatment_effects, function(est_treatment_effects) { mean((est_treatment_effects$pred - true_trt_effects$eff))})

  # TRUE if tree is exactly correct - CIT paper uses the correct number of splits, but shouldn't order effect the correct number??
  selected_correct_splits <- sapply(selected_trees, function(t) {list(row.names(t$splits)) %in% correct_splits})


  # TRUE if one of the trees in the generated sequence is correct
  correct_tree_in_sequence <- any(sapply(tree_list1, function(t) {list(row.names(t$splits)) %in% correct_splits}))

  # Number of noise variables selected
  numb_noise <- sapply(selected_trees, function(t) {sum(t$frame$var %in% noise_var)})

  return(selected_tree_size %>%
           mutate(selected_correct_splits = selected_correct_splits,
                  correct_tree_in_sequence = correct_tree_in_sequence,
                  mse = mse,
                  bias = bias,
                  numb_noise = numb_noise))
}

correct_splits <- list()
noise_var <- c('em1', 'em2', 'em3', 'em4')
true_trt_effect_func <- function(df) {
  df %>%
    mutate(eff = 10)
}


######################## Setting 1: 4 groups, same intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise_var <- c('em3', 'em4')
beta <- 0.3

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + beta * (0.5 * em1 - 0.8 * em2))) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)

}

#lambdas <- c(1.5, 1.75, 2, 3, 4)


selected_tree_data <- data.frame()
  for (i in 1:50) {

    ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise_var = noise_var,
                                                 num_exposure_cats = 20, em_spec = 1)
    selected_tree_data <- rbind(selected_tree_data, as.data.frame(ret))
  }



save(selected_tree_data, file = paste0(dir_out,'selected_tree_data.with.stopping.setting6.lambda3000-5000.50trials.RData'))

load(paste0(dir_out,'selected_tree_data.with.stopping.setting1.beta0.3.lambda3000-5000.50trials.RData'))

selected_tree_data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() +
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))


################## Setting 2: 3 groups, same intercept
correct_splits <- list(c('em2', 'em1'))
noise_var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)

}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected_tree_data <- data.frame()
for (i in 1:50) {

  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise_var = noise_var,
                                               num_exposure_cats = 20, em_spec = 2)
  selected_tree_data <- rbind(selected_tree_data, as.data.frame(ret))
}

save(selected_tree_data, file = paste0(dir_out,'selected_tree_data.with.stopping.setting2.sd10.lambda3000-5000.50trials.RData'))

load(paste0(dir_out,'no.stopping.RData'))


selected_tree_data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() +
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))



######################## Setting 3: 4 groups, different intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise_var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 - 0.8 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)

}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected_tree_data <- data.frame()
for (i in 1:50) {

  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise_var = noise_var,
                                               num_exposure_cats = 20, heterogenous_intercept = TRUE)
  selected_tree_data <- rbind(selected_tree_data, as.data.frame(ret))
}

save(selected_tree_data, file = paste0(dir_out,'selected_tree_data.with.stopping.setting3.lambda3000-5000.50trials.RData'))

load(paste0(dir_out,'no.stopping.RData'))

selected_tree_data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() +
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.with.stopping50trials.png'))


selected_tree_data %>%
  ggplot(aes(x = lambda, y = bias, group = lambda)) + geom_boxplot() +
  scale_x_continuous(breaks = lambdas)

selected_tree_data %>% summarise_at(vars(numb_noise, selected_correct_splits, bias, mse), mean)

################## Setting 4: 3 groups, different intercept
correct_splits <- list(c('em2', 'em1'))
noise_var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)

}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected_tree_data <- data.frame()
for (i in 1:50) {

  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise_var = noise_var,
                                               num_exposure_cats = 20, em_spec = 2, heterogenous_intercept = TRUE)
  selected_tree_data <- rbind(selected_tree_data, as.data.frame(ret))
}

save(selected_tree_data, file = paste0(dir_out,'selected_tree_data.with.stopping.setting4.lambda3000-5000.50trials.RData'))

selected_tree_data_combined <- data.frame()
for (i in 1:5){
  load(file = paste0(dir_out,'selected_tree_data.with.stopping.setting', i, '.lambda3000-5000.50trials.RData'))
  selected_tree_data_combined <- selected_tree_data_combined %>% rbind(selected_tree_data %>% cbind(setting = i))
}

setting.labs <- c('Setting 1', 'Setting 2', 'Setting 3', 'Setting 4', 'Setting 5')
names(setting.labs) <- c(1,2,3,4,5)

selected_tree_data_combined %>%
  filter(lambda == 3000) %>%
  ggplot(aes(y = mse, group = lambda)) + geom_boxplot() +
  facet_wrap(vars(setting), nrow = 1, labeller = labeller(setting = setting.labs)) +
  scale_x_continuous(breaks = lambdas)
ggsave(file = paste0(dir_out, 'mse.with.stopping2.png'))

selected_tree_data_combined %>%
  filter(lambda == 3000) %>%
  ggplot(aes(y = bias, group = lambda)) + geom_boxplot() +
  facet_wrap(vars(setting), nrow = 1, labeller = labeller(setting = setting.labs)) +
  scale_x_continuous(breaks = lambdas)
ggsave(file = paste0(dir_out, 'bias.with.stopping2.png'))

selected_tree_data_combined %>% group_by(setting) %>%
  summarise_at(vars(numb_noise, selected_correct_splits, bias, mse), mean)







load(paste0(dir_out,'no.stopping.RData'))


selected_tree_data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() +
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))







selected_tree_data %>%
  ggplot(aes(x = lambda, fill = selected_correct_splits, y = ..count../100)) + geom_histogram(position = 'dodge') +
  scale_x_continuous(breaks = lambdas) +
  labs(y = 'Proportion of Simulation trials')
ggsave(paste0(dir_out, 'Prop.correct.trees.no.stopping.png'))

selected_tree_data %>%
  ggplot(aes(x = lambda, y = numb_noise, group = lambda)) + geom_boxplot() +
  scale_x_continuous(breaks = lambdas) +
  labs(y = 'Proportion of Simulation trials')
ggsave(paste0(dir_out, 'Prop.correct.trees.no.stopping.png'))

selected_tree_data %>%
  mutate('Pruned to Correct Size' = (floor(tree_size) == 3)) %>%
  ggplot(aes(x = lambda, fill = `Pruned to Correct Size`, y = ..count../100)) + geom_histogram(position = 'dodge') +
  scale_x_continuous(breaks = lambdas) +
  labs(y = 'Proportion of Simulation trials')
ggsave(paste0(dir_out, 'Pruning_accuracy.no.stopping.png'))
ggsave(paste0(dir_out, 'Pruning_accuracy_bootstrap_n100000.png'))


selected_tree_data %>%
  ggplot(aes(x = lambda, y = floor(tree_size), group = lambda)) + geom_boxplot() +
  geom_hline(yintercept = 3) +
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruned_tree_size_crossfit.png'))
ggsave(paste0(dir_out, 'Pruned_tree_size_bootstrapn100000.png'))



selected_tree_data %>%
  ggplot(aes(x = lambda, y = selected_tree_size, group = lambda)) + geom_boxplot() +
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruned_tree_size.png'))


selected_tree_data %>%
  mutate('Pruned to Correct Size' = as.logical(selected_correct_splits)) %>%
  ggplot(aes(x = lambda, fill = `Pruned to Correct Size`)) + geom_histogram(position = 'dodge') +
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruning_accuracy.png'))


selected_tree_data %>%
  group_by(lambda) %>%
  summarise(prop_correct = sum(selected_correct_splits)/n(),
            avg_tree_size = mean(selected_tree_size),
            prop_correct_in_sequence = sum(correct_tree_in_sequence)/n(),
            prop_correct_size = sum(selected_tree_size == 3)/n())

ret <- evaluate_simulated_continuous_outcome(lambda = 1, correct_splits = correct_splits, num_exposure_cats = 20)

# complex_vals1 <- data.frame(complex_val = evaluate_sequence(tree_list1, matched_validation_sample, matched_exploration_sample, lambda),
#            pruning_step = 1:length(tree_list1),
#            tree_size = tree_sizes1)
#
# selected_tree_size1 <- complex_vals1[which(complex_vals1$complex_val == max(complex_vals1$complex_val)), 'tree_size']
#
# overall_effect2 <- lm(Y ~ w + 1, data = matched_validation_sample)$coefficients['w']
#
# parms.used2 <- list(w  = matched_validation_sample$w,
#                     # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
#                     Y   = matched_validation_sample$Y,
#                     overall_effect = overall_effect2)
#
# lists2 <- create_sequence(matched_validation_sample, ulist_used, parms.used2)
#
# tree.list2 <- lists2[[1]]
# g_h.list2 <- lists2[[2]]
#
# tree.sizes2 <- sapply(tree.list2, function(t) {nrow(t$splits)})
# complex.vals2 <- data.frame(complex_val = evaluate_sequence(tree.list2, matched_exploration_sample, matched_validation_sample, lambda),
#                             pruning_step = 1:length(tree.list2),
#                             tree_size = tree.sizes2)
#
# selected.tree.size2 <- complex.vals2[which(complex.vals2$complex_val == max(complex.vals2$complex_val)), 'tree_size']

#selected_tree_size <- mean(c(selected_tree_size1, selected.tree.size2))

# pps calculation doesn't complete
# pps = 1
# for(i in 1:(nrow(inference_sample)-1)){
#   for(j in (i+1):nrow(inference_sample)){
#     a=b=0
#     if(true_trt_effects$eff[i] == true_trt_effects$eff[j]){a = 1}
#
#     # Treat observations in the 1-observation terminal node be in the same node
#     if(is.na(est_treatment_effects$pred[i]) | is.na(est_treatment_effects$pred[j])) {
#       b = 1
#     } else if (est_treatment_effects$pred[i] == est_treatment_effects$pred[j]) {
#       b = 1
#     }
#     pps = pps - abs(a-b)/choose(nrow(inference_sample), 2)
#   }
# }

