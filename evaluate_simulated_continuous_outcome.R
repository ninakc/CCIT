library(dplyr)
library(truncnorm)
library(CausalGPS)
library(caret)
library(rpart)

dir_out = '/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/'

# wrapper function to simulate data, fit tree, and evaluate outcomes
evaluate_simulated_continuous_outcome <- function(lambdas, correct_splits, noise.var, num_exposure_cats, gps_spec = 1, em_spec = 1, 
                                                  sample_size = 20000, heterogenous_intercept = FALSE) {
  synth_data <- generate_syn_data_het(sample_size = 20000, outcome_type = 'continuous',
                                       gps_spec = 1, em_spec = em_spec, cova_spec = 1, heterogenous_intercept = heterogenous_intercept, 
                                       em_as_confounder = FALSE, outcome_sd = 1, beta = 0.3)
  
  a.vals <- seq(min(synth_data$treat), max(synth_data$treat), length.out = num_exposure_cats)
  
  delta_n <- a.vals[2] - a.vals[1]
  
  synth_data <- 
    synth_data %>%
    mutate(treat_level = cut(treat, breaks = a.vals))

  synth_data <- 
    split_dataset(data = synth_data, exposure_bins = synth_data$treat_level)
  
  exploration.sample <- 
    synth_data %>% filter(subsample == 'exploration')
  
  val.sample <- 
    synth_data %>% filter(subsample == 'validation')
  
  inference.sample <- 
    synth_data %>% filter(subsample == 'inference')
  
  matched.exploration.sample <- 
    stratified_GPS_matching(exploration.sample, delta_n) %>%
    mutate(id = row_number()) %>%
    #rename(w = treat) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
  matched.validation.sample <- 
    stratified_GPS_matching(val.sample, delta_n) %>%#
    mutate(id = row_number()) %>%
    #rename(w = treat) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
    
  overall_effect1 <- lm(Y ~ w + 1, data = matched.exploration.sample)$coefficients['w']
    
  parms.used1 <- list(w  = matched.exploration.sample$w,
                     # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                     Y   = matched.exploration.sample$Y,
                     overall_effect = overall_effect1)   
  
  lists1 <- create.sequence(matched.exploration.sample, ulist.used, parms.used1)

  tree.list1 <- lists1[[1]]
  g.h.list1 <- lists1[[2]]
  tree.sizes1 <- sapply(tree.list1, function(t) {ifelse(is.null(t$splits), 0,nrow(t$splits))})
  
  # tree size selected via the validation sample (no cross-fitting for now)
  selected.tree.sizes <- data.frame()
  for (b in 1:1) {
    boot_sample_index <- sample(1:nrow(matched.validation.sample), replace = TRUE)
    boot.matched.validation.sample <- matched.validation.sample[boot_sample_index, ]
    
    complex.vals1 <- evaluate.sequence(tree.list1, boot.matched.validation.sample, matched.exploration.sample, lambdas)
    
    selected.tree.size1 <- complex.vals1 %>% group_by(lambda) %>%
      filter(complex.val == max(complex.val)) %>%
      select(tree.size) %>%
      ungroup()
    selected.tree.sizes <- rbind(selected.tree.sizes, selected.tree.size1)
    
  }
  
  selected.tree.size <- selected.tree.sizes %>%
    group_by(lambda) %>%
    summarise_at(vars(tree.size), function(s) {round(mean(s))}) %>%
    ungroup()
  
  print(selected.tree.size$tree.size)
  print(tree.sizes1)
  selected.trees <- lapply(selected.tree.size$tree.size, function(s) {tree.list1[[which(tree.sizes1 == s)]]})
  names(selected.trees) <- lambdas
  est.treatment.effects <-  lapply(selected.trees, 
                                   function(t) {mutate(inference.sample, pred = predict(t, newdata = inference.sample))})
    
  
  true_trt_effects <- true_trt_effect_func(inference.sample)
  
  mse <- sapply(est.treatment.effects, function(est.treatment.effect) { mean((est.treatment.effect$pred - true_trt_effects$eff)^2)})
  
  bias <- sapply(est.treatment.effects, function(est.treatment.effect) { mean((est.treatment.effect$pred - true_trt_effects$eff))})
  
  # TRUE if tree is exactly correct - CIT paper uses the correct number of splits, but shouldn't order effect the correct number??
  selected.correct.splits <- sapply(selected.trees, function(t) {list(row.names(t$splits)) %in% correct_splits})
    
  
  # TRUE if one of the trees in the generated sequence is correct
  correct.tree.in.sequence <- any(sapply(tree.list1, function(t) {list(row.names(t$splits)) %in% correct_splits}))
  
  # Number of noise variables selected
  numb.noise <- sapply(selected.trees, function(t) {sum(t$frame$var %in% noise.var)})
  
  return(selected.tree.size %>%
           mutate(selected.correct.splits = selected.correct.splits,
                  correct.tree.in.sequence = correct.tree.in.sequence,
                  mse = mse, 
                  bias = bias,
                  numb.noise = numb.noise))
}

correct_splits <- list()
noise.var <- c('em1', 'em2', 'em3', 'em4')
true_trt_effect_func <- function(df) {
  df %>%
    mutate(eff = 10)
}


######################## Setting 1: 4 groups, same intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')
beta <- 0.3

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + beta * (0.5 * em1 - 0.8 * em2))) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}
  
#lambdas <- c(1.5, 1.75, 2, 3, 4)


selected.tree.data <- data.frame()
  for (i in 1:50) {
    
    ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise.var = noise.var, 
                                                 num_exposure_cats = 20, em_spec = 1)
    selected.tree.data <- rbind(selected.tree.data, as.data.frame(ret))
  }
  


save(selected.tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting6.lambda3000-5000.50trials.RData'))

load(paste0(dir_out,'selected.tree.data.with.stopping.setting1.beta0.3.lambda3000-5000.50trials.RData'))

selected.tree.data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))


################## Setting 2: 3 groups, same intercept
correct_splits <- list(c('em2', 'em1'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected.tree.data <- data.frame()
for (i in 1:50) {
  
  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise.var = noise.var, 
                                               num_exposure_cats = 20, em_spec = 2)
  selected.tree.data <- rbind(selected.tree.data, as.data.frame(ret))
}

save(selected.tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting2.sd10.lambda3000-5000.50trials.RData'))

load(paste0(dir_out,'no.stopping.RData'))


selected.tree.data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))



######################## Setting 3: 4 groups, different intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 - 0.8 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected.tree.data <- data.frame()
for (i in 1:50) {
  
  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise.var = noise.var, 
                                               num_exposure_cats = 20, heterogenous_intercept = TRUE)
  selected.tree.data <- rbind(selected.tree.data, as.data.frame(ret))
}

save(selected.tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting3.lambda3000-5000.50trials.RData'))

load(paste0(dir_out,'no.stopping.RData'))

selected.tree.data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.with.stopping50trials.png'))


selected.tree.data %>%
  ggplot(aes(x = lambda, y = bias, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)

selected.tree.data %>% summarise_at(vars(numb.noise, selected.correct.splits, bias, mse), mean)

################## Setting 4: 3 groups, different intercept
correct_splits <- list(c('em2', 'em1'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected.tree.data <- data.frame()
for (i in 1:50) {
  
  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise.var = noise.var, 
                                               num_exposure_cats = 20, em_spec = 2, heterogenous_intercept = TRUE)
  selected.tree.data <- rbind(selected.tree.data, as.data.frame(ret))
}

save(selected.tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting4.lambda3000-5000.50trials.RData'))

selected.tree.data.combined <- data.frame()
for (i in 1:5){
  load(file = paste0(dir_out,'selected.tree.data.with.stopping.setting', i, '.lambda3000-5000.50trials.RData'))
  selected.tree.data.combined <- selected.tree.data.combined %>% rbind(selected.tree.data %>% cbind(setting = i))
}

setting.labs <- c('Setting 1', 'Setting 2', 'Setting 3', 'Setting 4', 'Setting 5')
names(setting.labs) <- c(1,2,3,4,5)

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  ggplot(aes(y = mse, group = lambda)) + geom_boxplot() + 
  facet_wrap(vars(setting), nrow = 1, labeller = labeller(setting = setting.labs)) + 
  scale_x_continuous(breaks = lambdas)
ggsave(file = paste0(dir_out, 'mse.with.stopping2.png'))

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  ggplot(aes(y = bias, group = lambda)) + geom_boxplot() + 
  facet_wrap(vars(setting), nrow = 1, labeller = labeller(setting = setting.labs)) + 
  scale_x_continuous(breaks = lambdas)
ggsave(file = paste0(dir_out, 'bias.with.stopping2.png'))

selected.tree.data.combined %>% group_by(setting) %>%
  summarise_at(vars(numb.noise, selected.correct.splits, bias, mse), mean)







load(paste0(dir_out,'no.stopping.RData'))


selected.tree.data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))







selected.tree.data %>%
  ggplot(aes(x = lambda, fill = selected.correct.splits, y = ..count../100)) + geom_histogram(position = 'dodge') + 
  scale_x_continuous(breaks = lambdas) + 
  labs(y = 'Proportion of Simulation trials') 
ggsave(paste0(dir_out, 'Prop.correct.trees.no.stopping.png'))

selected.tree.data %>%
  ggplot(aes(x = lambda, y = numb.noise, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas) + 
  labs(y = 'Proportion of Simulation trials')
ggsave(paste0(dir_out, 'Prop.correct.trees.no.stopping.png'))

selected.tree.data %>%
  mutate('Pruned to Correct Size' = (floor(tree.size) == 3)) %>%
  ggplot(aes(x = lambda, fill = `Pruned to Correct Size`, y = ..count../100)) + geom_histogram(position = 'dodge') + 
  scale_x_continuous(breaks = lambdas) + 
  labs(y = 'Proportion of Simulation trials')
ggsave(paste0(dir_out, 'Pruning_accuracy.no.stopping.png'))
ggsave(paste0(dir_out, 'Pruning_accuracy_bootstrap_n100000.png'))


selected.tree.data %>%
  ggplot(aes(x = lambda, y = floor(tree.size), group = lambda)) + geom_boxplot() + 
  geom_hline(yintercept = 3) + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruned_tree_size_crossfit.png'))
ggsave(paste0(dir_out, 'Pruned_tree_size_bootstrapn100000.png'))



selected.tree.data %>%
  ggplot(aes(x = lambda, y = selected.tree.size, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruned_tree_size.png'))


selected.tree.data %>%
  mutate('Pruned to Correct Size' = as.logical(selected.correct.splits)) %>%
  ggplot(aes(x = lambda, fill = `Pruned to Correct Size`)) + geom_histogram(position = 'dodge') + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruning_accuracy.png'))


selected.tree.data %>%
  group_by(lambda) %>%
  summarise(prop_correct = sum(selected.correct.splits)/n(), 
            avg_tree_size = mean(selected.tree.size), 
            prop_correct_in_sequence = sum(correct.tree.in.sequence)/n(),
            prop_correct_size = sum(selected.tree.size == 3)/n())

ret <- evaluate_simulated_continuous_outcome(lambda = 1, correct_splits = correct_splits, num_exposure_cats = 20)

# complex.vals1 <- data.frame(complex.val = evaluate.sequence(tree.list1, matched.validation.sample, matched.exploration.sample, lambda), 
#            pruning_step = 1:length(tree.list1),
#            tree.size = tree.sizes1)
# 
# selected.tree.size1 <- complex.vals1[which(complex.vals1$complex.val == max(complex.vals1$complex.val)), 'tree.size']
# 
# overall_effect2 <- lm(Y ~ w + 1, data = matched.validation.sample)$coefficients['w']
# 
# parms.used2 <- list(w  = matched.validation.sample$w,
#                     # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
#                     Y   = matched.validation.sample$Y,
#                     overall_effect = overall_effect2)   
# 
# lists2 <- create.sequence(matched.validation.sample, ulist.used, parms.used2)
# 
# tree.list2 <- lists2[[1]]
# g.h.list2 <- lists2[[2]]
# 
# tree.sizes2 <- sapply(tree.list2, function(t) {nrow(t$splits)})
# complex.vals2 <- data.frame(complex.val = evaluate.sequence(tree.list2, matched.exploration.sample, matched.validation.sample, lambda), 
#                             pruning_step = 1:length(tree.list2),
#                             tree.size = tree.sizes2)
# 
# selected.tree.size2 <- complex.vals2[which(complex.vals2$complex.val == max(complex.vals2$complex.val)), 'tree.size']

#selected.tree.size <- mean(c(selected.tree.size1, selected.tree.size2))

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

