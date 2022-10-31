# Load in relevant scripts
# all R scripts listed below must be lcoated in source_dir
source_dir <- '/Users/ninakatz-christy/Documents/CCIT/Code/CCIT'
dir_out <- '/Users/ninakatz-christy/Documents/CCIT/Code/CCIT/simulated_data_09122022/'
source(paste0(source_dir, '/evaluate.sequence.R'))
source(paste0(source_dir, '/run_simu.R'))
source(paste0(source_dir, '/create.sequence.R'))
source(paste0(source_dir, '/rpart_funcs.R'))
source(paste0(source_dir, '/stratified_GPS_matching.R'))
source(paste0(source_dir, '/split_dataset.R'))
source(paste0(source_dir, '/generate_synthetic_data_covs.R'))
source(paste0(source_dir, '/generate_synthetic_data_outcome.R'))
source(paste0(source_dir, '/CCIT'))
library("devtools")
install_github("fasrc/CausalGPS", ref="master")
library(CausalGPS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(truncnorm)
library(caret)
library(rpart)
library(rpart.plot)
library(scales)
library(data.table)
library(truncnorm)
library(caret)

######################## Run Simulation without generating outcome

# list of different regularization parameters to use. This will only impact how much the full tree is pruned
lambdas <- c(1, 2, 3)

# generate covariates with no confounding, split into subsamples, and perform gps matching
matched <- run_simu(num_exposure_cats = 20, 
                    gps_spec = 1, 
                    sample_size = 20000, n_trials = 0)
save(matched, file=paste0(source_dir, '/simulated_data_09122022/matched.RData'))
load(paste0(source_dir, '/simulated_data/matched.RData'))

# generate covariates with confounding
matched.c <- run_simu(num_exposure_cats = 20, 
                      gps_spec = 2, sample_size = 20000, n_trials = 0)
save(matched.c, file=paste0(source_dir, '/simulated_data_09122022/matched.c.RData'))
load(paste0(source_dir, '/simulated_data/matched.c.RData'))

# plot overall covariate balance
unmatched.no.c.explor <- matched$exploration.sample_covs %>% 
  as.data.table
matched.no.c.explor <- matched$matched.exploration.sample %>% 
  as.data.table
unmatched.c.explor <- matched.c$exploration.sample_covs %>% 
  as.data.table
matched.c.explor <- matched.c$matched.exploration.sample %>% 
  as.data.table

unmatched.no.c.explor <- matched$exploration.sample_covs %>% 
  as.data.table
matched.no.c.explor <- matched$matched.exploration.sample %>% 
  as.data.table
unmatched.c.explor <- matched.c$exploration.sample_covs %>% 
  as.data.table
matched.c.explor <- matched.c$matched.exploration.sample %>% 
  as.data.table

unmatched_cor.no.c <- 
  absolute_corr_fun(unmatched.no.c.explor[,'treat'], unmatched.no.c.explor[,1:6])
matched_cor.no.c <- 
  absolute_corr_fun(matched.no.c.explor[,'w'], matched.no.c.explor[,6:11])
unmatched_cor.c <- 
  absolute_corr_fun(unmatched.c.explor[,'treat'], unmatched.c.explor[,1:6])
matched_cor.c <- 
  absolute_corr_fun(matched.c.explor[,'w'], matched.c.explor[,6:11])

abs_cor <- 
  rbind(data.frame(cov = c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                   #names(synth_data)[3:length(names(synth_data))],
                   cor = unmatched_cor.no.c$absolute_corr, matched = FALSE, confounding = 'No Confounding'), 
        data.frame(cov = c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                   #names(synth_data)[3:length(names(synth_data))],
                   cor = matched_cor.no.c$absolute_corr, matched = TRUE, confounding = 'No Confounding'), 
        data.frame(cov = c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                   #names(synth_data)[3:length(names(synth_data))],
                   cor = unmatched_cor.c$absolute_corr, matched = FALSE, confounding = 'Confounding'), 
        data.frame(cov = c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                   #names(synth_data)[3:length(names(synth_data))],
                   cor = matched_cor.c$absolute_corr, matched = TRUE, confounding = 'Confounding'))
ggplot(abs_cor, aes(x = cov, y = cor, color = matched, group = matched)) + geom_point() + geom_line() + 
  facet_grid(rows = vars(confounding)) + 
  theme(axis.text.x = element_text(angle = 45))
ggsave(paste0(source_dir, '/cov_balance_overall.png'), width = 10)

# plot covariate balance within effect modifier strata
abs_cor <- data.frame()
for(i in c(0,1)) {
  for (j in c(0,1)) {
    for (k in c(0,1)) {
      for (l in c(0,1)) {
        unmatched.no.c.explor <- matched$exploration.sample_covs %>% 
          filter(em1 == i, em2 == j, em3 == k, em4 == l) %>% 
          as.data.table
        matched.no.c.explor <- matched$matched.exploration.sample %>% 
          filter(em1 == i, em2 == j, em3 == k, em4 == l) %>% 
          as.data.table
        unmatched.c.explor <- matched.c$exploration.sample_covs %>% 
          filter(em1 == i, em2 == j, em3 == k, em4 == l) %>% 
          as.data.table
        matched.c.explor <- matched.c$matched.exploration.sample %>% 
          filter(em1 == i, em2 == j, em3 == k, em4 == l) %>% 
          as.data.table
        
        unmatched_cor.no.c <- 
          absolute_corr_fun(unmatched.no.c.explor[,'treat'], unmatched.no.c.explor[,1:6])
        matched_cor.no.c <- 
          absolute_corr_fun(matched.no.c.explor[,'w'], matched.no.c.explor[,6:11])
        unmatched_cor.c <- 
          absolute_corr_fun(unmatched.c.explor[,'treat'], unmatched.c.explor[,1:6])
        matched_cor.c <- 
          absolute_corr_fun(matched.c.explor[,'w'], matched.c.explor[,6:11])
        
        abs_cor <- 
          rbind(abs_cor,
                data.frame(cov = c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                           #names(synth_data)[3:length(names(synth_data))],
                           cor = unmatched_cor.no.c$absolute_corr, matched = FALSE, confounding = 'No Confounding', em1 = i, em2 = j, em3 = k, em4 = l), 
                data.frame(cov = c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                           #names(synth_data)[3:length(names(synth_data))],
                           cor = matched_cor.no.c$absolute_corr, matched = TRUE, confounding = 'No Confounding', em1 = i, em2 = j, em3 = k, em4 = l), 
                data.frame(cov = c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                           #names(synth_data)[3:length(names(synth_data))],
                           cor = unmatched_cor.c$absolute_corr, matched = FALSE, confounding = 'Confounding', em1 = i, em2 = j, em3 = k, em4 = l), 
                data.frame(cov = c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                           #names(synth_data)[3:length(names(synth_data))],
                           cor = matched_cor.c$absolute_corr, matched = TRUE, confounding = 'Confounding', em1 = i, em2 = j, em3 = k, em4 = l))
        
      }}}}
ggplot(abs_cor, aes(x = cov, y = cor, color = matched, group = matched)) + geom_point() + geom_line() + 
  facet_grid(rows = vars(confounding), cols = vars(em1, em2, em3, em4)) + 
  theme(axis.text.x = element_text(angle = 45))
ggsave(paste0(source_dir, '/cov_balance_by_strata.png'), width = 10)

######################## Setting 1: 4 groups, same intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')
beta <- 1

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + beta * (0.5 * em1 - 0.8 * em2))) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, 
                    gps_spec = 2, em_spec = 1,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = matched.c$exploration.sample_covs, val.sample_covs = matched.c$val.sample_covs, 
                    inference.sample_covs = matched.c$inference.sample_covs,
                    matched.exploration.sample = matched.c$matched.exploration.sample, matched.validation.sample = matched.c$matched.validation.sample,
                    matched.inference.sample = matched.c$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting1.RData'))

################## Setting 2: 3 groups, same intercept
correct_splits <- list(c('em2', 'em1'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, 
                    gps_spec = 2, em_spec = 2,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = matched.c$exploration.sample_covs, val.sample_covs = matched.c$val.sample_covs, 
                    inference.sample_covs = matched.c$inference.sample_covs,
                    matched.exploration.sample = matched.c$matched.exploration.sample, matched.validation.sample = matched.c$matched.validation.sample,
                    matched.inference.sample = matched.c$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting2.RData'))

######################## Setting 3: 4 groups, different intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 - 0.8 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, 
                    gps_spec = 1, em_spec = 1,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, heterogenous_intercept = TRUE, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = matched.c$exploration.sample_covs, val.sample_covs = matched.c$val.sample_covs, 
                    inference.sample_covs = matched.c$inference.sample_covs,
                    matched.exploration.sample = matched.c$matched.exploration.sample, matched.validation.sample = matched.c$matched.validation.sample,
                    matched.inference.sample = matched.c$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting3.RData'))

################## Setting 4: 3 groups, different intercept
correct_splits <- list(c('em2', 'em1'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 2,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, heterogenous_intercept = TRUE, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = matched.c$exploration.sample_covs, val.sample_covs = matched.c$val.sample_covs, 
                    inference.sample_covs = matched.c$inference.sample_covs,
                    matched.exploration.sample = matched.c$matched.exploration.sample, matched.validation.sample = matched.c$matched.validation.sample,
                    matched.inference.sample = matched.c$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting4.RData'))


########################## Setting 5: no effect modification
correct_splits <- list(c())
noise.var <- c('em1', 'em2', 'em3', 'em4')
true_trt_effect_func <- function(df) {
  df %>%
    mutate(eff = 10)
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 0,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, heterogenous_intercept = FALSE, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = matched.c$exploration.sample_covs, val.sample_covs = matched.c$val.sample_covs, 
                    inference.sample_covs = matched.c$inference.sample_covs,
                    matched.exploration.sample = matched.c$matched.exploration.sample, matched.validation.sample = matched.c$matched.validation.sample,
                    matched.inference.sample = matched.c$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting5.RData'))

###################### Combining results for different settings
selected.tree.data.combined <- data.frame()
#load(file = paste0(dir_out,'selected.tree.data.with.stopping.setting', 1, '.lambda3000-5000.100trials.resampled.outcome2.RData'))
#selected.tree.data.combined <- selected.tree.data.combined %>% rbind(tree.data %>% cbind(setting = 1))

for (i in 1:5){
  load(file = paste0(dir_out,'selected.tree.data.with.stopping.setting', i, '.RData'))
  selected.tree.data.combined <- selected.tree.data.combined %>% rbind(tree.data %>% cbind(setting = i))
}

setting.labs <- c('Setting 1', 'Setting 2', 'Setting 3', 'Setting 4', 'Setting 5')
names(setting.labs) <- c(1,2,3,4,5)

selected.tree.data.combined %>%
  filter(lambda == 1) %>%
  tidyr::gather("metric", "value", c(bias, mse)) %>%
  ggplot(aes(y = value, group = lambda)) + geom_boxplot() + 
  facet_grid(rows = vars(metric), cols = vars(setting), labeller = labeller(setting = setting.labs), scales = 'free_y') + 
  labs(y=NULL) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file = paste0(dir_out, 'mse.biase.with.stopping.no.c.png'))

selected.tree.data.combined %>%
  group_by(beta) %>%
  summarise_at(vars(selected.correct.splits), mean)

############################ Setting 1 with Different Beta Values
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')

betas <- c(0.1, 0.15,0.2,0.3,1)
for (beta in betas) {
  true_trt_effect_func <- function(df) {
    df %>%
      mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
      mutate(eff = 10*(1 + beta * (0.5 * em1 - 0.8 * em2))) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
  }
  
  results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 1,
                      true_trt_effect_func = true_trt_effect_func, beta = beta,
                      sample_size = 2000, heterogenous_intercept = FALSE, stopping.rule = FALSE, n_trials = 50, 
                      exploration.sample_covs = matched.c$exploration.sample_covs, val.sample_covs = matched.c$val.sample_covs, 
                      inference.sample_covs = matched.c$inference.sample_covs,
                      matched.exploration.sample = matched.c$matched.exploration.sample, matched.validation.sample = matched.c$matched.validation.sample,
                      matched.inference.sample = matched.c$matched.inference.sample)
  tree.data <- results$results
  save(tree.data, file = paste0(dir_out,'selected.tree.data.no.stopping.setting1.lambda3000-5000.100trials.resampled.outcome.beta', beta, '.c.RData'))
}


## results for setting 1 with different beta
selected.tree.data.combined <- data.frame()
for (beta in betas){
  load(file = paste0(dir_out,'selected.tree.data.with.stopping.setting1.lambda3000-5000.100trials.resampled.outcome.beta', beta, '.no.c.RData'))
  selected.tree.data.combined <- selected.tree.data.combined %>% rbind(tree.data %>% cbind(beta = beta, stopping = TRUE))
}
for (beta in betas){
  load(file = paste0(dir_out,'selected.tree.data.no.stopping.setting1.lambda3000-5000.100trials.resampled.outcome.beta', beta, '.c.RData'))
  selected.tree.data.combined <- selected.tree.data.combined %>% rbind(tree.data %>% cbind(beta = beta, stopping = FALSE))
}

beta.labs <- sapply(betas, function(x) {paste0('Beta = ', x)})
names(beta.labs) <- as.character(betas)

selected.tree.data.combined %>%
  filter(lambda == 1) %>%
  tidyr::gather("metric", "value", c(bias, mse)) %>%
  ggplot(aes(y = value, group = lambda)) + geom_boxplot() + 
  facet_grid(rows = vars(metric), cols = vars(beta), labeller = labeller(beta = beta.labs), scales = 'free_y') + 
  labs(y=NULL) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file = paste0(dir_out, 'mse.bias.no.stopping.c.beta.png'))

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  group_by(beta, stopping) %>%
  summarise(percent.correct.tree = mean(selected.correct.splits), n = n())

