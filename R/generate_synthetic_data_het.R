#' @title
#' Generate Synthetic Data for CausalGPS Package - Modified by Nina
#'
#' @description
#' Generates synthetic data set based on different GPS models and covariates.
#'
#' @param sample_size Number of data samples.
#' 
#' @param outcome_type continuous or binary outcome
#' 
#' @param outcome_sd Standard deviation used to generate the outcome in the synthetic data set
#' 
#' @param gps_spec A numerical value (1-7) that indicates the GPS model
#' used to generate synthetic data. See the code for more details.

#' @param cova_spec A numerical value (1-2) to modify the covariates. If 
#' cova_spec = 1, covariates are mutually independent. See the
#' code for more details.
#' 
#' @param heterogenous_intercept Allow intercepts to differ by effect modifier level
#' 
#' @param em_spec A numerical value (1-2) that indicates the effect modification model.
#' em_spec = 1 -> separate main effects of effect modifiers on treatment effect
#' (i.e. no interaction between effect modifiers)
#' em_spec = 2 -> interaction between effect modifiers
#' 
#' @param em_as_confounder Effect modifier 1 has linear association with treatment 
#' (i.e. is a confounder)
#'
#' @return
#' \code{synthetic_data}: The function returns a data.frame saved the
#'  constructed synthetic data.
#'
#' @export
#'
#' @examples
#'
#' set.seed(298)
#' s_data <- generate_syn_data(sample_size=100,
#'                                   outcome_sd = 10, gps_spec = 1,
#'                                   cova_spec = 1)
#'
#'Required Libraries: truncnorm, dplyr
generate_syn_data_het <- function(sample_size=10000, outcome_type = 'continuous', outcome_sd = 1,
                              gps_spec = 1, cova_spec = 1, em_spec = 1, heterogenous_intercept = FALSE, em_as_confounder = FALSE, 
                              beta = 1) {

  if (sample_size < 0 || !is.numeric(sample_size)){
    stop("'sample_size' should be a positive ineteger numer.")
   }

  #TODO: Check other input arguments.

  size <- sample_size

  #pre-treatment variables (confounders)
  cf  <- MASS::mvrnorm(n = size,
                       mu = c(0,0,0,0),
                       Sigma = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
                       ncol=4))
  colnames(cf) <- c('cf1', 'cf2', 'cf3', 'cf4')

  cf5 <- sample(c((-2):2), size, replace = TRUE)
  cf6 <- stats::runif(size, min=-3, max=3)

  # true effect modifiers
  em1 <- rbinom(size, 1, 0.5)
  em2 <- rbinom(size, 1, 0.4)
  
  # null effect modifiers
  em3 <- rbinom(size, 1, 0.3)
  em4 <- rbinom(size, 1, 0.2)
  
  if (gps_spec == 1) {
    # modified Xiao's propensity function to ensure treatment >= 0. Note that mu is not the expected propensity.
    # mu <- (- 0.8 + 0.1 * cf[ ,1] + 0.1 * cf[ ,2] - 0.1 * cf[ ,3]
    #        + 0.2 * cf[ ,4] + 0.1 * cf5 + 0.1 * cf6) * 9 + 17
    mu <- 3
    treat <- rtruncnorm(size,a=0,b=Inf,mean=mu,sd=5)
    # currently not using other gps_spec values
  } else if (gps_spec == 2) {
    treat <- ((- 0.8 + 0.1 * cf[ ,1] + 0.1 * cf[ ,2] - 0.1 * cf[ ,3]
              + 0.2 * cf[ ,4] + 0.1 * cf5 + 0.1 * cf6) * 15 + 22 + stats::rt(size,2))

    treat[which(treat < (-5))] <- (-5)
    treat[which(treat > (25))] <- (25)

  } else if (gps_spec == 3) {

    treat <- ((- 0.8 + 0.1 * cf[ , 1] + 0.1 * cf[ , 2]- 0.1 *cf[ ,3] + 0.2 * cf [ , 4]
               + 0.1 * cf5 + 0.1 * cf6) * 9
               + 1.5 * cf[ , 3] ^ 2 + stats::rnorm(size, mean = 0, 5) + 15)

  } else if (gps_spec == 4) {

    treat <- (49 * exp((-0.8 + 0.1 * cf[ ,1] + 0.1 * cf[ , 2] - 0.1 * cf[ , 3]
            + 0.2 * cf[ , 4] + 0.1 * cf5 + 0.1 * cf6))
            / (1 + exp((-0.8 + 0.1 * cf[,1] + 0.1 * cf[ , 2] - 0.1 * cf[ , 3]
            + 0.2 * cf[ , 4] + 0.1 * cf5 + 0.1 * cf6))) - 6 + stats::rnorm(size, sd=5))

  } else if (gps_spec == 5) {

    treat <- (42 / (1 + exp((-0.8 + 0.1 * cf[ , 1] + 0.1 * cf[ , 2]- 0.1 * cf[ , 3]
           + 0.2 * cf[,4] + 0.1 * cf5 + 0.1 * cf6))) - 18 + stats::rnorm(size,sd=5))

  } else if (gps_spec == 6) {

    treat <- (log(abs(-0.8 + 0.1 * cf[ , 1] + 0.1 * cf[ , 2] - 0.1 * cf[ , 3]
             + 0.2 * cf[ , 4] + 0.1 * cf5 + 0.1 * cf6)) * 7 + 13 + stats::rnorm(size,sd=4))

  } else if (gps_spec == 7) {

    treat <- ((-0.8 + 0.1 * cf[,1] + 0.1 * cf[,2] - 0.1 * cf[,3] + 0.2 * cf[,4]
             + 0.1 * cf5 + 0.1 * cf6) * 15 + 22 + stats::rt(size,2)) #+ rcauchy(size)
  } else {

    stop(paste("gps_spec: ", gps_spec, ", is not a valid value."))

  }
  
  if (em_as_confounder) {
    treat <- treat + em1*3
  }

  
  # bin exposure level into 20 categories - only necessary for grouped outcomes in the binary case
  if (outcome_type == 'binary') {
    treat_df <- data.frame(treat = treat) %>%
      mutate(treat_level = cut(treat, breaks = 20)) %>%
      dplyr::group_by(treat_level) %>% 
      dplyr::mutate(avg_treat = mean(treat))
  }
  else if (outcome_type == 'continuous') {
    treat_df <- data.frame(treat = treat) %>%
      dplyr::mutate(avg_treat = treat)
  }
  
    
  # produce outcome Y
  mu <- as.numeric()
  Y <- as.numeric()
  
  # em2 = 0, em1 = 0 -> 10
  # em2 = 0, em1 = 1 -> 15
  # em2 = 1, em1 = 1 -> 7
  # em2 = 1, em2 = 0 -> 2
  for (i in 1:size) {
    if (em_spec == 0) {
      mu[i] <- -10 - sum(c(2,2,3,-1) * cf[i, ]) - 2 * (cf5[i] + cf6[i]) + 10 * treat_df$avg_treat[i]
    }
    if (em_spec == 1) {
      mu[i] <- -10 - sum(c(2, 2, 3, -1) * cf[i, ]) - 2 * (cf5[i] + cf6[i]) + 
        10 * (1 + beta * (0.5 * em1[i] - 0.8 * em2[i])) * treat_df$avg_treat[i]
    }
    else if (em_spec == 2) {
      mu[i] <- -10 - sum(c(2, 2, 3, -1) * cf[i, ]) - 2 * (cf5[i] + cf6[i]) + 
        10 * (1 + beta * (0.5 * em1[i] * em2[i] + 0.2 * em2[i])) * treat[i]
    }
    if (heterogenous_intercept) {
      mu[i] <- mu[i] + 20 * (em1[i] + 2*em2[i])
    }
    if (outcome_type == 'continuous') {
      Y[i] <- mu[i] + stats::rnorm(1, mean=0, sd=outcome_sd)
    }
    if (outcome_type == 'binary') {
      Y[i] <- rpois(1, mu[i])
    }
  }
  # if (outcome_type == 'binary') {
  #   simulated_data <- 
  #     as.data.frame(cbind(treat_level = treat_df$treat_level, avg_treat = treat_df$avg_treat, 
  #                         mu, em1, em2, em3, em4)) %>%
  #     mutate(treat_level = cut(treat, breaks = 20)) %>%
  #     group_by(em1, em2, em3, em4, treat_level, avg_treat, mu) %>%
  #     summarise(Y = rpois(1, n()*mu), 
  #                         n = n()) %>%
  #     mutate_at(vars(em1, em2, em3, em4), as.factor)
  #   #colnames(simulated_data)[2:12]<-c("cf1","cf2","cf3","cf4","cf5","cf6","em1","em2","em3","em4")
  # }

  # # why are the covariates changed after the propensity score model is used to define treatment??
  # if (cova_spec == 1) {
  #   cf = cf
  # } else if (cova_spec == 2) {
  # 
  #   cf[,1] <- exp(cf[ ,1] / 2)
  #   cf[,2] <- (cf[ ,2] / (1 + exp(cf[ ,1]))) + 10
  #   cf[,3] <- (cf[ ,1] * cf[ ,3]/25 + 0.6) ^ 3
  #   cf[,4] <- (cf[ ,2] + cf[ ,4] + 20) ^ 2
  # 
  # } else {
  #   stop(paste("cova_spec: ", cova_spec, ", is not a valid value."))
  # }
  
  if (outcome_type == 'continuous') {
    simulated_data<-data.frame(cbind(Y,treat,cf, cf5, cf6, em1, em2, em3, em4)) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    colnames(simulated_data)[3:12]<-c("cf1","cf2","cf3","cf4","cf5","cf6","em1","em2","em3","em4")
  }
  if (outcome_type == 'binary') {
    simulated_data<-data.frame(cbind(Y,treat,cf, cf5, cf6, em1, em2, em3, em4, mu)) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    colnames(simulated_data)[3:12]<-c("cf1","cf2","cf3","cf4","cf5","cf6","em1","em2","em3","em4", 'mu')
  }
  return(simulated_data)
}
 

