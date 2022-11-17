#' @title
#' Generate Synthetic Data for CausalGPS Package - Modified by Nina
#'
#' @description
#' Generates synthetic data set based on different GPS models and covariates.
#'
#'Required Libraries: truncnorm, dplyr
generate_syn_data_covs <- function(sample_size=10000, gps_spec = 1) {

  if (sample_size < 0 || !is.numeric(sample_size)){
    stop("'sample_size' should be a positive ineteger numer.")
  }

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
  if (gps_spec == 1) { # no confounding (treatment is independent of covariates)
    # modified Xiao's propensity function to ensure treatment >= 0. Note that mu is not the expected propensity.
    mu <- 3
    treat <- rtruncnorm(size,a=0,b=Inf,mean=mu,sd=5)
  } else if (gps_spec == 2) { # with confounding
    mu <- (- 0.8 + 0.1 * cf[ ,1] + 0.1 * cf[ ,2] - 0.1 * cf[ ,3]
           + 0.2 * cf[ ,4] + 0.1 * cf5 + 0.1 * cf6) * 9 + 17

    treat <- rtruncnorm(size,a=0,b=Inf,mean=mu,sd=5)

  } else {

    stop(paste("gps_spec: ", gps_spec, ", is not a valid value."))

  }

  simulated_covariates <- as.data.frame(cf) %>%
    mutate(cf5 = cf5, cf6 = cf6, em1 = em1, em2 = em2, em3 = em3, em4 = em4, treat = treat)

  return(simulated_covariates)
}


