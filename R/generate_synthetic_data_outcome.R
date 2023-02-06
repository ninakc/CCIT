#' @title
#' Generate outcome for synthetic data
#'
#' @description
#' Generates synthetic outcome variable given covariate, effect modifier, and
#' treatment variables.
#'
#' @param cf A data.frame of confounder variables
#' @param em A data.frame of effect modifier variables
#' @param treat: A vector of treatment variable
#' @param outcome_sd: standard deviation of outcome variable
#' @param em_spec: specification of effect modifier function
#' (possible values: 0, 1, 2)
#'   - em_spec = 0: no effect modification
#'   - em_spec = 1:independent effect modification of em1 and em2 (4 causal effect groups)
#'   - em_spec = 2 -> interactive effect modification of em1 and em2 (3 causal effect groups)
#' @param heterogenous_intercept: whether the intercept of the CRF should vary by effect modifier levels
#' @param beta: strength of effect modification (only effective if em_spec is 1 or 2)
generate_syn_data_outcome <- function(cf, em, treat,
                                      outcome_sd = 1,
                              em_spec = 1, heterogenous_intercept = FALSE,
                              beta = 1) {


  cf14 <- as.matrix(cf[, c("cf1", "cf2", "cf3", "cf4")])
  cf5 <- cf[["cf5"]]
  cf6 <- cf[["cf6"]]
  # produce outcome Y
  Y <- as.numeric()
  em1 <- as.numeric(as.character(em[["em1"]]))
  em2 <- as.numeric(as.character(em[["em2"]]))
  em3 <- as.numeric(as.character(em[["em3"]]))
  em4 <- as.numeric(as.character(em[["em4"]]))
  if (em_spec == 0) {
    mu <- -10 - sum(c(2,2,3,-1) * cf14) - 2 * (cf5 + cf6) + 10 * treat
  }
  if (em_spec == 1) {
    mu <- -10 - sum(c(2, 2, 3, -1) * cf14) - 2 * (cf5 + cf6) +
      10 * (1 + beta * (0.5 * em1 - 0.8 * em2)) * treat
  }
  else if (em_spec == 2) {
    mu <- -10 - sum(c(2, 2, 3, -1) * cf14) - 2 * (cf5 + cf6) +
      10 * (1 + beta * (0.5 * em1 * em2 + 0.2 * em2)) * treat
  }
  if (heterogenous_intercept) {
    mu <- mu + 20 * (em1 + 2*em2)
  }
  Y <- mu + stats::rnorm(length(treat), mean=0, sd=outcome_sd)

    simulated_data<-data.frame(cbind(Y,treat,cf14, cf5, cf6, em1, em2, em3, em4)) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    colnames(simulated_data)[3:12]<-c("cf1","cf2","cf3","cf4","cf5","cf6","em1","em2","em3","em4")

  return(simulated_data)
}


