#' @title
#' Stratify GPS matching
#'
#' @description
#' TBD
#'
#' @param data: A data frame with data to be matched.
#' @param delta_n: A bin width parameter to pass to the CausalGPS package.
#' @param exposure_name: The name of the exposure variable.
#' @param confounders_names: A vector of strings representing names of the
#' confounding variables.
#' @param names_of_strata_vars: A vector of strings representing names of
#' variables to stratify by.
#' @param outcome_name: The name of the outcome variable.

stratified_GPS_matching <- function(data, delta_n, exposure_name,
                                    confounders_names,
                                    names_of_strata_vars,
                                    outcome_name) {

  if (length(names_of_strata_vars) > 0) {
    strata_var <- names_of_strata_vars[1]
    data[,strata_var] <- as.factor(data[,strata_var])
    strata_lvls <- levels(data[,strata_var])
    combined_matched_data <- data.frame()
    for (strata_lvl in strata_lvls)
      {
      combined_matched_data <-
        rbind(combined_matched_data,
               stratified_GPS_matching(
                filter(data, get(strata_var) == strata_lvl),
                delta_n = delta_n,
                exposure_name = exposure_name,
                confounders_names = confounders_names,
                names_of_strata_vars = names_of_strata_vars[names_of_strata_vars != strata_var],
                outcome_name = outcome_name
              )
            )
    }
    return(combined_matched_data)
  }
  else { # only a single stratum left
    if(nrow(data) > 0) {
      sub_pop <- data %>%
        tibble::rowid_to_column('new_id') #row index of sub-pop
      pseudo_pop_fit <- CausalGPS::generate_pseudo_pop(ifelse(is.na(outcome_name), 1:nrow(sub_pop), sub_pop[,outcome_name]),  # why is the outcome required?
                                            sub_pop[,exposure_name],
                                            dplyr::select(sub_pop, confounders_names),
                                            ci_appr = "matching",
                                            pred_model = "sl",
                                            gps_model = "parametric",
                                            use_cov_transform = TRUE,
                                            transformers = list("pow2", "pow3"),
                                            # bin_seq = bin_seq,
                                            sl_lib = c("m_xgboost"),
                                            params = list('xgb_max_depth' = c(3,4,5),
                                                          'xgb_nrounds'=c(10,20,30,40,50,60),
                                                          #"xgb_nrounds"=50,
                                                          #"xgb_max_depth"=6,
                                                          "xgb_eta"=0.3,
                                                          "xgb_min_child_weight"=1),
                                            nthread=5, # number of cores, you can change,
                                            covar_bl_method = "absolute",
                                            covar_bl_trs = 0.1,
                                            covar_bl_trs_type = "mean",
                                            trim_quantiles = c(0.05,0.95), # trimed, you can change
                                            optimized_compile = FALSE, #created a column counter for how many times matched,
                                            max_attempt = 1,
                                            matching_fun = "matching_l1",
                                            delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                            scale = 1.0)
      weights <- count(pseudo_pop_fit$pseudo_pop, row_index)
      sub_pop <- left_join(sub_pop, weights, by = c('new_id' = 'row_index'))
      sub_pop <- sub_pop %>% mutate(n = ifelse(is.na(n), 0, n))
      sub_pop <- as.data.frame(lapply(sub_pop, rep, sub_pop$n))
      return(dplyr::select(sub_pop, -n))
    }
    else {
      return(data.frame())}
  }
}
