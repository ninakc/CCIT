#' @title
#' Split data set
#'
#' @description
#' Splits data set into training, validation, and inference, stratifying
#' by exposure level.
#'
#' Required Libraries: caret, dplyr
#'
#' @param data: A data.frame with variable "treat" corresponding to
#' exposure level
#' @param num_exposure_cats: the number of categories to bin the exposure
#' level into for stratification

split_dataset <- function(data, num_exposure_cats) {

  a.vals <- seq(min(data$treat), max(data$treat), length.out = num_exposure_cats)

  data <-
    data %>%
    mutate(treat_level = cut(treat, breaks = a.vals)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)

  # split data into subsamples, stratifying on exposure level bins
  subsample.index <- caret::createFolds(as.factor(data$treat_level), list = TRUE, k = 3)
  data <- data %>%
    mutate(subsample = ifelse(row_number() %in% subsample.index[[2]], 'exploration',
                              ifelse(row_number() %in% subsample.index[[3]], 'validation', 'inference')))

  return(dplyr::select(data, -treat_level))
}

