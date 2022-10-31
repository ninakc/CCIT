#'Required Libraries: caret, dplyr
# split dataset into training, validation, and inference 
split_dataset <- function(data, exposure_bins) {
  subsample.index <- createFolds(as.factor(exposure_bins), list = TRUE, k = 3)
  # data <- data %>%
  #   mutate(subsample = ifelse(row_number() %in% subsample.index[[2]], 'exploration', 
  #                             ifelse(row_number() %in% subsample.index[[3]], 'validation', 
  #                                    ifelse(row_number() %in% subsample.index[[4]], 'generation', 'inference'))))
  data <- data %>%
    mutate(subsample = ifelse(row_number() %in% subsample.index[[2]], 'exploration', 
                              ifelse(row_number() %in% subsample.index[[3]], 'validation', 'inference')))
  
  return(data)
}

