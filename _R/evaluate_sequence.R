#' @description
#' Evaluate list of trees returned by create_sequence function
#'
#' @param tree_list: list of trees returned by create_sequence function
#' @param val: validation dataframe with outcome Y, exposure variable treat
#' @param exploration_dat: exploration dataframe with outcome Y, exposure variable treat
#' @param lambdas: regularization values to use
#'
#' @returns dataframe with one row for each tree in sequence and lambda value, representing the ability of the tree to explain the effect heterogeneityß

evaluate_sequence <- function(tree_list, val, exploration_dat, lambdas) {
  val = as.data.frame(val)

  # Storage space for the complexity value associated with each candidate tree
  #complex_val = rep(NA, length(tree_list))
  complex_val <- data.frame()

  # Computing G on validation set for each tree in sequence
  for (m in 1:length(tree_list)){
    tree_used = tree_list[[m]]

    # If only root node there is no internal node
    if(nrow(tree_used$frame) == 1){
      goodness_test = 0
      numb_int = 0

    } else { # If at least one split

      is_leaf <- (tree_used$frame$var == "<leaf>")
      goodness_test = 0
      # Finding the test data falling in each terminal node (each non-terminal node??)
      numb_int = sum(!is_leaf)

      right_ind  <- 1
      last_right <- NULL

      for (h in 1:dim(tree_used$frame)[1]){
        #print(h)

        # Finding observations in validation sample falling in that node
        if (tree_used$frame$var[h] == "<leaf>") {
          next
        } else if (tree_used$frame$n[h] == dim(exploration_dat)[1]){ # if this node has all training obs
          val_sample_used <- val
        } else{
          if (!is.null(last_left)){
            val_sample_used <- last_left[[1]]
          } else {
            val_sample_used <- as.data.frame(last_right[[right_ind - 1]])
            if (length(last_right) > 2){
              last_right   <- last_right[1:(right_ind - 2)]
            } else if (length(last_right) == 2){
              last_right   <- list(last_right[[1]])
            }else {
              last_right <- NULL
            }
            right_ind <- right_ind - 1
          }
        }
        row_ind <- sum((tree_used$frame$var[1:h]) != "<leaf>")
        split_used <- tree_used$splits[row_ind, 4]
        var_used <- tree_used$frame$var[h] # splitting variable
        col_ind <- which(colnames(val_sample_used) == var_used)

        lvls <- levels(val[, col_ind])
        # lvls[tree_used$csplit[split_used, ] == 1] returns the value of the EM variable that goes to the left
        val_sample_left  <- val_sample_used[val_sample_used[, col_ind] %in% lvls[tree_used$csplit[split_used, ] == 1], ]
        val_sample_right <- val_sample_used[val_sample_used[, col_ind] %in% lvls[tree_used$csplit[split_used,] == 3], ]
        if (tree_used$frame$var[h+1] != "<leaf>"){
          last_left   <- list(val_sample_left)
        } else{
          last_left   <- NULL
        }
        # index of node to the right
        which_right <- as.numeric(rownames(tree_used$frame)[h+1]) + 1
        # if the next node is not a leaf, store current data as last_right
        if (tree_used$frame$var[as.numeric(rownames(tree_used$frame)) == which_right] != "<leaf>"){
          last_right[[right_ind]]   <- val_sample_right
          right_ind <- right_ind + 1
        }
        f <- as.formula(paste0('Y ~ ', var_used, '*treat'))
        lmod <- lm(f, data = val_sample_used)
        t <- abs(coef(summary(lmod))[4,"t value"])
        goodness_test <- goodness_test + t**2
      } # End h
    } # End if loop
    # Calculating complexity value

    complex_val <- rbind(complex_val, data.frame(complex_val = goodness_test - lambdas * numb_int,
                                        lambda = lambdas,
                                        pruning_step = m,
                                        tree_size = ifelse(is.null(tree_used$splits), 0, nrow(tree_used$splits))))
    #names(complex_val[[m]]) <- lambdas
  } # End m loop
  return (complex_val)
}

