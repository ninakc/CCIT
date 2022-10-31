evaluate.sequence <- function(tree.list, val, exploration.dat, lambdas) {
  val = as.data.frame(val)
  
  # Storage space for the complexity value associated with each candidate tree
  #complex.val = rep(NA, length(tree.list))
  complex.val <- data.frame()
  
  # Computing G on validation set for each tree in sequence
  for (m in 1:length(tree.list)){
    tree.used = tree.list[[m]] 
    
    # If only root node there is no internal node
    if(nrow(tree.used$frame) == 1){
      goodness.test = 0
      numb.int = 0
      
    } else { # If at least one split
      
      is.leaf <- (tree.used$frame$var == "<leaf>")
      goodness.test = 0
      # Finding the test data falling in each terminal node (each non-terminal node??)
      numb.int = sum(!is.leaf)
      
      right.ind  <- 1
      last.right <- NULL
      
      for (h in 1:dim(tree.used$frame)[1]){
        #print(h)

        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>") {
          next
        } else if (tree.used$frame$n[h] == dim(exploration.dat)[1]){ # if this node has all training obs
          val.sample.used <- val
        } else{
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right   <- last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right   <- list(last.right[[1]])
            }else {
              last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
        }
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used <- tree.used$frame$var[h] # splitting variable
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        lvls <- levels(val[, col.ind])
        # lvls[tree.used$csplit[split.used, ] == 1] returns the value of the EM variable that goes to the left
        val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
        val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left   <- list(val.sample.left)
        } else{
          last.left   <- NULL
        }
        # index of node to the right
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        # if the next node is not a leaf, store current data as last.right
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]]   <- val.sample.right
          right.ind <- right.ind + 1
        }
        f <- as.formula(paste0('Y ~ ', var.used, '*w'))
        lmod <- lm(f, data = val.sample.used)
        t <- abs(coef(summary(lmod))[4,"t value"])
        goodness.test <- goodness.test + t**2
      } # End h
    } # End if loop
    # Calculating complexity value

    complex.val <- rbind(complex.val, data.frame(complex.val = goodness.test - lambdas * numb.int, 
                                        lambda = lambdas, 
                                        pruning.step = m,
                                        tree.size = ifelse(is.null(tree.used$splits), 0, nrow(tree.used$splits))))
    #names(complex.val[[m]]) <- lambdas
  } # End m loop
  return (complex.val)
}

