#' @description
#' This function calculates the sequence of candidate trees
#'
#'
#' @param exploration_dat: exploration dataframe with outcome Y and exposure variable treat
#' @param ulist_used: user-specified functions to supply to rpart (return value of define_u_list)
#' @param parms_used: additional parameters to supply to rpart
#'
#'
#' @returns: list of nested decision trees and corresponding goodness metrics
create_sequence = function(exploration_dat, ulist_used, parms_used){

  exploration_dat <- exploration_dat %>% tibble::rowid_to_column('dummy')

  # Fit a large tree using the user written splitting functions
  a <- rpart::rpart(
    dummy ~ em1 + em2 + em3 + em4,
    # f,
    data = exploration_dat,
    method   = ulist_used,
    parms    = parms_used,
    control  = rpart.control(cp = -10, minbucket = 1, maxsurrogate = 0, maxcompete = 0))

  if (dim(a$frame)[1] == 1){ # Deal with the root only tree

    tree_list = list(a)
    g_h_list = list(Inf)
    return(list(tree_list = tree_list, g_h_list = g_h_list))

  } else {

    # Finding which variables are leaf nodes
    is_leaf <- (a$frame$var == "<leaf>")

    # A function adapted from the partykit package that identifies the rows of the frame
    # which correspond to child nodes of row i in frame matrix
    rpart_kids <- function(i, is_leaf) {
      if (is_leaf[i]) return(NULL)
      else return(c(i + 1L,
                    which((cumsum(!is_leaf[-(1L:i)]) + 1L) == cumsum(is_leaf[-(1L:i)]))[1L] + 1L + i))
    }

    # Finding goodness of the split
    a$frame$split_stat = 0
    a$frame$split_stat[!is_leaf] = a$splits[, 3]

    # Calculating the g(h) parameter for each non-terminal node
    g_h = rep(0, nrow(a$frame))
    for(i in 1:nrow(a$frame)){
      if(is_leaf[i]){g_h[i] = Inf} else{
        # Find all kids of node i
        kids_i = i
        stop_loop = FALSE
        while(stop_loop == FALSE){
          kids_old = kids_i
          for(j in 1:length(kids_i)){
            kids_i = unique(c(kids_i, rpart_kids(kids_i[j], is_leaf)))
          }
          if(length(kids_i) == length(kids_old)){stop_loop = TRUE}
          # Calculating g_h for node i
          g_h[i] = sum(a$frame$split_stat[kids_i])/sum(a$frame$split_stat[kids_i] != 0)
        }
      }
    }

    # Adding g_h to frame
    a$frame$g_h = g_h

    # Start pruning
    # First tree is the large tree
    tree_list = list(a)
    g_h_list = list(0)
    # pruned.node = list(NULL)
    stop.prune = FALSE
    k = 1

    while(stop.prune == FALSE){
      tree_used = tree_list[[k]]
      # Calculating the g(h) parameter for each non-terminal node
      tree_used$frame$g_h = rep(0, nrow(tree_used$frame))
      is_leaf.prune <- (tree_used$frame$var == "<leaf>")
      # Setting splitting statistics for new terminal nodes to 0
      tree_used$frame$split_stat[is_leaf.prune] = 0

      # Calculating the g(h) function for each non-terminal node
      for(i in 1:nrow(tree_used$frame)){
        if(is_leaf.prune[i]){tree_used$frame$g_h[i] = Inf} else{
          # Find all kids of node i
          kids_i = i
          stop_loop = FALSE
          while(stop_loop == FALSE){
            kids_old = kids_i
            for(j in 1:length(kids_i)){
              kids_i = unique(c(kids_i, rpart_kids(kids_i[j], is_leaf.prune)))
            }
            if(length(kids_i) == length(kids_old)){stop_loop = TRUE}
            tree_used$frame$g_h[i] = sum(tree_used$frame$split_stat[kids_i])/sum(tree_used$frame$split_stat[kids_i] != 0)
          }
        }
      }

      # Finding the value which minimizes g(h) (among internal nodes)
      to.prune = which.min(tree_used$frame$g_h)
      #print(tree_used$frame[to.prune,])
      # Finding the minimum g_h value
      g_h.min = min(tree_used$frame$g_h)

      # Find all kids of node to.prune
      kids_i = to.prune
      stop_loop = FALSE
      while(stop_loop == FALSE){
        kids_old = kids_i
        for(j in 1:length(kids_i)){
          kids_i = unique(c(kids_i, rpart_kids(kids_i[j], is_leaf.prune)))
        }
        if(length(kids_i) == length(kids_old)){stop_loop = TRUE}
      }


      # Finding number of splits to prune
      split.to.prune = length(kids_i[which(!is_leaf.prune[kids_i])])

      # Creating the new splits and frames for new tree
      splits.new = tree_used$splits[-c(sum(!is_leaf.prune[1:to.prune]):(sum(!is_leaf.prune[1:to.prune]) + split.to.prune - 1)), ]
      frame.new = tree_used$frame[-setdiff(kids_i, to.prune), ]

      # Changing all nodes that were internal nodes and are now terminal node to terminal nodes
      frame.new$var[to.prune] =  "<leaf>"

      tree.new = tree_used
      tree.new$frame = frame.new
      print(class(splits.new))
      if("matrix" %in% class(splits.new)){tree.new$splits = splits.new}
      if("numeric" %in% class(splits.new)){
        tree.new$splits = matrix(splits.new, nrow = 1)
        colnames(tree.new$splits) = colnames(tree_used$splits)
      }

      # Changing the terminal node for $where in rpart object
      tree.new$where = tree_used$where
      tree.new$where[tree.new$where %in% kids_i] = to.prune
      tree.new$where[tree.new$where > max(kids_i)] = tree.new$where[tree.new$where > max(kids_i)] - length(kids_i) + 1
      tree.new$where = as.integer(tree.new$where)

      k = k+1
      # Add tree and lambda to the list
      tree_list[[k]] <- tree.new
      g_h_list[[k]] <- g_h.min
      # pruned.node[[k]] <-

      if(sum(tree.new$frame$var == "<leaf>") == 1){stop.prune = TRUE}
    }
    return(list(tree_list = tree_list, g_h_list = g_h_list))
  }
}
