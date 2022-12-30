#' @description
#' Perform CCIT algorithm to estimate conditional average treatment effect
#'
#'
#' @param matched_exploration_sample_outcomes: exploration dataframe with outcome Y, exposure variable treat
#' @param matched_validation_sample_outcomes: validation dataframe with outcome Y, exposure variable treat
#' @param matched_inference_sample_outcomes: inference dataframe with outcome Y, exposure variable treat
#' @param lambdas: vector of values to use as the regularization parameter in the CCIT algorithm
#' @param stopping_rule: boolean to indicate whether the tree-splitting algorithm should stop when the estimated interaction effect is <1/10 of the overall effect
#'
#' @returns list of:
#'       est_treatment_effects: For each lambda value, dataframe with CATE estimates for all observations
#'       selected_trees: For each lambda value, selected decision tree
#'       tree_list: all decisions trees in the sequence
#'       selected_tree_size: For each lambda value, size of selected decision tree
CCIT <- function(matched_exploration_sample_outcomes = NULL, matched_validation_sample_outcomes = NULL, matched_inference_sample_outcomes = NULL,
                 lambdas, stopping_rule) {

  overall_effect <- lm(Y ~ treat + 1, data = matched_exploration_sample_outcomes)$coefficients["treat"]

  parms_used <- list(
    treat = matched_exploration_sample_outcomes$treat,
    Y = matched_exploration_sample_outcomes$Y,
    overall_effect = overall_effect
  )
  ulist_used <- define_u_list(stopping_rule = stopping_rule)
  lists <- create_sequence(matched_exploration_sample_outcomes, ulist_used, parms_used)

  tree_list <- lists[[1]]
  g_h_list <- lists[[2]]
  tree_sizes <- sapply(tree_list, function(t) {
    ifelse(is.null(t$splits), 0, nrow(t$splits))
  })

  complex_vals <- evaluate_sequence(tree_list, matched_validation_sample_outcomes, matched_exploration_sample_outcomes, lambdas)

  selected_tree_size <- complex_vals %>%
    group_by(lambda) %>%
    filter(complex_val == max(complex_val)) %>%
    select(tree_size) %>%
    ungroup()

  selected_trees <- lapply(selected_tree_size$tree_size, function(s) {
    tree_list[[which(tree_sizes == s)]]
  })
  names(selected_trees) <- lambdas

  # determine which selected subgroup each inference observation is in. There are a number of packages that do this, but the most recent versions
  # are not installed on the RCE
  matched_inference_sample_subgroups <- lapply(
    selected_trees,
    function(t) {
      matched_inference_sample_new <- matched_inference_sample_outcomes
      matched_inference_sample_new$subgroup <- 0
      leaf_nodes <- t$frame %>%
        tibble::rownames_to_column() %>%
        filter(var == "<leaf>") %>%
        .$rowname %>%
        as.numeric()
      for (n in leaf_nodes) {
        rule <- path.rpart(t, n)
        if (nrow(t$frame) == 1) {
          ind <- 1:nrow(matched_inference_sample_outcomes)
        } else {
          rule_2 <- sapply(rule[[1]][-1], function(x) strsplit(x, "(?<=[><=])(?=[^><=])|(?<=[^><=])(?=[><=])", perl = TRUE))
          splitting_vars <- sapply(rule_2, function(rule) {
            rule[1]
          })
          comparison_operators <- sapply(rule_2, function(rule) {
            ifelse(rule[2] == "=", "==", rule[2])
          })
          comparison_vals <- sapply(rule_2, function(rule) {
            rule[3]
          })
          ind <- sapply(seq(length(rule_2)), function(i) {
            get(comparison_operators[i])(matched_inference_sample_outcomes[, splitting_vars[i]], comparison_vals[i])
          }) %>% apply(1, all)

        }
        matched_inference_sample_new[ind, ]$subgroup <- n
      }
      matched_inference_sample_new <- matched_inference_sample_new %>%
        mutate_at(vars(subgroup), as.factor)
      return(matched_inference_sample_new)
    }
  )

  # estimate treatment effects on the inference subsample by fitting linear models within each selected subgroup
  est_treatment_effects <- lapply(
    matched_inference_sample_subgroups,
    function(df) {
      df %>%
        group_by(subgroup) %>%
        do(model = lm(Y ~ treat + 1, data = .)) %>%
        mutate(est = summary(model)$coefficients["treat", 1]) %>%
        right_join(df, by = "subgroup")
    }
  )
  return(list(est_treatment_effects = est_treatment_effects, selected_trees = selected_trees, tree_list = tree_list, selected_tree_size = selected_tree_size))
}
