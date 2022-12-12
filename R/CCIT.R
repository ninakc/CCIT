#' @title
#' Estimate conditional average treatment effect
#'
#' @description
#' Performs CCIT algorithm to estimate conditional average treatment effect.
#'
#'
#' @param matched.exploration.sample.outcomes: The exploration dataframe with
#' outcome Y, exposure variable treat.
#' @param matched.validation.sample.outcomes: The validation dataframe
#' with outcome Y, exposure variable treat.
#' @param matched.inference.sample.outcomes: The inference dataframe
#' with outcome Y, exposure variable treat.
#' @param lambdas: A vector of values to use as the regularization parameter in
#' the CCIT algorithm.
#' @param stopping.rule: A boolean value to indicate whether the tree-splitting
#'  algorithm should stop when the estimated interaction effect is <1/10 of
#'  the overall effect.
#'
#' @returns list of:
#'       - est.treatment.effects: For each lambda value, dataframe with CATE
#'         estimates for all observations
#'       - selected.trees: For each lambda value, selected decision tree
#'       - tree.list: all decisions trees in the sequence
#'       - selected.tree.size: For each lambda value, size of selected decision
#'         tree
#'
CCIT <- function(matched.exploration.sample.outcomes,
                 matched.validation.sample.outcomes,
                 matched.inference.sample.outcomes,
                 lambdas,
                 stopping.rule) {

  overall_effect <- lm(Y ~ treat + 1,
                       data = matched.exploration.sample.outcomes)$coefficients["treat"]

  parms.used <- list(
    treat = matched.exploration.sample.outcomes$treat,
    Y = matched.exploration.sample.outcomes$Y,
    overall_effect = overall_effect
  )
  ulist.used <- define_u.list(stopping.rule = stopping.rule)
  lists <- create.sequence(matched.exploration.sample.outcomes,
                           ulist.used,
                           parms.used)

  tree.list <- lists[[1]]
  g.h.list <- lists[[2]]
  tree.sizes <- sapply(tree.list, function(t) {
    ifelse(is.null(t$splits), 0, nrow(t$splits))
  })

  complex.vals <- evaluate.sequence(tree.list, matched.validation.sample.outcomes,
                                    matched.exploration.sample.outcomes, lambdas)

  selected.tree.size <- complex.vals %>%
    group_by(lambda) %>%
    filter(complex.val == max(complex.val)) %>%
    select(tree.size) %>%
    ungroup()

  selected.trees <- lapply(selected.tree.size$tree.size, function(s) {
    tree.list[[which(tree.sizes == s)]]
  })
  names(selected.trees) <- lambdas

  # determine which selected subgroup each inference observation is in.
  # There are a number of packages that do this, but the most recent versions
  # are not installed on the RCE
  matched.inference.sample.subgroups <- lapply(
    selected.trees,
    function(t) {
      matched.inference.sample.new <- matched.inference.sample.outcomes
      matched.inference.sample.new$subgroup <- 0
      leaf_nodes <- t$frame %>%
        tibble::rownames_to_column() %>%
        filter(var == "<leaf>") %>%
        .$rowname %>%
        as.numeric()
      for (n in leaf_nodes) {
        rule <- path.rpart(t, n)
        if (nrow(t$frame) == 1) {
          ind <- 1:nrow(matched.inference.sample.outcomes)
        } else {
          rule_2 <- sapply(rule[[1]][-1],
                           function(x) strsplit(x, "(?<=[><=])(?=[^><=])|(?<=[^><=])(?=[><=])",
                                                perl = TRUE))
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
            get(comparison_operators[i])(matched.inference.sample.outcomes[, splitting_vars[i]],
                                         comparison_vals[i])
          }) %>% apply(1, all)

        }
        matched.inference.sample.new[ind, ]$subgroup <- n
      }
      matched.inference.sample.new <- matched.inference.sample.new %>%
        mutate_at(vars(subgroup), as.factor)
      return(matched.inference.sample.new)
    }
  )

  # estimate treatment effects on the inference subsample by fitting linear
  # models within each selected subgroup
  est.treatment.effects <- lapply(
    matched.inference.sample.subgroups,
    function(df) {
      df %>%
        group_by(subgroup) %>%
        do(model = lm(Y ~ treat + 1, data = .)) %>%
        mutate(est = summary(model)$coefficients["treat", 1]) %>%
        right_join(df, by = "subgroup")
    }
  )
  return(list(est.treatment.effects = est.treatment.effects,
              selected.trees = selected.trees,
              tree.list = tree.list,
              selected.tree.size = selected.tree.size))
}
