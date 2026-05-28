# combiners.R
# Embedded and adapted from jointest (https://github.com/livioivil/jointest)
# Original author: Livio Finos

#--------------------------------------------
# Internal combiner function (not exported)
#--------------------------------------------
combine <- function(x, combine_with = "Fisher") {
  if (is.function(combine_with)) {
    return(combine_with(x))
  }
  
  combine_with <- match.arg(combine_with, 
                            c("Fisher", "Liptak", "Tippett", "Mahalanobis"))
  
  switch(combine_with,
    Fisher = {
      T_obs <- -2 * sum(log(x))
      list(
        T_obs   = T_obs,
        p.value = stats::pchisq(T_obs, df = 2 * length(x), lower.tail = FALSE)
      )
    },
    Liptak = {
      T_obs <- sum(stats::qnorm(1 - x))
      list(
        T_obs   = T_obs,
        p.value = stats::pnorm(T_obs, 
                               mean = 0, 
                               sd   = sqrt(length(x)), 
                               lower.tail = FALSE)
      )
    },
    Tippett = {
      T_obs <- min(x)
      list(
        T_obs   = T_obs,
        p.value = 1 - (1 - T_obs)^length(x)
      )
    },
    Mahalanobis = {
      stop("Mahalanobis combining function requires the full permutation 
            distribution. Use a custom function instead.")
    }
  )
}

#--------------------------------------------
#' Combine Tests
#'
#' Combines p-values from multiple tests using a combining function.
#'
#' @param x A \code{joint_flipscores} object or a numeric vector of p-values.
#' @param combine_with A string specifying the combining method 
#'   (\code{"Fisher"}, \code{"Liptak"}, \code{"Tippett"}) or a custom function.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with combined test statistic and p-value.
#' @export
combine_tests <- function(x, combine_with = "Fisher", ...) {
  UseMethod("combine_tests")
}

#' @export
combine_tests.default <- function(x, combine_with = "Fisher", ...) {
  if (!is.numeric(x)) stop("`x` must be a numeric vector of p-values.")
  result <- combine(x, combine_with = combine_with)
  data.frame(
    combine_with = if (is.function(combine_with)) "custom" else combine_with,
    T_obs        = result$T_obs,
    p.value      = result$p.value
  )
}

#' @export
combine_tests.joint_flipscores <- function(x, combine_with = "Fisher", ...) {
  # Extract p-values from each model's flipscores result
  pvals <- sapply(x$results, function(r) r$p.value)
  combine_tests.default(pvals, combine_with = combine_with, ...)
}

#--------------------------------------------
#' Combine Contrasts
#'
#' Combines p-values across contrasts (tested terms) from a 
#' \code{joint_flipscores} object.
#'
#' @param x A \code{joint_flipscores} object.
#' @param combine_with A string specifying the combining method or a function.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with combined statistics per contrast.
#' @export
combine_contrasts <- function(x, combine_with = "Fisher", ...) {
  UseMethod("combine_contrasts")
}

#' @export
combine_contrasts.joint_flipscores <- function(x, combine_with = "Fisher", ...) {
  # Gather all summary tables from each flipscores result
  tables <- lapply(x$results, function(r) {
    s <- r$score_table
    if (is.null(s)) stop("No score_table found in flipscores result.")
    s
  })
  
  # Get all unique tested terms across models
  all_terms <- unique(unlist(lapply(tables, rownames)))
  
  results <- lapply(all_terms, function(term) {
    # Extract p-values for this term across all models that have it
    pvals <- sapply(tables, function(tbl) {
      if (term %in% rownames(tbl)) tbl[term, "p.value"] else NA
    })
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) == 0) return(NULL)
    
    res        <- combine(pvals, combine_with = combine_with)
    data.frame(
      term         = term,
      combine_with = if (is.function(combine_with)) "custom" else combine_with,
      n_models     = length(pvals),
      T_obs        = res$T_obs,
      p.value      = res$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, results)
}