# join_flipscores.R
# Embedded and adapted from jointest (https://github.com/livioivil/jointest)
# Original author: Livio Finos

#' Join Flip Scores
#'
#' Runs \code{flipscores()} on multiple models and combines the results
#' into a \code{joint_flipscores} object.
#'
#' @param models A list of \code{glm} objects or formulas.
#' @param combine_with Combining function: \code{"Fisher"} (default), 
#'   \code{"Liptak"}, \code{"Tippett"}, or a custom function.
#' @param score_type Type of score. See \code{\link{flipscores}}.
#' @param tested_terms Tested terms. See \code{\link{flipscores}}.
#' @param n_flips Number of flips. See \code{\link{flipscores}}.
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to \code{\link{flipscores}}.
#'
#' @return A \code{joint_flipscores} object.
#' @keywords internal
join_flipscores <- function(models,
                            combine_with  = "Fisher",
                            score_type    = "standardized",
                            tested_terms  = NULL,
                            n_flips       = 1000,
                            seed          = NULL,
                            ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.list(models)) {
    stop("`models` must be a list of glm objects or formulas.")
  }
  
  # Run flipscores on each model
  results <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    
    if (!inherits(model, "glm")) {
      stop(sprintf("Element %d of `models` is not a glm object.", i))
    }
    
    flipscores(
      formula      = model,
      score_type   = score_type,
      tested_terms = tested_terms,
      n_flips      = n_flips,
      ...
    )
  })
  
  # Name results if models are named
  if (!is.null(names(models))) {
    names(results) <- names(models)
  } else {
    names(results) <- paste0("model_", seq_along(results))
  }
  
  # Combine p-values
  combined <- combine_tests.joint_flipscores(
    structure(list(results = results), class = "joint_flipscores"),
    combine_with = combine_with
  )
  
  structure(
    list(
      results      = results,
      combined     = combined,
      combine_with = combine_with,
      call         = match.call()
    ),
    class = "joint_flipscores"
  )
}