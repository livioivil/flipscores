# utils_jointest.R
# Internal utilities embedded from jointest
# Original author: Livio Finos

#--------------------------------------------
# Rename jointest class to joint_flipscores
# Used as a post-processing step if any jointest
# internals return "jointest" class objects
#--------------------------------------------
.rename_jointest_class <- function(x) {
  if (inherits(x, "jointest")) {
    class(x) <- gsub("jointest", "joint_flipscores", class(x))
  }
  x
}


#--------------------------------------------
# Validate that a list contains only glm objects
#--------------------------------------------
.validate_glm_list <- function(models) {
  if (!is.list(models)) {
    stop("`models` must be a list.")
  }
  not_glm <- which(!sapply(models, inherits, "glm"))
  if (length(not_glm) > 0) {
    stop(sprintf(
      "Elements at positions %s are not glm objects.",
      paste(not_glm, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#--------------------------------------------
# Validate that a list contains only formulas
#--------------------------------------------
.validate_formula_list <- function(formulas) {
  if (!is.list(formulas)) {
    stop("`formulas` must be a list.")
  }
  not_formula <- which(!sapply(formulas, inherits, "formula"))
  if (length(not_formula) > 0) {
    stop(sprintf(
      "Elements at positions %s are not formula objects.",
      paste(not_formula, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#--------------------------------------------
# Extract score table from a flipscores object
# Handles different possible structures robustly
#--------------------------------------------
.get_score_table <- function(fs_obj) {
  # Try common slot names used in flipscores objects
  if (!is.null(fs_obj$score_table))   return(fs_obj$score_table)
  if (!is.null(fs_obj$table))         return(fs_obj$table)
  if (!is.null(fs_obj$coefficients))  return(fs_obj$coefficients)

  # Try summary
  s <- tryCatch(summary(fs_obj)$score_table, error = function(e) NULL)
  if (!is.null(s)) return(s)

  stop("Cannot extract score table from flipscores object.")
}

#--------------------------------------------
# Extract p-values from a flipscores object
#--------------------------------------------
.get_pvalues <- function(fs_obj) {
  tbl <- .get_score_table(fs_obj)

  if ("p.value" %in% colnames(tbl)) {
    return(tbl[, "p.value"])
  }

  stop("Cannot find p.value column in score table.")
}

#--------------------------------------------
# Extract tested terms from a flipscores object
#--------------------------------------------
.get_tested_terms <- function(fs_obj) {
  tbl <- .get_score_table(fs_obj)
  rownames(tbl)
}

#--------------------------------------------
# Safe model name extractor
# Returns names or generates default ones
#--------------------------------------------
.get_model_names <- function(models) {
  nms <- names(models)
  if (is.null(nms) || any(nms == "")) {
    nms <- paste0("model_", seq_along(models))
  }
  nms
}

#--------------------------------------------
# Build a summary data frame from a
# joint_flipscores object (used internally
# by print/summary/plot methods)
#--------------------------------------------
.build_summary_df <- function(x) {
  do.call(rbind, lapply(names(x$results), function(nm) {
    tbl <- tryCatch(
      .get_score_table(x$results[[nm]]),
      error = function(e) NULL
    )
    if (is.null(tbl)) return(NULL)

    df <- as.data.frame(tbl)
    df$model <- nm
    df$term  <- rownames(tbl)
    rownames(df) <- NULL

    # Reorder columns: model, term first
    cols <- c("model", "term", setdiff(colnames(df), c("model", "term")))
    df[, cols, drop = FALSE]
  }))
}
