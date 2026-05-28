# joint_flipscores.R
# S3 methods for the "joint_flipscores" class

#--------------------------------------------
#' Print method for joint_flipscores
#' @param x A \code{joint_flipscores} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.joint_flipscores <- function(x, ...) {
  cat("\n=== Joint Flip Scores ===\n")
  cat("Number of models:", length(x$results), "\n")
  cat("Combining method:", x$combine_with, "\n\n")

  cat("--- Per-model results ---\n")
  for (nm in names(x$results)) {
    cat(sprintf("\n[%s]\n", nm))
    print(x$results[[nm]])
  }

  cat("\n--- Combined test ---\n")
  print(x$combined)

  invisible(x)
}

#--------------------------------------------
#' Summary method for joint_flipscores
#' @param object A \code{joint_flipscores} object.
#' @param ... Additional arguments (currently unused).
#' @export
summary.joint_flipscores <- function(object, ...) {
  cat("\n=== Summary: Joint Flip Scores ===\n")
  cat("Number of models  :", length(object$results), "\n")
  cat("Combining method  :", object$combine_with, "\n\n")

  # Per-model summary table
  cat("--- Per-model p-values ---\n")
  per_model <- do.call(rbind, lapply(names(object$results), function(nm) {
    r <- object$results[[nm]]
    # Extract overall p-value from flipscores object
    pval <- tryCatch(r$score_table[, "p.value"], error = function(e) NA)
    data.frame(
      model        = nm,
      p.value      = if (length(pval) > 1) NA else pval,
      stringsAsFactors = FALSE
    )
  }))
  print(per_model, row.names = FALSE)

  cat("\n--- Combined test ---\n")
  print(object$combined)

  invisible(object)
}

#--------------------------------------------
#' Plot method for joint_flipscores
#'
#' Plots the score distributions for each model.
#'
#' @param x A \code{joint_flipscores} object.
#' @param ... Additional arguments passed to \code{plot}.
#' @export
plot.joint_flipscores <- function(x, ...) {
  n_models <- length(x$results)

  # Set up grid layout
  old_par <- par(mfrow = c(ceiling(n_models / 2), min(n_models, 2)),
                 mar   = c(4, 4, 3, 1))
  on.exit(par(old_par))

  for (nm in names(x$results)) {
    r <- x$results[[nm]]
    tryCatch(
      plot(r, main = nm, ...),
      error = function(e) {
        # Fallback: plot score distribution manually
        scores <- r$score_table[, "T_obs"]
        hist(scores,
             main  = nm,
             xlab  = "Score",
             col   = "steelblue",
             border = "white",
             ...)
      }
    )
  }
  invisible(x)
}

#--------------------------------------------
#' p.adjust method for joint_flipscores
#'
#' Adjusts p-values from a \code{joint_flipscores} object for
#' multiple comparisons.
#'
#' @param p A \code{joint_flipscores} object.
#' @param method Correction method. See \code{\link[stats]{p.adjust.methods}}.
#' @param n Number of comparisons. See \code{\link[stats]{p.adjust}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The \code{joint_flipscores} object with adjusted p-values.
#' @export
p.adjust.joint_flipscores <- function(p, method = "BH", n = length(p), ...) {
  # Adjust per-model score tables
  p$results <- lapply(p$results, function(r) {
    if (!is.null(r$score_table) && "p.value" %in% colnames(r$score_table)) {
      r$score_table[, "p.value.adj"] <- stats::p.adjust(
        r$score_table[, "p.value"],
        method = method
      )
    }
    r
  })

  # Adjust combined p-value if present
  if (!is.null(p$combined) && "p.value" %in% colnames(p$combined)) {
    p$combined[, "p.value.adj"] <- stats::p.adjust(
      p$combined[, "p.value"],
      method = method
    )
  }

  p
}

#--------------------------------------------
#' Coerce joint_flipscores to data frame
#'
#' Extracts all score tables from a \code{joint_flipscores} object
#' into a single tidy data frame.
#'
#' @param x A \code{joint_flipscores} object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with columns \code{model}, \code{term},
#'   and all columns from the score tables.
#' @export
as.data.frame.joint_flipscores <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  do.call(rbind, lapply(names(x$results), function(nm) {
    r   <- x$results[[nm]]
    tbl <- r$score_table
    if (is.null(tbl)) return(NULL)
    cbind(
      data.frame(model = nm, stringsAsFactors = FALSE),
      tbl
    )
  }))
}
