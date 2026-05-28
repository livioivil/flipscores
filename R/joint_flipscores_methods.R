# joint_flipscores.R
# S3 methods for the "joint_flipscores" class

#--------------------------------------------
#' Print method for joint_flipscores
#'
#' @param x A \code{joint_flipscores} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.joint_flipscores <- function(x, ...) {
  cat("\n=== Joint Flip Scores ===\n")
  cat("Number of models:", length(x$mods), "\n\n")

  cat("--- Summary table ---\n")
  print(x$summary_table, row.names = FALSE)

  invisible(x)
}


#--------------------------------------------
#' Summary method for joint_flipscores
#'
#' @param object A \code{joint_flipscores} object.
#' @param ... Additional arguments (currently unused).
#' @export
summary.joint_flipscores <- function(object, ...) {
  cat("\n=== Summary: Joint Flip Scores ===\n")
  cat("Number of models:", length(object$mods), "\n\n")

  cat("--- Per-model results ---\n")
  for (nm in names(object$mods)) {
    cat(sprintf("\n[%s]\n", nm))
    tab <- object$mods[[nm]]$summary_table
    if (!is.null(tab)) print(tab, row.names = FALSE)
  }

  cat("\n--- Overall summary table ---\n")
  print(object$summary_table, row.names = FALSE)

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
  n_models <- length(x$mods)
  old_par  <- par(mfrow = c(ceiling(n_models / 2), min(n_models, 2)),
                  mar   = c(4, 4, 3, 1))
  on.exit(par(old_par))

  for (nm in names(x$mods)) {
    tryCatch(
      plot(x$mods[[nm]], main = nm, ...),
      error = function(e) {
        # fallback: plot Tspace distribution of first tested coefficient
        ts <- x$mods[[nm]]$Tspace[, 1]
        hist(ts,
             main   = nm,
             xlab   = "Score",
             col    = "steelblue",
             border = "white",
             ...)
      }
    )
  }
  invisible(x)
}


#--------------------------------------------
#' Coerce joint_flipscores to data frame
#'
#' @param x A \code{joint_flipscores} object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with all summary tables stacked, with a
#'   \code{Model} column prepended.
#' @export
as.data.frame.joint_flipscores <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  x$summary_table
}
