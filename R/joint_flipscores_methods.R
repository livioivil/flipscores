# joint_flipscores.R
# S3 methods for the "joint_flipscores" class

#--------------------------------------------
# internal helper: print first and last n rows
# adapted from jointest
#--------------------------------------------
.trim <- function(x, n = 2, ...) {
  nr <- nrow(x)
  if (is.null(nr) || nr <= 2 * n) {
    print(x, ...)
  } else {
    print(head(x, n), ...)
    cat(sprintf("... %s rows ...\n", nr - 2 * n))
    print(tail(x, n), ...)
  }
  invisible(x)
}

#--------------------------------------------
#' Print method for joint_flipscores
#'
#' @param x A \code{joint_flipscores} object.
#' @param n Number of rows to show at head and tail of summary table.
#' @param ... Additional arguments (currently unused).
#' @export
print.joint_flipscores <- function(x, n = 2, ...) {
  msg <- "== Joining n = %s models"
  cat(sprintf(msg, length(unique(x$summary_table$Model))))
  cat("\n\n")
  .trim(x$summary_table, n = n)
  invisible(x)
}

#--------------------------------------------
#' Summary method for joint_flipscores
#'
#' @param object A \code{joint_flipscores} object.
#' @param digits Number of digits to print. Default \code{4}.
#' @param ... Additional arguments (currently unused).
#' @export
summary.joint_flipscores <- function(object, digits = 4, ...) {
  tab <- object$summary_table
  tab$.assign <- NULL   # remove internal column if present
  print(tab, digits = digits)
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
             xlab   = "score",
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
