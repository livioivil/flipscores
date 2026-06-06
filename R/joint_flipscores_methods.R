#' @title Methods for jfs objects
#' @name jfs-methods
#' @description
#' Collection of methods for objects of class \code{jfs}.
#' @param x an object of class \code{jfs}
#' @param ... additional arguments
NULL

# S3 methods for the "jfs" class

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
#' Print method for jfs
#'
#' @param x A \code{jfs} object.
#' @param n Number of rows to show at head and tail of summary table.
#' @param ... Additional arguments (currently unused).
#' @export
print.jfs <- function(x, n = 2, ...) {
  msg <- "== Joining n = %s objects"
  model_col <- if ("model" %in% names(x$summary_table)) "model" else "Model"
  cat(sprintf(msg, length(unique(x$summary_table[[model_col]]))))
  cat("\n\n")
  .trim(x$summary_table, n = n)
  invisible(x)
}

#--------------------------------------------
#' Summary method for jfs
#'
#' @param object A \code{jfs} object.
#' @param digits Number of digits to print. Default \code{4}.
#' @param ... Additional arguments (currently unused).
#' @export
summary.jfs <- function(object, digits = 4, ...) {
  tab <- object$summary_table
  tab$.assign <- NULL   # remove internal column if present
  print(tab, digits = digits)
  invisible(object)
}


#--------------------------------------------
#' Plot method for jfs
#'
#' Plots the score distributions for each model.
#'
#' @param x A \code{jfs} object.
#' @param ... Additional arguments passed to \code{plot}.
#' @noRd
plot.jfs <- function(x, ...) {
  objects <- .joint_objects(x)
  n_models <- length(objects)
  old_par  <- par(mfrow = c(ceiling(n_models / 2), min(n_models, 2)),
                  mar   = c(4, 4, 3, 1))
  on.exit(par(old_par))

  for (nm in names(objects)) {
    tryCatch(
      plot(objects[[nm]], main = nm, ...),
      error = function(e) {
        # fallback: plot Tspace distribution of first tested coefficient
        ts <- objects[[nm]]$Tspace[, 1]
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
#' Combine jfs objects
#'
#' @param ... Objects to combine. Supported objects are \code{jfs},
#'   \code{flipscores}, \code{glm}, \code{fs_contrasts}, and
#'   \code{fs_combined}.
#'
#' @return A \code{jfs} object.
#' @method c jfs
#' @export
c.jfs <- function(..., recursive = FALSE) {
  dots <- list(...)
  if (length(dots) == 0) {
    return(NULL)
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }

  pieces <- lapply(seq_along(dots), function(i) {
    .as_jfs_component(dots[[i]], dot_names[i])
  })

  nrows <- vapply(pieces, function(x) nrow(as.matrix(x$Tspace)), integer(1))
  if (length(unique(nrows)) > 1) {
    stop("Cannot combine Tspace matrices with different numbers of rows.",
         call. = FALSE)
  }

  objects <- unlist(lapply(pieces, .joint_objects), recursive = FALSE)
  if (is.null(names(objects)) || any(!nzchar(names(objects)))) {
    names(objects) <- paste0("object", seq_along(objects))
  }

  out <- list(
    Tspace = do.call(cbind, lapply(pieces, `[[`, "Tspace")),
    summary_table = .rbind_fill(lapply(pieces, `[[`, "summary_table")),
    objects = objects,
    call = match.call()
  )
  class(out) <- "jfs"
  out
}


#--------------------------------------------
#' Coerce jfs to data frame
#'
#' @param x A \code{jfs} object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with all summary tables stacked, with a
#'   \code{Model} column prepended.
#' @export
as.data.frame.jfs <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  x$summary_table
}
