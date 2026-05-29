# utils_jointest.R
# Embedded and adapted from jointest (https://github.com/livioivil/jointest)
# Original author: Livio Finos

#--------------------------------------------
# Set model names from their response variables
#--------------------------------------------
.set_mods_names <- function(mods) {
  nms <- names(mods)
  if (is.null(nms) || any(nms == "")) {
    nms <- sapply(mods, function(mod) {
      resp <- tryCatch(
        as.character(formula(mod)[[2]]),
        error = function(e) NULL
      )
      if (is.null(resp))
        paste0("mod", which(sapply(mods, identical, mod)))
      else
        resp
    })
  }
  nms
}

#--------------------------------------------
# Get all coefficient names for each model as a list
#--------------------------------------------
.get_all_coeff_names_list <- function(mods) {
  lapply(mods, function(mod) names(coefficients(mod)))
}

#--------------------------------------------
# Bind Tspace from all models column-wise
#--------------------------------------------
.get_all_Tspace <- function(mods) {
  Tspaces <- lapply(mods, function(mod) mod$Tspace)
  do.call(cbind, Tspaces)
}

#--------------------------------------------
# Build summary table from a single flipscores object
# matches jointest column structure:
# model, .assign, coefficient, estimate, score, se, z, pcor, p
#--------------------------------------------
.get_summary_table_from_flipscores <- function(x, model_name = NULL) {
  pvals  <- x$p.values
  if (is.null(pvals)) return(NULL)

  # observed test statistics (last row of Tspace = observed)
  Tobs <- sapply(seq_along(pvals), function(i) {
    x$Tspace[nrow(x$Tspace), i]
  })

  # coefficient estimates
  estimates <- tryCatch(x$coefficients[names(pvals)], error = function(e) rep(NA, length(pvals)))

  # standard errors and z values from scores attributes
  sd_scores <- attributes(x$scores)$sd
  if (is.null(sd_scores)) sd_scores <- rep(NA, length(pvals))

  z_vals <- Tobs / sd_scores

  # partial correlations
  nrm <- attributes(x$scores)$nrm
  if (is.null(nrm)) nrm <- rep(NA, length(pvals))
  pcor <- unlist(colSums(x$scores)) / nrm

  # .assign: term assignment from model matrix
  assign_vec <- tryCatch(
    attr(model.matrix(x), "assign"),
    error = function(e) rep(NA, length(pvals))
  )
  names(assign_vec) <- colnames(model.matrix(x))
  assign_vals <- assign_vec[names(pvals)]

  data.frame(
    model       = if (!is.null(model_name)) model_name else NA_character_,
    .assign     = assign_vals,
    coefficient = names(pvals),
    estimate    = estimates,
    score       = unlist(colSums(x$scores)),
    se          = sd_scores,
    z           = z_vals,
    pcor        = pcor,
    p           = pvals,
    row.names   = NULL,
    stringsAsFactors = FALSE
  )
}

#--------------------------------------------
# Bind summary tables from all models
#--------------------------------------------
.get_all_summary_table <- function(mods) {
  tabs <- lapply(names(mods), function(nm) {
    tab <- mods[[nm]]$summary_table
    if (is.null(tab)) return(NULL)
    tab$model <- nm
    tab
  })
  do.call(rbind, tabs)
}

#--------------------------------------------
# Bind summary tables from combined results
# (output of .npc2jointest list)
#--------------------------------------------
.get_all_summary_table_combined <- function(res) {
  tabs <- lapply(names(res), function(nm) {
    tab <- res[[nm]]$summary_table
    if (is.null(tab)) return(NULL)
    tab
  })
  do.call(rbind, tabs)
}
#' Print method for joint_flipscores
#'
#' @param x A \code{joint_flipscores} object.
#' @param n Number of rows to show at head and tail of summary table.
#' @param ... Additional arguments (currently unused).
#' @export
print.joint_flipscores <- function(x, n = 2, ...) {
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  msg <- "== Joining n = %s models"
  cat(sprintf(msg, length(unique(x$summary_table$model))))
  cat("\n\n")
  .trim(x$summary_table, n = n)
  invisible(x)
}

#' Summary method for joint_flipscores
#'
#' @param object A \code{joint_flipscores} object.
#' @param digits Number of digits to print. Default \code{4}.
#' @param ... Additional arguments (currently unused).
#' @export
summary.joint_flipscores <- function(object, digits = 4, ...) {
  cat("\nCall: ")
  print(object$call)
  cat("\n")
  tab <- object$summary_table
  tab$.assign <- NULL  # remove internal column
  print(tab, digits = digits)
  invisible(object)
}
