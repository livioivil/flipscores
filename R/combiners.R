# combiners.R
# Embedded and adapted from jointest (https://github.com/livioivil/jointest)
# Original author: Livio Finos

#--------------------------------------------
# Internal combiner - NOT exported
#--------------------------------------------
.combine <- function(x, combine_with = "Fisher") {
  if (is.function(combine_with)) {
    return(combine_with(x))
  }

  combine_with <- match.arg(combine_with,
                            c("Fisher", "Liptak", "Tippett"))

  switch(combine_with,
    Fisher = {
      T_obs <- -2 * sum(log(x))
      list(
        T_obs   = T_obs,
        p.value = stats::pchisq(T_obs, df = 2 * length(x),
                                lower.tail = FALSE)
      )
    },
    Liptak = {
      T_obs <- sum(stats::qnorm(1 - x))
      list(
        T_obs   = T_obs,
        p.value = stats::pnorm(T_obs,
                               mean       = 0,
                               sd         = sqrt(length(x)),
                               lower.tail = FALSE)
      )
    },
    Tippett = {
      T_obs <- min(x)
      list(
        T_obs   = T_obs,
        p.value = 1 - (1 - T_obs)^length(x)
      )
    }
  )
}


#--------------------------------------------
#' Combine Tests
#'
#' Combines p-values from a \code{joint_flipscores} object or a numeric
#' vector using a combining function.
#'
#' @param x A \code{joint_flipscores} object or a numeric vector of p-values.
#' @param combine_with A string specifying the combining method:
#'   \code{"Fisher"} (default), \code{"Liptak"}, \code{"Tippett"},
#'   or a custom function.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with the combined test statistic and p-value.
#' @export
combine_tests <- function(x, combine_with = "Fisher", ...) {
  UseMethod("combine_tests")
}

#' @export
combine_tests.default <- function(x, combine_with = "Fisher", ...) {
  if (!is.numeric(x)) stop("`x` must be a numeric vector of p-values.")
  result <- .combine(x, combine_with = combine_with)
  data.frame(
    combine_with = if (is.function(combine_with)) "custom" else combine_with,
    T_obs        = result$T_obs,
    p.value      = result$p.value,
    stringsAsFactors = FALSE
  )
}

#' @export
combine_tests.joint_flipscores <- function(x, combine_with = "Fisher", ...) {
  # extract overall p-value per model from summary_table
  pvals <- sapply(x$mods, function(mod) {
    tab <- mod$summary_table
    if (is.null(tab)) return(NA)
    # combine within-model p-values first using Fisher
    # then combine across models
    .combine(tab$p.value, combine_with = combine_with)$p.value
  })
  pvals <- pvals[!is.na(pvals)]
  combine_tests.default(pvals, combine_with = combine_with, ...)
}


#--------------------------------------------
#' Combine Contrasts
#'
#' Combines p-values across models for each tested coefficient
#' from a \code{joint_flipscores} object.
#'
#' @param x A \code{joint_flipscores} object.
#' @param combine_with A string specifying the combining method:
#'   \code{"Fisher"} (default), \code{"Liptak"}, \code{"Tippett"},
#'   or a custom function.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with one row per coefficient, with combined
#'   test statistic and p-value.
#' @export
combine_contrasts <- function(x, combine_with = "Fisher", ...) {
  UseMethod("combine_contrasts")
}

#' @export
combine_contrasts.joint_flipscores <- function(x, combine_with = "Fisher", ...) {
  # gather all summary tables
  tabs <- lapply(names(x$mods), function(nm) {
    tab <- x$mods[[nm]]$summary_table
    if (is.null(tab)) return(NULL)
    tab$Model <- nm
    tab
  })
  all_tabs <- do.call(rbind, tabs)

  if (is.null(all_tabs)) stop("No summary tables found in joint_flipscores object.")

  # get unique coefficients across all models
  all_coeffs <- unique(all_tabs$Coeff)

  results <- lapply(all_coeffs, function(coeff) {
    pvals <- all_tabs$p.value[all_tabs$Coeff == coeff]
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) == 0) return(NULL)

    res <- .combine(pvals, combine_with = combine_with)
    data.frame(
      Coeff        = coeff,
      combine_with = if (is.function(combine_with)) "custom" else combine_with,
      n_models     = length(pvals),
      T_obs        = res$T_obs,
      p.value      = res$p.value,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}


#--------------------------------------------
#' p.adjust method for joint_flipscores
#'
#' Adjusts p-values in a \code{joint_flipscores} object for
#' multiple comparisons.
#'
#' @param p A \code{joint_flipscores} object.
#' @param method Correction method. See \code{\link[stats]{p.adjust.methods}}.
#' @param n Number of comparisons. See \code{\link[stats]{p.adjust}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The \code{joint_flipscores} object with adjusted p-values
#'   added as \code{p.value.adj} column in each model's summary table
#'   and in the overall summary table.
#' @export
p.adjust.joint_flipscores <- function(p, method = "BH", n = length(p), ...) {
  # adjust per-model summary tables
  p$mods <- lapply(p$mods, function(mod) {
    if (!is.null(mod$summary_table) &&
        "p.value" %in% colnames(mod$summary_table)) {
      mod$summary_table$p.value.adj <- stats::p.adjust(
        mod$summary_table$p.value,
        method = method
      )
    }
    mod
  })

  # adjust overall summary table
  if (!is.null(p$summary_table) &&
      "p.value" %in% colnames(p$summary_table)) {
    p$summary_table$p.value.adj <- stats::p.adjust(
      p$summary_table$p.value,
      method = method
    )
  }

  p
}
