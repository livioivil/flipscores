# combiners.R
# Embedded and adapted from jointest (https://github.com/livioivil/jointest)
# Original author: Livio Finos

#--------------------------------------------
# internal NPC combining functions
#--------------------------------------------
.comb_funct_fisher <- function(p) {
  T_obs <- -2 * sum(log(p))
  list(
    T_obs   = T_obs,
    p.value = stats::pchisq(T_obs, df = 2 * length(p), lower.tail = FALSE)
  )
}

.comb_funct_liptak <- function(p) {
  T_obs <- sum(stats::qnorm(1 - p))
  list(
    T_obs   = T_obs,
    p.value = stats::pnorm(T_obs, mean = 0,
                           sd = sqrt(length(p)), lower.tail = FALSE)
  )
}

.comb_funct_tippett <- function(p) {
  T_obs <- min(p)
  list(
    T_obs   = T_obs,
    p.value = 1 - (1 - T_obs)^length(p)
  )
}

.comb_funct_maxT <- function(Tspace) {
  # Tspace: matrix of permutation distribution (rows=perms, cols=tests)
  # combine by taking the maximum
  apply(Tspace, 1, max)
}

.comb_funct_mahalanobis <- function(Tspace) {
  # Tspace: matrix of permutation distribution
  mahalanobis_npc(Tspace)
}

#--------------------------------------------
# internal: run NPC on a subset of columns of Tspace
#--------------------------------------------
.npc2jointest <- function(id, mods, combined, tail = 0,
                          comb_funct = "maxT") {
  ids   <- combined[[id]]
  Tsub  <- mods$Tspace[, ids, drop = FALSE]

  if (is.character(comb_funct)) {
    comb_funct <- match.arg(comb_funct,
                            c("maxT", "Mahalanobis", "Fisher",
                              "Liptak", "Tippett"))
    Tcomb <- switch(comb_funct,
                    maxT        = apply(Tsub, 1, max),
                    Mahalanobis = mahalanobis_npc(Tsub),
                    Fisher      = apply(Tsub, 1, function(r) -2 * sum(log(.T2p(r, tail)))),
                    Liptak      = apply(Tsub, 1, function(r)
                      sum(stats::qnorm(1 - .T2p(r, tail)))),
                    Tippett     = apply(Tsub, 1, function(r) min(.T2p(r, tail)))
    )
  } else if (is.function(comb_funct)) {
    Tcomb <- apply(Tsub, 1, comb_funct)
  } else {
    stop("comb_funct must be a string or a function.")
  }

  # p-value: proportion of permutations >= observed (first row)
  T_obs <- Tcomb[1]
  p_val <- mean(Tcomb >= T_obs)

  # build summary row
  #smr <- mods$summary_table[ids[1], , drop = FALSE]
  #smr$score   <- T_obs
  #smr$p       <- p_val

    comb_name = names(combined)[id]
    if (is.null(comb_name))
      comb_name = "combined"
    colnames(Tcomb) = comb_name


    Coeff=mods$summary_table[combined[[id]],"coefficient"]
    Coeff=unique(Coeff)
    if(length(Coeff) > 1) Coeff="many"
    smr = data.frame(coefficient = Coeff,
                     stat = comb_funct,
                     ntests = max(1, length(combined[[id]])),
                     S = Tspace[1],
                     p = .t2p_only_first(Tspace, tail = 1))

  list(
    Tspace        = matrix(Tcomb, ncol = 1,
                           dimnames = list(NULL, names(combined)[id])),
    summary_table = smr
  )
}

# internal helper: T to p-value
.T2p <- function(T_obs, tail = 0) {
  if (tail == 0)  return(2 * pmin(mean(T_obs >= T_obs[1]),
                                  mean(T_obs <= T_obs[1])))
  if (tail > 0)   return(mean(T_obs >= T_obs[1]))
  if (tail < 0)   return(mean(T_obs <= T_obs[1]))
}

#--------------------------------------------
#' Combine Tests
#'
#' Combines tests from a \code{joint_flipscores} object using a
#' non-parametric combination (NPC) method.
#'
#' @param mods A \code{joint_flipscores} object.
#' @param comb_funct Combining function: \code{"maxT"} (default),
#'   \code{"Mahalanobis"}, \code{"Fisher"}, \code{"Liptak"},
#'   \code{"Tippett"}, or a custom function.
#' @param by Character vector of column names in \code{summary_table}
#'   to group by before combining. If \code{NULL}, all tests are combined.
#' @param by_list A named list of integer vectors specifying custom groupings
#'   of column indices in \code{Tspace}. Overrides \code{by} if provided.
#' @param tail Tail direction: \code{0} = two-sided (default),
#'   \code{1} = right, \code{-1} = left.
#' @return A \code{joint_flipscores} object with combined results.
#' @export
#' @examples
#' # flipscores jointly on all models and all coefficients
#' set.seed(1)
#' n <- 50
#' df <- data.frame(
#'   y1  = rnorm(n),
#'   y2  = rnorm(n),
#'   yb  = rbinom(n, 1, 0.5),
#'   y1p = rpois(n, lambda = 3),
#'   y2p = rpois(n, lambda = 5),
#'   x1  = rnorm(n),
#'   x2  = rnorm(n)
#' )
#'
#' res=flipscores(
#'   list(y1 ~ x1 + x2, y2 ~ x1 + x2),
#'   family = gaussian(),
#'   data   = df
#' )
#' summary(res)
#' summary(combine_tests(res))
#' summary(combine_tests(res, by="model"))
#' summary(combine_tests(res, by="coefficient"))
#' res2=combine_contrasts(res)
#' summary(res2)
#' #custom combinations:
#' coeffs=c("(Intercept)","X","Z1","Z2")
#' coeffs_ids=lapply(coeffs,grep,res2$summary_table$Coeff)
#' names(coeffs_ids)=coeffs
#' summary(combine_tests(res2,by_list =   coeffs_ids))
combine_tests <- function(mods, comb_funct = "maxT",
                          by = NULL, by_list = NULL, tail = 0) {
  UseMethod("combine_tests")
}

#' @export
combine_tests.joint_flipscores <- function(mods, comb_funct = "maxT",
                                           by = NULL, by_list = NULL,
                                           tail = 0) {
  if (!is.null(by_list)) {
    combined <- by_list
  } else {
    if (is.null(by)) {
      # combine all tests overall
      combined <- list(overall = seq_len(ncol(mods$Tspace)))
    } else {
      # group by columns in summary_table
      smr      <- apply(mods$summary_table[, by, drop = FALSE],
                        1, paste, collapse = ".")
      uniq_nm  <- unique(smr)
      combined <- lapply(uniq_nm, function(nm) which(smr == nm))
      names(combined) <- uniq_nm
    }
  }

  res <- lapply(seq_along(combined), .npc2jointest,
                mods = mods, combined = combined,
                tail = tail, comb_funct = comb_funct)
  names(res) <- names(combined)

  out <- list(
    Tspace        = .get_all_Tspace(res),
    summary_table = .get_all_summary_table_combined(res)
  )
  class(out) <- c("joint_flipscores", class(out))
  out
}

#--------------------------------------------
#' Combine Contrasts
#'
#' Combines tests across contrasts of factor variables to obtain
#' a global test per factor (analogous to ANOVA).
#'
#' @param mods A \code{joint_flipscores} object.
#' @param comb_funct Combining function. Default \code{"Mahalanobis"}.
#' @param tail Tail direction. Default \code{0} (two-sided).
#' @return A \code{joint_flipscores} object with combined results.
#' @export
combine_contrasts <- function(mods, comb_funct = "Mahalanobis", tail = 0) {
  UseMethod("combine_contrasts")
}
