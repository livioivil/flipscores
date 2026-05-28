# flipscores.R
# Extended version of flipscores() with dispatch for:
#   1. formula        -> original behavior
#   2. list of glms   -> join_flipscores()
#   3. list of formulas -> convert to glms, then join_flipscores()
#   4. matrix response formula -> convert to list of formulas, then case 3)

#' Flip Scores Test
#'
#' Extended version of \code{flipscores()} supporting joint testing via
#' \code{join_flipscores()} when \code{formula} is a list of models/formulas
#' or a formula with matrix response.
#'
#' @param formula A \code{formula}, a list of \code{glm} objects, a list of
#'   \code{formula} objects, or a \code{formula} with a matrix response.
#' @param family A description of the error distribution. Used when
#'   \code{formula} is a single formula or a list of formulas.
#' @param data A data frame. Used when \code{formula} is a single formula
#'   or a list of formulas.
#' @param score_type Type of score to use. One of \code{"standardized"}
#'   (default), \code{"sign"}, \code{"raw"}, \code{"Williams"}.
#' @param tested_terms A character vector of terms to test. If \code{NULL},
#'   all terms are tested.
#' @param n_flips Number of random flips (sign changes) for the permutation
#'   distribution. Default is \code{1000}.
#' @param combine_with Combining method for joint tests: \code{"Fisher"}
#'   (default), \code{"Liptak"}, \code{"Tippett"}, or a custom function.
#'   Only used in cases 2), 3), 4).
#' @param seed Random seed for reproducibility. Only used in cases 2), 3), 4).
#' @param ... Additional arguments passed to \code{glm()} or the internal
#'   \code{flipscores} engine.
#'
#' @return A \code{flipscores} object (case 1) or a \code{joint_flipscores}
#'   object (cases 2, 3, 4).
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' df <- data.frame(
#'   y1 = rnorm(n),
#'   y2 = rnorm(n),
#'   y3 = rpois(n, lambda = 2),
#'   yb = rbinom(n, 1, 0.5),
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   grp = factor(rep(c("A", "B"), n / 2))
#' )
#'
#' # --- Case 2: list of glm objects ---
#' # Models can have different families and responses
#' models <- list(
#'   glm(y1 ~ x1 + x2, data = df, family = gaussian()),
#'   glm(y2 ~ x1 + x2, data = df, family = gaussian()),
#'   glm(yb ~ x1 + x2, data = df, family = binomial())
#' )
#' flipscores(models)
#'
#' # With different tested_terms
#' flipscores(models, tested_terms = "x1")
#'
#' # --- Case 3: list of formulas, same family ---
#' formulas_gauss <- list(
#'   y1 ~ x1 + x2,
#'   y2 ~ x1 + x2
#' )
#' flipscores(formulas_gauss, family = gaussian(), data = df)
#'
#' # List of formulas, different families not directly supported in case 3)
#' # use case 2) instead (pre-fit glm objects with desired families)
#'
#' # --- Case 4: matrix response, gaussian ---
#' flipscores(cbind(y1, y2) ~ x1 + x2, family = gaussian(), data = df)
#'
#' # Case 4: matrix response, poisson
#' # (note: y3 is count data, y1/y2 are rounded to non-negative integers)
#' df$y1p <- rpois(n, lambda = 3)
#' df$y2p <- rpois(n, lambda = 5)
#' flipscores(cbind(y1p, y2p) ~ x1 + x2, family = poisson(), data = df)
flipscores <- function(formula,
                       family       = gaussian(),
                       data         = NULL,
                       score_type   = "standardized",
                       tested_terms = NULL,
                       n_flips      = 1000,
                       combine_with = "Fisher",
                       seed         = NULL,
                       ...) {

  ##############################################################
  # CASE 2: list of glm objects
  ##############################################################
  if (is.list(formula) && all(sapply(formula, inherits, "glm"))) {

    message("flipscores: list of glm objects detected -> joint test")
    return(
      join_flipscores(
        models       = formula,
        combine_with = combine_with,
        score_type   = score_type,
        tested_terms = tested_terms,
        n_flips      = n_flips,
        seed         = seed,
        ...
      )
    )
  }

  ########################
  ##############################################################
  # CASE 3: list of formulas
  ##############################################################
  if (is.list(formula) && all(sapply(formula, inherits, "formula"))) {

    message("flipscores: list of formulas detected -> converting to glms")

    models <- lapply(formula, function(f) {
      # Explicitly pass the formula as a string to avoid scoping issues
      stats::glm(
        formula = stats::as.formula(deparse(f)),
        family  = family,
        data    = data,
        ...
      )
    })

    return(
      join_flipscores(
        models       = models,
        combine_with = combine_with,
        score_type   = score_type,
        tested_terms = tested_terms,
        n_flips      = n_flips,
        seed         = seed
      )
    )
  }

  ##############################################################
  # CASE 4: formula with matrix response -> list of formulas
  ##############################################################
  if (inherits(formula, "formula")) {
    lhs <- formula[[2]]
    rhs <- formula[[3]]

    is_matrix_response <- .is_matrix_lhs(lhs, data = data,
                                         parent_env = parent.frame())

    if (is_matrix_response) {
      message("flipscores: matrix response detected -> converting to list of formulas")

      resp_names <- .extract_matrix_response_names(lhs)

      # Fix: build formulas as strings and reconstruct them
      # to avoid scoping issues in downstream lapply calls
      rhs_str <- deparse(rhs)
      formulas <- lapply(resp_names, function(yn) {
        stats::as.formula(
          paste0(yn, " ~ ", rhs_str)
          # note: no env argument, formula will be evaluated
          # in glm() with the `data` argument
        )
      })

      # Recurse as case 3) - which now also uses deparse() fix
      return(
        flipscores(
          formula      = formulas,
          family       = family,
          data         = data,
          score_type   = score_type,
          tested_terms = tested_terms,
          n_flips      = n_flips,
          combine_with = combine_with,
          seed         = seed,
          ...
        )
      )
    }
  }

  ##############################################################
  # CASE 1: standard formula -> original flipscores engine
  ##############################################################
  .flipscores_engine(
    formula      = formula,
    family       = family,
    data         = data,
    score_type   = score_type,
    tested_terms = tested_terms,
    n_flips      = n_flips,
    ...
  )
}


#--------------------------------------------
# Helper: detect matrix LHS in formula
#--------------------------------------------
.is_matrix_lhs <- function(lhs, data = NULL, parent_env = parent.frame()) {
  # Case: cbind(y1, y2, ...)
  if (is.call(lhs) && deparse(lhs[[1]]) == "cbind") return(TRUE)

  # Case: bare symbol that resolves to a matrix
  if (is.symbol(lhs)) {
    obj <- tryCatch(
      eval(lhs, envir = data, enclos = parent_env),
      error = function(e) NULL
    )
    if (!is.null(obj) && is.matrix(obj) && ncol(obj) > 1) return(TRUE)
  }

  FALSE
}

#--------------------------------------------
# Helper: extract response variable names from matrix LHS
#--------------------------------------------
.extract_matrix_response_names <- function(lhs) {
  # Case: cbind(y1, y2, ...)
  if (is.call(lhs) && deparse(lhs[[1]]) == "cbind") {
    # Fix: use deparse() on each element to get clean variable name strings
    # instead of returning language objects that cause scoping issues
    return(sapply(as.list(lhs)[-1], function(x) deparse(x)))
  }

  stop(
    "Matrix response must be specified as cbind(y1, y2, ...) in the formula. ",
    "Bare matrix variables are not supported."
  )
}
