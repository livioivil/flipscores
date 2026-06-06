#' Flip Scores Test
#'
#' Extended version of \code{flipscores()} supporting joint testing when
#' \code{formula} is a list of glm objects, a list of formulas, or a formula
#' with a matrix response.
#'
#' @param formula One of:
#'   \enumerate{
#'     \item A \code{formula} object (standard use, returns \code{flipscores} object)
#'     \item A list of fitted \code{glm} objects (returns \code{jfs} object)
#'     \item A list of \code{formula} objects (returns \code{jfs} object)
#'     \item A \code{formula} with matrix response e.g. \code{cbind(y1,y2) ~ x}
#'       (returns \code{jfs} object)
#'   }
#' @param family Error distribution for \code{glm}. Used in cases 1, 3, 4.
#' @param data A data frame. Used in cases 1, 3, 4.
#' @param score_type Type of score. One of \code{"standardized"} (default),
#'   \code{"orthogonalized"}, \code{"effective"}, \code{"basic"}.
#' @param n_flips Number of random flips. Default \code{5000}.
#' @param alternative Alternative hypothesis. Default \code{"two.sided"}.
#' @param id Optional cluster id vector.
#' @param seed Random seed for reproducibility.
#' @param to_be_tested Character vector or list of coefficient names to test.
#'   If a character vector, the same terms are intersected with each model's
#'   coefficients. If a list, each element applies to the corresponding model.
#'   If \code{NULL} all coefficients are tested.
#' @param tested_coeffs Alias for \code{to_be_tested} for compatibility.
#' @param flips Optional matrix of pre-computed flips.
#' @param precompute_flips Logical, default \code{TRUE}.
#' @param ... Additional arguments passed to \code{glm()} or the internal
#'   flipscores engine.
#'
#' @return A \code{flipscores} object (case 1) or a \code{jfs}
#'   object (cases 2, 3, 4).
#'
#' @examples
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
#' # Case 1: single formula
#' res=flipscores(y1 ~ x1 + x2, family = gaussian(), data = df)
#' summary(res)
#'
#' # Case 2: list of glm objects
#' models <- list(
#'   glm(y1 ~ x1 + x2, data = df, family = gaussian()),
#'   glm(y2 ~ x1 + x2, data = df, family = gaussian()),
#'   glm(yb ~ x1 + x2, data = df, family = binomial())
#' )
#' res=flipscores(models)
#' summary(res)
#' res=flipscores(models, to_be_tested = "x1")
#' summary(res)
#'
#' # Case 3: list of formulas, same family
#' res=flipscores(
#'   list(y1 ~ x1 + x2, y2 ~ x1 + x2),
#'   family = gaussian(),
#'   data   = df
#' )
#' summary(res)
#'
#' # Case 4: matrix response, gaussian
#' res=flipscores(cbind(y1, y2) ~ x1 + x2, family = gaussian(), data = df)
#' summary(res)
#'
#' # Case 4: matrix response, poisson
#' res=flipscores(cbind(y1p, y2p) ~ x1 + x2, family = poisson(), data = df)
#' summary(res)
#'
#' # more examples
#'n=20
#'DD=data.frame(X=rnorm(n),Z1=rnorm(n),Z2=rnorm(n))
#'DD$Y=DD$Z1+DD$X+rnorm(n)
#'# Run four glms abd combine it in a list
#'mod1=glm(Y~X+Z1+Z2,data=DD)
#'mod2=glm(Y~X+poly(Z1,2)+Z2,data=DD)
#'mod3=glm(Y~X+poly(Z1,2)+poly(Z2,2),data=DD)
#'mod4=glm(Y~X+Z1+poly(Z2,2),data=DD)
#'mods=list(mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)
#'res=flipscores(mods, to_be_tested = "X")
#'summary(res)
#'
#' @export
flipscores <- function(formula,
                       family           = gaussian(),
                       data             = NULL,
                       score_type       = "standardized",
                       n_flips          = 5000,
                       alternative      = "two.sided",
                       id               = NULL,
                       seed             = NULL,
                       to_be_tested     = NULL,
                       tested_coeffs    = NULL,
                       flips            = NULL,
                       precompute_flips = TRUE,
                       ...) {

  # tested_coeffs is an alias for to_be_tested
  if (!is.null(tested_coeffs) && is.null(to_be_tested))
    to_be_tested <- tested_coeffs

  ##############################################################
  # CASE 2: list of glm objects
  ##############################################################
  if (is.list(formula) && all(sapply(formula, inherits, "glm"))) {
    #message("flipscores: list of glm objects detected -> joint test")

    for(i in 1:length(formula))
      formula[[i]]$call$data=eval(formula[[i]]$call$data, parent.frame())

    out <- .join_flipscores(
      mods         = formula,
      tested_coeffs = to_be_tested,
      n_flips      = n_flips,
      flips        = flips,
      score_type   = score_type,
      seed         = seed,
      ...
    )
    out$call <- match.call()
    return(out)
  }

  ##############################################################
  # CASE 3: list of formulas
  ##############################################################
  if (is.list(formula) && all(sapply(formula, inherits, "formula"))) {
    #message("flipscores: list of formulas detected -> converting to glms")


    models <- lapply(formula, function(f) {
      do.call(
        stats::glm,
        list(
          formula = stats::as.formula(paste(deparse(f), collapse = " ")),
          family  = family,
          data    = data
        )
      )
    })

    out <- .join_flipscores(
      mods         = models,
      tested_coeffs = to_be_tested,
      n_flips      = n_flips,
      flips        = flips,
      score_type   = score_type,
      seed         = seed,
      ...
    )
    out$call <- match.call()
    return(out)
  }

  ##############################################################
  # CASE 4: formula with matrix response -> list of formulas
  ##############################################################
  if (inherits(formula, "formula")) {
    lhs <- formula[[2]]
    rhs <- formula[[3]]

    if (.is_matrix_lhs(lhs, data = data, parent_env = parent.frame())) {
     # message("flipscores: matrix response detected -> converting to list of formulas")

      original_call <- match.call()
      resp_names    <- .extract_matrix_response_names(lhs)
      rhs_str       <- paste(deparse(rhs), collapse = " ")

      mods <- lapply(resp_names, function(yn) {
        do.call(
          stats::glm,
          list(
            formula = stats::as.formula(paste0(yn, " ~ ", rhs_str)),
            family  = family,
            data    = data
          )
        )
      })

      out <- .join_flipscores(
        mods          = mods,
        score_type       = score_type,
        n_flips          = n_flips,
        alternative      = alternative,
        id               = id,
        seed             = seed,
        tested_coeffs     = to_be_tested,
        flips            = flips,
        precompute_flips = precompute_flips,
        ...
      )
      out$call <- original_call
      return(out)
    }
  }

  ##############################################################
  # CASE 1: standard formula -> original flipscores engine
  ##############################################################
  .flipscores_engine(
    formula          = formula,
    family           = family,
    data             = data,
    score_type       = score_type,
    n_flips          = n_flips,
    alternative      = alternative,
    id               = id,
    seed             = seed,
    to_be_tested     = to_be_tested,
    flips            = flips,
    precompute_flips = precompute_flips,
    .user_call       = match.call(),
    ...
  )
}


#--------------------------------------------
# Helper: detect matrix LHS in formula
#--------------------------------------------
.is_matrix_lhs <- function(lhs, data = NULL, parent_env = parent.frame()) {
  if (is.call(lhs) && deparse(lhs[[1]]) == "cbind") return(TRUE)
  if (is.symbol(lhs)) {
    obj <- tryCatch(
      eval(lhs, envir = as.list(data), enclos = parent_env),
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
  if (is.call(lhs) && deparse(lhs[[1]]) == "cbind") {
    return(sapply(as.list(lhs)[-1], deparse))
  }
  stop(
    "Matrix response must be specified as cbind(y1, y2, ...) in the formula. ",
    "Bare matrix variables are not supported."
  )
}


#--------------------------------------------
# .flipscores_engine: original flipscores() body
#--------------------------------------------
.flipscores_engine <- function(formula,
                               family,
                               data,
                               score_type       = "standardized",
                               n_flips          = 5000,
                               alternative      = "two.sided",
                               id               = NULL,
                               seed             = NULL,
                               to_be_tested     = NULL,
                               flips            = NULL,
                               precompute_flips = TRUE,
                               .user_call       = NULL,
                               ...) {
  mf <- match.call()
  mf$.user_call <- NULL
  fs_call <- if (is.null(.user_call)) mf else .user_call

  # force evaluation of all arguments immediately
  # to prevent match.call() storing unevaluated symbols
  # that fail when accessed as logicals later
  precompute_flips <- force(precompute_flips)
  alternative      <- force(alternative)
  score_type       <- force(score_type)
  n_flips          <- force(n_flips)
  seed             <- force(seed)
  id               <- force(id)

  # save original formula, data symbol and data value
  original_formula   <- formula
  original_formula_call <- fs_call$formula
  original_data      <- data
  original_data_call <- fs_call$data  # unevaluated symbol e.g. `df`

  score_type <- match.arg(score_type,
                          c("orthogonalized", "standardized",
                            "effective", "basic", "my_lab"))

  m <- match(c("score_type", "n_flips", "alternative", "id",
               "seed", "flips", "precompute_flips"), names(mf), 0L)
  m <- m[m > 0]
  flip_param_call       <- mf[c(1L, m)]
  flip_param_call[[1L]] <- .flip_test

  # evaluate all flip parameters explicitly
  flip_param_call$id               <- eval(flip_param_call$id,
                                           parent.frame())
  flip_param_call$alternative      <- eval(flip_param_call$alternative,
                                           parent.frame())
  flip_param_call$flips            <- eval(flip_param_call$flips,
                                           parent.frame())
  flip_param_call$precompute_flips <- eval(flip_param_call$precompute_flips,
                                           parent.frame())
  flip_param_call$seed             <- eval(flip_param_call$seed,
                                           parent.frame())
  flip_param_call$score_type       <- eval(flip_param_call$score_type,
                                           parent.frame())
  flip_param_call$n_flips          <- eval(flip_param_call$n_flips,
                                           parent.frame())
  flip_param_call$nobservations    <- eval(mf$nobservations,
                                           parent.frame())

  # set defaults for any NULL values
  if (is.null(flip_param_call$precompute_flips))
    flip_param_call$precompute_flips <- TRUE
  if (is.null(flip_param_call$score_type))
    flip_param_call$score_type <- "standardized"
  if (is.null(flip_param_call$n_flips))
    flip_param_call$n_flips <- 5000

  mf$nobservations <- NULL

  if (!is.null(list(...)$parms_DV)) {
    flip_param_call$parms_DV <- list(...)$parms_DV
    mf$parms_DV <- NULL
  }

  # handle to_be_tested
  m2 <- match("to_be_tested", names(mf), 0L)
  if (m2 == 0)
    to_be_tested <- NULL
  else {
    m <- c(m, m2)
    to_be_tested <- eval(mf[[m2]], parent.frame())
  }

  if (!is.null(flip_param_call$id) && (score_type == "orthogonalized")) {
    warning("Use of id is not possible with score_type=='orthogonalized', yet. Nothing done.")
    return(NULL)
  }

  # remove flip-specific params from mf, keep only glm params
  if (length(m) > 0) mf <- mf[-m]
  mf$offset <- eval(mf$offset, parent.frame())

  # use formula directly (already evaluated as function argument)
  model <- formula
  mf$formula <- model

  if (inherits(model, "formula")) {
    if (is.character(family) && family == "negbinom") {
      mf[[1L]] <- quote(MASS::glm.nb)
      mf$family <- NULL
    } else {
      mf[[1L]] <- quote(glm)
      mf$family <- family
    }
    param_x_ORIGINAL <- mf$x
    mf$x    <- TRUE
    model   <- eval(mf, parent.frame())
    if (!is.character(family) || family != "negbinom")
      model$call$family <- family
  } else if (inherits(model, "glm")) {
    param_x_ORIGINAL <- TRUE
    # refit using do.call to avoid any symbol lookup issues
    # extract all arguments from the existing model
    # refit_args        <- list(
    #   formula = formula(model),
    #   family  = model$family,
    #   data    = model$model,  # use stored model frame directly
    #   x       = TRUE
    # )
    # # preserve any extra arguments (offset, weights, subset etc)
    # extra_args <- model$call[!names(model$call) %in%
    #                            c("", "formula", "family", "data", "x")]
    # extra_args[[1]] <- NULL  # remove function name
    # if (length(extra_args) > 0) {
    #   for (nm in names(extra_args))
    #     refit_args[[nm]] <- eval(extra_args[[nm]],
    #                              envir = model$model)
    # }
    # model <- do.call(stats::glm, refit_args)
    # model$call$family <- model$family
    model <- update(model,x=TRUE)

  } else {
    stop("'formula' must be a formula or a glm object.")
  }

  if (is.null(model$y)) model$y <- model$model[, 1]

  # resolve to_be_tested
  if (is.null(to_be_tested))
    to_be_tested <- colnames(model[["x"]])
  else {
    if (is.numeric(to_be_tested))
      to_be_tested <- colnames(model[["x"]])[to_be_tested]
    to_be_tested <- eval(to_be_tested, parent.frame())
  }

  # handle flips
  if (!is.null(flip_param_call$flips)) {
    flip_param_call$precompute_flips <- FALSE
    flip_param_call$n_flips <- nrow(flip_param_call$flips)
  } else if (flip_param_call$precompute_flips) {
    set.seed(seed)
    flip_param_call$flips <- .make_flips(
      max(nrow(model$model),
          ifelse(is.null(flip_param_call$nobservations),
                 0L,
                 flip_param_call$nobservations)),
      flip_param_call$n_flips,
      flip_param_call$id
    )
  }

  # compute scores for each tested coefficient
  results <- lapply(to_be_tested,
                    socket_compute_scores_and_flip,
                    model,
                    flip_param_call = flip_param_call)

  model$scores   <- data.frame(lapply(results, function(x) x[[1]]$scores))
  nrm            <- sapply(results,
                           function(x) attributes(x[[1]]$scores)$scale_objects$nrm)
  std_dev        <- sapply(results,
                           function(x) attributes(x[[1]]$scores)$sd)
  model$Tspace   <- data.frame(lapply(results, function(x) x[[1]]$Tspace))
  model$p.values <- sapply(results, function(x) x[[1]]$p.values)

  flip_param_call$flips <- NULL

  attr(model$scores, "nrm")       <- nrm
  attr(model$scores, "sd")        <- std_dev
  attr(model$scores, "resid_std") <- data.frame(
    lapply(results, function(x) attr(x[[1]]$scores, "resid_std"))
  )

  names(attributes(model$scores)$resid_std) <-
    names(nrm)            <-
    names(std_dev)        <-
    names(model$scores)   <-
    names(model$Tspace)   <-
    names(model$p.values) <- to_be_tested

  # build a clean glm call for update() - only valid glm() arguments
  glm_call          <- call("glm")
  glm_call$formula  <- original_formula
  glm_call$family   <- family
  glm_call$data     <- original_data      # actual data frame
  if (!is.null(mf$offset))  glm_call$offset  <- mf$offset
  if (!is.null(mf$weights)) glm_call$weights <- mf$weights
  if (!is.null(mf$subset))  glm_call$subset  <- mf$subset

  # build the user-facing flipscores call for print/summary/update
  fs_call[[1L]]         <- quote(flipscores)
  fs_call$formula       <- original_formula_call
  fs_call$data          <- original_data_call
  fs_call$nobservations <- NULL
  fs_call$parms_DV      <- NULL

  model$call            <- glm_call
  model$flipscores_call <- fs_call
  model$flip_param_call <- flip_param_call
  model$score_type      <- score_type
  model$summary_table=.get_summary_table_from_flipscores(model)

  if (is.null(param_x_ORIGINAL) || (!param_x_ORIGINAL)) model$x <- NULL

  class(model) <- c("flipscores", class(model))
  return(model)
}
