#' Flip Scores Test
#'
#' Extended version of \code{flipscores()} supporting joint testing when
#' \code{formula} is a list of glm objects, a list of formulas, or a formula
#' with a matrix response.
#'
#' @param formula One of:
#'   \enumerate{
#'     \item A \code{formula} object (standard use, returns \code{flipscores} object)
#'     \item A list of fitted \code{glm} objects (returns \code{joint_flipscores} object)
#'     \item A list of \code{formula} objects (returns \code{joint_flipscores} object)
#'     \item A \code{formula} with matrix response e.g. \code{cbind(y1,y2) ~ x}
#'       (returns \code{joint_flipscores} object)
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
#' @param tested_coeffs an alias for  \code{to_be_tested}. \code{tested_coeffs} overcomes \code{to_be_tested} if it is not \code{NULL}.
#' @param flips Optional matrix of pre-computed flips. If provided, shared
#'   across all models in cases 2, 3, 4.
#' @param precompute_flips Logical, default \code{TRUE}.
#' @param ... Additional arguments passed to \code{glm()} or the internal
#'   flipscores engine.
#'
#' @return A \code{flipscores} object (case 1) or a \code{joint_flipscores}
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
#'
#'# more examples
#'n=20
#'D=data.frame(X=rnorm(n),Z1=rnorm(n),Z2=rnorm(n))
#'D$Y=D$Z1+D$X+rnorm(n)
#'# Run four glms abd combine it in a list
#'mod1=glm(Y~X+Z1+Z2,data=D)
#'mod2=glm(Y~X+poly(Z1,2)+Z2,data=D)
#'mod3=glm(Y~X+poly(Z1,2)+poly(Z2,2),data=D)
#'mod4=glm(Y~X+Z1+poly(Z2,2),data=D)
#'mods=list(mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)
#'res=flipscores(mods, to_be_tested = "X")
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
                       tested_coeffs    = NULL,  # alias for to_be_tested
                       flips            = NULL,
                       precompute_flips = TRUE,
                       ...) {

  # handle tested_coeffs as alias for to_be_tested
  if (!is.null(tested_coeffs) )
    to_be_tested <- tested_coeffs

  # ... rest unchanged ...

  # capture the unevaluated data expression HERE
  # before it gets renamed to `data` inside .flipscores_engine()
  # this preserves the original symbol name (e.g. `df`) for display
  data_call <- match.call()$data

  ##############################################################
  # CASE 2: list of glm objects
  ##############################################################
  if (is.list(formula) && all(sapply(formula, inherits, "glm"))) {
   # message("flipscores: list of glm objects detected -> joint test")

    # capture caller environment BEFORE entering .join_flipscores
    # so that data references in glm calls are resolved correctly
    caller_env <- parent.frame()

    out <- .join_flipscores(
      mods         = formula,
      to_be_tested = to_be_tested,
      n_flips      = n_flips,
      flips        = flips,
      score_type   = score_type,
      seed         = seed,
      caller_env   = caller_env,
      ...
    )
    out$call <- match.call()
    return(out)
  }

  ##############################################################
  # CASE 3: list of formulas
  ##############################################################
  if (is.list(formula) && all(sapply(formula, inherits, "formula"))) {
    message("flipscores: list of formulas detected -> converting to glms")

    caller_env <- parent.frame()
    .data   <- data
    .family <- family

    models <- lapply(formula, function(f) {
      do.call(stats::glm, list(
        formula = stats::as.formula(paste(deparse(f), collapse = " ")),
        family  = .family,
        data    = .data
      ))
    })

    out <- .join_flipscores(
      mods         = models,
      to_be_tested = to_be_tested,
      n_flips      = n_flips,
      flips        = flips,
      score_type   = score_type,
      seed         = seed,
      caller_env   = caller_env,
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
      #message("flipscores: matrix response detected -> converting to list of formulas")

      # capture original call BEFORE recursing
      original_call <- match.call()

      resp_names <- .extract_matrix_response_names(lhs)
      rhs_str    <- paste(deparse(rhs), collapse = " ")
      formulas   <- lapply(resp_names, function(yn) {
        stats::as.formula(paste0(yn, " ~ ", rhs_str))
      })

      out <- flipscores(
        formula          = formulas,
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
        ...
      )
      out$call <- original_call  # overwrite with the true original call
      return(out)
    }
  }

  ##############################################################
  # CASE 1: standard formula -> original flipscores engine
  ##############################################################
  result <- .flipscores_engine(
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
    ...
  )

  # patch the data slot in flipscores_call with the original
  # unevaluated expression (e.g. `df` instead of the full data frame)
  if (!is.null(result$flipscores_call) && !is.null(data_call)) {
    result$flipscores_call$data <- data_call
  }

  result
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
# renamed to avoid infinite dispatch loop
#--------------------------------------------
.join_flipscores <- function(mods,
                             to_be_tested  = NULL,
                             n_flips       = 5000,
                             flips         = NULL,
                             score_type    = "standardized",
                             statistics    = "t",
                             seed          = NULL,
                             output_models = TRUE,
                             caller_env    = parent.frame(),
                             ...) {

  if (!is.null(seed)) set.seed(seed)

  # DO NOT modify mods[[i]]$call$data
  # instead pass data directly when calling .flipscores_engine
  mods <- lapply(seq_along(mods), function(i) {
    # extract data directly from the fitted model
    model_data <- mods[[i]]$model

    temp <- .flipscores_engine(
      formula       = mods[[i]],
      family        = mods[[i]]$family,
      data          = model_data,    # pass extracted data explicitly
      score_type    = score_type,
      flips         = FLIPS,
      to_be_tested  = to_be_tested[[i]],
      nobservations = n_obs,
      ...
    )
    if (statistics %in% "t") {
      temp$summary_table <- .get_summary_table_from_flipscores(temp,
                                                               model_name = names(mods)[i])
    }
    temp
  })
  names(mods) <- .set_mods_names(mods)

  # ... rest of function unchanged ...

  fs_call <- mf <- match.call()

  # save BEFORE anything modifies formula or data
  original_formula <- formula
  original_data    <- data
  # save the UNEVALUATED data argument for display purposes
  # this preserves the symbol name (e.g. `df`) if user passed a named object
  original_data_call <- fs_call$data  # still unevaluated at this point


  score_type <- match.arg(score_type,
                          c("orthogonalized", "standardized",
                            "effective", "basic", "my_lab"))

  # identify flip-specific parameters in the call
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
    model <- update(model, x = TRUE)
    if (is.null(model$call$family))
      model$call$family <- model$family
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
  # this is used by update() inside gcor(), anova(), confint(), etc.
  glm_call          <- call("glm")
  glm_call$formula  <- original_formula
  glm_call$family   <- family
  glm_call$data     <- original_data
  if (!is.null(mf$offset))  glm_call$offset  <- mf$offset
  if (!is.null(mf$weights)) glm_call$weights <- mf$weights
  if (!is.null(mf$subset))  glm_call$subset  <- mf$subset

  # build the full flipscores call for print/summary display
  # this is stored in model$flipscores_call
  fs_call[[1L]]            <- quote(flipscores)
  fs_call$formula          <- original_formula
  fs_call$data             <- original_data_call
  fs_call$family           <- family
  fs_call$score_type       <- score_type
  fs_call$n_flips          <- flip_param_call$n_flips
  fs_call$alternative      <- alternative
  fs_call$id               <- id
  fs_call$seed             <- seed
  fs_call$to_be_tested     <- to_be_tested
  fs_call$precompute_flips <- precompute_flips
  fs_call$flips            <- NULL
  fs_call$nobservations    <- NULL
  fs_call$parms_DV         <- NULL

  # model$call       -> clean glm call, used by update(), gcor(), anova()
  # model$flipscores_call -> full flipscores call, used by print/summary
  model$call            <- glm_call
  model$flipscores_call <- fs_call
  model$flip_param_call <- flip_param_call
  model$score_type      <- score_type

  if (is.null(param_x_ORIGINAL) || (!param_x_ORIGINAL)) model$x <- NULL

  class(model) <- c("flipscores", class(model))
  return(model)
}
