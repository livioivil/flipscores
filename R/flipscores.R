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
#' @param id Optional cluster identifiers, also used to align score rows across
#'   models. Use a vector with one value per score row, a column name, or a
#'   one-sided formula such as \code{~ids}. Column names and formulas are
#'   evaluated in \code{data}, or separately in each fitted model's original
#'   data for a list of models. With \code{NULL}, rows are matched by row names.
#'   Repeated IDs are aligned by occurrence, with the maximum multiplicity
#'   across models, and all occurrences share the same random sign.
#'   A supplied vector is applied unchanged to every model and must match each
#'   model's number of score rows, after \code{subset} and missing-value removal.
#'   A length mismatch is an error; values are never recycled or truncated.
#'   For models with different fitted row counts, consider \code{id = ~ids} or
#'   \code{id = "ids"}, with an ID column present in every model's data.
#' @param seed Random seed for reproducibility.
#' @param to_be_tested Character vector or list of coefficient names to test.
#'   If a character vector, the same terms are intersected with each model's
#'   coefficients. If a list, each element applies to the corresponding model.
#'   If \code{NULL} all coefficients are tested.
#' @param tested_coeffs Alias for \code{to_be_tested} for compatibility.
#' @param flips Optional matrix of pre-computed flips. Columns may be named by
#'   observations or, when \code{id} is supplied, by cluster identifiers. All
#'   observations in the same cluster must receive the same sign.
#'   For joint tests with \code{id}, observation columns follow the common
#'   rows in \code{alignment$rows}; alternatively supply columns named by ID.
#' @param precompute_flips Logical, default \code{TRUE}. If \code{FALSE}, GLM
#'   flips are generated on demand using the same observation or cluster mapping.
#' @param ... Additional arguments passed to \code{glm()} or the internal
#'   flipscores engine.
#'
#' @return A \code{flipscores} object (case 1) or a \code{jfs}
#'   object (cases 2, 3, 4). Joint objects retain local score matrices in
#'   \code{objects}, with resolved identifiers in \code{attr(scores, "id")},
#'   and an \code{alignment} list containing the common \code{rows} and local
#'   row-to-common-row \code{index} maps.
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
#' # Different dataset sizes: resolve IDs separately in each model's data
#' df$ids <- rep(seq_len(n / 2), each = 2)
#' df_short <- df[-c(2, 5, 9), ]
#' models_by_id <- list(
#'   full = glm(y1 ~ x1 + x2, data = df),
#'   shorter = glm(y2 ~ x1 + x2, data = df_short)
#' )
#' # A shared vector of length 50 cannot match the shorter model's 47 rows:
#' tryCatch(
#'   flipscores(models_by_id, id = df$ids, to_be_tested = "x1", n_flips = 99),
#'   error = function(e) message(conditionMessage(e))
#' )
#' # These forms support different row counts and cluster flips by the same IDs:
#' res_by_id <- flipscores(models_by_id, id = ~ids, to_be_tested = "x1",
#'                         n_flips = 99, seed = 1)
#' res_by_name <- flipscores(models_by_id, id = "ids", to_be_tested = "x1",
#'                           n_flips = 99, seed = 1)
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
  if (is.list(formula) && all(sapply(formula, inherits, c("glm","jfs","flipscores","fs_contrasts","fs_lm")))) {
    #message("flipscores: list of glm objects detected -> joint test")

    out <- .join_flipscores(
      mods         = formula,
      tested_coeffs = to_be_tested,
      n_flips      = n_flips,
      flips        = flips,
      score_type   = score_type,
      seed         = seed,
      id           = id,
      alternative  = alternative,
      precompute_flips = precompute_flips,
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
      id           = id,
      alternative  = alternative,
      precompute_flips = precompute_flips,
      ...
    )
    out$call <- match.call()
    return(out)
  }

  ##############################################################
  # CASE 4: formula with matrix response -> list of formulas
  ##############################################################
  if (inherits(formula, "formula")) {

    ## TODO: se family gaussian e link=identity fare funzionina per lm simile
    # a clip e poi direttamente as.jfs
    #


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
  if(inherits(formula, c("fs_lm","fs_contrasts"))){
    return(do.call(update,match.call()))
  } else {
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
                            "effective", "basic"))

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
  flip_param_call$obs_names    <- eval(mf$obs_names,
                                           parent.frame())
  if(is.null(  flip_param_call$obs_names)){
    flip_param_call$nobservations    <-   eval(mf$nobservations,parent.frame()) }
  else
      flip_param_call$nobservations=length(flip_param_call$obs_names)
  # set defaults for any NULL values
  if (is.null(flip_param_call$precompute_flips))
    flip_param_call$precompute_flips <- TRUE
  if (is.null(flip_param_call$score_type))
    flip_param_call$score_type <- "standardized"
  if (is.null(flip_param_call$n_flips))
    flip_param_call$n_flips <- 5000

  mf$nobservations <- NULL
  mf$obs_names <- NULL

  if (!is.null(list(...)$parms_DV)) {
    flip_param_call$parms_DV <- list(...)$parms_DV
    mf$parms_DV <- NULL
  }
  for (nm in c("cluster_index", "n_clusters")) {
    value <- list(...)[[nm]]
    if (!is.null(value)) flip_param_call[[nm]] <- value
    if (nm %in% names(mf)) mf[[nm]] <- NULL
  }
  mf$output_flips <- NULL

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
    # Keep the fitted rows and their order; no refit of the full model is needed.
    if (is.null(model[["x"]])) model$x <- stats::model.matrix(model)

  } else {
    stop("'formula' must be a formula or a glm object.")
  }

  if (is.null(model$y)) model$y <- model$model[, 1]
  resolved_id <- .fs_resolve_id(id, model,
                                if (inherits(original_formula, "formula")) data else NULL)
  flip_param_call$id <- resolved_id
  # Do not pad local score objects to the union of observations across models.
  flip_param_call$nobservations <- NULL
  flip_param_call$obs_names <- NULL

  # resolve to_be_tested
  if (is.null(to_be_tested))
    to_be_tested <- colnames(model[["x"]])
  else {
    if (is.numeric(to_be_tested))
      to_be_tested <- colnames(model[["x"]])[to_be_tested]
    to_be_tested <- eval(to_be_tested, parent.frame())
  }

  # handle flips
  if (!is.null(flips)) {
    flip_param_call$flips <- .fs_prepare_flips(flips, rownames(model$x), resolved_id)
    flip_param_call$precompute_flips <- FALSE
    flip_param_call$n_flips <- nrow(flip_param_call$flips)

  } else if (flip_param_call$precompute_flips) {
    set.seed(seed)
    flip_param_call$flips <- .make_flips(
      nrow(model$x),
      flip_param_call$n_flips,
      flip_param_call$id
    )
  }

  if (!is.null(flip_param_call$flips)) {
    colnames(flip_param_call$flips) <- rownames(model$x)
  } else if (is.null(seed)) {
    flip_param_call$seed <- sample.int(.Machine$integer.max, 1L)
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
  attr(model$scores, "id")        <- resolved_id
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
  glm_call <- model$call
  glm_call$formula <- stats::formula(model)
  if (!is.character(family) || family != "negbinom")
    glm_call$family <- model$family
  model_data <- .fs_model_data(model,
                               if (inherits(original_formula, "formula")) original_data else NULL)
  if (!is.null(model_data)) glm_call$data <- model_data

  # build the user-facing flipscores call for print/summary/update
  fs_call[[1L]]         <- quote(flipscores)
  fs_call$formula       <- original_formula_call
  fs_call$data          <- original_data_call
  fs_call$nobservations <- NULL
  fs_call$parms_DV      <- NULL
  fs_call$cluster_index <- NULL
  fs_call$n_clusters <- NULL
  if (inherits(original_formula, "glm") && !is.null(model_data)) {
    # Calls created while joining models must remain usable after this frame exits.
    fs_call$formula <- stats::formula(model)
    fs_call$data <- model_data
    fs_call$family <- model$family
    fs_call$id <- id
    for (nm in setdiff(names(glm_call), c("", "formula", "data", "family")))
      if (!nm %in% names(fs_call)) fs_call[[nm]] <- glm_call[[nm]]
  }

  model$call            <- glm_call
  model$flipscores_call <- fs_call
  model$flip_param_call <- flip_param_call
  model$score_type      <- score_type
  model$summary_table=.get_summary_table_from_flipscores(model)

  if (is.null(param_x_ORIGINAL) || (!param_x_ORIGINAL)) model$x <- NULL

  class(model) <- c("flipscores", class(model))
  return(model)
}
