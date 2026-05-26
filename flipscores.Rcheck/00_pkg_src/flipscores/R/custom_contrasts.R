#' Build custom contrasts for flip-score tests
#'
#' @description
#' Creates emmeans-like contrasts for factors in a fitted model. The supported
#' syntax is intentionally close to common \code{emmeans} usage:
#' \code{~ factor}, \code{pairwise ~ factor}, \code{~ factor1:factor2}, and
#' \code{pairwise ~ factor1:factor2 | by_factor}.
#'
#' @param model a fitted \code{glm} or \code{flipscores} object.
#' @param specs a one-sided or two-sided formula. Put \code{pairwise} on the
#' left-hand side to request all pairwise comparisons; otherwise simple
#' treatment contrasts against the first level combination are generated.
#' @return an object of class \code{custom_contrasts}.
#' @export
custom_contrasts <- function(model, specs) {
  if (!inherits(specs, "formula")) {
    stop("specs must be a formula.", call. = FALSE)
  }

  model_data <- .fs_model_frame(model)
  parsed <- .fs_parse_contrast_specs(specs)
  contrast_vars <- parsed$contrast_vars
  by_vars <- parsed$by_vars
  all_vars <- c(contrast_vars, by_vars)

  missing_vars <- setdiff(all_vars, names(model_data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in the model data: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  factor_vars <- names(model_data)[vapply(model_data, .fs_is_factor_like, logical(1))]
  not_factor <- setdiff(all_vars, factor_vars)
  if (length(not_factor) > 0) {
    stop("Only factor or character variables can be used in custom contrasts. ",
         "Non-factor variables: ", paste(not_factor, collapse = ", "),
         call. = FALSE)
  }

  levels_list <- lapply(all_vars, .fs_factor_levels, data = model_data)
  names(levels_list) <- all_vars

  contrast_levels <- levels_list[contrast_vars]
  cell_grid <- expand.grid(contrast_levels, KEEP.OUT.ATTRS = FALSE,
                           stringsAsFactors = FALSE)
  cell_names <- .fs_interaction_names(cell_grid)

  if (length(by_vars) == 0) {
    grid <- cell_grid
    grid_names <- cell_names
    pieces <- list(.fs_contrast_block(cell_names, parsed$is_pairwise))
  } else {
    by_levels <- levels_list[by_vars]
    by_grid <- expand.grid(by_levels, KEEP.OUT.ATTRS = FALSE,
                           stringsAsFactors = FALSE)
    grid <- .fs_cross_grid(cell_grid, by_grid)
    grid_names <- .fs_interaction_names(grid[, all_vars, drop = FALSE])
    pieces <- vector("list", nrow(by_grid))

    for (i in seq_len(nrow(by_grid))) {
      by_label <- .fs_interaction_names(by_grid[i, , drop = FALSE])
      block <- .fs_contrast_block(cell_names, parsed$is_pairwise)
      block$matrix <- .fs_embed_block(block$matrix, nrow(by_grid), i)
      block$names <- paste0(block$names, " | ", by_label)
      pieces[[i]] <- block
    }
  }

  contrast_matrix <- do.call(rbind, lapply(pieces, `[[`, "matrix"))
  contrast_names <- unlist(lapply(pieces, `[[`, "names"), use.names = FALSE)
  colnames(contrast_matrix) <- grid_names
  rownames(contrast_matrix) <- contrast_names

  out <- list(
    matrix = contrast_matrix,
    names = contrast_names,
    grid = grid[, all_vars, drop = FALSE],
    combinations = grid_names,
    levels = levels_list,
    contrast_vars = contrast_vars,
    by_vars = by_vars,
    variables = all_vars,
    is_pairwise = parsed$is_pairwise,
    contrast_type = if (parsed$is_pairwise) "pairwise" else "simple",
    formula = specs
  )
  class(out) <- "custom_contrasts"
  out
}

#' Apply custom contrasts
#'
#' @param contrasts a \code{custom_contrasts} object.
#' @param model optional fitted model used only to attach metadata.
#' @return a contrast matrix.
#' @export
apply_custom_contrasts <- function(contrasts, model = NULL) {
  if (!inherits(contrasts, "custom_contrasts")) {
    stop("contrasts must be a custom_contrasts object.", call. = FALSE)
  }

  out <- contrasts$matrix
  rownames(out) <- contrasts$names
  attr(out, "grid") <- contrasts$grid
  attr(out, "variables") <- contrasts$variables
  attr(out, "contrast_vars") <- contrasts$contrast_vars
  attr(out, "by_vars") <- contrasts$by_vars
  if (!is.null(model)) {
    attr(out, "model_formula") <- stats::formula(model)
  }
  class(out) <- c("custom_contrast_matrix", "matrix")
  out
}

#' Flip-score tests for custom contrasts
#'
#' @description
#' Defines custom factor contrasts and applies the flip-score test to each
#' contrast. This is the convenience function that connects
#' \code{custom_contrasts()} to the flip-score machinery.
#'
#' @param model a fitted \code{glm} or \code{flipscores} object.
#' @param specs contrast formula; see \code{custom_contrasts()}.
#' @param score_type type of score, as in \code{flipscores()}.
#' @param n_flips number of sign flips.
#' @param alternative \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' @param id optional cluster id.
#' @param seed optional random seed.
#' @param flips optional precomputed flip matrix.
#' @param precompute_flips whether to precompute the flip matrix.
#' @param ... currently unused.
#' @return an object of class \code{flipscores_contrasts}.
#' @export
flipscores_contrasts <- function(model, specs,
                                 score_type = "standardized",
                                 n_flips = 5000,
                                 alternative = "two.sided",
                                 id = NULL,
                                 seed = NULL,
                                 flips = NULL,
                                 precompute_flips = TRUE,
                                 ...) {
  score_type <- match.arg(score_type,
                          c("orthogonalized", "standardized",
                            "effective", "basic"))
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

  if (!inherits(model, "glm")) {
    stop("model must be a glm or flipscores object.", call. = FALSE)
  }

  model <- .fs_model_with_x(model)
  contrasts <- custom_contrasts(model, specs)
  contrast_matrix <- apply_custom_contrasts(contrasts, model)
  ref_grid <- .fs_reference_grid(model, contrasts)
  grid_x <- .fs_model_matrix(model, ref_grid)

  if (ncol(grid_x) != length(stats::coef(model))) {
    stop("Could not align the reference grid with the model coefficients.",
         call. = FALSE)
  }

  if (anyNA(stats::coef(model))) {
    active <- !is.na(stats::coef(model))
    if (any(abs(grid_x[, !active, drop = FALSE]) > .Machine$double.eps)) {
      stop("Custom contrasts involving aliased coefficients are not supported.",
           call. = FALSE)
    }
    model$x <- model$x[, active, drop = FALSE]
    grid_x <- grid_x[, active, drop = FALSE]
  }

  if (!is.null(flips)) {
    n_flips <- nrow(flips)
    precompute_flips <- FALSE
  } else if (precompute_flips) {
    set.seed(seed)
    flips <- .make_flips(nrow(model$x), n_flips, id)
  }

  tests <- vector("list", nrow(contrast_matrix))
  estimates <- numeric(nrow(contrast_matrix))
  names(estimates) <- rownames(contrast_matrix)

  for (i in seq_len(nrow(contrast_matrix))) {
    coef_contrast <- as.vector(contrast_matrix[i, , drop = FALSE] %*% grid_x)
    estimates[i] <- sum(coef_contrast * stats::coef(model)[seq_along(coef_contrast)])
    tests[[i]] <- .fs_one_contrast_test(
      model = model,
      coef_contrast = coef_contrast,
      score_type = score_type,
      n_flips = n_flips,
      alternative = alternative,
      seed = seed,
      flips = flips,
      precompute_flips = precompute_flips
    )
  }

  table <- data.frame(
    contrast = rownames(contrast_matrix),
    estimate = estimates,
    Score = vapply(tests, function(x) x$Tspace[1], numeric(1)),
    p.value = vapply(tests, function(x) x$p.values[[1]], numeric(1)),
    row.names = NULL,
    check.names = FALSE
  )

  out <- list(
    call = match.call(),
    model = model,
    contrasts = contrasts,
    contrast_matrix = contrast_matrix,
    reference_grid = ref_grid,
    Tspace = do.call(cbind, lapply(tests, `[[`, "Tspace")),
    p.values = table$p.value,
    table = table,
    score_type = score_type,
    n_flips = n_flips,
    alternative = alternative
  )
  colnames(out$Tspace) <- table$contrast
  names(out$p.values) <- table$contrast
  class(out) <- "flipscores_contrasts"
  out
}

#' @export
print.custom_contrasts <- function(x, ...) {
  cat("Custom contrasts\n")
  cat("Formula: ", deparse1(x$formula), "\n", sep = "")
  cat("Type: ", x$contrast_type, "\n", sep = "")
  cat("Variables: ", paste(x$variables, collapse = ", "), "\n", sep = "")
  cat("Level combinations: ", length(x$combinations), "\n", sep = "")
  cat("Contrasts: ", nrow(x$matrix), "\n\n", sep = "")
  if (nrow(x$matrix) > 0) {
    print(utils::head(x$matrix, min(5, nrow(x$matrix))))
  }
  invisible(x)
}

#' @export
print.custom_contrast_matrix <- function(x, ...) {
  cat("Custom contrast matrix\n")
  cat("Dimensions: ", nrow(x), " contrasts x ", ncol(x), " cells\n\n", sep = "")
  print(utils::head(unclass(x), min(5, nrow(x))))
  invisible(x)
}

#' @export
print.flipscores_contrasts <- function(x, ...) {
  cat("Flip-score custom contrasts\n")
  cat("score_type = ", x$score_type, ", n_flips = ", x$n_flips,
      ", alternative = ", x$alternative, "\n\n", sep = "")
  print(x$table, row.names = FALSE)
  invisible(x)
}

.fs_one_contrast_test <- function(model, coef_contrast, score_type, n_flips,
                                  alternative, seed, flips,
                                  precompute_flips) {
  if (all(abs(coef_contrast) < .Machine$double.eps)) {
    stop("A contrast maps to a zero coefficient contrast.", call. = FALSE)
  }

  basis <- .fs_contrast_basis(coef_contrast)
  z0 <- model$x %*% basis$nuisance
  z1 <- model$x %*% basis$tested
  colnames(z0) <- paste0("N", seq_len(ncol(z0)))
  colnames(z1) <- "contrast"

  weights <- model$prior.weights
  if (is.null(weights)) {
    weights <- rep(1, NROW(model$y))
  }
  offset <- model$offset
  if (is.null(offset)) {
    offset <- rep(0, NROW(model$y))
  }

  model0 <- stats::glm.fit(
    x = z0,
    y = model$y,
    weights = weights,
    offset = offset,
    family = model$family
  )
  model0$x <- z0
  model0$y <- model$y
  model0$model <- data.frame(model$y, z0, check.names = FALSE)
  model0$prior.weights <- weights
  model0$offset <- offset
  class(model0) <- c("glm", "lm")

  scores <- compute_scores(model0, z1, score_type = score_type)
  flip_call <- as.call(list(
    .flip_test,
    n_flips = n_flips,
    alternative = alternative,
    flips = flips,
    seed = seed,
    precompute_flips = precompute_flips
  ))
  socket_compute_flip(scores, flip_call)[[1]]
}

.fs_contrast_basis <- function(coef_contrast) {
  coef_contrast <- as.numeric(coef_contrast)
  p <- length(coef_contrast)
  tested <- coef_contrast / sum(coef_contrast^2)

  if (p == 1) {
    nuisance <- matrix(numeric(0), nrow = 1, ncol = 0)
  } else {
    decomp <- svd(matrix(coef_contrast, nrow = 1), nv = p)
    nuisance <- decomp$v[, -1, drop = FALSE]
  }

  list(nuisance = nuisance, tested = matrix(tested, ncol = 1))
}

.fs_parse_contrast_specs <- function(specs) {
  spec_text <- paste(deparse(specs), collapse = "")
  spec_text <- gsub("[[:space:]]+", "", spec_text)
  parts <- strsplit(spec_text, "~", fixed = TRUE)[[1]]
  if (length(parts) == 1) {
    lhs <- ""
    rhs <- parts[1]
  } else {
    lhs <- parts[1]
    rhs <- parts[2]
  }

  is_pairwise <- identical(tolower(lhs), "pairwise")
  by_split <- strsplit(rhs, "|", fixed = TRUE)[[1]]
  contrast_expr <- by_split[1]
  by_expr <- if (length(by_split) > 1) by_split[2] else ""

  contrast_vars <- all.vars(stats::as.formula(paste("~", contrast_expr)))
  by_vars <- if (nzchar(by_expr)) {
    all.vars(stats::as.formula(paste("~", by_expr)))
  } else {
    character(0)
  }

  if (length(contrast_vars) == 0) {
    stop("No contrast variables found in specs.", call. = FALSE)
  }

  list(
    contrast_vars = contrast_vars,
    by_vars = by_vars,
    is_pairwise = is_pairwise
  )
}

.fs_model_frame <- function(model) {
  mf <- stats::model.frame(model)
  if (is.null(mf)) {
    stop("Could not recover the model data.", call. = FALSE)
  }
  mf
}

.fs_model_with_x <- function(model) {
  if (is.null(model$x)) {
    model <- stats::update(model, x = TRUE)
  }
  if (is.null(model$y)) {
    model$y <- stats::model.response(stats::model.frame(model))
  }
  model
}

.fs_is_factor_like <- function(x) {
  is.factor(x) || is.character(x)
}

.fs_factor_levels <- function(var, data) {
  if (is.factor(data[[var]])) {
    levels(data[[var]])
  } else {
    unique(as.character(data[[var]]))
  }
}

.fs_interaction_names <- function(grid) {
  if (ncol(grid) == 1) {
    as.character(grid[[1]])
  } else {
    apply(grid, 1, paste, collapse = ":")
  }
}

.fs_contrast_block <- function(cell_names, pairwise) {
  n_cells <- length(cell_names)
  if (n_cells <= 1) {
    return(list(matrix = matrix(numeric(0), nrow = 0, ncol = n_cells),
                names = character(0)))
  }

  if (pairwise) {
    pairs <- utils::combn(n_cells, 2)
    out <- matrix(0, nrow = ncol(pairs), ncol = n_cells)
    names <- character(ncol(pairs))
    for (i in seq_len(ncol(pairs))) {
      out[i, pairs[1, i]] <- 1
      out[i, pairs[2, i]] <- -1
      names[i] <- paste(cell_names[pairs[1, i]], "-", cell_names[pairs[2, i]])
    }
  } else {
    out <- matrix(0, nrow = n_cells - 1, ncol = n_cells)
    for (i in 2:n_cells) {
      out[i - 1, 1] <- -1
      out[i - 1, i] <- 1
    }
    names <- paste(cell_names[-1], "-", cell_names[1])
  }
  colnames(out) <- cell_names
  rownames(out) <- names
  list(matrix = out, names = names)
}

.fs_embed_block <- function(block, n_blocks, block_id) {
  out <- matrix(0, nrow = nrow(block), ncol = ncol(block) * n_blocks)
  cols <- ((block_id - 1) * ncol(block) + 1):(block_id * ncol(block))
  out[, cols] <- block
  out
}

.fs_cross_grid <- function(cell_grid, by_grid) {
  out <- vector("list", nrow(by_grid))
  for (i in seq_len(nrow(by_grid))) {
    out[[i]] <- cbind(cell_grid, by_grid[rep(i, nrow(cell_grid)), , drop = FALSE])
  }
  do.call(rbind, out)
}

.fs_reference_grid <- function(model, contrasts) {
  model_data <- .fs_model_frame(model)
  response <- all.vars(stats::formula(model))[1]
  predictors <- setdiff(names(model_data), response)
  requested_grid <- contrasts$grid

  base <- as.data.frame(lapply(model_data[predictors], .fs_reference_value),
                        stringsAsFactors = FALSE)
  grid <- base[rep(1, nrow(requested_grid)), , drop = FALSE]

  for (var in names(requested_grid)) {
    grid[[var]] <- requested_grid[[var]]
  }

  for (var in predictors) {
    if (is.factor(model_data[[var]])) {
      grid[[var]] <- factor(grid[[var]], levels = levels(model_data[[var]]))
    }
  }
  grid
}

.fs_reference_value <- function(x) {
  if (is.factor(x)) {
    levels(x)[1]
  } else if (is.character(x)) {
    unique(x)[1]
  } else if (is.numeric(x) || is.integer(x)) {
    mean(x, na.rm = TRUE)
  } else if (is.logical(x)) {
    FALSE
  } else {
    x[1]
  }
}

.fs_model_matrix <- function(model, newdata) {
  tt <- stats::delete.response(stats::terms(model))
  mf <- stats::model.frame(tt, data = newdata, xlev = model$xlevels)
  stats::model.matrix(tt, data = mf, contrasts.arg = model$contrasts)
}
