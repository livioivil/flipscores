# Resolve identifiers before the score engine replaces the original model data.
.fs_model_data <- function(model, data = NULL, envir = parent.frame()) {
  if (!is.null(data)) return(as.data.frame(data))
  model_env <- environment(stats::formula(model))
  if (is.null(model_env)) model_env <- envir
  value <- tryCatch(eval(model$call$data, model_env), error = function(e) NULL)
  if (is.null(value))
    value <- tryCatch(eval(model$call$data, envir), error = function(e) NULL)
  if (is.null(value)) value <- model$model
  if (!is.null(value)) as.data.frame(value) else NULL
}

.fs_validate_id <- function(id, n) {
  if (is.null(id)) return(NULL)
  if (!is.atomic(id) || !is.null(dim(id)) || length(id) != n)
    stop("'id' must be a vector with one value per score row (", n, ").",
         call. = FALSE)
  if (anyNA(id)) stop("'id' must not contain missing values.", call. = FALSE)
  unname(id)
}

.fs_model_rows <- function(model) {
  if (inherits(model, "fs_contrasts")) return(.fs_model_rows(model$model))
  if (inherits(model, "fs_lm")) return(rownames(model$info$D$X))
  if (inherits(model, "lm")) {
    if (is.matrix(model[["x"]])) return(rownames(model[["x"]]))
    return(rownames(stats::model.matrix(model)))
  }
  rownames(model$scores)
}

.fs_resolve_id <- function(id, model, data = NULL, envir = parent.frame()) {
  if (is.null(id)) return(NULL)
  rows <- .fs_model_rows(model)
  if (!inherits(id, "formula") && !(is.character(id) && length(id) == 1L))
    return(.fs_validate_id(id, length(rows)))

  if (inherits(model, "fs_contrasts")) model <- model$model
  if (inherits(model, "fs_lm")) {
    if (is.null(data)) {
      model_env <- environment(model$info$formula)
      if (is.null(model_env)) model_env <- envir
      data <- tryCatch(eval(model$call$data, model_env), error = function(e) NULL)
    }
  } else {
    data <- .fs_model_data(model, data, envir)
  }
  if (is.null(data))
    stop("The model's original data are required to evaluate 'id'.", call. = FALSE)
  data <- as.data.frame(data)

  if (inherits(id, "formula")) {
    if (length(id) != 2L)
      stop("'id' must be a one-sided formula, for example ~ids.", call. = FALSE)
    missing_vars <- setdiff(all.vars(id), names(data))
    if (length(missing_vars))
      stop("'id' column(s) missing from model data: ",
           paste(missing_vars, collapse = ", "), ".", call. = FALSE)
    value <- eval(id[[2L]], data, environment(id))
  } else {
    if (!id %in% names(data))
      stop("'id' column missing from model data: ", id, ".", call. = FALSE)
    value <- data[[id]]
  }
  if (!is.atomic(value) || !is.null(dim(value)) || length(value) != nrow(data))
    stop("The 'id' expression must return one value per data row.", call. = FALSE)
  index <- match(rows, rownames(data))
  if (anyNA(index))
    stop("Cannot match the model's score rows to its original data for 'id'.",
         call. = FALSE)
  .fs_validate_id(value[index], length(rows))
}

# Flatten only for bookkeeping; the returned jfs retains its object hierarchy.
.fs_joint_models <- function(models) {
  out <- list()
  collect <- function(x, path) {
    if (inherits(x, "jfs") && !inherits(x, "fs_contrasts")) {
      children <- x$objects
      names(children) <- .set_mods_names(children)
      for (i in seq_along(children))
        collect(children[[i]], paste(path, names(children)[i], sep = "."))
    } else {
      out[[length(out) + 1L]] <<- x
      names(out)[length(out)] <<- path
    }
  }
  names(models) <- .set_mods_names(models)
  for (i in seq_along(models)) collect(models[[i]], names(models)[i])
  names(out) <- make.unique(names(out))
  out
}

.fs_score_alignment <- function(rows, ids) {
  use_id <- !all(vapply(ids, is.null, logical(1)))
  if (!use_id) {
    common <- unique(unlist(rows, use.names = FALSE))
    index <- lapply(rows, match, table = common)
    aligned_rows <- data.frame(id = common, occurrence = 1L, row.names = common)
  } else {
    ids <- lapply(ids, as.character)
    common <- unique(unlist(ids, use.names = FALSE))
    codes <- lapply(ids, match, table = common)
    counts <- lapply(codes, tabulate, nbins = length(common))
    maxima <- Reduce(pmax, counts)
    offsets <- cumsum(c(0L, utils::head(maxima, -1L)))
    index <- lapply(codes, function(code) {
      occurrence <- as.integer(stats::ave(seq_along(code), code, FUN = seq_along))
      as.integer(offsets[code] + occurrence)
    })
    aligned_rows <- data.frame(
      id = rep(common, maxima), occurrence = sequence(maxima),
      row.names = paste0(".fs_row_", seq_len(sum(maxima))))
  }
  list(by = if (use_id) "id" else "rownames", rows = aligned_rows, index = index)
}

# Materialize structural NA values only when the aligned view is requested.
# Score metadata and the test calculations always use the local score objects.
.fs_aligned_scores <- function(object) {
  models <- .fs_joint_models(object$objects)
  alignment <- object$alignment
  if (is.null(alignment)) stop("The object has no score alignment.", call. = FALSE)
  blocks <- lapply(seq_along(models), function(i) {
    scores <- as.matrix(models[[i]]$scores)
    block <- matrix(NA_real_, nrow(alignment$rows), ncol(scores),
                    dimnames = list(rownames(alignment$rows),
                                    paste(names(models)[i], colnames(scores), sep = ".")))
    block[alignment$index[[i]], ] <- scores
    block
  })
  out <- do.call(cbind, blocks)
  if (alignment$by == "id") attr(out, "id") <- alignment$rows$id
  out
}

# Accept observation columns (possibly a superset), or columns named by cluster.
.fs_prepare_flips <- function(flips, rows, id = NULL) {
  if (!is.matrix(flips) || !is.numeric(flips) || nrow(flips) < 1L ||
      anyNA(flips) || any(!flips %in% c(-1, 1)))
    stop("'flips' must be a matrix of -1 and +1 values.", call. = FALSE)
  if (is.null(colnames(flips))) {
    if (ncol(flips) != length(rows))
      stop("'flips' must have one column per score row when unnamed.", call. = FALSE)
    colnames(flips) <- rows
  } else {
    if (anyDuplicated(colnames(flips)))
      stop("'flips' column names must be unique.", call. = FALSE)
    if (all(rows %in% colnames(flips))) {
      flips <- flips[, rows, drop = FALSE]
    } else if (!is.null(id) && all(as.character(id) %in% colnames(flips))) {
      flips <- flips[, as.character(id), drop = FALSE]
      colnames(flips) <- rows
    } else {
      stop("'flips' has wrong observation or cluster names.", call. = FALSE)
    }
  }
  if (!is.null(id) && any(flips != flips[, match(id, id), drop = FALSE]))
    stop("'flips' must use the same sign for all rows with the same 'id'.",
         call. = FALSE)
  flips
}
