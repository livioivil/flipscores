# join_flipscores.R
# Embedded and adapted from jointest (https://github.com/livioivil/jointest)
# Original author: Livio Finos
#' Coerce to jfs objects
#'
#' Generic function to coerce objects to class 'jfs'.
#'
#' @param x An object.
#' @param ... Additional arguments passed to methods.
#' @export
as.jfs <- function(x, ...) {
  UseMethod("as.jfs")
}

#' Coerce to jfs objects from contrasts
#'
#' Coerce contrast objects to class 'jfs'.
#'
#' @param x An object (typically from contrast methods).
#' @param ... Additional arguments passed to methods.
#' @export
as.jfs.fs_contrasts <- function(x, ...) {
  # Implementazione della coercizione per i contrasti di flipscores
  if (!inherits(x, "fs_contrasts")) {
    stop("The object must be of class 'fs_contrasts'")
  }

  # Aggiungi qui la logica specifica o restituisci l'oggetto con classe aggiornata
  class(x) <- unique(c("jfs", class(x)))
  return(x)
}

#' @keywords internal
.join_flipscores <- function(mods,
                             tested_coeffs  = NULL,
                             n_flips       = 5000,
                             flips         = NULL,
                             score_type    = "standardized",
                             statistics    = "t",
                             seed          = NULL,
                             output_models = TRUE,
                             id            = NULL,
                             alternative   = "two.sided",
                             precompute_flips = TRUE,
                             ...) {
  if (!length(mods)) stop("The model list must not be empty.", call. = FALSE)
  if (!is.null(id) && score_type == "orthogonalized")
    stop("'id' is not supported for orthogonalized scores.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  names(mods) <- .set_mods_names(mods)
  mods_names <- names(mods)
  leaves <- .fs_joint_models(mods)
  rows <- lapply(leaves, .fs_model_rows)
  # Resolve once on each original model, before constructing any null model.
  ids <- lapply(leaves, function(model) .fs_resolve_id(id, model))
  names(ids) <- names(rows) <- names(leaves)
  alignment <- .fs_score_alignment(rows, ids)
  cluster_names <- unique(alignment$rows$id)

  if (!is.null(flips)) {
    if (!is.null(id) && !is.null(colnames(flips)) &&
        !all(rownames(alignment$rows) %in% colnames(flips)) &&
        all(cluster_names %in% colnames(flips))) {
      FLIPS <- .fs_prepare_flips(flips, cluster_names)
    } else {
      expanded <- .fs_prepare_flips(flips, rownames(alignment$rows),
                                    if (!is.null(id)) alignment$rows$id else NULL)
      FLIPS <- expanded[, match(cluster_names, alignment$rows$id), drop = FALSE]
    }
    colnames(FLIPS) <- cluster_names
    n_flips <- nrow(FLIPS)
  } else if (precompute_flips) {
    FLIPS <- make_flips(length(cluster_names), n_flips, obs_names = cluster_names)
  } else {
    FLIPS <- NULL
    # Each coefficient regenerates the same global sequence of cluster signs.
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  }

  if (is.list(tested_coeffs) && length(tested_coeffs) != length(mods))
    stop("'tested_coeffs' must have one element per model.", call. = FALSE)
  position <- 0L
  run <- function(model, tested) {
    if (inherits(model, "jfs") && !inherits(model, "fs_contrasts")) {
      children <- lapply(model$objects, run, tested = tested)
      model$objects <- children
      model$Tspace <- .get_all_Tspace(children)
      model$summary_table <- .get_all_summary_table(children)
      child_leaves <- .fs_joint_models(children)
      model$alignment <- .fs_score_alignment(
        lapply(child_leaves, .fs_model_rows),
        lapply(child_leaves, function(x) attr(x$scores, "id")))
      return(model)
    }
    position <<- position + 1L
    i <- position
    cluster_index <- match(if (is.null(id)) rows[[i]] else as.character(ids[[i]]),
                           cluster_names)
    local_flips <- if (!is.null(FLIPS)) FLIPS[, cluster_index, drop = FALSE] else NULL
    if (!is.null(local_flips)) colnames(local_flips) <- rows[[i]]
    if (inherits(model, "glm")) {
      if (is.null(tested)) tested <- names(stats::coef(model))
      if (is.character(tested)) tested <- intersect(tested, names(stats::coef(model)))
      if (!length(tested)) stop("No requested coefficients are present in model ",
                                names(leaves)[i], ".", call. = FALSE)
      result <- flipscores(model, score_type = score_type, n_flips = n_flips,
                            alternative = alternative, to_be_tested = tested,
                            flips = local_flips, id = ids[[i]], seed = seed,
                            precompute_flips = precompute_flips,
                            cluster_index = cluster_index,
                            n_clusters = length(cluster_names), ...)
      # Do not retain calls referencing this temporary frame or generated flip
      # matrices. Preserve the original ID specification for later updates.
      result$flipscores_call$id <- id
      result$flipscores_call$flips <- if (!is.null(flips)) local_flips else NULL
      result$flipscores_call$to_be_tested <- tested
      result$flipscores_call$score_type <- score_type
      result$flipscores_call$n_flips <- n_flips
      result$flipscores_call$alternative <- alternative
      result$flipscores_call$seed <- seed
      result$flipscores_call$precompute_flips <- precompute_flips
      return(result)
    }
    # These engines need a matrix, even when GLM flips are generated on demand.
    if (is.null(local_flips)) {
      set.seed(seed)
      local_flips <- .mk_rnd_flips(length(cluster_names), n_flips)[,
                                                                 cluster_index, drop = FALSE]
      colnames(local_flips) <- rows[[i]]
    }
    if (inherits(model, "fs_contrasts")) {
      return(flipscores_contrasts(model$model, linfct = model$linfct,
                                  score_type = score_type, n_flips = n_flips,
                                  alternative = alternative, flips = local_flips,
                                  id = ids[[i]], seed = seed))
    }
    if (inherits(model, "fs_lm")) {
      if (is.null(tested)) tested <- unique(model$summary_table$coefficient)
      result <- .fs_lm(model$info$D, flips = local_flips,
                       alternative = alternative, tested_coeffs = tested)
      result$call <- model$call
      result$info <- model$info
      result$n_flips <- n_flips
      attr(result$scores, "id") <- ids[[i]]
      class(result) <- class(model)
      return(result)
    }
    stop("Unsupported model class in 'formula'.", call. = FALSE)
  }
  mods <- lapply(seq_along(mods), function(i) {
    tested <- if (is.list(tested_coeffs)) tested_coeffs[[i]] else tested_coeffs
    run(mods[[i]], tested)
  })
  names(mods) <- mods_names
  out <- list(Tspace = .get_all_Tspace(mods),
              summary_table = .get_all_summary_table(mods),
              objects = mods, alignment = alignment, call = match.call())
  class(out) <- c("jfs", "list")
  out
}

#' @title Coerce fs_contrasts to jfs
#' @name as.jfs
#' @description Coerces an \code{fs_contrasts} object to class \code{jfs}.
#' @param x an object of class \code{fs_contrasts}
#' @param ... further arguments (currently unused)
#' @return an object of class \code{jfs}
#' @export as.jfs

as.jfs.default <- function(object, ...) {
  Tspace=object$Tspace
  summary_table=object$summary_table
  call=object$call

  object$call <- object$Tspace <- object$summary_table <- NULL
  jfs_obj <- list(
    call      = call,
    summary_table = summary_table,
    Tspace = Tspace,
    objects      = object,
  )
  class(jfs_obj) <- "jfs"
  jfs_obj
}
