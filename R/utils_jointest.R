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
.get_summary_table_from_flipscores <- function(object){
  tab = as.data.frame(summary.flipscores(object)$coefficients)
  tab = tab[!is.na(tab[, "score"]), ]

  names(tab) <- c("estimate", "score", "se", "z", "pcor", "p")
  # colnames(tab)[ncol(tab)]="p"
  mm=model.matrix(object)
  .assign=attr(mm,"assign")
  .assign=.assign[dimnames(mm)[[2]]%in%rownames(tab)]

  tab = cbind( .assign=.assign,
               coefficient = rownames(tab),
               tab)
}

.get_summary_table_from_fs_contrasts <- function(object) {
  tab <- object$table
  if (is.null(tab)) {
    stop("A fs_contrasts object must contain a table.", call. = FALSE)
  }

  out <- data.frame(
    .assign = NA_integer_,
    coefficient = tab$contrast,
    estimate = tab$estimate,
    score = tab$Score,
    p = tab$p.value,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  extra <- tab[setdiff(names(tab), c("contrast", "estimate", "Score", "p.value"))]
  data.frame(out, extra, check.names = FALSE)
}

.get_summary_table_from_object <- function(object) {
  if (inherits(object, "flipscores")) {
    return(.get_summary_table_from_flipscores(object))
  }
  if (inherits(object, c("fs_contrasts", "contrast_flipscores",
                         "contrasts_flipscores"))) {
    return(.get_summary_table_from_fs_contrasts(object))
  }
  if (!is.null(object$summary_table)) {
    return(object$summary_table)
  }
  stop("Cannot extract a summary table from this object.", call. = FALSE)
}

.joint_objects <- function(object) {
  objects <- object$objects
  if (is.null(objects)) {
    objects <- object$mods
  }
  if (is.null(objects)) {
    objects <- list(object)
  }
  objects
}

.as_jfs_component <- function(object, name = NULL) {
  if (inherits(object, "jfs")) {
    out <- object
    if (is.null(out$objects)) {
      out$objects <- .joint_objects(out)
    }
    return(out)
  }

  if (inherits(object, "glm") && !inherits(object, "flipscores")) {
    object <- flipscores(object)
  }

  if (!inherits(object, c("flipscores", "fs_contrasts",
                          "contrast_flipscores", "contrasts_flipscores"))) {
    stop("Cannot combine object of class ",
         paste(class(object), collapse = ", "), ".", call. = FALSE)
  }

  objects <- list(object)
  if (!is.null(name) && nzchar(name)) {
    names(objects) <- name
  }
  if (is.null(names(objects)) || !nzchar(names(objects)[1])) {
    names(objects) <- "object1"
  }

  out <- list(
    Tspace = object$Tspace,
    summary_table = .get_all_summary_table(objects),
    objects = objects,
    call = object$call
  )
  class(out) <- unique(c("jfs", class(out)))
  out
}

.rbind_fill <- function(tabs) {
  tabs <- Filter(Negate(is.null), tabs)
  if (length(tabs) == 0) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(tabs, names), use.names = FALSE))
  tabs <- lapply(tabs, function(tab) {
    missing <- setdiff(all_names, names(tab))
    for (nm in missing) {
      tab[[nm]] <- NA
    }
    tab[all_names]
  })

  out <- do.call(rbind, tabs)
  rownames(out) <- NULL
  out
}


#--------------------------------------------
# Bind summary tables from all models
#--------------------------------------------
.get_all_summary_table <- function(mods,mods_name=NULL){
  if(is.null(mods_name)) mods_name=names(mods)
  if (is.null(mods_name)) mods_name <- rep("", length(mods))
  res=lapply(1:length(mods), function(i) {
    tab <- .get_summary_table_from_object(mods[[i]])
    label <- mods_name[i]
    if (is.na(label) || identical(label, "")) {
      label <- paste0("object", i)
    }
    if (!"model" %in% names(tab)) {
      tab <- data.frame(model = rep(label, nrow(tab)), tab,
                        check.names = FALSE)
    } else if (all(is.na(tab$model) | tab$model == "")) {
      tab$model <- label
    }
    tab
  })
  .rbind_fill(res)
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
