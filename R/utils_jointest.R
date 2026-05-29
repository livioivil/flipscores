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
  tab = as.data.frame(summary(object)$coefficients)
  tab = tab[!is.na(tab[, "Score"]), ]

  names(tab) <- c("estimate", "score", "se", "z", "pcor", "p")
  # colnames(tab)[ncol(tab)]="p"
  mm=model.matrix(object)
  .assign=attr(mm,"assign")
  .assign=.assign[dimnames(mm)[[2]]%in%rownames(tab)]

  tab = cbind( .assign=.assign,
               coefficient = rownames(tab),
               tab)
}


#--------------------------------------------
# Bind summary tables from all models
#--------------------------------------------
.get_all_summary_table <- function(mods,mods_name=NULL){
  if(is.null(mods_name)) mods_name=names(mods)
  res=lapply(1:length(mods), function(i) {
    cbind(model=names(mods)[i],
          mods[[i]]$summary_table)
  })
  res=do.call(rbind,res)
  rownames(res)=NULL
  res
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
#' Print method for joint_flipscores
#'
#' @param x A \code{joint_flipscores} object.
#' @param n Number of rows to show at head and tail of summary table.
#' @param ... Additional arguments (currently unused).
#' @export
print.joint_flipscores <- function(x, n = 2, ...) {
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  msg <- "== Joining n = %s models"
  cat(sprintf(msg, length(unique(x$summary_table$model))))
  cat("\n\n")
  .trim(x$summary_table, n = n)
  invisible(x)
}

#' Summary method for joint_flipscores
#'
#' @param object A \code{joint_flipscores} object.
#' @param digits Number of digits to print. Default \code{4}.
#' @param ... Additional arguments (currently unused).
#' @export
summary.joint_flipscores <- function(object, digits = 4, ...) {
  cat("\nCall: ")
  print(object$call)
  cat("\n")
  tab <- object$summary_table
  tab$.assign <- NULL  # remove internal column
  print(tab, digits = digits)
  invisible(object)
}
