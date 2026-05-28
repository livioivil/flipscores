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
# Bind summary tables from all models,
# adding a Model column
#--------------------------------------------
.get_all_summary_table <- function(mods) {
  tabs <- lapply(names(mods), function(nm) {
    tab <- mods[[nm]]$summary_table
    if (is.null(tab)) return(NULL)
    cbind(Model = nm, tab)
  })
  do.call(rbind, tabs)
}

#--------------------------------------------
# Build summary table from a single
# flipscores object
#--------------------------------------------
.get_summary_table_from_flipscores <- function(x) {
  pvals <- x$p.values
  if (is.null(pvals)) return(NULL)

  Tobs <- sapply(seq_along(pvals), function(i) {
    ts <- x$Tspace[, i]
    ts[nrow(x$Tspace)]
  })

  data.frame(
    Coeff             = names(pvals),
    T_obs             = Tobs,
    p.value           = pvals,
    row.names         = NULL,
    stringsAsFactors  = FALSE
  )
}
