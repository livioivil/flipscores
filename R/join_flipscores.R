# join_flipscores.R
# Embedded and adapted from jointest (https://github.com/livioivil/jointest)
# Original author: Livio Finos

#' @keywords internal
.join_flipscores <- function(mods,
                             to_be_tested  = NULL,
                             n_flips       = 5000,
                             flips         = NULL,
                             score_type    = "standardized",
                             statistics    = "t",
                             seed          = NULL,
                             output_models = TRUE,
                             ...) {

  if (!is.null(seed)) set.seed(seed)

  # embed actual data frame directly from model$model into each model's call
  # do NOT use eval() on the data symbol - it may resolve to
  # stats::D or other functions that shadow the user's variable name
  # embed actual data frame directly from model$model into each model's call
  # embed actual data frame directly from model$model into each model's call
  for (i in seq_along(mods))
    mods[[i]]$call$data <- mods[[i]]$model

  names(mods) <- .set_mods_names(mods)

  # handle to_be_tested
  if (is.null(to_be_tested)) {
    to_be_tested <- .get_all_coeff_names_list(mods)
  } else if (!is.list(to_be_tested)) {
    temp         <- .get_all_coeff_names_list(mods)
    to_be_tested <- gsub(" ", "", to_be_tested)
    to_be_tested <- lapply(temp, function(nms)
      intersect(to_be_tested, gsub(" ", "", nms)))
  }

  # compute max n_obs across models
  # use mod$model directly instead of model.matrix()
  # to avoid symbol lookup issues with data
  n_obs_rn <- sapply(mods, function(mod)
    max(as.numeric(rownames(mod$model))))
  n_obs_rn <- max(n_obs_rn)
  n_obs    <- sapply(mods, function(mod) length(mod$y))
  n_obs    <- max(n_obs, n_obs_rn)

  mods_names <- names(mods)

  # generate shared flips across all models
  if (is.null(flips)) {
    FLIPS <- make_flips(n_obs = n_obs, n_flips = n_flips)
  } else {
    FLIPS <- flips
  }

  # run .flipscores_engine on each model
  # pass model$model as data to avoid symbol lookup issues
  mods <- lapply(seq_along(mods), function(i) {
    temp <- .flipscores_engine(
      formula       = mods[[i]],
      family        = mods[[i]]$family,
      data          = mods[[i]]$model,
      score_type    = score_type,
      flips         = FLIPS,
      to_be_tested  = to_be_tested[[i]],
      nobservations = n_obs,
      ...
    )
    if (statistics %in% "t") {
      temp$summary_table <- .get_summary_table_from_flipscores(
        temp,
        model_name = mods_names[i]
      )
    }
    temp
  })

  if (is.null(mods_names)) {
        names(mods) <- paste0("mod", seq_along(mods))
      } else {
        names(mods) <- mods_names
      }

      out <- list(
        Tspace        = .get_all_Tspace(mods),
        summary_table = .get_all_summary_table(mods),
        mods          = mods,
        call          = match.call()
      )

      class(out) <- c("joint_flipscores", class(out))
      out
}
