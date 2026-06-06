# join_flipscores.R
# Embedded and adapted from jointest (https://github.com/livioivil/jointest)
# Original author: Livio Finos

#' @keywords internal
.join_flipscores <- function(mods,
                             tested_coeffs  = NULL,
                             n_flips       = 5000,
                             flips         = NULL,
                             score_type    = "standardized",
                             statistics    = "t",
                             seed          = NULL,
                             output_models = TRUE,
                             ...) {
    if(!is.null(seed)) set.seed(seed)

    for(i in 1:length(mods))
      mods[[i]]$call$data=eval(mods[[i]]$call$data, parent.frame())

    names(mods) = .set_mods_names(mods)
    if (is.null(tested_coeffs)) {
      tested_coeffs = .get_all_coeff_names_list(mods)
    }
    if (!is.list(tested_coeffs)) {
      temp = .get_all_coeff_names_list(mods)

      tested_coeffs = gsub(" ", "", tested_coeffs)
      tested_coeffs = lapply(temp, function(nms) intersect(tested_coeffs,
                                                           gsub(" ", "", nms)))
    }

    n_obs_rn = sapply(mods, function(mod) max(as.numeric(rownames(model.matrix(mod)))))
    n_obs_rn = max(n_obs_rn)
    n_obs=sapply(mods, function(mod) length(mod$y))
    n_obs=max(n_obs,n_obs_rn)


    mods_names=names(mods)
    if(is.null(flips)){
      FLIPS=make_flips(n_obs=n_obs,n_flips=n_flips)
    }else{
      FLIPS = flips
    }
    mods = lapply(1:length(mods), function(i) {

      temp = flipscores(formula = eval(mods[[i]],parent.frame()), score_type = score_type,
                        flips = eval(FLIPS), to_be_tested = tested_coeffs[[i]],
                        output_flips=FALSE,nobservations=n_obs,...
      )
      if (statistics %in% c("t")) {
        temp$summary_table=.get_summary_table_from_flipscores(temp)
        temp
      }})

    if(is.null(mods_names)){
      names(mods)=paste0("mod",1:length(mods))
    } else
      names(mods) = mods_names

    out=list(Tspace=.get_all_Tspace(mods),
             summary_table=.get_all_summary_table(mods),
             objects=mods,
             call = match.call())
    class(out) <- unique(c("jfs", class(out)))
    out
  }
