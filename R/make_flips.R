#################
#' It creates a \code{n_flips}x\code{n_obs} matrix of random +1 and -1. 
#' The first row is made by ones (i.e. the observed test statistic is computed)
#' @param n_obs number of observations
#' @param n_flips number of flips
#' @export
#' 
make_flips <- function(n_obs,n_flips){
  .make_flips(n_obs,n_flips)
}


.make_flips <- function(n_obs,n_flips,id=NULL){
  if(is.null(id)){  
    flips = .mk_rnd_flips(n_obs,n_flips)
  } else {
    unique_id=unique(id)
    n_id=length(unique_id)
    temp = .mk_rnd_flips(n_id,n_flips)
    flips=matrix(NA,n_flips,n_obs)
    for(i in 1:n_id)
      flips[,id==unique_id[i]]=temp[,i]
    # for(i in unique(id))
    #   for(j in which(id==i)) 
    #     flips[,j]=temp[,i]
  }
  flips
}

.mk_rnd_flips <- function(n_obs,n_flips){
  matrix(c(rep(1,n_obs),sample(c(-1L,+1L),(n_flips-1)*n_obs,replace=TRUE)),
         n_flips,n_obs,byrow = TRUE)
}
