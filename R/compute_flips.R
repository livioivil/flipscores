compute_flips<- function(scores,alternative="two.sided",
                         n_flips=5000,
                         flips=NULL,
                         #  seed=NULL,
                         statTest="sum",
                         output_flips=FALSE,
                         ...)
{
  NULL
  # DA FARE. chiama socket_compute_flip()
}
#################
.flip_test<- function(scores,ftail,
                      flips=NULL,
                      n_flips=NULL,
                      .score_fun,
                      output_flips=FALSE,
                      seed=NULL,
                      precompute_flips=TRUE,
                      ...){
  
  
  ############### only for internal experiments
  #######START
  if(attributes(scores)$score_type=="my_lab") {
    .score_fun <- function(flp,Y,Xt) Xt%*%flp%*%Y
    n_obs=nrow(scores)
    Xt=attributes(scores)$scale_objects$Xt
    
    Tobs=  .score_fun(diag(n_obs),scores,Xt=Xt)
    #        set.seed(seed)
    Tspace=as.vector(c(Tobs,replicate(n_flips-1,{
      flp<-flip::rom(n_obs)
      .score_fun(flp,scores,Xt)
    })))
    #       set.seed(NULL)
    # Tspace=.sum2t(Tspace,
    #                sumY2 = sum(Y^2,na.rm = TRUE),
    #               n=sum(!is.na(Y)))
    
    p.values=.t2p(ftail(unlist(Tspace)))
    # named vector?
    
    out=list(Tspace=Tspace,p.values=p.values)
    names(out$p.values)=names(scores)
    return(out)  
  }
  ######### END
  ##########################################
  
  #      browser()
  n_obs=nrow(scores)
  Tobs=  .score_fun(rep(1,n_obs),scores)
  #      set.seed(seed)
  if(precompute_flips){
    #  browser()
    Tspace=as.vector(c(Tobs,
                       sapply(1:(n_flips-1),
                              function(i).score_fun(flips[i,],scores))))
    
  } else {
    set.seed(seed)
    Tspace=as.vector(c(Tobs,replicate(n_flips-1,{
      .score_fun(sample(c(-1,1),n_obs, replace = T),scores)
    })))
  }
  # TODO: decidere se meglio standardizzare così o con fisher stimata (credo questa seconda)
  # if(score_type=="effective"||score_type=="orthogonalized") 
  #  Tspace=.sum2t(Tspace,
  #                sumY2 = sum(Y^2,na.rm = TRUE),
  #                n=attributes(Y)$scale_objects$df.residual)
  
  p.values=.t2p(ftail(unlist(Tspace)))
  # named vector?
  
  out=list(Tspace=Tspace,p.values=p.values)
  names(out$p.values)=colnames(scores)
  return(out)
}