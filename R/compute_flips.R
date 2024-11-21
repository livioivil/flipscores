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
  
  
  Tspace=.flip_test_no_pval(scores,flips=flips,
                                       n_flips=n_flips,
                                       .score_fun=.score_fun,
                                       output_flips=output_flips,
                                       seed=seed,
                                       precompute_flips=precompute_flips)
  p.values=.t2p(ftail(unlist(Tspace)))
  # named vector?
  
  out=list(Tspace=Tspace,p.values=p.values)
  names(out$p.values)=colnames(scores)
  return(out)
}

.flip_test_no_pval<- function(scores,
                      flips=NULL,
                      n_flips=NULL,
                      .score_fun,
                      output_flips=FALSE,
                      seed=NULL,
                      precompute_flips=TRUE,
                      ...){
  
  
  ##########################################
  
  #      browser()
  n_obs=nrow(scores)
  #      set.seed(seed)
  if(!is.null(flips)){
    #  browser()
    Tspace=as.vector(c(sapply(1:(n_flips),
                              function(i).score_fun(flips[i,],scores))))
    
  } else {
    set.seed(seed)
    Tspace=as.vector(c(.score_fun(rep(1,n_obs),scores),replicate(n_flips-1,{
      .score_fun(sample(c(-1,1),n_obs, replace = T),scores)
    })))
  }

  return(Tspace)
}
