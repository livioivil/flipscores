compute_flips<- function(scores,alternative="two.sided",
                         n_flips=5000,
                         flips=NULL,
                         #  seed=NULL,
                         statTest="sum",
                         #output_flips=FALSE,
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
                      #output_flips=FALSE,
                      seed=NULL,
                      precompute_flips=TRUE,
                      cluster_index=NULL,
                      n_clusters=NULL,
                      ...){


  Tspace=.flip_test_no_pval(scores,flips=flips,
                            n_flips=n_flips,
                            .score_fun=.score_fun,
                            #output_flips=output_flips,
                            seed=seed,
                            precompute_flips=precompute_flips,
                            cluster_index=cluster_index,n_clusters=n_clusters)
  p.values=.t2p_only_first(ftail(unlist(Tspace)))
  # named vector?

  out=list(Tspace=Tspace,p.values=p.values)
  names(out$p.values)=colnames(scores)
  return(out)
}

.flip_test_no_pval<- function(scores,
                              flips=NULL,
                              n_flips=NULL,
                              .score_fun,
                              #output_flips=FALSE,
                              seed=NULL,
                              precompute_flips=TRUE,
                              cluster_index=NULL,
                              n_clusters=NULL,
                              ...){


  ##########################################

  #      browser()

  #      set.seed(seed)
  if(!is.null(flips)){
    #  browser()
    n_flips=nrow(flips)
    flips=flips[,rownames(scores),drop=FALSE]
    Tspace <- vapply(seq_len(n_flips),
                     function(i) as.numeric(.score_fun(flips[i,],scores)), numeric(1))

  } else {
    n_obs=nrow(scores)
    if (is.null(cluster_index)) {
      id <- attr(scores,"id")
      cluster_index <- if (is.null(id)) seq_len(n_obs) else match(id,unique(id))
    }
    if (is.null(n_clusters)) n_clusters <- max(cluster_index)
    set.seed(seed)
    Tspace <- numeric(n_flips)
    Tspace[1L] <- as.numeric(.score_fun(rep(1,n_obs),scores))
    if (n_flips > 1L) {
      for (i in seq.int(2L,n_flips)) {
        signs <- sample(c(-1L,1L),n_clusters,replace=TRUE)
        Tspace[i] <- as.numeric(.score_fun(signs[cluster_index],scores))
      }
    }
  }

  return(Tspace)
}



.flip_test_lm<- function(scores,
                         flips=NULL,
                         n_flips=NULL,
                         seed=NULL,
                         ...){


  Tspace=sapply(1:nrow(flips),
                function(i).score_std(flips[i,],scores))
  if(is.matrix(Tspace)) Tspace=t(Tspace)
  if(is.vector(Tspace)) Tspace=matrix(Tspace)

  score_sd0=sqrt(attr(scores,"Xnorm2"))*
    attr(scores,"scale_objects")$resid_sd
  if(any(!is.finite(score_sd0)|score_sd0<=0))
    stop("The null-model score scale must be finite and positive.")
  Tspace=sweep(Tspace,2L,score_sd0,FUN="/")

  return(Tspace)
}
