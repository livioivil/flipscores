
###########
get_head_flip_out <- function(x){
  n_flips <- x$n_flips
  if (is.null(n_flips) && !is.null(x$flip_param_call)) {
    n_flips <- x$flip_param_call$n_flips
  }
  if(!is.null(dim(n_flips)))
    n_flips=nrow(n_flips)
  if(length(grep("Negative Binomial",x$family$family))==0)
  {paste("Flip Score Test:
         score_type =",x$score_type,
         ", n_flips =",n_flips,"\n")}
  else
    paste("Flip Score Test:
          score_type =",x$score_type,
          ", n_flips =",n_flips,
          ", theta =",round(x$theta,digits=5),"\n")
}

.intercept_in_Z_and_count_family <- function(mod){
  has_intercept=TRUE
  # check if family is count data
  if(mod$family$family%in% c("poisson", "quasipoisson", "Negative Binomial"))
    {
    Z=model.matrix(mod)
    if(ncol(Z)==0) {
      has_intercept=FALSE
      } else {
        P=Z%*%solve(t(Z)%*%Z)%*%t(Z)
        ones=matrix(1,nrow(Z),1)
        has_intercept=!(sum(abs(P%*%ones-ones))>1E-10)
      }
    if(!has_intercept) { #if intercept is the space spanned by Z
    warning("The Normalized Generalized Partial Correlation (Determination) Coefficient for Count families without interncept in the null model has not implemented, yet. NA will be returned.")
    }
  }
  has_intercept
}


.t2p <- function(Tspace,tail=0){
  Tspace=.set_tail(Tspace,tail=tail)
  if(is.vector(Tspace)){
    P = rank(-Tspace , ties.method = "max", na.last = "keep")/length(Tspace)
  } else {
    P = apply(-Tspace, 2, rank, ties.method = "max", na.last = "keep")/nrow(Tspace)
    P = as.matrix(P)
  }
  P
}

.t2p_only_first <- function(Tspace,tail=0){
  Tspace=.set_tail(Tspace,tail=tail)
  if(is.vector(Tspace)){
    P = mean(Tspace>=Tspace[1])
  } else if(ncol(Tspace)==1){
    P = mean(Tspace[,]>=Tspace[1,1])
  } else    {
    P = apply(Tspace,2, function(Tsp)mean(Tspace>=Tspace[1]))
  }
  P
}


.set_tail <- function(Tspace,tail=0){
  if((tail==0)|(tail=="two.sided"))
    Tspace=abs(Tspace) else
      if((tail<0)|(tail=="less"))
        Tspace=-Tspace

      Tspace
}
.find_common_pattern <- function(vettore) {
  if (length(vettore) < 2) {
    return(vettore)
  }

  if(length(grep(":",vettore[1]))>0){
    vettore_splt=strsplit(vettore,":")
    pttrns=sapply(1:length(vettore_splt[[1]]),function(i){
      .find_common_pattern (sapply(vettore_splt,function(x) x[i]) )
    })
    return(paste(pttrns,collapse=":"))
  }

  chars=sapply(vettore,strsplit,"")
  #matrix of all chracters
  charsMat=suppressWarnings(do.call(cbind,chars))
  # ask row-wise if they are all equals
  all_eqs=apply(charsMat,1,function(x)length(unique(x))==1)
  # this is the first different character
  if(all(all_eqs)) {
    common_pattern=vettore[1]
  } else {
    common_pattern=substr(vettore[1],1,max(which.min(all_eqs)-1,1))
  }


  return(common_pattern)
}
