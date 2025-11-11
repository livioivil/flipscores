# for standardized:
.score_std=function(flp,scr_eff) {
  # scr_eff # un vettore
  numerator=crossprod(flp,scr_eff) #t(scr_eff)%*%flp
  # A<- attributes(scr_eff)$scale_objects$A
  # m<- attributes(scr_eff)$scale_objects$m
  if (all(sign(flp)==1)|(all(sign(flp)==-1))){
    denominator = 1
  } else {
    denominator = 1 - sum((colSums(attributes(scr_eff)$scale_objects$A[flp==1,,drop=FALSE])
                           -colSums(attributes(scr_eff)$scale_objects$A[flp==-1,,drop=FALSE]))^2)
  }
  numerator/(denominator**0.5)
}

#for effective and others:
.score <- function(flp,Y) flp%*%Y

get_std_dev_score <- function(fit,x2){
  # see also statmod::glm.scoretest()
  w <- fit$weights
  r <- fit$residuals
  if (any(w <= 0)) {
    r <- r[w > 0]
    x2 <- x2[w > 0]
    w <- w[w > 0]
  }
  fixed.dispersion <- (fit$family$family %in% c("poisson",
                                                "binomial"))
  if (fixed.dispersion)
    dispersion <- 1
  else if (fit$df.residual > 0) {
    dispersion <- sum(w * fit$residuals^2)/fit$df.residual
  }

  ws <- sqrt(w)
  if(!is.null(fit$qr)) {x2.1w <- qr.resid(fit$qr, ws * x2)
  } else  x2.1w = ws * x2
  sqrt(colSums(as.matrix(x2.1w * x2.1w))*dispersion)
}

#transform sum stat into t stat
.sum2t <- function(stat,sumY2,n){
  # sumY2=sum(Y^2,na.rm = TRUE)
  # n=sum(!is.na(Y))
  # print(sumY2)
  # print(stat)
  # stat0=stat
  sumY2=sumY2*(n**0.5)
  # if(any((sumY2-(stat^2)/n)*(n/(n-1))<0)) browser()
  stat=stat/sqrt((sumY2-(stat^2)/n)*(n/(n-1)))
  # print(stat)
# if(any(is.na(stat))) browser()
    stat
}
# .sum2t <- function(Tvector, sumY2, n){
#   dev_resid=(sumY2-(Tvector^2/n))
#   Tvector/sqrt(dev_resid/(n-1)*n)
# }
#### compute p-value
.t2p  <- function(pvls){
  mean(as.vector(pvls)>=pvls[1])
}

#################################
mahalanobis_npc <- function(permT){
  if(ncol(permT)==0) return(rep(0,nrow(permT)))
  permT=as.matrix(permT)
  dimnames(permT)=NULL

  if(ncol(permT)==1) {
    permT=as.vector(permT)
    return(abs(permT)/(sum(permT^2)^.5))}
  #else
  ei=eigen(t(permT)%*%permT)
  ei$vectors=ei$vectors[,ei$values>1E-12,drop=FALSE]
  ei$values=ei$values[ei$values>1E-12]
  dst=rowSums((permT%*%ei$vect%*%diag(ei$val^-.5))^2)
  dst=dst*nrow(permT)
  return(dst)
}

#performs mahalanobis_npc() on selected columns of permT
mahalanobis_npc_multi <- function(ids_list,permT){

  ff=function(ids,permT) mahalanobis_npc(permT[,ids,drop=FALSE])
  out=sapply(ids_list,ff,permT)
  (out)
}


# i and exclude are indices of the columns of model.frame x
socket_compute_scores <- function(i,model,score_type,nobservations=NULL,parms_DV=NULL){
  if(is.numeric(i)) {
    i=colnames(model$x)[i]
  }
  if(is.character(i)) {
    #check if it is present
    if(!any(i%in%colnames(model$x))) warning("Variable(s) ",paste(sep=", ",setdiff(i,colnames(model$x)))," is(are) not present in the model")
  }

  i_id=sapply(i,function(ii)which(colnames(model$x)==ii))

  #if(is.character(i)) i=which(colnames(model$x)==i)
  #to avoid re-run a flipscores everytime:
  attributes(model)$class= attributes(model)$class[attributes(model)$class!="flipscores"]
  tested_X=model[["x"]][, i_id, drop = FALSE]
  # model$x=model$x[,-c(i,exclude),drop=FALSE]
  if(ncol(model[["x"]])>0)
    colnames(model[["x"]])=paste0("V",1:ncol(model$x))

  model$call$data=data.frame(model[["y"]],model[["x"]][,-i_id,drop=FALSE])
  yname=all.vars( model$call$formula)[[1]]
  names(model$call$data)[1]=yname

  # frml=update(as.formula(model$call$formula), formula(paste(yname,"~0+",paste(colnames(model[["x"]]),collapse =" + "))))
  # model$call$formula=as.formula(frml)
  # model$call$formula=update( model$call$formula,formula(paste("~.",paste("-",colnames(model[["x"]])[i],collapse=""))))
  # frml=as.formula(paste(yname,"~0+."))
  model$call$formula=as.formula(paste(yname,"~0+."))
  # if(!is.null(model$offset)){
  #   offs<-model$offset
  #   model$call$formula=update( model$call$formula,formula(paste("~+offset(offs)")))
  # }
  model$call$score_type=NULL
  model$call$n_flips = NULL
  model$call$flips = NULL

  if(length(grep("Negative Binomial",model$family$family))==1)
    model$call[[1]]=quote(glm.nb) else
      model$call[[1]]=quote(glm)
  model_i <-update(model)
  # print(flip_param_call$score_type)
  # browser()
    scores=compute_scores(model0 = model_i,model1 = tested_X,score_type=score_type,nobservations=nobservations,parms_DV=parms_DV)
}


#####################
# i and exclude are indices of the columns of model.frame x
socket_compute_flip <- function(scores,flip_param_call){

#  flip_param_call$score_type=attributes(scores)$score_type


  # scores=as.matrix(unlist(scores[,]))
  if(is.null(flip_param_call$alternative)) flip_param_call$alternative = "two.sided"
  if(flip_param_call$alternative=="two.sided") flip_param_call$ftail <- function(Tspace) abs(Tspace) else
    if(flip_param_call$alternative=="less") flip_param_call$ftail <- function(Tspace) -Tspace else
      if(flip_param_call$alternative=="greater") flip_param_call$ftail <- function(Tspace) Tspace
      flip_param_call$alternative=NULL

      score_type=attributes(scores)$score_type
      score_type=match.arg(score_type,c("orthogonalized","standardized","effective","basic"))
      if(score_type=="standardized") flip_param_call$.score_fun <- .score_std else
        flip_param_call$.score_fun <- .score

      # if(flip_param_call$precompute_flips){
      #   set.seed(seed)
      #   flip_param_call$flips=.make_flips(nrow(scores),flip_param_call$n_flips)
      # }
  results=lapply(1:ncol(scores), function(id_col){
    score1=scores[,id_col,drop=FALSE]
    attributes(score1)$scale_objects=attributes(scores)$scale_objects[[id_col]]
    attributes(score1)$score_type=attributes(scores)$score_type
    attributes(score1)$sd=attributes(scores)$sd
    attributes(score1)$resid_std=attributes(scores)$resid_std
    flip_param_call$scores=score1
    res=eval(flip_param_call, parent.frame())
    res$scores=score1
    res
  })


  # flip_param_call$Y=scores
  #
  # results=eval(flip_param_call, parent.frame())
  # results=.flip_test(Y=scores,score_type=score_type,
  # alternative=alternative,n_flips=n_flips,
  # seed=seed,
  # statTest="sum")
  # results$scores=scores
  results
}

############################################
# i and exclude are indices of the columns of model.frame x
socket_compute_scores_and_flip <- function(i,model,exclude=NULL,
                                           flip_param_call){
  scores  <- socket_compute_scores(i,model,score_type=flip_param_call$score_type,
                                   nobservations=flip_param_call$nobservations,
                                   parms_DV = flip_param_call$parms_DV)
  results <- socket_compute_flip (scores,flip_param_call)
}

##########################
get_X <- function(model0,model1){
  if(is(model1,"glm")){
    mm=model.matrix(model1)
    vars1=colnames(mm)
    vars0=colnames(model.matrix(model0))
    model1=mm[,setdiff(vars1,vars0),drop=FALSE]
  } else if(!is.matrix(model1))
    model1=as.matrix(model1)
  # else is already a matrix
  model1
}
######### get the matrices D, V of glm model (for compute_scores) and also dispersion parameter

get_par_expo_fam <- function(model0){
  if(("lm"%in%class(model0))&(!("glm"%in%class(model0)))){
    Dhat<- Vhat <- rep(1,length(model0$fitted.values))
    return(list(D=Dhat, V=Vhat))
  } else if(("glm"%in%class(model0))){ # negbinom is comprised in this condition

    #The following lines are taken from the mdscore package on CRAN
    #by Antonio Hermes M. da Silva-Junior, Damiao N. da Silva and Silvia L. P. Ferrari

    # in response scale:
    mu.est <- model0$fitted.values
    eta.est <- model0$family$linkfun(mu.est)
    V <- if(model0$family[[1]] == "gaussian") quote(1)
    else as.list(model0$family$variance)[[2]]

    if(model0$family[[2]] %in% c("log", "cloglog", "logit")){
      mu <- switch(model0$family[[2]],
                   log     = quote(exp(eta)),
                   cloglog = quote(1 - exp(-exp(eta))),
                   logit   = quote(exp(eta)/(1 + exp(eta))))
    } else mu <- as.list(model0$family$linkinv)[[2]]

    if(model0$family[[2]] %in% c("probit","cauchit")){
      Dmu <- switch(model0$family[[2]],
                    probit  = quote(exp(-eta^2/2)/sqrt(2*pi)),
                    cauchit = quote(1/(pi*(1+eta^2))))
    } else{
      Dmu <- D(mu,"eta")
    }
    Dhat<-as.vector(eval(Dmu, list(eta= eta.est)))
    Vhat<-as.vector(eval(V, list(mu= mu.est,.Theta=model0$theta)))
    return(list(D=Dhat, V=Vhat))
  } else {
    warning("Class of the model not detected, homoscedasticity and canonical link are assumed.")
    Dhat<-Vhat<-rep(1,length(model0$fitted.values))
    return(list(D=Dhat, V=Vhat))}
}

