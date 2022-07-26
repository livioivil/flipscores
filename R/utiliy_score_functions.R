# for standardized:
.score_std=function(scr_eff,flp) {
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
  if(is.na(numerator/(denominator**0.5))) browser()
  numerator/(denominator**0.5)
}

#for effective and others:
.score <- function(Y,flp) flp%*%Y

#transform sum stat into t stat
.sum2t <- function(stat,sumY2,n){
  # sumY2=sum(Y^2,na.rm = TRUE)
  # n=sum(!is.na(Y))
  stat=stat/sqrt((sumY2-(stat^2)/n)*(n/(n-1)))
  stat
}

#### compute p-value
.t2p  <- function(pvls){
  mean(as.vector(pvls)>=pvls[1])
}
#################
.flip_test<- function(Y,score_type="standardized",alternative="two.sided",
                      n_flips=5000,
                      seed=NULL,
                      statTest="sum",...){
  
  if(alternative=="two.sided") ff <- function(Tspace) abs(Tspace) else
    if(alternative=="less") ff <- function(Tspace) -Tspace else
      if(alternative=="greater") ff <- function(Tspace) Tspace
      
      score_type=match.arg(score_type,c("orthogonalized","standardized","effective","basic"))
      if(score_type=="standardized") .score_fun <- .score_std else
        .score_fun <- .score
      
      nobs=nrow(Y)
      Tobs=  .score_fun(Y,rep(1,nobs))
      set.seed(seed)
      Tspace=data.frame(as.vector(c(Tobs,replicate(n_flips-1,{
        flp<-sample(c(-1,1),nobs, replace = T)
        .score_fun(Y,flp)
      }))))
      set.seed(NULL)
      # if(score_type=="effective"||score_type=="orthogonalized") 
      #   Tspace=.sum2t(Tspace,
      #                 sumY2 = sum(Y^2,na.rm = TRUE),
      #                 n=sum(!is.na(Y)))
      # 
      p.values=flipscores:::.t2p(ff(unlist(Tspace)))
      # named vector?
      
      out=list(Tspace=Tspace,p.values=p.values)
      names(out$p.values)=names(Y)
      return(out)
}

#################################
mahalanobis_npc <- function(permT){
  if(ncol(permT)==0) return(rep(0,nrow(permT)))
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
socket_compute_scores <- function(i,model,score_type){
  if(is.numeric(i)) {
    i=colnames(model$x)[i]
  }
  if(is.character(i)) {
    #check if it is present
    if(!all(i%in%colnames(model$x))) warning("Variable(s) ",paste(sep=", ",setdiff(i,colnames(model$x)))," is(are) not present in the model")
  }
  
  i_id=sapply(i,grep,colnames(model$x))
  
  #if(is.character(i)) i=which(colnames(model$x)==i)
  #to avoid re-run a flipscores everytime:
  attributes(model)$class= attributes(model)$class[attributes(model)$class!="flipscores"]
  tested_X=model[["x"]][, i_id, drop = FALSE]
  # model$x=model$x[,-c(i,exclude),drop=FALSE]
  if(ncol(model[["x"]])>0)
    colnames(model[["x"]])=paste0("V",1:ncol(model$x))
  
  model$call$data=data.frame(model[["y"]],model[["x"]][,-i_id,drop=FALSE])
  yname=as.character(model$call$formula[[2]])
  names(model$call$data)[1]=yname
  
  # frml=update(as.formula(model$call$formula), formula(paste(yname,"~0+",paste(colnames(model[["x"]]),collapse =" + "))))
  # model$call$formula=as.formula(frml)
  # model$call$formula=update( model$call$formula,formula(paste("~.",paste("-",colnames(model[["x"]])[i],collapse=""))))
  frml=as.formula(paste(yname,"~0+."))
  model$call$formula=as.formula(paste(yname,"~0+."))
  if(!is.null(model$offset)){
    offs<-model$offset
    model$call$formula=update( model$call$formula,formula(paste("~.+offset(offs)")))
  }
  model$call$score_type=NULL
  model$call$perms = NULL
  model$call[[1]]=quote(glm)
  model_i <-update(model)
  # print(flip_param_call$score_type)
  # browser()
  scores=compute_scores(model0 = model_i,model1 = tested_X,score_type=score_type)
}


#####################
# i and exclude are indices of the columns of model.frame x
socket_compute_flip <- function(scores,flip_param_call,score_type){
  ############### fit the H1 model and append the scores (refitted under H0s)
  
  ###############################
  ## compute flips
  
  ### TODO RENDERE PI AGILE INPUT DI id (es formula se possibile?) 
  # + quality check
  if(!is.null(flip_param_call$id)&&
     (!(flip_param_call$score_type%in%c("orthogonalized"))))
    scores=lapply(scores,rowsum,id)
  # scores=as.matrix(unlist(scores[,]))

  results=lapply(1:ncol(scores), function(id_col){
    score1=scores[,id_col,drop=FALSE]
    attributes(score1)$scale_objects=attributes(scores)$scale_objects[[id_col]]
    flip_param_call$Y=score1
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
  scores  <- socket_compute_scores(i,model,score_type=flip_param_call$score_type)
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
    Dhat<- Vhat <-1
    return(list(D=Dhat, V=Vhat))
  } else if(("glm"%in%class(model0))){
    
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
    
    # } else if(("glm"%in%class(model0))&&("negbin"%in%class(model0))){
    #   
    #   mu.est <- model0$fitted.values
    #   eta.est <- model0$family$linkfun(mu.est)
    #   V <- as.list(model0$family$variance)[[2]]
    #   #TODO mettere le altre 2 link
    #   mu <-quote(exp(eta))
    #   Dmu <- D(mu,"eta")
    #   a <- eval(V, list(mu= mu.est,.Theta=model0$theta)) / eval(Dmu, list(eta= eta.est))
    #   return(a)
    #   
  } else { 
    warning("Class of the model not detected, homoscedasticity and canonical link are assumed.")
    Dhat<-Vhat<-1
    return(list(D=Dhat, V=Vhat))}
}