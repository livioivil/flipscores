# for standardized:
.score_std=function(scr_eff,flp) {
  # scr_eff # un vettore
  numeratore=t(scr_eff)%*%flp
  aflp= attributes(scr_eff)$scale_objects$a * flp
  denominatore= sqrt((t(aflp)  %*% attributes(scr_eff)$scale_objects$B %*% aflp)/length(scr_eff))
  numeratore/denominatore
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
t2p  <- function(pvls){
  mean(as.vector(pvls)>=pvls[1])
}
#################
.flip_test<- function(Y,score_type="standardized",alternative="two.sided",
                      n_flips=5000,seed=NULL,
                      statTest="sum"){
  
  if(alternative=="two.sided") ff <- function(Tspace) abs(Tspace) else
    if(alternative=="less") ff <- function(Tspace) -Tspace else
      if(alternative=="greater") ff <- function(Tspace) Tspace
      
      score_type=match.arg(score_type,c("orthogonalized","standardized","effective","basic"))
      if(score_type=="standardized") .score_fun <- .score_std else
        .score_fun <- .score
      
      n=nrow(Y)
      Tobs=  .score_fun(Y,rep(1,n))
      set.seed(seed)
      Tspace=data.frame(as.vector(c(Tobs,replicate(n_flips-1,{
        flp<-sample(c(-1,1),n, replace = T)
        .score_fun(Y,flp)
      }))))
      set.seed(NULL)
      if(score_type=="effective"||score_type=="orthogonalized") 
        Tspace=.sum2t(Tspace,
                      sumY2 = sum(Y^2,na.rm = TRUE),
                      n=sum(!is.na(Y)))
      
      p.values=t2p(ff(unlist(Tspace)))
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
  dst=abs(rowSums(permT%*%ei$vect%*%diag(ei$val^-.5)))
  return(dst)
  #manca ^2 e *nrow(permT). omesso per semplificare i calcoli
}

#performs mahalanobis_npc() on selected columns of permT
mahalanobis_npc_multi <- function(ids_list,permT){
  
  ff=function(ids,permT) mahalanobis_npc(permT[,ids,drop=FALSE])
  out=laply(ids_list,ff,permT)
  t(out)
}


# i and exclude are indices of the columns of model.frame x
socket_compute_scores_and_flip <- function(i,model,exclude=NULL,flip_param_call#score_type,id,n_flips,
                                           # alternative="two-sided",seed=NULL
                                            ){
  model$x=model.matrix(model)
  if(is.character(i)) i=which(colnames(model$x)==i)
  #to avoid re-run a flipscores everytime:
  attributes(model)$class= attributes(model)$class[attributes(model)$class!="flipscores"]
  tested_X=model[["x"]][, i, drop = FALSE]
  # model$x=model$x[,-c(i,exclude),drop=FALSE]
  if(ncol(model[["x"]])>0)
    colnames(model[["x"]])=paste0("V",1:ncol(model$x))
  
  model$call$data=data.frame(model$y,model$x)
  yname=as.character(model$call$formula[[2]])
  names(model$call$data)[1]=yname
  
  frml=update(as.formula(model$call$formula), formula(paste(yname,"~0+",paste(colnames(model[["x"]]),collapse =" + "))))
  model$call$formula=as.formula(frml)
  model$call$formula=update( model$call$formula,formula(paste("~.",paste("-",colnames(model[["x"]])[i],collapse=""))))
  if(!is.null(model$offset)){
    offs<-model$offset
    model$call$formula=update( model$call$formula,formula(paste("~.+offset(offs)")))
  }
  model_i <-update(model)
  # print(flip_param_call$score_type)
  # browser()
  scores=compute_scores(model0 = model_i,model1 = tested_X,score_type=flip_param_call$score_type)
  
  ############### fit the H1 model and append the scores (refitted under H0s)
  
  ###############################
  ## compute flips
  
  ### TODO RENDERE PI AGILE INPUT DI id (es formula se possibile?) 
  # + quality check
  if(!is.null(flip_param_call$id)&&(score_type%in%c("effective"))) #"standardized"
    scores=lapply(scores,rowsum,id)
  # scores=as.matrix(unlist(scores[,]))
  
  #  call to flip::flip()
  flip_param_call$Y=scores
  
  results=eval(flip_param_call, parent.frame())
  # results=.flip_test(Y=scores,score_type=score_type,
  # alternative=alternative,n_flips=n_flips,
  # seed=seed,
  # statTest="sum")
  results$scores=scores
  results
}

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
    Dhat<-diag(1, length(model0$y))
    Vhat<-diag(1, length(model0$y))
    a <- c(1, length(model0$y))
    return(list(a=a, D=Dhat, V=Vhat))
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
    }else mu <- as.list(model0$family$linkinv)[[2]]
    
    if(model0$family[[2]] %in% c("probit","cauchit")){
      Dmu <- switch(model0$family[[2]],
                    probit  = quote(exp(-eta^2/2)/sqrt(2*pi)),
                    cauchit = quote(1/(pi*(1+eta^2)))) 
    } else{
      Dmu <- D(mu,"eta")
    }
    Dhat<-diag(eval(Dmu, list(eta= eta.est)), length(mu.est))
    Vhat<-diag(eval(V, list(mu= mu.est,.Theta=model0$theta)), length(mu.est))
    a <- diag(Vhat) / diag(Dhat)
    return(list(a=a, D=Dhat, V=Vhat))
    
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
  } else { warning("Class of the model not detected, canonical link is assumed.")
    Dhat<-diag(1, length(model0$y))
    Vhat<-diag(1, length(model0$y))
    a <- c(1, length(model0$y))
    return(list(a=a, D=Dhat, V=Vhat))}
}
