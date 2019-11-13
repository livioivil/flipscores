# i and exclude are indices of the columns of model.frame x
socket_compute_scores <- function(i,model,exclude=NULL,score_type){
  #to avoid re-run a flipscores everytime:
  attributes(model)$class= attributes(model)$class[attributes(model)$class!="flipscores"]

  model$call$data=data.frame(model$y,model$x[,-c(i,exclude),drop=FALSE])
  yname=as.character(model$call$formula[[2]])
  model$call$formula=as.formula(paste(yname,"~0+."))
  names(model$call$data)[1]=yname
  model_i <-update(model)
  compute_scores(model0 = model_i,model1 = model$x[,i,drop=FALSE],score_type=score_type)
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
######### get the a scaling value of glm model (for compute_scores)

get_a_expo_fam <- function(model0){
  if(("lm"%in%class(model0))&(!("glm"%in%class(model0)))){
    return(1)
  } else if(("glm"%in%class(model0))){
    
    #The following lines for obtaining 'a' are taken from the mdscore package on CRAN
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
    
    Dmu <- D(mu,"eta")
    a <- eval(V, list(mu= mu.est,.Theta=model0$theta)) / eval(Dmu, list(eta= eta.est))
    return(a)
    
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
  } else { warning("Class of the model not detected, scale parameter is set to 1.")
    return(1)}
}