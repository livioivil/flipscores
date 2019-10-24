### in flipscores_glm

compute_scores_glm <- function(mf,score_type){
  mf$x=TRUE
  model <- eval(mf, parent.frame())
  
  # compute scores
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  singular.ok <- model.extract(mf, "singular.ok")
  if(is.null(singular.ok)) singular.ok = TRUE
  start <- model.extract(mf, "start")
  
  
  socket_compute_scores <- function(i,model){
    fit <- eval(call(model$method,
                     x = model$x[,-i,drop=FALSE], y = model$y, weights = as.vector(model.weights(mf)), start = start, etastart = etastart,
                     mustart = model$mustart, offset = model$offset, family = model$family,
                     control = model$control, intercept = FALSE, singular.ok = singular.ok))
    fit$x=model$x[,-i,drop=FALSE]
    compute_scores(fit0 = fit,X = model$x[,i,drop=FALSE],score_type=score_type)
  }
  
  model$scores=sapply(1:ncol(model$x),socket_compute_scores,model)
  colnames(model$scores)=colnames(model$x)
  model
}


######### get the a scaling value of glm model (for compute_scores)

get_a_expo_fam <- function(model0){
  if(class(model0)[1]=="lm"){
    a <- 1
  } else{

    #The following lines for obtaining 'a' are taken from the mdscore package on CRAN
    #by Antonio Hermes M. da Silva-Junior, Damiao N. da Silva and Silvia L. P. Ferrari

    # in response scale:
    mu.est <- model0$fitted.values
    # mu.est <- model0$fitted.values
    eta.est <- model0$family$linkfun(mu.est)

    V <- if(model0$family[[1]] == "gaussian") quote(1) else
      as.list(model0$family$variance)[[2]]


    if(model0$family[[2]] %in% c("log", "cloglog", "logit")){
      mu <- switch(model0$family[[2]],
                   log     = quote(exp(eta)),
                   cloglog = quote(1 - exp(-exp(eta))),
                   logit   = quote(exp(eta)/(1 + exp(eta))))
    }else mu <- as.list(model0$family$linkinv)[[2]]

    if(stringr::word(model0$family$family,1)=="Negative"){.Theta=model0$theta}  
    
    Dmu <- D(mu,"eta")
    a <- eval(V, list(mu= mu.est)) / eval(Dmu, list(eta= eta.est))
  }
  a
}
