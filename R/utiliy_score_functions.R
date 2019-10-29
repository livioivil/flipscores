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