#' compute_scores
#' @param model0 a \code{glm} object with the model under the null hypothesis (i.e. the covariates, the nuisance parameters).
#' @param model1 a \code{glm} or a \code{matrix} (or \code{vector}). If it is a \code{glm} object, it has the model under the alternative hypothesis. The variables in \code{model1} are the same variables in \code{model0} plus one or more variables to be tested.  Alternatively, if
#' \code{model1} is a \code{matrix}, it contains the tested variables column-wise.
#' @param score_type The type of score that is computed. It is "orthogonalized", "effective" or "basic".
#' "effective" and "orthogonalized" take into account the nuisance estimation.
#' @description 
#' Same usage as \code{anova.glm}. 
#' The parameter \code{id}  is used too, 
#' if present in \code{model0} (with priority) or in \code{model1}.
#'
#' @author Jesse Hemerik, Riccardo De Santis, Vittorio Giatti, Jelle Goeman and Livio Finos
#' @examples
#' set.seed(1)
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rpois(n=20,lambda=exp(Z+X))
#' mod0=glm(Y~Z,family="poisson")
#' (scr0=compute_scores(model0 = mod0, model1 = X, score_type = "effective"))
#' @export

compute_scores <- function(model0, model1,score_type){
  score_type=match.arg(score_type,c("orthogonalized","standardized","effective","basic","maybe"))
  if(missing(score_type))
    stop("test type is not specified or recognized")
  X=get_X(model0,model1)
  if(!is.null(model0))  
  {id_model1=if(is.list(model1)) eval(model1$flip_param_call$id) 
  else 
    NULL
  } else id_model1=NULL
  
  rm(model1)
  
  ###############
  if(is.null(model0$x)||is.null(ncol(model0$x))) {
    call <- getCall(model0)
    call$x=TRUE
    term <- terms(model0)
    env <- attr(term, ".Environment")
    model0=eval(call, env, parent.frame())
  }
  
  
  if(ncol(model0$x)==0) return(model0$y)
  if(is.null(model0$y)) model0$y=model0$model[,1]
  Z=model0$x
  residuals=(model0$y-model0$fitted.values)
  if(is.null(model0$weights))  sqrtW=rep(1,length(residuals)) else sqrtW=(as.numeric(model0$weights)**0.5)
  pars_score<-get_par_expo_fam(model0)
  D_vect<-pars_score$D
  V_vect<-pars_score$V
  sqrtinvV_vect<-V_vect**(-0.5)
  

  #BASIC SCORE
  if (score_type == "basic") {
    B=t(t(X * D_vect) %*% (diag(sqrtinvV_vect**2, nrow = length(sqrtW))))
    scores = B * (sqrtinvV_vect*residuals) * (1/length(model0$y)^0.5)
    scale_objects=list(nrm=sum(B^2)*sum((sqrtinvV_vect*residuals)^2))
  } else
    ##  EFFECTIVE SCORE 
  if(score_type=="effective"){
    A<-(sqrtW)*Z
    B<-X*(sqrtW)-t(crossprod(crossprod(A,X*(sqrtW)),solve(crossprod(A),t(A))))
    scores=B*sqrtinvV_vect*residuals/(length(model0$y)**0.5)
    scale_objects=list(nrm=sum(B^2)*sum((sqrtinvV_vect*residuals)^2))
  } else
    ##  STANDARDIZED SCORE
    if(score_type=="standardized"){
      A<-(sqrtW)*Z
      B<-X*(sqrtW)-t(crossprod(crossprod(A,X*(sqrtW)),solve(crossprod(A),t(A))))
      scores=B*sqrtinvV_vect*residuals
      AA = crossprod(A)
      AB = crossprod(A,B)
      m = sum(B * B) - crossprod(AB, solve(AA, AB))
      scale_objects=list(AA=AA,A=A, B=B,m=m,nrm=sum(B^2)*sum((sqrtinvV_vect*residuals)^2))
    } else  
      if(score_type=="maybe"){
        U=svd((sqrtW*Z),nv=0)$u
        B=((diag(nrow(Z))-tcrossprod(U))%*%(X*(sqrtW)))
        m = sum(B^2)
        scores=B*sqrtinvV_vect*residuals
        nrm=sum((sqrtinvV_vect*residuals)^2)*sum(B^2)
        scale_objects=list(U=U,B=B,m=m,nrm=nrm)
      } else  
        #ORTHO EFFECTIVE SCORE
  if(score_type=="orthogonalized"){
        OneMinusH = diag(nrow(Z)) - ((sqrtW)* Z) %*% solve(t(Z) %*% ((sqrtW**2) * Z)) %*% t(Z * (sqrtW))
        deco=svd((V_vect^0.5)*OneMinusH,nv = 0)
        deco$d[deco$d<1E-12]=0
        B=(t(X*sqrtW)%*%OneMinusH*(invV_vect**0.5))
        scores=t(B%*%deco$u)*(t(deco$u)%*%(sqrtinvV_vect*residuals))[,]*(1/length(model0$y)**0.5)
        nrm=sum(B^2)*sum((sqrtinvV_vect*residuals)^2)
        scale_objects=list(U=U,B=B,m=m,nrm=nrm)
  }
  attr(scores,"scale_objects")=scale_objects
  rownames(scores)=names(model0$fitted.values)
  return(scores)
}
