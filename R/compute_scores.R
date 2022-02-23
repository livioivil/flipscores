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
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rpois(n=20,lambda=exp(Z+X))
#' mod0=glm(Y~Z,family="poisson")
#' (scr0=compute_scores(model0 = mod0, model1 = X, score_type = "effective"))
#' @export

compute_scores <- function(model0, model1,score_type){
  score_type=match.arg(score_type,c("orthogonalized","standardized","effective","basic"))
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
  if(is.null(model0$weights))  W=rep(1,length(residuals)) else W=(as.numeric(model0$weights))
  pars_score<-get_par_expo_fam(model0)
  D<-pars_score$D
  V<-pars_score$V
  invV_vect<-diag(V**(-1))
  
  #BASIC SCORE
  if(score_type=="basic"){
    scores=t(t(X)%*%D%*%(diag(invV_vect)))*residuals*(1/length(model0$y)**0.5)
  } else
    ##  EFFECTIVE SCORE OR "standardized"
  if(score_type%in%c("standardized","effective")){
      OneMinusH=diag(nrow(Z))-((W**0.5)*Z)%*%solve(t(Z*W)%*%Z)%*%t(Z*(W**0.5))
      scores=t(t(X*(W**0.5))%*%(OneMinusH))*((invV_vect**0.5)*residuals)*(1/length(model0$y)**0.5)
    } else
        #ORTHO EFFECTIVE SCORE
  if(score_type=="orthogonalized"){
        sqrtW=diag(sqrt(diag(W)))
        OneMinusH=(diag(nrow(Z))-(W**0.5)%*%Z%*%solve(t(Z)%*%W%*%Z)%*%t(Z)%*%(W**0.5))
        deco=svd((V^0.5)%*%OneMinusH,nv = 0)
        deco$d[deco$d<1E-12]=0
        scores=t(t(X)%*%sqrtW%*%OneMinusH%*%(diag(invV_vect)**0.5)%*%deco$u)*(t(deco$u)%*%residuals)[,]*(1/length(model0$y)**0.5)
  }
  
  if(score_type=="standardized"){
    a=OneMinusH%*%((W**0.5)*X)
    B=OneMinusH
    scale_objects=list(a=a, B=B)
    attr(scores,"scale_objects")=scale_objects
  }
  rownames(scores)=names(model0$fitted.values)
  return(scores)
}
