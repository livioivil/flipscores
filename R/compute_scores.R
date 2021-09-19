#' compute_scores
#' @param model0 a \code{glm} object with the model under the null hypothesis (i.e. the covariates, the nuisance parameters).
#' @param model1 a \code{glm} or a \code{matrix} (or \code{vector}). If it is a \code{glm} object, it has the model under the alternative hypothesis. The variables in \code{model1} are the same variables in \code{model0} plus one or more variables to be tested.  Alternatively, if
#' \code{model1} is a \code{matrix}, it contains the tested variables column-wise.
#' @param score_type The type of score that is computed. It is "orthogonalized", "effective" or "basic".
#' Default is "orthogonalized". "effective" and "orthogonalized" take into account the nuisance estimation.
#' @description 
#' Same usage as \code{anova.glm}. 
#' The parameter \code{id}  is used too, 
#' if present in \code{model0} (with priority) or in \code{model1}.
#'
#' @author Jesse Hemerik, Vittorio Giatti, Jelle Goeman and Livio Finos
#' @examples
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rpois(n=20,lambda=exp(Z+X))
#' mod0=glm(Y~Z,family="poisson")
#' (scr0=compute_scores(model0 = mod0, model1 = X, score_type = "effective"))
#' @export

compute_scores <- function(model0, model1,score_type="effective"){
  score_type=match.arg(score_type,c("orthogonalized","standardized","effective","basic"))
  X=get_X(model0,model1)
  if(!is.null(model0))
  {id_model1=if(is.list(model1)) eval(model1$id) 
    else 
      NULL
  } else id_model1=NULL

  rm(model1)
  a <- get_a_expo_fam(model0)
  ###############
  if(is.null(model0$x)||is.null(ncol(model0$x))) {
    # model0=update(model0,x=TRUE)
    call <- getCall(model0)
    call$x=TRUE
    term <- terms(model0)
    env <- attr(term, ".Environment")
    model0=eval(call, env, parent.frame())
  }

  
  if(ncol(model0$x)==0) return(model0$y)
  Z=model0$x
  residuals=(model0$y-model0$fitted.values)/a
  W=diag(as.numeric(model0$weights))


  #BASIC SCORE
  if(score_type=="basic"){
    scores=(X*(residuals))
  } else
    ##  EFFECTIVE SCORE
    if(score_type=="effective"){
      OneMinusH=(diag(nrow(Z))-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
      scores=(OneMinusH%*%X)*residuals
      
    } else
      ##   SCORE standardized
      if(score_type=="standardized"){
        OneMinusH=(diag(nrow(Z))-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
        a=t(OneMinusH)%*%X
        B=OneMinusH%*%W%*%t(OneMinusH)
        scores=a*(residuals)
        var_obs=(t(a)%*%W%*%a)[,]
        e_var_flp=(t(a)%*%diag(diag(B))%*%a)[,]
        scale_obs=sqrt(var_obs/e_var_flp)
        attr(scores,"scale_obs")=scale_obs
        scale_objects=list(a=a, B=B)
        attr(scores,"scale_objects")=scale_objects
      } else
        #ORTHO EFFECTIVE SCORE
      if(score_type=="orthogonalized"){
        sqrtW=diag(sqrt(diag(W)))
        # OneMinusHtilde_0=(diag(n)-sqrtW%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z)%*%sqrtW)
        OneMinusH=(diag(nrow(Z))-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
        deco=svd(OneMinusH* matrix(diag(sqrtW),dim(OneMinusH)[1],dim(OneMinusH)[1],byrow = TRUE),
                 nv = 0)
        
        deco$d[deco$d<1E-12]=0
        scores=
          (t(t(X)%*%OneMinusH%*%deco$u)*
          (t(deco$u)%*%(residuals)))
        
      }
  
  # if(!is.null(model0$id)){
  #   scores=rowsum(scores,eval(model0$id)) 
  #   } else if(!is.null(id_model1)){
  #     scores=rowsum(scores,id_model1)      
  #   }
  
  # mean_adjustment
  # median_adjustment
  # 
  return(scores)
}
