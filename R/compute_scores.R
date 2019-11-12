#' compute_scores
#' @param model0 a \code{glm} object with the model under the null hypothesis (i.e. the covariates, the nuisances parameters).
#' @param model1 a \code{glm} or a \code{matrix} (or \code{vector}). If it is a \code{glm} object, it has the model under the alternative hypothesis. The variables in \code{model1} are the same variables in \code{model0} plus one or more variables to be tested.  Alternatively, if
#' \code{model1} is a \code{matrix}, it contains the tested variables column-wise.
#' @param score_type The type of score that is computed, either "orthogonalized", "effective" or "basic". 
#' By default is "orthogonalized". "effective" and "orthogonalized" takes into account nuisance estimation.
#'
#' @examples
#'
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rpois(n=20,lambda=exp(Z+X))
#' mod0=glm(Y~Z,family="poisson")
#' scr0=compute_scores(model0 = mod0, X = X, score_type = "basic")
#' flip:::flip(scr0)
#'
#' mod1=glm(Y~Z+X,family="poisson")
#' scr1=compute_scores(model0 = mod1, X = X, score_type = "basic")
#' flip:::flip(scr1)
#'
#' @export

compute_scores <- function(model0, model1,score_type="orthogonalized"){
  X=get_X(model0,model1)
  rm(model1)
  a <- get_a_expo_fam(model0)
  ###############
  if(is.null(model0$x)||(length(model0$x)==0)) model0=update(model0,x=TRUE)
  Z=model0$x
  residuals=(model0$y-model0$fitted.values)/a
  W=diag(as.numeric(model0$weights))


  #BASIC SCORE
  if(score_type=="basic"){
    return(as.vector(X*(residuals)))
  } else
    ##  EFFECTIVE SCORE
    if(score_type=="effective"){
      OneMinusH=(diag(nrow(Z))-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
      return(as.vector(
        (t(X)%*%OneMinusH*(residuals))))
    } else
      #ORTHO EFFECTIVE SCORE
      if(score_type=="orthogonalized"){
        sqrtW=diag(sqrt(diag(W)))
        # OneMinusHtilde_0=(diag(n)-sqrtW%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z)%*%sqrtW)
        OneMinusH=(diag(nrow(Z))-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
        OneMinusHtilde=(sqrt(1/diag(W))) * OneMinusH * matrix(diag(sqrtW),dim(OneMinusH)[1],dim(OneMinusH)[1],byrow = TRUE)
        #equivalent to:
        # sqrtWinv=diag(sqrt(1/diag(W)))
        # OneMinusHtilde=sqrtWinv%*% OneMinusH%*%sqrtW
        OneMinusHtilde=(OneMinusHtilde+t(OneMinusHtilde))/2
        
        deco=svd(OneMinusHtilde* matrix(diag(sqrtW),dim(OneMinusH)[1],dim(OneMinusH)[1],byrow = TRUE))
        # equivalent to:
        # deco=svd(OneMinusHtilde%*%sqrtW)
        
        deco$d[deco$d<1E-12]=0
        return(
          t(t(X)%*%(deco$v)%*%diag(deco$d))*
            as.vector(diag(deco$d)%*%t(deco$v)%*%diag(diag(W)^-1)%*%(residuals))
          )
        
      }
}
# OLD: 
# compute_scores <- function(model0, X,score_type="orthogonalized"){
#   
#   a <- get_a_expo_fam(model0)
#   ###############
#   Z=model0$x
#   residuals=(model0$y-model0$fitted.values)/a
#   W=diag(as.numeric(model0$weights))
#   
#   
#   #BASIC SCORE
#   if(score_type=="basic"){
#     return(as.vector(X*(residuals)))
#   } else
#     ##  EFFECTIVE SCORE
#     if(score_type=="effective"){
#       OneMinusH=(nrow(Z)-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
#       return(as.vector(
#         (t(X)%*%OneMinusH*(residuals))))
#     } else
#       #ORTHO EFFECTIVE SCORE
#       if(score_type=="orthogonalized"){
#         sqrtW=diag(sqrt(diag(W)))
#         # OneMinusHtilde_0=(diag(n)-sqrtW%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z)%*%sqrtW)
#         OneMinusH=(diag(nrow(Z))-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
#         OneMinusHtilde=(sqrt(1/diag(W))) * OneMinusH * matrix(diag(sqrtW),dim(OneMinusH)[1],dim(OneMinusH)[1],byrow = TRUE)
#         #equivalent to:
#         # sqrtWinv=diag(sqrt(1/diag(W)))
#         # OneMinusHtilde=sqrtWinv%*% OneMinusH%*%sqrtW
#         OneMinusHtilde=(OneMinusHtilde+t(OneMinusHtilde))/2
#         
#         deco=svd(OneMinusHtilde* matrix(diag(sqrtW),dim(OneMinusH)[1],dim(OneMinusH)[1],byrow = TRUE))
#         # equivalent to:
#         # deco=svd(OneMinusHtilde%*%sqrtW)
#         
#         deco$d[deco$d<1E-12]=0
#         return(as.vector(
#           (t(X)%*%(deco$v)%*%diag(deco$d))*
#             t(diag(deco$d)%*%t(deco$v)%*%diag(diag(W)^-1)%*%(residuals))))
#       }
# }
