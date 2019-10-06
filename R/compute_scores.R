#' compute_scores
#' @examples
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rpois(n=20,lambda=exp(Z+X))
#' mod0=glm(Y~Z,family="poisson")
#' mod1=glm(Y~Z+X,family="poisson")
#' (scr=compute_scores(fit0 = mod0,))
#' flip:::flip(scr)
#' @export

compute_scores <- function(fit0, X,score_type="orthogonalized"){

  a <- get_a_expo_fam(fit0)
  ###############
  Z=fit0$x
  residuals=(fit0$y-fit0$fitted.values)/a
  W=diag(as.numeric(fit0$weights))


  #BASIC SCORE
  if(score_type=="basic"){
    return(as.vector(X*(residuals)))
  } else
    ##  EFFECTIVE SCORE
    if(score_type=="effective"){
      OneMinusH=(nrow(Z)-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
      return(as.vector(
        (t(X)%*%OneMinusH*(residuals))))
    } else
      #ORTHO EFFECTIVE SCORE
      if(score_type=="orthogonalized"){
        sqrtW=diag(sqrt(diag(W)))
        sqrtWinv=diag(sqrt(1/diag(W)))
        # OneMinusHtilde_0=(diag(n)-sqrtW%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z)%*%sqrtW)
        OneMinusH=(diag(nrow(Z))-W%*%Z%*%solve(t(Z*diag(W))%*%Z)%*%t(Z))
        OneMinusHtilde=sqrtWinv%*% OneMinusH%*%sqrtW
        #equivalent:
        # print(OneMinusHtilde_0-OneMinusHtilde)
        OneMinusHtilde=(OneMinusHtilde+t(OneMinusHtilde))/2
        deco=svd(OneMinusHtilde%*%sqrtW)
        deco$d[deco$d<1E-12]=0
        return(as.vector(
          (t(X)%*%(deco$v)%*%diag(deco$d))*
            t(diag(deco$d)%*%t(deco$v)%*%diag(diag(W)^-1)%*%(residuals))))
      }
}
