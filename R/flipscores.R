#' @title Robust testing in GLMs, by sign-flipping score contributions
#'
#' @description Provides robust tests for testing in GLMs, by sign-flipping score contributions. The tests are often robust against overdispersion, heteroscedasticity and, in some cases, ignored nuisance variables.
#' @param score_type The type of score that is computed, either "orthogonalized", "effective" or "basic". 
#' By default is "orthogonalized". 
#' Both "orthogonalized effective score" and "effective score" take account of nuisance estimation and they provide the same test statistic. In case of small samples "effective score" might be anti-conservative. "orthogonalized effective score" gives a solution for this issue. 
#' Note that in case of a big model matrix, the "orthogonalized effective score" requires several time.
#'
#' @param n_flips The number of random flips of the scores contributions
#' When \code{n_flips} is equal or larger than the maximum number of possible flips (i.e. n^2), all possible flips are performed. 
#' Default is 5000.
#'
#' @param id a \code{vector} identifying the clustered observations. If \code{NULL} (default) observations are assumed to be independent. 
#' @param alternative Should be "greater", "less" or "two.sided" (default)
#'
#' @usage flipscores(formula, family, data, score_type = "orthogonalized", n_flips=1000, ...)
#'
#' @return glm class object with sign-flip score test.
#' See also the related functions (\code{summary.flipscores}, \code{anova.flipscores}, \code{print.flipscores}). 
#'
#' @details \code{flipscores} borrow the same parameters from function \code{glm} (and \code{glm.nb}). See these helps for more details about parameters such as \code{formula},
#' \code{data}, \code{family}. Note: in order to use Negative Binomial family, \code{family} reference must have quotes (i.e. \code{family="negbinom"}). 
#'
#' @author Livio Finos, Vittorio Giatti, Jesse Hemerik and Jelle Goeman
#'
#' @seealso \code{\link{anova.flipscores}}, \code{\link{summary.flipscores}}, \code{\link[flip]{flip}}
#'
#' @name flipscores
#' 
#' @references "Robust testing in generalized linear models by sign-flipping score contributions" by J.Hemerik, J.Goeman and L.Finos.
#' 
#' @examples
#' set.seed(1)
#' dt=data.frame(X=rnorm(20),
#'    Z=factor(rep(LETTERS[1:3],length.out=20)))
#' dt$Y=rpois(n=20,lambda=exp(Z))
#' mod=flipscores(Y~Z+X,data=dt,family="poisson")
#' summary(mod)
#' 
#'
#' @export

flipscores<-function(formula, family, data,
                         score_type = "orthogonalized",
                         n_flips=5000, 
                         id = NULL, 
                         ...){
  # catturo la call,
  fs_call <- mf <- match.call()

  if(match(c("alternative"), names(call), 0L)){
    if (alternative == "less" | alternative == "smaller") {alternative = -1}
    if (alternative == "two.sided") {alternative = 0}
    if (alternative == "greater" | alternative == "larger") {alternative = 1}}
  else alternative=0

  score_type=match.arg(score_type,c("orthogonalized","effective","basic"))
  if(missing(score_type))
    stop("test type is not specified or recognized")

  # individuo i parametri specifici di flip score
  m <- match(c("score_type","n_flips","alternative","id"), names(mf), 0L)
  m <- m[m>0]
  flip_param_call= mf[c(1L,m)]
  #rinomino la funzione da chiamare:
  flip_param_call[[1L]]=quote(flip::flip)
  names(flip_param_call)[names(flip_param_call)=="alternative"]="tail"

  # mi tengo solo quelli buoni per glm
  if(length(m)>0) mf <- mf[-m]
  
  param_x_ORIGINAL=mf$x
  #set the model to fit
  if(!is.null(mf$family)&&(mf$family=="negbinom")){
    mf[[1L]]=quote(MASS::glm.nb)
    mf$family=NULL
    } else{
      mf[[1L]]=quote(glm)
    }
  
  #compute H1 model
  mf$x=TRUE
  model <- eval(mf, parent.frame())
  
  
  #compute H0s models
  model$scores=sapply(1:ncol(model$x),socket_compute_scores,model,score_type=score_type)
  colnames(model$scores)=colnames(model$x)
  
  ############### fit the H1 model and append the scores (refitted under H0s)
  
  ###############################
  ## compute flips
  
  ### TODO RENDERE PIù AGILE INPUT DI id (es formula se possibile?) 
  # + quality check
  if(!is.null(flip_param_call$id))
    model$scores=rowsum(model$scores,eval(flip_param_call$id))
  
  #  call to flip::flip()
  flip_param_call$Y=model$scores
  flip_param_call$statTest = "sum"
  results=eval(flip_param_call, parent.frame())
  
  ### output
  model$call=fs_call
  model$id=flip_param_call$id
  model$Tspace=results@permT
  model$p.values=flip:::p.value(results)
  model$score_type=score_type
  model$n_flips=n_flips

  
  if(is.null(param_x_ORIGINAL)||(!param_x_ORIGINAL)) model$x=NULL
  # class(model) <- 
  class(model) <- c("flipscores", class(model))
  return(model)
}
