#' @title Performing robust testing with GLM's estimation
#'
#' @description Used to both fit generalized linear models and perform a sign-flip score test on coefficients.
#'
#' @param formula an object of class "formula" of the model to be fitted.
#' 
#' @param data data frame used. A list or environment (or object coercible by as.data.frame to a data frame)
#' containing the variables in the model.
#' 
#' @param family distribution of the response variable. Note: in order to use Negative Binomial family, family reference must have quotes (family="negbinom"). 
#'  
#' @param cluster identifies subgrouped observations. Note: must be a vector.
#' 
#' @param score_type The type of score that is computed, either "basic", "effective" or "orthogonalized". By default is "orthogonalized". 
#' #' Using "effective" and "orthogonalized" takes into account nuisance estimation.
#' 
#' @param alternative Should be "greater", "less" or "two.sided". By default is "two.sided"
#'
#' @param n_flips The number of random flips of the scores contributions
#' When \code{n_flips} is equal or larger than the maximum number of possible flips (i.e. n^2), all possible flips are performed. 
#' Default is 5000.
#'
#' @usage flipscores_glm(formula, family, data, score_type = "orthogonalized", n_flips=1000, ...)
#'
#' @return glm class object with sign-flip score test.
#' See also the related functions (summary.flipscores, plot.flipscores, print.flipscores). 
#'
#' @examples
#' data(iris)
#' data=iris[iris$Species!="setosa",]
#' data$Species=factor(data$Species)
#' m1 = flipscores_glm(Species~.+Petal.Width*Petal.Length, family=binomial(link = "logit"), data=data, score_type="orthogonalized", n_flips=1000)
#' summary(m1)
#' 
#' data = as.data.frame(Titanic)
#' m1 = flipscores_glm(Freq~., family=poisson(link = "log"), data=data, score_type="orthogonalized", n_flips=1000)
#' summary(m1)
#'
#' @docType package
#'
#' @author Livio Finos and Vittorio Giatti
#'
#' @seealso flip
#'
#' @name flipscores_glm
#' 
#' @references "Robust testing in generalized linear models by sign-flipping score contributions" by J.Hemerik, J.Goeman and L.Finos.
#'
#' @export

flipscores_glm<-function(formula, family, data,
                         score_type = "orthogonalized",
                         n_flips=5000, 
                         cluster = NULL, 
                         ...){
  # catturo la call,
  mf <- match.call()

  if(match(c("alternative"), names(call), 0L)){
    if (alternative == "less" | alternative == "smaller") {alternative = -1}
    if (alternative == "two.sided") {alternative = 0}
    if (alternative == "greater" | alternative == "larger") {alternative = 1}}
  else alternative=0

  score_type=match.arg(score_type,c("orthogonalized","effective","basic"))
  if(missing(score_type))
    stop("test type is not specified or recognized")

  # individuo i parametri specifici di flip score
  m <- match(c("score_type","n_flips","alternative","cluster"), names(mf), 0L)
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
  
  yname=as.character(model$call$formula[[2]])
  
  #compute H0s models
  socket_compute_scores <- function(i,model){
    model$call$data=data.frame(model$y,model$x[,-i])
    model$call$formula=as.formula(paste(yname,"~0+."))
    names(model$call$data)[1]=yname
    model_i <-update(model)
    compute_scores(fit0 = model_i,X = model$x[,i,drop=FALSE],score_type=score_type)
  }
  model$scores=sapply(1:ncol(model$x),socket_compute_scores,model)
  colnames(model$scores)=colnames(model$x)
  
  ############### fit the H1 model and append the scores (refitted under H0s)
  
  ###############################
  ## compute flips
  
  ### TODO RENDERE PIù AGILE INPUT DI cluster (es formula se possibile?) 
  # + quality check
  if(!is.null(flip_param_call$cluster))
    model$scores=rowsum(model$scores,eval(flip_param_call$cluster))
  
  #  call to flip::flip()
  flip_param_call$Y=model$scores
  flip_param_call$statTest = "sum"
  results=eval(flip_param_call, parent.frame())
  
  ### output
  model$Tspace=results@permT/nrow(model$scores)
  model$p.values=flip:::p.value(results)
  model$score_type=score_type
  model$n_flips=n_flips

  
  if(is.null(param_x_ORIGINAL)||(!param_x_ORIGINAL)) model$x=NULL
  # class(model) <- 
  class(model) <- c("flipscores", class(model))
  return(model)
}
