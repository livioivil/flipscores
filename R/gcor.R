#' Compute Generalized Partial Correlations for GLM Variables
#'
#' This function computes the generalized partial correlation coefficient \eqn{r} 
#' for each specified variable in a generalized linear model. For each variable,
#' it refits the null model excluding that variable and computes the cosine similarity
#' between the residualized predictor and standardized residuals.
#'
#' @param full_glm A fitted GLM object of class `glm`.
#' @param variables Character vector of variable names (referred to the 
#'   model.matrix, that is pay attention for factors) for which to compute 
#'   generalized partial correlations. If `NULL` (default), computes for all 
#'   non-intercept variables in the model.
#'
#' @return A data frame with two columns:
#'   \item{variable}{The variable name}
#'   \item{r}{The generalized partial correlation coefficient}
#'
#' @details
#' The generalized partial correlation \eqn{\r} measures the association between
#' a predictor and response after adjusting for all other variables in the model.
#' It is defined as the cosine similarity between the residualized predictor 
#' \eqn{X_r} and standardized residuals \eqn{Y_r}:
#' 
#' \deqn{\r = \frac{X_r^\top Y_r}{\|X_r\| \|Y_r\|}}
#' 
#' where:
#' \itemize{
#'   \item \eqn{X_r = (I - H)W^{1/2}X} is the residualized predictor
#'   \item \eqn{Y_r = V^{-1/2}(Y - \hat{\mu})} is the standardized residual vector
#'   \item \eqn{H} is the hat matrix for the nuisance covariates
#'   \item \eqn{W = DV^{-1}D} is the weight matrix
#'   \item \eqn{V} is the variance matrix and \eqn{D} is the derivative matrix
#' }
#' 
#' The function uses `flipscores:::get_par_expo_fam()` to compute \eqn{V} and \eqn{D}
#' consistently with the flipscores package methodology.
#'
#' @examples
#' set.seed(1)
#' dt=data.frame(X=rnorm(20),
#'    Z=factor(rep(LETTERS[1:3],length.out=20)))
#' dt$Y=rpois(n=20,lambda=exp(dt$Z=="C"))
#' mod=flipscores(Y~Z+X,data=dt,family="poisson",n_flips=1000)
#' summary(mod)
#'  
#' # Compute generalized partial correlations for all variables
#' (results <- gcor(mod))
#' 
#' # Compute for specific variables only
#' gcor(mod, variables = c("X", "ZC"))
#' 
#' @export

#library(flipscores)
#sapply(dir("./to_flipscores/",pattern = ".R$",full.names = TRUE),source)


gcor <- function(full_glm, variables = NULL,intercept_too=FALSE) {
  
  # Extract model components
  Y <- full_glm$y
  X_full <- model.matrix(full_glm)
  family <- full_glm$family
  n <- length(Y)
  
  # Get all variable names (excluding intercept)
  all_vars <- colnames(X_full)
  if(!intercept_too) all_vars <- all_vars[all_vars != "(Intercept)"]
  
  # If variables not specified, use all non-intercept variables
  if (is.null(variables)) {
    variables <- all_vars
  } else {
    # Validate specified variables
    missing_vars <- setdiff(variables, all_vars)
    if (length(missing_vars) > 0) {
      stop("Variables not found in model: ", paste(missing_vars, collapse = ", "))
    }
  }
  # 

  results=lapply(variables,socket_compute_gcor,
                 full_glm)
  results <- data.frame(
    variable = as.character(variables),
    r = as.numeric(results),
    stringsAsFactors = FALSE
  )

  
  return(results)
}



#################
#' compute_gcor
#' @param model0 a \code{glm} object with the model under the null hypothesis (i.e. the covariates, the nuisance parameters).
#' @param model1 a \code{glm} or a \code{matrix} (or \code{vector}). If it is a \code{glm} object, it has the model under the alternative hypothesis. The variables in \code{model1} are the same variables in \code{model0} plus one or more variables to be tested.  Alternatively, if
#' \code{model1} is a \code{matrix}, it contains the tested variables column-wise.
#' @param score_type The type of score that is computed. It is "orthogonalized", "effective" or "basic".
#' "effective" and "orthogonalized" take into account the nuisance estimation.
#' @param ... other arguments.
#' @description
#' Same usage as \code{anova.glm}.
#' The parameter \code{id}  is used too,
#' if present in \code{model0} (with priority) or in \code{model1}.
#'
#' @author Paolo Girardi and Livio Finos
#' @examples
#' library(rsq)
#' set.seed(1)
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rpois(n=20,lambda=exp(Z+.5*X))
#' model0=glm(Y~Z+1,family="poisson",x=TRUE)
#' model1=glm(Y~Z+X,family="poisson")
#' X=data.frame(X=X)
#' scr0=compute_gcor(model0 = model0,  X)
#' scr0$part_cor^2
#' sc=part_cor_new(model0 = model0,  model1)
#' sapply(c('v','kl','sse','lr','n'),function(type) rsq.partial(model1,model0,type=type)$partial.rsq)
#' plot(scr0$IHX,scr0$IHY)
#'
#' set.seed(1)
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rbinom(n=20,1,prob=plogis((-.5*Z+2*X)))
#' model0=glm(Y~Z+1,family="binomial")
#' model1=glm(Y~Z+X,family="binomial")
#' X=data.frame(X=X)
#' scr0=compute_gcor(model0 = model0,  X)
#' scr0$part_cor^2
#' sapply(c('v','kl','sse','lr','n'),function(type) rsq.partial(model1,model0,type=type)$partial.rsq)
#' plot(scr0$IHX,scr0$IHY)
#'
#' set.seed(1)
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rnorm(n=20,mean=(Z+2*X))
#' model0=lm(Y~Z+1)
#' model1=lm(Y~Z+X)
#' scr0=compute_gcor(model0 = model0,  X)
#' scr0$part_cor^2
#' cor(resid(modelX),resid(model0))^2
#' sapply(c('v','kl','sse','lr','n'),function(type) rsq.partial(model1,model0,type=type)$partial.rsq)
#' plot(scr0$IHX,scr0$IHY)
#' library(sensemakr)
#' partial_r2(model1)
#'
#'
#' 
#' library(pima)
#' data("pimads")
#' summary(pimads)
#' pimads=na.omit(pimads)
#' mod=glm(diabetes ~ . -npreg,data=pimads,family=binomial,x=TRUE)
#' scr0=compute_gcor(model0 = mod,  pimads$npreg)
#' scr0$part_cor^2
#' plot(scr0$IHX,scr0$IHY)
#' cor(scr0$IHX,scr0$IHY)^2 
#' # come mai sono diversi?
#' # perchè i residui non hanno media 0. perchè?
#' mean(scr0$IHX)
#' #[1] -0.0004750917
#' mean(scr0$IHY)
#' #[1] -2.90525e-15
#' # negli lm invece mean(scr0$IHX) è sempre 0!
#' 
#' 
#' #esempio 50% 1 in binomiale. vedi NAGELKERKE 1991
#' X= Y = rep(0:1,2)
#' model0=glm(Y~1,family="binomial")
#' model1=glm(Y~1+X,family="binomial")
#' X=data.frame(X=X)
#' scr0=compute_gcor(model0 = model0,  X)
#' scr0$part_cor^2
#' sapply(c('v','kl','sse','lr','n'),function(type) rsq.partial(model1,model0,type=type)$partial.rsq)
#' plot(scr0$IHX,scr0$IHY)

# @export

compute_gcor <- function(model0, X, ...){
  X=get_X(model0,X)
  Z=model.matrix(model0)
  
  
  # no variables in the null model
  if(ncol(Z)==0){
    IHY <- matrix(model0$y)        
    IHY <- matrix(IHY)/sqrt(sum(IHY^2))
    IHX = rep(1,nrow(IHY))
    IHX <- IHX/sqrt(nrow(IHX))
  } else { # at least one covariate
    if(is.null(model0$y)) model0$y=model0$model[,1]
    Z=Z
    ##  EFFECTIVE SCORE
    residuals=(model0$y-model0$fitted.values)
    if(is.null(model0$weights))  sqrtW=rep(1,length(residuals)) else sqrtW=(as.numeric(model0$weights)**0.5)
    if(is.null(list(...)$parms_DV))
      parms_DV<-flipscores:::get_par_expo_fam(model0) else parms_DV=list(...)$parms_DV
    D_vect<-rep(parms_DV$D,length.out=length(residuals))
    V_vect<-rep(parms_DV$V,length.out=length(residuals))
    # W = DV^−1D
    sqrtinvV_vect<-V_vect**(-0.5)
    sqrtW=diag(D_vect*sqrtinvV_vect)
    residuals=sqrtinvV_vect*residuals
    IHY <- matrix(residuals)
    IH=diag(nrow(IHY))-sqrtW%*%Z%*%solve(t(Z)%*%sqrtW^2%*%Z)%*%t(Z)%*%sqrtW 
    IHX <- t(t(X)%*%sqrtW%*%IH)
    H=sqrtW%*%Z%*%solve(t(Z)%*%(sqrtW^2)%*%Z)%*%t(Z)%*%sqrtW
    HsqrtinvV <- t(diag(sqrtinvV_vect)%*%H)
    IHsqrtinvV <- t(diag(sqrtinvV_vect)%*%(diag(nrow(IHY))-H))
    # ##  Alternativa
    # IHY <- residuals*(sqrtinvV_vect)
    # IHY <- matrix(IHY / sqrt(sum(IHY^2)))
    # IH=diag(nrow(IHY))-sqrtW%*%Z%*%solve(t(Z)%*%diag(D_vect*V_vect^-1*D_vect)%*%Z)%*%t(Z)%*%sqrtW 
    # IHX <- t(t(X)%*%sqrtW%*%IH)
    
  }
  
  part_cor=as.vector(t(IHX)%*%IHY/ sqrt(sum(residuals^2)))
  if(ncol(IHX)>1)
    part_cor <-part_cor*(colSums(IHX^2)^-.5)
  else
    part_cor <-part_cor/sqrt(sum(IHX^2))
  names(part_cor)=colnames(X)
  return(part_cor)
  # return(list(IHX=IHX, IHY=IHY,part_cor=part_cor,
  #             HsqrtinvV=HsqrtinvV,
  #             IHsqrtinvV=IHsqrtinvV,
  #             sqrtV_vect=V_vect**(0.5),
  #             H=H,
  #             sqrtW=sqrtW))
}
