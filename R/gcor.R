#' Compute Generalized Partial Correlations for GLM terms
#'
#' This function computes the generalized partial correlation coefficient \eqn{r}
#' for each specified variable in a generalized linear model. For each variable,
#' it refits the null model excluding that variable and computes the cosine similarity
#' between the residualized predictor and standardized residuals.
#'
#' @param full_glm A fitted GLM object of class `glm`.
#' @param terms Character vector of variable names (referred to the
#'   model.matrix, that is pay attention for factors) for which to compute
#'   generalized partial correlations. If `NULL` (default), computes for all
#'   non-intercept terms in the model.
#' @param normalize FALSE by default.
#' @param intercept_too Logical indicating whether to include the intercept
#'   as a variable. Default is FALSE.
#' @param algorith Only used if \code{normalize} is \code{TRUE}. `"auto"` by default. It choose between `"intercept_only"`, `"brute_force"` and `"multi_start"`
#' @param algorithm.control Only used if \code{normalize} is \code{TRUE}. `list` of control parameters:
#' `n_exact` Integer specifying the sample size threshold for using exact
#'   methods (brute force). Default is 15.
#' `thresholds` Numeric vector of threshold values for multi-start initialization.
#' `n_random` Integer number of random starts for multi-start optimization.
#' `max_iter` Integer maximum number of iterations per start.
#' `topK` Integer number of top candidates to consider at each iteration.
#' `tol` Numeric tolerance for convergence.
#' `patience` Integer number of iterations without improvement before stopping.
#'
#' @return if \code{normalize} is \code{FALSE}, a data frame with five columns:
#'   \item{variable}{The variable name}
#'   \item{r}{The generalized partial correlation coefficient}
#'   while if \code{normalize} is \code{TRUE}
#'   \item{terms}{The variable name}
#'   \item{r}{The generalized partial correlation coefficient}
#'   \item{r_n}{The normalized generalized partial correlation coefficient}
#'   \item{null_model}{The null model used to compute the generalized (partial) correlation}
#'   \item{algorithm}{The algorithm used to compute the upper/lower bounds of the generalized partial correlation coefficient (to compute its normalized version)}
#'   \item{exact}{logical}
#'
#' @details
#' The generalized partial correlation \eqn{\r} measures the association between
#' a predictor and response after adjusting for all other terms in the model.
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
#' The normalized generalized partial correlation is computed as:
#' \deqn{
#' r_n = \begin{cases}
#' +r / r_+ & \text{if } r > 0 \\
#' -r / r_- & \text{if } r < 0
#' \end{cases}
#' }
#' where \eqn{r_+} is the maximum possible correlation and \eqn{r_-} is the minimum.
#'
#' When the (full) model has only intercept and only one predictor X, the
#' generalized (non partial) correlation is computed and the normalization factor
#' for X is exact.
#' In the more general case with more predictors, for sample sizes \eqn{n \leq n_{\text{exact}}}, brute force search is used to find
#' the exact extrema. For larger sample sizes, a greedy multi-start algorithm is employed:
#' \itemize{
#'   \item Multiple starting points are generated using thresholding and random sampling
#'   \item From each start, coordinates are greedily flipped to improve the correlation
#'   \item Early stopping is used when no improvements are found for several iterations
#'   \item The best solution across all starts is returned
#' }
#' This approach provides a good trade-off between computational efficiency and solution
#' quality for large problems where brute force is infeasible.
#'
#'
#' @author Livio Finos and Paolo Girardi
#' @examples
#' set.seed(1)
#' dt=data.frame(X=rnorm(20),
#'    Z=factor(rep(LETTERS[1:3],length.out=20)))
#' dt$Y=rpois(n=20,lambda=exp(dt$Z=="C"))
#' mod=flipscores(Y~Z+X,data=dt,family="poisson",n_flips=1000)
#' summary(mod)
#'
#' # Compute generalized partial correlations for all terms
#' (results <- gcor(mod))
#'
#' # Compute for specific terms only
#' gcor(mod, terms = c("X", "ZC"))
#'
#' gcor(mod, terms = c("X", "ZC"),normalize=TRUE)
#'
#'
#' gcor(mod, intercept_too=TRUE, normalize=TRUE)
#' set.seed(123)
#' dt=data.frame(X=rnorm(20),
#'    Z=factor(rep(LETTERS[1:3],length.out=20)))
#' dt$Y=rbinom(n=20,prob=plogis((dt$Z=="C")*2),size=1)
#' mod=flipscores(Y~Z+X,data=dt,family="binomial",n_flips=1000)
#' summary(mod)
#'
#' (results <- gcor(mod,normalize=TRUE))
#' # Compute for specific terms only
#' gcor(mod, terms = c("X", "ZC"),normalize=TRUE)
#'
#'
#' @export

#library(flipscores)
#sapply(dir("./to_flipscores/",pattern = ".R$",full.names = TRUE),source)


gcor <- function(full_glm, terms = NULL,
                 normalize=FALSE,
                 intercept_too=FALSE,
                 algorithm = "auto",
                 algorithm.control=list(n_exact = 15,
                                        thresholds = c(-.1, 0, .1),
                                        n_random = max(1,13+log(1/nrow(model.matrix(full_glm)))),
                                        max_iter = 1000,
                                        topK = max(10,min(100,length(nrow(model.matrix(full_glm)))/10)),
                                        tol = 1e-12,
                                        patience = 10)

                 ) {
  # Extract model components
  Y <- full_glm$y
  X_full <- model.matrix(full_glm)
  family <- full_glm$family
  #n <- length(Y)

  # Check if intercept is present in full model
  has_intercept <- attr(terms(full_glm$formula), "intercept") == 1

  if (has_intercept) {
    null_frml <- "~1"
  } else {
    null_frml <- "~0"
  }

  # Get all variable names (excluding intercept)
  all_vars <- colnames(X_full)
  if(!intercept_too) all_vars <- all_vars[all_vars != "(Intercept)"]

  # If terms not specified, use all non-intercept terms
  if (is.null(terms)) {
    terms <- all_vars
  } else {
    # Validate specified terms
    missing_vars <- setdiff(terms, all_vars)
    if (length(missing_vars) > 0) {
      stop("terms not found in model: ", paste(missing_vars, collapse = ", "))
    }
  }
  #

  if(normalize){
    # If terms not specified, use all non-intercept terms
    if (is.null(terms)) {
      terms <- all_vars
    } else {
      # Validate specified terms
      missing_vars <- setdiff(terms, all_vars)
      if (length(missing_vars) > 0) {
        stop("terms not found in model: ", paste(missing_vars, collapse = ", "))
      }
    }
    #

    results=lapply(terms,socket_compute_gcor_normalized_binom,
                   full_glm,algorithm=algorithm,
                   algorithm.control=algorithm.control)
    results=do.call(rbind,results)


    all_vars= c(null_frml,all_vars)

    if("(Intercept)"%in%all_vars)
      all_vars=all_vars["(Intercept)"!=all_vars]

    results <- data.frame(
      terms = paste0("~",terms),
      r=results$r,
      r_n=results$r_n,
      null_model = sapply(terms,function(i) paste0(setdiff(all_vars,i),collapse   ="+")),
      algorithm=results$algorithm,
      is.exact=results$is.exact)
    return(results)
  } else {
    results=sapply(terms,socket_compute_gcor,
                   full_glm)

    all_vars= c(null_frml,all_vars)
    results <- data.frame(
      terms = paste0("~",terms),
      gcor = results,
      null_model = sapply(terms,function(i) paste0(setdiff(all_vars,i),collapse   ="+")),
      stringsAsFactors = FALSE
    )

    rownames(results)=NULL
    return(results)
  }
}



#################
#' compute_gcor
#' @param model0 a \code{glm} object with the model under the null hypothesis (i.e. the covariates, the nuisance parameters).
#' @param model1 a \code{glm} or a \code{matrix} (or \code{vector}). If it is a \code{glm} object, it has the model under the alternative hypothesis. The terms in \code{model1} are the same terms in \code{model0} plus one or more terms to be tested.  Alternatively, if
#' \code{model1} is a \code{matrix}, it contains the tested terms column-wise.
#' @param score_type The type of score that is computed. It is "orthogonalized", "effective" or "basic".
#' "effective" and "orthogonalized" take into account the nuisance estimation.
#' @param ... other arguments.
#' @description
#' Same usage as \code{anova.glm}.
#' The parameter \code{id}  is used too,
#' if present in \code{model0} (with priority) or in \code{model1}.
#'
#' @author Paolo Girardi and Livio Finos
#' @noRd
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



compute_gcor <- function(model0, X, compute_gR2=FALSE,...){
  X=get_X(model0,X)
  Z=model.matrix(model0)
  if(is.null(model0$y)) Y=model0$model[,1] else
    Y=model0$y

  valid_idx <- complete.cases(Y, Z,X)
  Z <- Z[valid_idx, , drop = FALSE]
  X <- X[valid_idx, , drop = FALSE]
  Y <- Y[valid_idx]

  # no terms in the null model
  if(ncol(Z)==0){
    IHY <- matrix(model0$y)
    IHY <- matrix(IHY)/sqrt(sum(IHY^2))
    IHX = X
    IHY = .nrmz(IHY)
  } else { # at least one covariate
    ##  EFFECTIVE SCORE
    residuals=(Y-model0$fitted.values)
    if(is.null(model0$weights))  sqrtW=rep(1,length(residuals)) else sqrtW=(as.numeric(model0$weights)**0.5)
    if(is.null(list(...)$parms_DV))
      parms_DV<-flipscores:::get_par_expo_fam(model0) else parms_DV=list(...)$parms_DV
    D_vect<-rep(parms_DV$D,length.out=length(residuals))
    V_vect<-rep(parms_DV$V,length.out=length(residuals))
    # W = DV^−1D
    sqrtinvV_vect<-V_vect**(-0.5)
    sqrtW=(D_vect*sqrtinvV_vect)
    residuals=sqrtinvV_vect*residuals
    IHY <- matrix(residuals)
    IHY=.nrmz(IHY)
    Z_weighted <- Z * sqrtW  # Equivalent to diag(w_sqrt) %*% Z
    IH <- .get_IH(Z_weighted)
    IHX <- t(t(X*sqrtW)%*%IH)

  }

  if(compute_gR2){
    gR2 <- as.numeric(t(IHY) %*% IHX  %*% solve(t(IHX)%*%IHX) %*% t(IHX) %*% IHY)
    return(gR2)
  } else {
    part_cor=as.vector(t(IHX)%*%IHY)
    if(ncol(IHX)>1)
      part_cor <-part_cor*(colSums(IHX^2)^-.5)
    else
      part_cor <-part_cor/sqrt(sum(IHX^2))
    names(part_cor)=colnames(X)
    return(part_cor)
  }
  # return(list(IHX=IHX, IHY=IHY,part_cor=part_cor,
  #             HsqrtinvV=HsqrtinvV,
  #             IHsqrtinvV=IHsqrtinvV,
  #             sqrtV_vect=V_vect**(0.5),
  #             H=H,
  #             sqrtW=sqrtW))
}
