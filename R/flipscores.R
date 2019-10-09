#' @title Robust testing in GLMs
#'
#' @description Provides two robust tests for testing in GLMs, by sign-flipping score contributions. The tests are often robust against overdispersion, heteroscedasticity and, in some cases, ignored nuisance variables.
#'
#' @param formula A glm or lm object representing the fit under H0.
#' @param data A glm or lm object representing the fit under H1. It contains one covariate in addition to those in model0. If that covariate is a factor, it should have at most 2 levels. Either model1 or X1 should be specified.
#' @param family Vector with the values of the covariate of interest. If X1 is a factor, it should have at most 2 levels. Either model1 or X1 should be specified.
#' @param score_type The type of score that is computed, either "basic", "effective" or "orthogonalized". Using "effective" takes into account nuisance estimation. ORTHO?? Default is "orthogonalized".
#' @param alternative Should be "greater", "less" or "two.sided". DEFAULT IS??
#' @param w The number of times that the scores are randomly sign-flipped.Note: the maximum number of possible permutation is n!^p.
#' Where n indicates the number of observations (rows) and p indicates the number of
#' covariates (columns). R typing: factorial(n)^p. Default is 1000.
#' @param scoretype The type of score that is computed, either "basic" or "effective". Using "effective" takes into account nuisance estimation.
#'
#' @usage flipscores(model0, model1, X1, alternative = "two.sided",  w=1E5, scoretype="basic")
#'
#' @return A p-value.
#' 
#' @examples
#' set.seed(8153)
#' n <-50
#' beta <- 2  #coeffient of interest (H0: beta=0)
#' gamma <- 1   #nuisance coefficient
#' 
#' dataset <- data.frame(x=NULL,z=NULL,y=NULL)
#' dataset[n,] <- NA
#' dataset[,1:2] <- data.frame(matrix(rnorm(n*2),n,2))  #generate covariates
#' colnames(dataset) <- c("x","z")
#' dataset$y <- rnbinom(n, mu = exp(dataset[,1]*beta + dataset[,2]*gamma), size=1 )
#' 
#' #Y has a negative binomial distribution but we assume a Poisson model:
#' modelz <- glm(y~z, family=poisson,data=dataset,x=TRUE)
#' 
#' 
#' ## Basic test:
#' pv1 <- flipscores(model0=modelz, X1=dataset[,1], alternative = "two.sided", scoretype="basic")
#' pv1  #p-value
#' 
#' ## Test that takes into account nuisance estimation:
#' pv2 <- flipscores(model0=modelz, X1=dataset[,1], alternative = "two.sided", scoretype="eff")
#' pv2  #p-value
#'
#' @docType package
#'
#' @name flipscores
#'
#' @export

flipscores=