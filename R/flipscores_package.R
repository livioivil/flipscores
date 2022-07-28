#' Robust Score Testing in GLMs, by Sign-Flip Contributions 
#'
#' @description It provides robust tests for testing in GLMs, by sign-flipping score contributions. The tests are often robust against overdispersion, heteroscedasticity and, in some cases, ignored nuisance variables.
#' @importFrom MASS glm.nb
#' @importFrom methods is
#' @importFrom stats D as.formula model.matrix sd summary.glm update formula getCall terms
#' @examples
#' set.seed(1)
#' dt=data.frame(X=rnorm(20),
#'    Z=factor(rep(LETTERS[1:3],length.out=20)))
#' dt$Y=rpois(n=20,lambda=exp(dt$X))
#' mod=flipscores(Y~Z+X,data=dt,family="poisson",score_type = "effective",x=TRUE,)
#' summary(mod)
#' 
#' # Anova test
#' anova(mod)
#' # or
#' mod0=flipscores(Y~Z,data=dt,family="poisson",score_type = "standardized",x=TRUE)
#' anova(mod0,mod)
#' # and
#' mod0=flipscores(Y~X,data=dt,family="poisson",score_type = "standardized")
#' anova(mod0,mod)
#' @docType package
#'
#' @author Livio Finos, Jelle Goeman and Jesse Hemerik, with contribution of Vittorio Giatti.
#' @name flipscores-package
NULL