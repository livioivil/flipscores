#' Robust tests for testing in GLMs, by sign-flipping score contributions
#'
#' @description Provides robust tests for testing in GLMs, by sign-flipping score contributions. The tests are often robust against overdispersion, heteroscedasticity and, in some cases, ignored nuisance variables.
#' @examples
#' set.seed(1)
#' dt=data.frame(X=rnorm(20),
#'    Z=factor(rep(LETTERS[1:3],length.out=20)))
#' dt$Y=rpois(n=20,lambda=exp(Z))
#' mod=flipscores(Y~Z+X,data=dt,family="poisson")
#' summary(mod)
#' 
#' # Anova test
#' anova(mod)
#' # or
#' mod0=flipscores(Y~Z,data=dt,family="poisson")
#' anova(mod,mod0)
#' 
#' @docType package
#'
#' @author Livio Finos, Jelle Goeman and Jesse Hemerik, with contribution of Vittorio Giatti.
#' @name flipscores-package
NULL