#' Robust score testing for linear models via sign-flipping
#'
#' Performs robust score tests for one or more coefficients in a linear model
#' (ordinary least squares) by sign‑flipping score contributions. The test is
#' valid under heteroscedasticity, overdispersion and, in some cases, ignored
#' nuisance variables. It is an efficient alternative to \code{flipscores()}
#' when the model is Gaussian and the link is identity, the \code{score_type}
#' implemented is the \code{"standardized"}.
#'
#' @param formula A formula describing the linear model (e.g. \code{y ~ x1 + x2}).
#'   The response can be a numeric vector or a matrix.
#' @param data An optional data frame containing the variables in the model.
#'   If \code{NULL}, variables are taken from the environment of the formula.
#' @param n_flips Number of random sign flips. Defaults to \code{5000}.
#'   When \code{n_flips >= n^2} (with \eqn{n} the number of observations)
#'   all possible flips are performed, giving an exact permutation test.
#' @param alternative Alternative hypothesis. Must be one of \code{"two.sided"}
#'   (default), \code{"greater"} or \code{"less"}.
#' @param seed Optional integer seed for reproducibility of the random flips.
#' @param tested_coeffs Character vector naming the coefficients to test.
#'   If \code{NULL} (default), all coefficients in the model are tested.
#' @param flips Optional matrix of pre‑computed sign flips. It must have
#'   \code{n_flips} rows and \eqn{n} columns, where each entry is \eqn{+1} or
#'   \eqn{-1}. If supplied, \code{n_flips} is ignored.
#' @param ... Additional arguments passed to internal functions (currently not used).
#'
#' @details
#' The function works by first projecting the response and the tested covariate
#' onto the orthogonal complement of the nuisance design matrix. The resulting
#' score contributions are then randomly sign‑flipped (with fixed signs across
#' observations) to generate a permutation distribution of the test statistic.
#'
#' The observed score statistic for a coefficient is
#' \eqn{X_r^\top Y_r / \sqrt{\sum X_r^2}}, where \eqn{X_r} and \eqn{Y_r} are the
#' residualised covariate and residualised response, respectively. Under the
#' null hypothesis the sign of each observation’s contribution is arbitrary,
#' so repeatedly flipping the signs yields an exact conditional test.
#'
#' The method is robust to misspecification of the variance structure and does
#' not require the errors to be normally distributed or homoscedastic.
#'
#' @return An object of class \code{c("fs_lm", ...)} with the following
#' components:
#' \describe{
#'   \item{\code{Tspace}}{A matrix of dimension \code{n_flips} by \code{k},
#'     where \code{k = length(tested_coeffs)}. The first row contains the
#'     observed test statistics; the remaining rows contain the test statistics
#'     obtained after random sign flips.}
#'   \item{\code{summary_table}}{A data frame with one row per tested coefficient,
#'     containing the observed score, the two‑sided or one‑sided p‑value,
#'     and optionally other summary statistics (e.g. standardised score).}
#'   \item{\code{info}}{A list with the original formula, the names of the
#'     design matrix columns, and the names of the response variables.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @references
#' Hemerik, J., Goeman, J. J., & Finos, L. (2020). Robust testing in generalized
#' linear models by sign‑flipping score contributions. \emph{Journal of the Royal
#' Statistical Society Series B: Statistical Methodology}, 82(3), 841‑864.
#' \doi{10.1111/rssb.12369}
#'
#' De Santis, R., Goeman, J. J., Hemerik, J., Davenport, S., & Finos, L. (2025).
#' Inference in Generalized Linear Models with Robustness to Misspecified Variances.
#' \emph{Journal of the American Statistical Association}, 1‑10.
#' \doi{10.1080/01621459.2025.2491775}
#'
#' @seealso
#' \code{\link{flipscores}} for the general GLM version;
#' \code{\link{summary.fs_lm}} for summarising results;
#' \code{\link{anova.fs_lm}} for joint tests across coefficients.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' data <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   y  = 2 * rnorm(n, mean = 0, sd = 1)   # pure noise (null holds)
#' )
#' data$y2  = data$x1+  rnorm(n, mean = 0, sd = 1)   # pure noise (null holds)
#'
#' # Test the coefficient of x1
#' res <- fs_lm(cbind(y,y2) ~ x1 + x2, data = data, n_flips = 999, alternative = "two.sided")
#' flipscores:::summary.fs_lm(res)
#'
#' # Test only x2
#' res2 <- fs_lm(y ~ x1 + x2, data = data, tested_coeffs = "x2")
#' flipscores:::summary.fs_lm(res2)
#'
#' @export
#'
flipscores_lm <- function(formula, data = NULL, n_flips = 5000,
                  alternative = "two.sided", seed = NULL, tested_coeffs = NULL,
                  flips = NULL, ...) {
original_call <- match.call()


if(!is.null(seed)) set.seed(seed)

D <- formula_to_matrices(formula, data = data)

n_obs=max(ifelse(is.null(nrow(D$X)), length(D$X), nrow(D$X)),
          ifelse(is.null(nrow(D$Y)), length(D$Y), nrow(D$Y)))

if(!exists("obs_names")) obs_names=rownames(D$X)
if(is.null(flips)){
  flips=make_flips(n_obs=n_obs,n_flips=n_flips,obs_names=obs_names)
}

flips=flips[,obs_names]

out=.fs_lm(D, flips=flips,alternative=alternative,
          tested_coeffs=tested_coeffs)

out$call <- original_call
class(out) <- c("fs_lm", class(out))
return(out)
}


# ###################################
formula_to_matrices <- function(formula, data) {

  # Build the model frame (handles NA, subset, etc.)
  mf <- model.frame(formula, data = data)

  # Right-hand side: design matrix (X), includes intercept by default
  X <- model.matrix(formula, data = mf)

  # Left-hand side: response matrix (Y)
  Y <- model.response(mf)
  if(is.vector(Y)) {
    Y= as.matrix(Y)
    colnames(Y)=as.character(formula[[2]])
  }

  list(Y = Y, X = X)
}


.fs_lm <- function(D,flips, alternative,
                  tested_coeffs=NULL){
  names_X=colnames(D$X)
  if(is.null(tested_coeffs)) tested_coeffs=names_X

  scores=lapply(tested_coeffs,function(i).get_scores(X=D$X[,i,drop=FALSE],
                                                     Y=D$Y,
                                                     Z=D$X[,setdiff(names_X,i),drop=FALSE]))

  Tspace=lapply(scores,.flip_test_lm,
                flips=flips)
  names(Tspace) <- names(scores) <- tested_coeffs


  summary_table=lapply(names(scores),function(i){
    temp=list(scores=scores[[i]],
              Tspace=Tspace[[i]],
              alternative=alternative)
    temp=.get_summary_table_from_fs_lm(temp)
    cbind(temp[,1:2,drop=FALSE],coefficient=i,temp[,-(1:2),drop=FALSE])
  })

  Tspace=do.call(cbind,Tspace)
  summary_table=do.call(rbind,summary_table)
  rownames(summary_table)=NULL
  list(Tspace=Tspace,
       summary_table=summary_table,
       info=list(formula=formula,
                D=D))
}



# for standardized (see in flipscores):
.score_std=function(flp,scores_objs) {
  # scr_eff # un vettore
  numerator=crossprod(flp,scores_objs$scores) #t(scr_eff)%*%flp
  if (all(sign(flp)==1)|(all(sign(flp)==-1))){
    denominator = 1
  } else {
    denominator = 1 - sum((colSums(scores_objs$vars_objs$A[flp==1,,drop=FALSE])
                           -colSums(scores_objs$vars_objs$A[flp==-1,,drop=FALSE]))^2)
  }
  as.vector(numerator/((denominator)**0.5))
}

######################
#X solo colonna
.get_scores<- function(X,Y,Z){
  Yr <- .get_IH(Z)%*%Y
  Q=qr.Q(qr(Z))
  Xr=crossprod(diag(nrow(Z))-tcrossprod(Q),X)
  m = sum(Xr^2)
  # we divide it by sqrt(m) which is the sd scaling factor of the observed test stat (i.e. effective and standardized have the same observed test stat)
  A=Xr[,]*Q/sqrt(m)
  scores=Xr[,]*Yr

  attr(scores,"scale_objects")=list(A=A)
  scores
}

###################
