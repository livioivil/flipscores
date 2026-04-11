#' Normalized Generalized Partial Correlation for Binomial GLMs
#'
#' Computes the normalized generalized partial correlation coefficient for binomial GLMs.
#' The normalization scales the correlation by its maximum possible absolute value.
#'
#' @param full_glm A fitted GLM object of class `glm` with binomial family.
#' @param terms Character vector of variable names for which to compute
#'   normalized correlations. If `NULL` (default), computes for all non-intercept
#'   terms in the model.
#' @param intercept_too Logical indicating whether to include the intercept
#'   as a variable. Default is FALSE.
#'
#' @return A data frame with five columns:
#'   \item{terms}{The variable name}
#'   \item{r}{The generalized partial correlation coefficient}
#'   \item{r_n}{The normalized generalized partial correlation coefficient}
#'   \item{null_model}{The null model used to compute the generalized (partial) correlation}
#'
#' @details
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
#' @examples
#' set.seed(123)
#' dt=data.frame(X=rnorm(20),
#'    Z=factor(rep(LETTERS[1:3],length.out=20)))
#' dt$Y=rbinom(n=20,prob=plogis((dt$Z=="C")*2),size=1)
#' mod=glm(Y~Z+X,data=dt,family="binomial")
#' summary(mod)
#'
#' (results <- gcor_normalized_binom(mod))
#' # Compute for specific terms only
#' gcor_normalized_binom(mod, terms = c("X", "ZC"))
#'
#'


#' @keywords internal
#' @noRd
# Robust version that avoids formula issues
compute_gcor_normalized_conditional <- function(model0, X, ...){
  X=get_X(model0,X)
  Z=model.matrix(model0)
  r=compute_gcor(model0, X)

  #########################
  # ATTENZIONE: TUTTA DA RIVEDERE!!!
  # ATTENZIONE: IMPLEMENTARE QUI
  # Family	  Default Link	Other Available Links
  # binomial	"logit"	"probit", "cloglog", "cauchit", "log"
  # gaussian	"identity"	"log", "inverse"
  # Gamma	"inverse"	"identity", "log"
  # inverse.gaussian	"1/mu^2"	"inverse", "identity", "log"
  # poisson	"log"	"identity", "sqrt"
  # quasi	"identity"	"log", "inverse", "sqrt"
  # quasibinomial	"logit"	"probit", "cloglog", "cauchit", "log"
  # quasipoisson	"log"	"identity", "sqrt"
  #
  if (model0$family$link =="identity"){
    #  verificare per quali famiglie questo vale:
    # credo che famiglie bounded non sempre valga, ad esempio binomial
    #  (double bounded) credo valga il solito limite, per poisson con identity
    #  link immagino valga per beta +inf ma non per beta -inf
    bound = sign(r$part_cor)
  } else if(model0$family$family %in% c("binomial","quasibinomial")){
    # r$IHX = .nrmz(r$IHX) #residuals(fit_X)
    # # # Optimal Y based on sign of X_r (for maximum positive correlation)
    # # if(r$part_cor>0){
    # #   Y_opt = X > 0
    # # } else {
    # #   Y_opt = X < 0
    # # }
    # # Optimal Y based on sign of X_r (for maximum positive correlation)
    # if(r$part_cor>0){
    #   Y_opt = r$IHX > 0
    # } else {
    #   Y_opt = r$IHX < 0
    # }
    # Y_opt=.nrmz(r$sqrtinvV_vect*(r$IH%*%Y_opt))
    #
    #   bound = t(.nrmz(r$IHX))%*%.nrmz(Y_opt)

      ######################################
      #X = X
      Z = model.matrix(model0)
      w <- p0 <- predict(model0, type="response")

      # weighted mean projection helper
      WZ = function(A) crossprod(Z, w * A)
      solveZ = solve(crossprod(Z, w * Z))

      project = function(A){
        Z %*% (solveZ %*% WZ(A))
      }

      Xr = r$IHX#X - project(X)

      mx = if(r$part_cor < 0) min(Xr) else max(Xr)
      S = which(abs(Xr - mx) < 1e-8)

      v = rep(0, length(X))
      v[S] = 1

      vr = r$IH%*%v

      bound = sum(Xr * vr) / sqrt(sum(Xr^2) * sum(vr^2))
  } else if(model0$family$link == "log"){
    # ATTENZIONE: pensato per poisson e gamma, ma cosa succedere per altre family?

        # X=round(X,5)
        # if(r$part_cor<0) mx=min(X) else mx=max(X)
        # # Find set S (maximum/minimum X)
        # S = which(X == mx)
        # Y_r_opt=rep(0,length(X))
        # Y_r_opt[S]=sqrt(1/predict(model0,type = "resp")[S])

        # r$IHX=round(r$IHX,5)
        # if(r$part_cor<0) mx=min(r$IHX) else mx=max(r$IHX)
        # # Find set S (maximum/minimum X)
        # S = which(r$IHX == mx)
        # Y_r_opt=rep(0,length(r$IHX))
        # Y_r_opt[S]=predict(model0,type = "response")[S]

        # Maximum theoretical correlation
        # bound = sum( Y_r_opt[S] * r$IHX[S]) / sqrt(sum(Y_r_opt[S]^2)*sum(r$IHX^2))

    # Residualized X
    Xr = r$IHX

    # Fitted means under H0
    w <- mu0 <- predict(model0, type = "response")

    # Select extremal set
    tol = 1e-8
    mx = if(r$part_cor < 0) min(Xr) else max(Xr)
    S = which(abs(Xr - mx) < tol)

    # Build v*
    v = rep(0, length(Xr))
    v[S] = mu0[S]

    # Weighted centering (intercept-only projection)

    v_bar_w = sum(w * v) / sum(w)
    v_r = v - v_bar_w

    # Compute bound
    bound = sum(Xr * v_r) / sqrt(sum(Xr^2) * sum(v_r^2))

  } else if(model0$family$link == "inverse"){
    # ATTENZIONE: pensato per  gamma, ma cosa succedere per altre family?
    X=round(X,5)
    if(r$part_cor<0) mx=max(X) else mx=min(X)
    # Find set S (maximum/minimum X)
    S = which(X == mx)
    Y_r_opt=rep(0,length(X))
    Y_r_opt[S]=sqrt(1/predict(model0,type = "resp")[S])
    # Maximum theoretical correlation
    bound = sum( Y_r_opt[S] * r$IHX[S]) / sqrt(sum(Y_r_opt^2)*sum(r$IHX^2))

  }
  #############################

  # Compute normalized correlation
  if (r$part_cor > 0) {
    normalized_r <- r$part_cor / bound
  } else if (r$part_cor < 0) {
    normalized_r <- -r$part_cor / bound
  } else {
    normalized_r <- 0
  }

  out=data.frame(r=r$part_cor,r_n=as.vector(normalized_r))
  return(out)
}

