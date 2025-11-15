#' Normalized Generalized Partial Correlation for Binomial GLMs
#'
#' Computes the normalized generalized partial correlation coefficient for binomial GLMs.
#' The normalization scales the correlation by its maximum possible absolute value.
#'
#' @param full_glm A fitted GLM object of class `glm` with binomial family.
#' @param variables Character vector of variable names for which to compute
#'   normalized correlations. If `NULL` (default), computes for all non-intercept
#'   variables in the model.
#' @param intercept_too Logical indicating whether to include the intercept
#'   as a variable. Default is FALSE.
#' @param algorith `"auto"` by default. It choose between `"intercept_only"`, `"brute_force"` and `"multi_start"`
#' @param algorithm.control `list` of coltrol parameters:
#' `n_exact` Integer specifying the sample size threshold for using exact
#'   methods (brute force). Default is 15.
#' `thresholds` Numeric vector of threshold values for multi-start initialization.
#' `n_random` Integer number of random starts for multi-start optimization.
#' `max_iter` Integer maximum number of iterations per start.
#' `topK` Integer number of top candidates to consider at each iteration.
#' `tol` Numeric tolerance for convergence.
#' `patience` Integer number of iterations without improvement before stopping.
#'
#' @return A data frame with two columns:
#'   \item{variable}{The variable name}
#'   \item{r}{The generalized partial correlation coefficient}
#'   \item{r_n}{The normalized generalized partial correlation coefficient}
#'   \item{algorithm}{The algorithm used to compute the upper/lower bounds of the generalized partial correlation coefficient (to compute its normalized version)}
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
#' mod=flipscores(Y~Z+X,data=dt,family="binomial",n_flips=1000)
#' summary(mod)
#'
#' # Compute generalized partial correlations for all variables
#' (results <- gcor(mod))
#'
#' # Compute for specific variables only
#' gcor(mod, variables = c("X", "ZC"))
#'
#' # Compute generalized partial correlations for all variables
#' (results <- gcor_normalized_binom(mod))
#'
#' # Compute for specific variables only
#' gcor_normalized_binom(mod, variables = c("X", "ZC"))
#'
#'
#' mod=glm(Y~1+X,data=dt,family="binomial")
#' summary(mod)
#' # Compute generalized partial correlations for all variables
#' (results <- gcor(mod))
#' # Compute generalized partial correlations for all variables
#' (results <- gcor_normalized_binom(mod))
#'

#' @export
gcor_normalized_binom <- function(full_glm, variables = NULL,
                                  intercept_too = FALSE,
                                  algorithm = "auto",
                                  algorithm.control=list(n_exact = 15, thresholds = c(-.1, 0, .1),
                                  n_random = 10, max_iter = 1000, topK = 10,
                                  tol = 1e-12, patience = 10)) {

  verbose=FALSE
  # Check if model is binomial
  if (full_glm$family$family != "binomial") {
    stop("Model must be from binomial family")
  }
  # Extract control parameters with defaults
  control <- list(
    n_exact = 15,
    thresholds = c(-.1, 0, .1),
    n_random = 10,
    max_iter = 1000,
    topK = 10,
    tol = 1e-12,
    patience = 10
  )
  control[names(algorithm.control)] <- algorithm.control



  # Extract model components
  Y <- full_glm$y
  X_full <- model.matrix(full_glm)
  family <- full_glm$family
  #n <- length(Y)

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

  results=lapply(variables,socket_compute_gcor_normalized_binom,
                 full_glm,algorithm=algorithm,
                 algorithm.control=algorithm.control)
  results = do.call(rbind,results)
  results <- cbind(
    variable = as.character(variables),
    results)
  rownames(results)=NULL
  return(results)
}






# Robust version that avoids formula issues
compute_gcor_normalized_binom <- function(model0, X,
                                          algorithm="auto",
                                          algorithm.control=
                                            list(n_exact = 15,
                                                 thresholds = c(-.1, 0, .1),
                                          n_random = 10, max_iter = 1000, topK = 10,
                                          tol = 1e-12, patience = 10), ...){
  verbose=FALSE
  X=get_X(model0,X)
  Z=model.matrix(model0)
  # Check if model is binomial
  if (model0$family$family != "binomial") {
    stop("Model must be from binomial family")
  }
  r=compute_gcor(model0, X)
  # # Check if link function is odd-symmetric
  link_fun <- model0$family$link
  # odd_symmetric_links <- c("logit", "probit", "cauchit", "identity")
  # is_odd_symmetric <- link_fun %in% odd_symmetric_links


  # Check if null model has intercept only
  intercept_only <- (ncol(Z) == 1) && (length(unique(Z))==1)

  # Determine maximum and minimum correlations
  if ((algorithm=="intercept_only")||intercept_only) {
    # Use exact method for empty/null model (r_minus = -r_plus always)
    r_plus <- exact_threshold_search_intercept_only(as.numeric(X))$r
    r_minus <- -r_plus
    algorithm="intercept_only"
  } else {
    # For non-empty null models
    if(is.null(algorithm.control$n_exact)) n_exact=15 else n_exact=algorithm.control$n_exact
    if ((algorithm=="brute_force")||(length(X) <= n_exact)) {
      bound <- bruteforce_r(Z, X, mode = ifelse(r>0,"max","min"),link = link_fun)$r
      algorithm="brute_force"
    } else { #then "multi_start"
      bound <- multi_start_r(Z, X, Y_user = model0$y, mode = ifelse(r>0,"max","min"), link = link_fun,
                                  thresholds = algorithm.control$thresholds,
                                  n_random = algorithm.control$n_random,
                                  max_iter = algorithm.control$max_iter,
                                  topK = algorithm.control$topK,
                                  tol = algorithm.control$tol,
                                  patience = algorithm.control$patience,
                                  verbose = verbose)$r
      algorithm="multi_start"
    }
  }

  # Compute normalized correlation
  if (r > 0) {
    normalized_r <- r / bound
  } else if (r < 0) {
    normalized_r <- -r / bound
  } else {
    normalized_r <- 0
  }
  out=data.frame(r=r,r_n=normalized_r,algorithm=algorithm)
  return(out)
}
