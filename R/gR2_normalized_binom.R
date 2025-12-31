#' Normalized Generalized R-squared for Binomial GLMs
#'
#' Computes the normalized generalized R-squared coefficient for binomial GLMs.
#' The normalization scales the R-squared by its maximum possible value for
#' the specified set of terms.
#'
#' @param full_glm A fitted GLM object of class `glm` with binomial family.
#' @param null_glm A fitted GLM object of class `glm` (the null model). If `NULL` (default),
#'   uses an empty model (intercept-only if intercept is present, otherwise no predictors).
#' @param terms Character vector of variable names for which to compute
#'   normalized R-squared. If `NULL` (default), computes for all non-intercept
#'   terms in the model.
#' @param algorithm `"auto"` by default. It chooses between `"brute_force"` and `"multi_start"`
#' @param algorithm.control `list` of control parameters:
#'   \itemize{
#'     \item `n_exact` Integer specifying the sample size threshold for using exact
#'       methods (brute force). Default is 15.
#'     \item `thresholds` Numeric vector of threshold values for multi-start initialization.
#'     \item `n_random` Integer number of random starts for multi-start optimization.
#'     \item `max_iter` Integer maximum number of iterations per start.
#'     \item `topK` Integer number of top candidates to consider at each iteration.
#'     \item `tol` Numeric tolerance for convergence.
#'     \item `patience` Integer number of iterations without improvement before stopping.
#'   }
#'
#' @return A list with components:
#'   \item{R2}{The generalized R-squared coefficient for the set of terms}
#'   \item{R2_n}{The normalized generalized R-squared coefficient}
#'   \item{algorithm}{The algorithm used to compute the maximum R-squared}
#'   \item{terms_tested}{The names of the terms included in the test}
#'
#' @details
#' The normalized generalized R-squared is computed as:
#' \deqn{
#' R^2_n = \frac{R^2}{R^2_{\max}}
#' }
#' where \eqn{R^2_{\max}} is the maximum possible R-squared value for the specified
#' set of terms.
#'
#' Different algorithms are used based on sample size:
#' \itemize{
#'   \item For small samples (\eqn{n \leq n_{\text{exact}}}), brute force search finds the exact maximum
#'   \item For larger samples, a greedy multi-start algorithm finds approximate maximum
#' }
#'
#' @examples
#' set.seed(123)
#' dt <- data.frame(X = rnorm(20),
#'                  Z = factor(rep(LETTERS[1:3], length.out = 20)))
#' dt$Y <- rbinom(n = 20, prob = plogis((dt$Z == "C") * 2), size = 1)
#' mod <- glm(Y ~ Z + X, data = dt, family = binomial)
#'
#' # Compute generalized partial correlations for all variables
#' (results <-  gR2_normalized_binom(mod))
#' # equivalent to
#' mod0=glm(Y~1,data=dt,family=binomial)
#' (results <-  gR2_normalized_binom(mod, mod0))
#'
#' # Compute for specific variables only
#' (results <-  gR2_normalized_binom(mod,terms = c("X","Z")))
#'
#'
#'
#' @export
# gR2_normalized_binom <- function(full_glm, null_glm = NULL, terms = NULL,
#                                  algorithm = "auto",
#                                  algorithm.control = list(n_exact = 15,
#                                                           thresholds = c(-.1, 0, .1),
#                                                           n_random = 10,
#                                                           max_iter = 1000,
#                                                           topK = 10,
#                                                           tol = 1e-12,
#                                                           patience = 10)) {
#
#
#   # Extract control parameters with defaults
#   control <- list(
#     n_exact = 15,
#     thresholds = c(-.1, 0, .1),
#     n_random = 10,
#     max_iter = 1000,
#     topK = 10,
#     tol = 1e-12,
#     patience = 10
#   )
#   control[names(algorithm.control)] <- algorithm.control
#
#
#   if(length(terms)>1){
#     temp=lapply(terms,.socket_compute_gR2_n,full_glm, null_glm)
#     return(do.call(rbind,temp))
#   } else
#     return(.socket_compute_gR2_n(terms,full_glm, null_glm))
#
# }

.socket_compute_gR2_n <- function(terms,full_glm, null_glm = NULL,algorithm,control){
  temp=.prepare_for_gR2(full_glm,null_glm,terms)
  Z=model.matrix(temp$null_glm)
  Y=temp$null_glm$y
  X=temp$X
  n=length(Y)
  # Compute generalized R-squared
  gR2=compute_gR2(temp$null_glm, X)

  if (full_glm$family$family != "binomial") {
    #warning("When normalize==TRUE, Model must be from binomial family")

    has_NO_intercept=!.intercept_in_Z_and_count_family(temp$null_glm)
    data.frame(
      terms = paste0("~ ",paste(colnames(temp$X),collapse = " + ")),
      gR2 = gR2,
      gR2_n = ifelse(has_NO_intercept,NA,gR2),
      algorithm = ifelse(has_NO_intercept,NA,"from theory"),
      exact=ifelse(has_NO_intercept,NA,TRUE),
      null_model = deparse(temp$null_glm$formula),
      stringsAsFactors = FALSE
    )
  } else { #is binomial family

    # Determine algorithm to use
    if (algorithm == "auto") {
      # Check if null model has intercept only
      intercept_only <- (ncol(Z) == 1) && (length(unique(Z))==1)
      algorithm_used="intercept_only"
      if (n <= control$n_exact) {
        algorithm_used <- "brute_force"
      } else {
        algorithm_used <- "multi_start"
      }
    } else {
      algorithm_used <- algorithm
    }

    # Compute maximum R-squared for the set of terms
    if (algorithm_used == "brute_force") {
      if (n > 20) {
        warning("Brute force with n > 20 may be computationally expensive")
      }
      max_result <- bruteforce_R2(Z, X)
      gR2_max <- max_result$R2
    } else if (algorithm_used == "brute_force") {
      if (n > 20) {
        warning("Brute force with n > 20 may be computationally expensive")
      }
      max_result <- exact_threshold_search_intercept_only(as.numeric(X))
      gR2_max <- max_result$R2
    } else {
      algorithm=="multi_start"
      max_result <- multi_start_R2(Z, X, Y_user = Y,
                                   link=full_glm$family$link,
                                   thresholds = control$thresholds,
                                   n_random = control$n_random,
                                   max_iter = control$max_iter,
                                   topK = control$topK,
                                   tol = control$tol,
                                   patience = control$patience,
                                   verbose = FALSE)
      gR2_max <- max_result$R2
    }

    # Compute normalized R-squared
    if (gR2_max > 0) {
      gR2_n <- gR2 / gR2_max
    } else {
      gR2_n <- 0
    }


    data.frame(
      terms = paste0("~ ",paste(colnames(temp$X),collapse = " + ")),
      gR2 = gR2,
      gR2_n = gR2_n,
      algorithm = algorithm_used,
      exact=ifelse(algorithm_used=="multi_start",FALSE,TRUE),
      null_model = deparse(temp$null_glm$formula),
      stringsAsFactors = FALSE
    )
  }
}
