#' Normalized Generalized R-squared for Binomial GLMs
#'
#' Computes the normalized generalized R-squared coefficient for binomial GLMs.
#' The normalization scales the R-squared by its maximum possible value for
#' the specified set of extra_terms.
#'
#' @param full_glm A fitted GLM object of class `glm` with binomial family.
#' @param null_glm A fitted GLM object of class `glm` (the null model). If `NULL` (default),
#'   uses an empty model (intercept-only if intercept is present, otherwise no predictors).
#' @param extra_terms Character vector of variable names for which to compute
#'   normalized R-squared. If `NULL` (default), computes for all non-intercept
#'   extra_terms in the model.
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
#'   \item{R2}{The generalized R-squared coefficient for the set of extra_terms}
#'   \item{R2_n}{The normalized generalized R-squared coefficient}
#'   \item{algorithm}{The algorithm used to compute the maximum R-squared}
#'   \item{extra_terms_tested}{The names of the extra_terms included in the test}
#'
#' @details
#' The normalized generalized R-squared is computed as:
#' \deqn{
#' R^2_n = \frac{R^2}{R^2_{\max}}
#' }
#' where \eqn{R^2_{\max}} is the maximum possible R-squared value for the specified
#' set of extra_terms.
#'
#' Different algorithms are used based on sample size:
#' \itemize{
#'   \item For small samples (\eqn{n \leq n_{\text{exact}}}), brute force search finds the exact maximum
#'   \item For larger samples, a greedy multi-start algorithm finds approximate maximum
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' dt <- data.frame(X = rnorm(20),
#'                  Z = factor(rep(LETTERS[1:3], length.out = 20)))
#' dt$Y <- rbinom(n = 20, prob = plogis((dt$Z == "C") * 2), size = 1)
#' mod <- glm(Y ~ Z + X, data = dt, family = binomial)
#'
#' # Compute normalized generalized R-squared for all extra_terms together
#' results <- gR2_normalized_binom(mod)
#' print(results)
#' results <- gR2_normalized_binom(mod,"X")
#'
#' #################################
#' # Compute for specific set of extra_terms
#' gR2_normalized_binom(mod, extra_terms = c("X", "ZC"))
#'
#' #' # Example with logistic regression
#' n <- 50
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' Z <- rnorm(n)
#' Y <- rbinom(n, 1, plogis(0.5 + 0.3*X1 + 0.2*X2 + 0.1*Z))
#' full_glm <- glm(Y ~ X1 + X2 + Z, family = binomial)
#'
#' # Test extra_terms X1 and X2
#' result <- gR2_normalized_binom(full_glm, extra_terms = c("X1", "X2"))
#' }
#'
#' @export
gR2_normalized_binom <- function(full_glm, null_glm = NULL, extra_terms = NULL,
                                 algorithm = "auto",
                                 algorithm.control = list(n_exact = 15,
                                                          thresholds = c(-.1, 0, .1),
                                                          n_random = 10,
                                                          max_iter = 1000,
                                                          topK = 10,
                                                          tol = 1e-12,
                                                          patience = 10)) {

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
  n <- length(Y)

  # Get all variable names
  all_vars <- colnames(X_full)

  # If extra_terms not specified, use all non-intercept extra_terms
  if (is.null(extra_terms)) {
    extra_terms <- all_vars[all_vars != "(Intercept)"]
  } else {
    # Validate specified extra_terms
    missing_vars <- setdiff(extra_terms, all_vars)
    if (length(missing_vars) > 0) {
      stop("extra_terms not found in model: ", paste(missing_vars, collapse = ", "))
    }
  }

  # Create null model if not provided
  if (is.null(null_glm)) {
    null_glm <- create_empty_model(full_glm)
  }

  # Validate null model
  if (!inherits(null_glm, "glm")) {
    stop("null_glm must be a glm object")
  }

  # Check model compatibility
  check_model_compatibility(full_glm, null_glm)

  # Extract design matrices
  Z <- model.matrix(null_glm)
  X <- X_full[, extra_terms, drop = FALSE]

  # Compute observed R-squared for the set of extra_terms
  R2_obs <- gR2(full_glm, null_glm = null_glm)

  # Determine algorithm to use
  if (algorithm == "auto") {
    if (n <= control$n_exact) {
      algorithm_used <- "brute_force"
    } else {
      algorithm_used <- "multi_start"
    }
  } else {
    algorithm_used <- algorithm
  }

  # Compute maximum R-squared for the set of extra_terms
  if (algorithm_used == "brute_force") {
    if (n > 20) {
      warning("Brute force with n > 20 may be computationally expensive")
    }
    max_result <- bruteforce_R2(Z, X)
    R2_max <- max_result$R2
  } else {
    max_result <- multi_start_R2(Z, X, Y_user = Y,
                                 thresholds = control$thresholds,
                                 n_random = control$n_random,
                                 max_iter = control$max_iter,
                                 topK = control$topK,
                                 tol = control$tol,
                                 patience = control$patience,
                                 verbose = FALSE)
    R2_max <- max_result$R2
  }

  # Compute normalized R-squared
  if (R2_max > 0) {
    R2_n <- R2_obs / R2_max
  } else {
    R2_n <- 0
  }

  data.frame(
    R2 = R2_obs,
    R2_n = R2_n,
    algorithm = algorithm_used,
    null_model = deparse(null_glm$formula),
    extra_terms_tested = paste0("~ ",paste(extra_terms,collapse = " + "))
  )
}
