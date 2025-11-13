#' Generalized R-squared for GLM Models
#'
#' Computes the generalized R-squared measure for nested generalized linear models.
#' The generalized R-squared measures the proportion of "variance" explained by
#' the additional predictors in the full model compared to the null model.
#'
#' @param full_glm A fitted GLM object of class `glm` (the full model).
#' @param null_glm A fitted GLM object of class `glm` (the null model). If `NULL` (default),
#'   uses an empty model (intercept-only if intercept is present, otherwise no predictors).
#' @param extra_terms A character vector of variable names or a formula specifying
#'   the additional terms in the full model compared to the null. If provided,
#'   this overrides `null_glm` and the null model is refitted excluding these terms.
#'
#' @return A numeric value representing the generalized R-squared measure.
#'
#' @details
#' The generalized R-squared is computed as:
#'
#' \deqn{gR^2 = \frac{Y_r^\top X_r (X_r^\top X_r)^{-1} X_r^\top Y_r}{Y_r^\top Y_r}}
#'
#' where:
#' \itemize{
#'   \item \eqn{Y_r = V^{-1/2}(Y - \hat{\mu}_0)} is the standardized residual vector from the null model
#'   \item \eqn{X_r = (I - H)W^{1/2}X} is the residualized additional predictors matrix
#'   \item \eqn{H} is the hat matrix for the null model
#'   \item \eqn{W = DV^{-1}D} is the weight matrix
#' }
#'
#' This measures the proportion of the standardized residual sum of squares explained
#' by the additional predictors in the full model.
#'
#' @examples
#' set.seed(1)
#' dt=data.frame(X=rnorm(20),
#'    Z=factor(rep(LETTERS[1:3],length.out=20)))
#' dt$Y=rpois(n=20,lambda=exp(dt$Z=="C"))
#' mod=glm(Y~Z+X,data=dt,family="poisson")
#' summary(mod)
#'
#' # Compute generalized partial correlations for all variables
#' (results <- gR2(mod))
#' # equivalente to
#' gR2(mod, extra_terms = c( "Z"))

#'
#' # Compute for specific variables only
#' gR2(mod, extra_terms = c( "Z"))
#'
#' @export
gR2 <- function(full_glm, null_glm = NULL, extra_terms = NULL) {

  # Validate full model
  if (!inherits(full_glm, "glm")) {
    stop("full_glm must be a glm object")
  }

  full_data <- full_glm$data
  full_family <- full_glm$family
  full_formula <- formula(full_glm)

  # Case 1: extra_terms provided - refit null model
  if (!is.null(extra_terms)) {
    null_glm <- create_null_from_extra_terms(full_glm, extra_terms)
  }

  # Case 2: null_glm not provided - use empty model (intercept-only if intercept present)
  if (is.null(null_glm)) {
    null_glm <- create_empty_model(full_glm)
  }

  # Validate null model
  if (!inherits(null_glm, "glm")) {
    stop("null_glm must be a glm object")
  }

  # Check model compatibility
  check_model_compatibility(full_glm, null_glm)

  # Compute generalized R-squared
  compute_gR2_from_models(full_glm, null_glm)
}

#' Create empty model (intercept-only if intercept present)
#'
#' @param full_glm The full GLM model
#' @return A fitted null GLM model with no predictors (or intercept-only)
#' @noRd
create_empty_model <- function(full_glm) {
  full_data <- full_glm$data
  full_family <- full_glm$family
  full_formula <- formula(full_glm)

  # Extract response variable
  response_var <- all.vars(full_formula)[1]

  # Check if intercept is present in full model
  has_intercept <- attr(terms(full_formula), "intercept") == 1

  if (has_intercept) {
    null_formula <- as.formula(paste(response_var, "~ 1"))
  } else {
    null_formula <- as.formula(paste(response_var, "~ 0"))
  }

  null_glm <- glm(null_formula, data = full_data, family = full_family)
  return(null_glm)
}

#' Create null model from extra terms specification
#'
#' @param full_glm The full GLM model
#' @param extra_terms Character vector or formula specifying extra terms
#' @return A fitted null GLM model
#' @noRd
create_null_from_extra_terms <- function(full_glm, extra_terms) {
  full_data <- full_glm$data
  full_family <- full_glm$family
  full_formula <- formula(full_glm)

  # Extract response variable
  response_var <- all.vars(full_formula)[1]

  # Get all terms from full model
  full_terms <- attr(terms(full_formula), "term.labels")

  # Convert extra_terms to character vector if it's a formula
  if (inherits(extra_terms, "formula")) {
    extra_terms <- all.vars(extra_terms)[-1]  # Remove response if present
  }

  # Remove terms that contain extra_terms (to handle factor variables)
  null_terms <- full_terms
  for (term in extra_terms) {
    # Remove terms that start with the variable name (handles factors like "factor_varlevel")
    null_terms <- null_terms[!grepl(paste0("^", term), null_terms)]
  }

  # Check if intercept should be included
  has_intercept <- attr(terms(full_formula), "intercept") == 1

  # Create null formula
  if (length(null_terms) == 0) {
    if (has_intercept) {
      null_formula <- as.formula(paste(response_var, "~ 1"))
    } else {
      null_formula <- as.formula(paste(response_var, "~ 0"))
    }
  } else {
    null_formula <- as.formula(paste(response_var, "~", paste(null_terms, collapse = " + ")))
  }

  # Fit null model
  null_glm <- glm(null_formula, data = full_data, family = full_family)
  return(null_glm)
}

#' Check compatibility between full and null models
#'
#' @param full_glm Full GLM model
#' @param null_glm Null GLM model
#' @return NULL, throws error if models are incompatible
#' @noRd
check_model_compatibility <- function(full_glm, null_glm) {
  # Check families
  if (!identical(full_glm$family$family, null_glm$family$family)) {
    stop("Full and null models must have the same family")
  }

  if (!identical(full_glm$family$link, null_glm$family$link)) {
    stop("Full and null models must have the same link function")
  }

  # Check if models are nested
  full_terms <- attr(terms(full_glm), "term.labels")
  null_terms <- attr(terms(null_glm), "term.labels")

  if (!all(null_terms %in% full_terms)) {
    stop("Null model is not nested within full model")
  }

  # Check response variables
  full_response <- all.vars(formula(full_glm))[1]
  null_response <- all.vars(formula(null_glm))[1]

  if (!identical(full_response, null_response)) {
    stop("Full and null models must have the same response variable")
  }

  invisible(NULL)
}

#' Compute generalized R-squared from full and null models
#'
#' @param full_glm Full GLM model
#' @param null_glm Null GLM model
#' @return Generalized R-squared value
#' @noRd
compute_gR2_from_models <- function(full_glm, null_glm) {
  # Extract components from both models
  Y <- full_glm$y
  X_full <- model.matrix(full_glm)
  X_null <- model.matrix(null_glm)
  family <- full_glm$family

  # Identify additional terms in full model
  full_terms <- colnames(X_full)
  null_terms <- colnames(X_null)
  extra_terms <- setdiff(full_terms, null_terms)

  if (length(extra_terms) == 0) {
    warning("No additional terms in full model compared to null model. Returning 0.")
    return(0)
  }

  # Get null model fitted values and components
  mu_null <- fitted(null_glm)
  eta_null <- predict(null_glm, type = "link")

  # Use flipscores internal function to get V and D under null model
  par_list <- get_par_expo_fam(null_glm)
  V <- par_list$V
  D <- par_list$D

  # Compute weights without creating diagonal matrices
  w_vals <- D^2 / V  # Vector of weights
  w_sqrt <- sqrt(w_vals)  # Vector of square roots

  # Standardized residual vector from null model
  Y_r <- (Y - mu_null) / sqrt(V)

  # Weighted design matrices using element-wise multiplication
  if (ncol(X_null) == 0) {
    H <- matrix(0, nrow = length(Y), ncol = length(Y))
  } else {
    Z_weighted <- X_null * w_sqrt  # Equivalent to diag(w_sqrt) %*% X_null
    H <- Z_weighted %*% solve(t(Z_weighted) %*% Z_weighted) %*% t(Z_weighted)
  }

  # Residualized additional predictors matrix
  X_extra <- X_full[, extra_terms, drop = FALSE]
  X_weighted <- X_extra * w_sqrt  # Equivalent to diag(w_sqrt) %*% X_extra
  X_r <- (diag(length(Y)) - H) %*% X_weighted

  # Remove any NA/NaN values
  valid_idx <- complete.cases(X_r, Y_r)
  X_r_clean <- X_r[valid_idx, , drop = FALSE]
  Y_r_clean <- Y_r[valid_idx]

  # Compute generalized R-squared using the correct formula
  if (length(Y_r_clean) == 0 || nrow(X_r_clean) == 0) {
    return(NA)
  }

  # gR^2 = Y_r^T X_r (X_r^T X_r)^{-1} X_r^T Y_r / (Y_r^T Y_r)
  YtX <- t(Y_r_clean) %*% X_r_clean
  XtX <- t(X_r_clean) %*% X_r_clean
  XtY <- t(X_r_clean) %*% Y_r_clean
  YtY <- t(Y_r_clean) %*% Y_r_clean

  # Handle case where XtX is singular
  if (rcond(XtX) < .Machine$double.eps) {
    warning("Design matrix of additional terms is singular. Generalized R-squared may be unreliable.")
  }

  gR2 <- as.numeric((YtX %*% solve(XtX) %*% XtY) / YtY)

  return(gR2)
}
