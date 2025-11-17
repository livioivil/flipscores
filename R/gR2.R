#' Generalized R-squared for GLM Models
#'
#' Computes the generalized R-squared measure for nested generalized linear models.
#' The generalized R-squared measures the proportion of "variance" explained by
#' the additional predictors in the full model compared to the null model.
#'
#' @param full_glm A fitted GLM object of class `glm` (the full model).
#' @param null_glm A fitted GLM object of class `glm` (the null model). If `NULL` (default),
#'   uses an empty model (intercept-only if intercept is present, otherwise no predictors).
#' @param terms A character vector of variable names or a formula specifying
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
#' @author Livio Finos and Paolo Girardi
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
#' # equivalent to
#' mod0=glm(Y~1,data=dt,family="poisson")
#' (results <- gR2(mod, mod0))
#'
#' # Compute for specific variables only
#' (results <- gR2(mod,terms = c("X","Z")))
#'
#' @export

gR2 <- function(full_glm, null_glm = NULL, terms = NULL) {
  .socket_compute_gR2 <- function(terms,full_glm, null_glm = NULL){
    temp=.prepare_for_gR2(full_glm,null_glm,terms)
    # Compute generalized R-squared
    gR2=compute_gR2(temp$null_glm, temp$X)
    data.frame(
      terms = paste0("~ ",paste(temp$terms,collapse = " + ")),
      gR2 = gR2,
      null_model = deparse(temp$null_glm$formula),
      stringsAsFactors = FALSE
    )
  }

    if(length(terms)>1){
   temp=lapply(terms,.socket_compute_gR2,full_glm, null_glm)
   return(do.call(rbind,temp))
  } else
  return(.socket_compute_gR2(terms,full_glm, null_glm))
}


.prepare_for_gR2 <- function(full_glm,null_glm,terms){
  # Validate full model
  if (!inherits(full_glm, "glm")) {
    stop("full_glm must be a glm object")
  }

  full_data <- full_glm$data
  full_family <- full_glm$family
  full_formula <- formula(full_glm)

  # Case 1: terms provided - refit null model
  if (!is.null(terms)) {
    null_glm <- create_null_from_terms(full_glm, terms)
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


  # Identify additional terms in full model

  temp=get_X(null_glm,full_glm,return_extra_terms=TRUE)

  return(list(null_glm=null_glm, X=temp$X,terms=temp$terms))
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
#' @param terms Character vector or formula specifying extra terms
#' @return A fitted null GLM model
#' @noRd
create_null_from_terms <- function(full_glm, terms) {
  full_data <- full_glm$data
  full_family <- full_glm$family
  full_formula <- formula(full_glm)

  # Extract response variable
  response_var <- all.vars(full_formula)[1]

  # Get all terms from full model
  full_terms <- attr(terms(full_formula), "term.labels")

  # Convert terms to character vector if it's a formula
  if (inherits(terms, "formula")) {
    terms <- all.vars(terms)[-1]  # Remove response if present
  }

  # Remove terms that contain terms (to handle factor variables)
  null_terms <- full_terms
  for (term in terms) {
    # Remove terms that start with the variable name (handles factors like "factor_varlevel")
    null_terms <- null_terms[!grepl(paste0("^", term), null_terms)]
  }

  # Check if intercept should be included
  has_intercept <- attr(terms(full_formula), "intercept") == 1

  # Create null formula
  if (length(null_terms) == 0) {
    if (has_intercept) {
      null_formula <- as.formula(paste(response_var, "~ 1"))
      null_formula <- as.formula(paste(response_var, "~ 0"))
    }
  } else {
    null_formula <- as.formula(paste(response_var, "~", paste(null_terms, collapse = " + ")))
  }

  # Avoid to re-run the flipscores if it is a flipscores-object
  attributes(full_glm)$class= attributes(full_glm)$class[attributes(full_glm)$class!="flipscores"]
  full_glm$call$score_type=NULL
  full_glm$call$n_flips = NULL
  full_glm$call$flips = NULL

  if(length(grep("Negative Binomial",full_glm$family$family))==1)
    full_glm$call[[1]]=quote(glm.nb) else
      full_glm$call[[1]]=quote(glm)
  # Fit null model
  null_glm <-update(full_glm,formula. = null_formula)
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
compute_gR2 <- function(null_glm,X) {
  # Extract components from both models
  Y <- null_glm$y
  Z <- model.matrix(null_glm)
  family <- null_glm$family


  # Remove any NA/NaN values
  valid_idx <- complete.cases(Y, Z,X)
  Z <- Z[valid_idx, , drop = FALSE]
  X <- X[valid_idx, , drop = FALSE]
  Y <- Y[valid_idx]


  # Get null model fitted values and components
  mu_null <- fitted(null_glm)

  # Use flipscores internal function to get V and D under null model
  par_list <- get_par_expo_fam(null_glm)
  V <- par_list$V
  D <- par_list$D

  # Compute weights without creating diagonal matrices
  w_sqrt <- D / V^.5  # Vector of square roots

  # Standardized residual vector from null model
  Y_r <- (Y - mu_null) / sqrt(V)

  # Weighted design matrices using element-wise multiplication
  if (ncol(Z) == 0) {
    IH <- diag(length(Y))
  } else {
    Z_weighted <- Z * w_sqrt  # Equivalent to diag(w_sqrt) %*% Z
    IH <- .get_IH(Z_weighted)
  }

  # Residualized additional predictors matrix
  X_weighted <- X * w_sqrt  # Equivalent to diag(w_sqrt) %*% X
  X_r <- IH %*% X_weighted


  # Compute generalized R-squared using the correct formula
  if (length(Y_r) == 0 || nrow(X_r) == 0) {
    return(NA)
  }

  # gR^2 = Y_r^T X_r (X_r^T X_r)^{-1} X_r^T Y_r / (Y_r^T Y_r)
  # YtX <- t(Y_r_clean) %*% X_r_clean
  X_r <- IH%*%X
  XtX <- t(X_r)%*%X_r
  # XtY <- t(X_r_clean) %*% Y_r_clean
  YtY <- t(Y_r) %*% Y_r

  # Handle case where XtX is singular
  if (rcond(XtX) < .Machine$double.eps) {
    warning("Design matrix of additional terms is singular. Generalized R-squared may be unreliable.")
  }

  gR2 <- as.numeric(t(Y_r) %*% X_r  %*% solve(XtX) %*% t(X_r) %*% Y_r / (YtY))

  return(gR2)

}
