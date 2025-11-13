
# helpers.R
# Shared utilities for rho and R2 greedy optimization

rho_from_pieces <- function(N, A, S) {
  if (A <= 0 || S <= 0) return(0)
  N / (A * sqrt(S))
}

compute_rho <- function(Y, mu, Vd, Xr) {
  Yr <- (Y - mu) / sqrt(Vd)
  sum(.nrmz(Yr) * .nrmz(Xr))
}

# corrected compute_rho_current: scalar X only
compute_rho_current <- function(Y, mu, Vd, Xr) {
  # Y: 0/1 vector
  # mu: fitted probabilities (length n)
  # Vd: variances mu*(1-mu) (length n)
  # Xr: residualized predictor; must be a vector or n x 1 matrix
  
  # standardize Yr
  Yr <- (Y - mu) / sqrt(Vd)
  
  # Ensure Xr is a vector (scalar predictor)
  if (is.matrix(Xr)) {
    if (ncol(Xr) == 1) {
      Xr <- as.numeric(Xr[,1])
    } else {
      stop("compute_rho_current: X_r has more than one column. rho is defined for a single predictor. Use R2 routines for multivariate X.")
    }
  } else {
    Xr <- as.numeric(Xr)
  }
  
  rho <- sum(.nrmz(Yr)*.nrmz(Xr))
  # Numerical safety: clip to [-1, 1]
  if(rho>1) browser()
  rho <- max(-1, min(1, rho))
  
  rho
}

compute_delta_rho_vec <- function(Y, mu, Vd, Xr) {
  Yr <- (Y - mu) / sqrt(Vd)
  N <- sum(Xr * Yr)
  S <- sum(Yr^2)
  A <- sqrt(sum(Xr^2))
  n <- length(Y)
  if (A == 0 || S == 0) return(rep(-Inf, n))
  deltas <- numeric(n)
  rho_old <- N / (A * sqrt(S))
  for (i in seq_len(n)) {
    delta_i <- (1 - 2 * Y[i]) / sqrt(Vd[i])
    Snew <- S + 2 * Yr[i] * delta_i + delta_i^2
    Nnew <- N + Xr[i] * delta_i
    rho_new <- rho_from_pieces(Nnew, A, Snew)
    deltas[i] <- rho_new - rho_old
  }
  deltas
}

threshold_init_Y_Xmat <- function(Z, X, tau = 0) {
  Xr <- .get_IH(Z) %*% X
  if(length(X)==1) v1=.nrmz(Xr) else{
    eigP <- svd(Xr,nv = 0,nu=1)
    v1 <- eigP$u[, 1]
  }
  as.integer(v1 >= tau)
}

threshold_init_Y <- function(nrmzXr, tau = 0) {
  # Xr <- .get_IH(Z)%*%X
  as.integer(nrmzXr >= tau)
}

.nrmz <- function(Xr){
  Xr/sqrt(sum(Xr^2))
}

.get_H <- function(Z){
  Z %*% solve(t(Z) %*% Z) %*% t(Z)
}

.get_IH <- function(Z){
  diag(nrow(Z))-.get_H(Z)
}

# i and exclude are indices of the columns of model.frame x
socket_compute_gcor <- function(i,model){
  mm=model.matrix(model)
  all_vars_names=colnames(mm)
  if(is.numeric(i)) {
    i=all_vars_names[i]
  }
  if(is.character(i)) {
    #check if it is present
    if(!any(i%in%all_vars_names)) 
      warning("Variable(s) ",paste(sep=", ",setdiff(i,all_vars_names))," is(are) not present in the model")
  }

  i_id=sapply(i,function(ii)which(all_vars_names==ii))

  #to avoid re-run a flipscores everytime:
  attributes(model)$class= attributes(model)$class[attributes(model)$class!="flipscores"]

  
  model$call$data=data.frame(model[["y"]],mm[,-i_id,drop=FALSE])
  yname=all.vars( model$call$formula)[[1]]
  names(model$call$data)[1]=yname

  model$call$formula=as.formula(paste(yname,"~0+."))
  model$call$score_type=NULL
  model$call$n_flips = NULL
  model$call$flips = NULL

  if(length(grep("Negative Binomial",model$family$family))==1)
    model$call[[1]]=quote(glm.nb) else
      model$call[[1]]=quote(glm)
  model_i <-update(model)
    gcor=compute_gcor(model0 = model_i,X = mm[,i_id,drop=FALSE])
}



# i and exclude are indices of the columns of model.frame x
socket_compute_gcor_normalized_binom <- function(i,model,algorithm="auto",
                                                 algorithm.control=list(
                                                   n_exact = 15, 
                                                   thresholds = c(-.1, 0, .1),
                                                   n_random = 10, 
                                                   max_iter = 1000, 
                                                   topK = 10,
                                                   tol = 1e-12, 
                                                   patience = 10), ...){
  mm=model.matrix(model)
  all_vars_names=colnames(mm)
  if(is.numeric(i)) {
    i=all_vars_names[i]
  }
  if(is.character(i)) {
    #check if it is present
    if(!any(i%in%all_vars_names)) 
      warning("Variable(s) ",paste(sep=", ",setdiff(i,all_vars_names))," is(are) not present in the model")
  }
  
  i_id=sapply(i,function(ii)which(all_vars_names==ii))
  
  #to avoid re-run a flipscores everytime:
  attributes(model)$class= attributes(model)$class[attributes(model)$class!="flipscores"]
  
  
  model$call$data=data.frame(model[["y"]],mm[,-i_id,drop=FALSE])
  yname=all.vars( model$call$formula)[[1]]
  names(model$call$data)[1]=yname
  
  model$call$formula=as.formula(paste(yname,"~0+."))
  model$call$score_type=NULL
  model$call$n_flips = NULL
  model$call$flips = NULL
  
  if(length(grep("Negative Binomial",model$family$family))==1)
    model$call[[1]]=quote(glm.nb) else
      model$call[[1]]=quote(glm)
  model_i <-update(model)
  gcor=compute_gcor_normalized_binom(model0 = model_i,X = mm[,i_id,drop=FALSE],
                                     algorithm=algorithm,
                                     algorithm.control=algorithm.control)
}



# fit_null_glm.R
# Implements fit_null_glm(Y, Z, X) expected by the helper and algorithm scripts.
# - Y: binary vector (0/1), length n
# - Z: n x q matrix of confounders (without intercept)
# - X: n x p matrix of predictors of interest (can be vector or matrix)
# Returns a list with:
#   mu      : fitted probabilities (length n)
#   V_diag  : vector of variances mu*(1-mu)
#   X_r     : residualized (weighted) design for X: X_r = (I - H) W^{1/2} X
#   H       : the weighted hat matrix for the null model (n x n)
#   fit     : the glm fit object (invisible use)
# fit_null_glm <- function(Y, Z, X) {
#   n <- length(Y)
#   if (is.null(Z)) {
#     Z <- matrix(nrow = n, ncol = 0)
#   }
#   # Ensure X is matrix
#   X <- as.matrix(X)
#   # Fit null logistic regression; use glm with binomial
#   fit <- glm(Y ~ Z-1, family = binomial())
#   mu <- fitted(fit)
#   # ensure mu in (epsilon, 1-epsilon)
#   eps <- 1e-12
#   mu <- pmin(pmax(mu, eps), 1 - eps)
#   Vd <- mu * (1 - mu)
#   # Weighted matrices
#   W <- Vd
#   Ws_half <- sqrt(W)
#   Wmat_sqrt <- diag(Ws_half, n, n)
#   # Compute H = W^{1/2} X (X^T W X)^{-1} X^T W^{1/2}
#   XtWX <- t(X) %*% (W * X)   # same as t(X) %*% diag(W) %*% X but faster
#   # regularize if near-singular
#   if (rcond(XtWX) < 1e-12) {
#     XtWX <- XtWX + diag(1e-8, ncol(XtWX))
#   }
#   XtWXinv <- solve(XtWX)
#   H <- Wmat_sqrt %*% X %*% XtWXinv %*% t(X) %*% Wmat_sqrt
#   # Compute X_r = (I - H) W^{1/2} X
#   WX <- Wmat_sqrt %*% X
#   IminusH <- diag(n) - H
#   X_r <- IminusH %*% WX
#   # Return list
#   list(mu = as.numeric(mu),
#        V_diag = as.numeric(Vd),
#        X_r = X_r,
#        H = H,
#        fit = fit)
# }

fit_null_glm <- function(Y, Z, X) {
  # Y: 0/1 vector length n
  # Z: n x q matrix or data.frame of confounders (DO NOT include intercept column here)
  # X: n-vector (variable of interest)
  n <- length(Y)
  # fit null model Y ~ Z (intercept included)
  fit <- glm(Y ~ Z-1, family = binomial(link = "logit"))
  mu <- fitted(fit)
  # numerical stability
  mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
  V_diag <- mu * (1 - mu)
  W_sqrt <- sqrt(V_diag)
  # Weighted Z (for H)
  Wz <- sweep(as.matrix(Z), 1, W_sqrt, "*")
  ZZ <- crossprod(Wz)   # Z^T W Z
  invZZ <- tryCatch(solve(ZZ), error = function(e) solve(ZZ))
  H <- tcrossprod(Wz %*% invZZ, Wz)  # n x n projection in W^{1/2}-space
  X_r <- (diag(n) - H) %*% (W_sqrt * X)
  list(mu = mu, V_diag = V_diag, X_r = X_r, fit = fit, H = H)
}