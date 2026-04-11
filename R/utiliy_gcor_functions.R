
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

  i_id=which(all_vars_names==i)

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
    gcor=compute_gcor(model0 = model_i,X = mm[,i,drop=FALSE])$part_cor
}



# i and exclude are indices of the columns of model.frame x
socket_compute_gcor_normalized_conditional <- function(i,model, ...){
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
  gcor=compute_gcor_normalized_conditional(model0 = model_i,X = mm[,i_id,drop=FALSE])
}


fit_null_glm <- function(Y, Z, X,link="logit") {
  # Y: 0/1 vector length n
  # Z: n x q matrix or data.frame of confounders (DO NOT include intercept column here)
  # X: n-vector (variable of interest)
  n <- length(Y)
  # fit null model Y ~ Z (intercept included)
  fit <- glm(Y ~ Z-1, family = binomial(link = link))
  mu <- fitted(fit)
  # numerical stability
  mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
  dv=get_par_expo_fam(fit)
  if(is.null(fit$weights))  sqrtW=rep(1,length(mu)) else
    sqrtW=(dv$D*dv$V^-0.5)

  sv=svd((sqrtW)* Z,nv=0)
  H= (sv$u) %*%  t(sv$u)
  OneMinusH = diag(nrow(Z)) - H

  X_r <-  OneMinusH %*% (sqrtW * X)
  list(mu = mu, V_diag = dv$V, X_r = X_r, fit = fit, H = H)
}
