# R2_algorithm.R
# Greedy multi-start algorithm for max R^2 with early stopping
# Brute force and umbrella function included
# Requires: helpers.R and a user-provided fit_null_glm()
#
# compute_R2 <- function(Y, Z, X) {
#   info <- fit_null_glm(Y, Z, X)
#   mu <- info$mu; Vd <- info$V_diag; Xr <- as.matrix(info$X_r)
#   Yr <- (Y - mu) / sqrt(Vd)
#   # P <- .get_H(Xr)
#   # num <- as.numeric(t(Yr) %*% P %*% Yr)
#   if(is.null(ncol(Xr))||(length(Xr)==0)){
#     rho2 <- as.numeric(sum(.nrmz(Yr)*.nrmz(Xr))^2)
#     return(rho2)
#   } else {
#     den <- sum(Yr^2)
#     if (den <= 0) return(0)
#     semiP <- svd(Xr,nv=0)$u
#     num <- as.numeric(sum((t(Yr) %*% semiP)^2))
#     return(num / den)
#   }
# }
#

compute_R2_and_info <- function(Y, Z, X,link) {
  info <- fit_null_glm(Y, Z, X,link=link)
  mu <- info$mu; Vd <- info$V_diag; Xr <- as.matrix(info$X_r)
  Yr <- (Y - mu) / sqrt(Vd)
  Yr <- .nrmz(Yr)
  if(is.null(ncol(Xr))||(length(Xr)==0)){
    semiP=matrix(.nrmz(Xr))
    R2 <- as.numeric(sum(Yr*semiP)^2)
    return(list(R2=R2, mu = mu, Vd = Vd, semiP = semiP))
  } else {
    semiP <- svd(Xr,nv=0)$u
    R2 <- as.numeric(sum((t(Yr) %*% semiP)^2))
    return(list(R2=R2, mu = mu, Vd = Vd, semiP = semiP))
  }
}

compute_R2_norefit <- function(Y,info) {
  Yr <- (Y - info$mu) / sqrt(info$Vd)
  Yr <- .nrmz(Yr)
  # P <- .get_H(Xr)
  # num <- as.numeric(t(Yr) %*% P %*% Yr)
  if(is.null(ncol(info$semiP))||(length(info$semiP)==0)){
    R2 <- as.numeric(sum(Yr*info$semiP)^2)
    return(R2)
  } else {
    R2 <- as.numeric(sum((t(Yr) %*% info$semiP)^2))
    return(R2)
  }
}

greedy_optimize_R2 <- function(Y_start, Z, X, link, max_iter = 1000,
                               topK = 5, tol = 1e-12, patience = 10,
                               verbose = TRUE) {
  Y <- as.integer(Y_start)
  n <- length(Y)
  iter <- 0
  no_improve <- 0
  fit_curr <- compute_R2_and_info(Y, Z, X,link=link)
  R2_curr <- fit_curr$R2

  while (iter < max_iter && no_improve < patience) {
    iter <- iter + 1
    if (verbose)
      cat(sprintf("Iter %d: initial, R2 = %.8f\n", iter, R2_curr))
    gains <- rep(-Inf, n)
    for (i in seq_len(n)) {
      Yt <- Y; Yt[i] <- 1 - Yt[i]
      R2_new <- compute_R2_norefit(Yt,fit_curr)
      gains[i] <- R2_new - R2_curr
    }
    best_idx <- order(gains, decreasing = TRUE)[1:min(topK, length(gains))]

    improved <- FALSE
    for (i in best_idx) {
      if (gains[i] <= tol) break
      Yt <- Y; Yt[i] <- 1 - Yt[i]
      fit_last <- compute_R2_and_info(Yt, Z, X,link=link)
      R2_new <- fit_last$R2
      if (R2_new > R2_curr + tol) {
        Y <- Yt
        R2_curr <- R2_new
        fit_curr <- fit_last
        improved <- TRUE
        no_improve <- 0
        if (verbose)
          cat(sprintf("Iter %d: flipped %d, R2 = %.8f\n", iter, i, R2_curr))
        break
      }
    }
    if (!improved) no_improve <- no_improve + 1
  }
  list(Y = Y, R2 = R2_curr, iterations = iter)
}

multi_start_R2 <- function(Z, X, Y_user = NULL,
                           link="logit",
                           thresholds = c(-.1,0,.1),
                           n_random = 10,
                           max_iter = 100,
                           topK = max(10,min(100,length(X)/10)),
                           tol = 1e-12,
                           patience = 10,
                           verbose = TRUE) {
  n <- nrow(Z)
  inits <- list()
  if (!is.null(Y_user)) inits[["user"]] <- as.integer(Y_user)

  Xr <- .get_IH(Z) %*% X
  if(length(X)==1) v1=.nrmz(Xr) else{
    eigP <- svd(Xr,nv = 0)
    v1 <- eigP$u
    coefs=.nrmz(rnorm(ncol(v1)))
    v1 <- v1%*%coefs
  }
  for (tau in thresholds) inits[[paste0("thr", tau)]] <-   as.integer(v1 >= tau)
  for (r in 1:n_random) inits[[paste0("rand", r)]] <- rbinom(n, 1, 0.5)

  best <- list(R2 = -Inf, Y = NULL)
  for (nm in names(inits)) {
    if (verbose)
      cat("\n evaluating", nm, "\n")
    res <- greedy_optimize_R2(inits[[nm]], Z, X,link=link,
                              max_iter = max_iter, topK = topK,
                              tol = tol, patience = patience,
                              verbose = verbose)
    if (res$R2 > best$R2) {
      best <- res
      if (verbose)
        cat("New best R2 =", round(best$R2, 6), "from", nm, "\n")
    }
  }
  best
}

# ----------- Brute force for small n -----------

bruteforce_R2 <- function(Z, X) {
  n <- nrow(Z)
  best_val <- -Inf
  best_Y <- NULL
  for (num in 0:(2^n - 1)) {
    Y <- as.integer(intToBits(num)[1:n])
    val <- compute_R2(Y, Z, X)
    if (val > best_val) {
      best_val <- val
      best_Y <- Y
    }
  }
  list(Y = best_Y, R2 = best_val)
}


