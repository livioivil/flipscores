# rho_algorithms.R
# Greedy multi-start algorithms for max and min r with early stopping
# Brute force and umbrella function included
# Requires: helpers.R and a user-provided fit_null_glm()

compute_delta_r_vec <- function(Y, mu, V_d, X_r) {
  # This function computes the approximate change in r from flipping each coordinate of Y

  # For each i = 1, ..., n:
  # 1. Let Y^{(i)} be Y with the i-th coordinate flipped
  # 2. Compute the change in standardized residual:
  #    ΔY_{r,i} = (1 - 2Y_i) / sqrt(V_{d,i})
  # 3. Approximate the change in r:
  #    δ_i = X_{r,i} · ΔY_{r,i}

  n <- length(Y)
  delta <- numeric(n)

  for (i in 1:n) {
    delta_Y <- (1 - 2 * Y[i]) / sqrt(V_d[i])
    delta[i] <- X_r[i] * delta_Y
  }

  return(delta)
}

fit_compute_r <- function(Y, Z, Xmat,link="logit") {
  info <- fit_null_glm(Y, Z, Xmat,link=link)
  compute_r(Y, info$mu, info$V_diag, info$X_r)
}

compute_r <- function(Y, mu, Vd, Xr) {
  Yr <- (Y - mu) / sqrt(Vd)
  num <- sum(.nrmz(Yr) * .nrmz(Xr))
}


compute_r_and_info <- function(Y, Z, X,link="logit") {
  info <- fit_null_glm(Y, Z, X,link=link)
  mu <- info$mu; Vd <- info$V_diag; Xr <- as.matrix(info$X_r)
  Yr <- (Y - mu) / sqrt(Vd)
  Yr <- .nrmz(Yr)
  Xr=.nrmz(Xr)
  r <- as.numeric(sum(Yr*Xr))
  return(list(r=r, mu = mu, Vd = Vd, Xr = Xr))
}

compute_r_norefit <- function(Y,info) {
  Yr <- (Y - info$mu) / sqrt(info$Vd)
  Yr <- .nrmz(Yr)
  r <- as.numeric(sum(Yr*info$Xr))
  return(r)
}


greedy_optimize_r <- function(Y_start=NULL, Z, X,link="logit",sign_mult = 1,
                                max_iter = 1000, topK = 5,
                                tol = 1e-12, patience = 10,
                                verbose = TRUE) {
  if(is.null(Y_start)) Y=rbinom(length(X),1,.5) else Y=Y_start
  Y <- as.integer(Y)
  iter <- 0
  no_improve <- 0
  improved <- TRUE

  fit_curr <- compute_r_and_info(Y, Z, X,link=link)
  r_curr <- fit_curr$r

  while (iter < max_iter && no_improve < patience) {
    iter <- iter + 1
#    if (verbose)
 #     cat(sprintf("\nIter %d: initial, r = %.8f\n", iter, r_curr))
    deltas <- compute_delta_r_vec(Y, fit_curr$mu, fit_curr$Vd, fit_curr$Xr) * sign_mult
    best_idx <- order(deltas, decreasing = TRUE)[1:min(topK, length(deltas))]

    improved <- FALSE
    for (i in best_idx) {
      if (deltas[i] <= tol) break
      Yt <- Y; Yt[i] <- 1 - Yt[i]
      fit_last <- compute_r_and_info(Yt, Z, X,link=link)
      r_new <- fit_last$r
      if (sign_mult * (r_new - r_curr) > tol) {
        Y <- Yt
        r_curr <- r_new
        fit_curr <- fit_last
        improved <- TRUE
        no_improve <- 0
      #  if (verbose)
      #    cat(sprintf("\nIter %d: flipped %d, r = %.8f;", iter, i, r_new))
        break

      }
    }
    if (!improved) no_improve <- no_improve + 1
  }

  final_r <- fit_compute_r(Y, Z, X,link=link)
  list(Y = Y, r = final_r, iterations = iter)
}

multi_start_r <- function(Z, X, Y_user = NULL,
                          link="logit",
                            thresholds = c(-.1,0,.1),
                            n_random = 10,
                            mode = "max",
                            max_iter = 1000,
                            topK = 10,
                            tol = 1e-12,
                            patience = 10,
                            verbose = TRUE) {
  n <- nrow(Z)
  topK=min(n,topK)
  inits <- list()
  sign_mult <- if (mode == "max") 1 else -1
  if (!is.null(Y_user)) inits[["user"]] <- as.integer(Y_user)
  nrmzXr=.nrmz(.get_IH(Z)%*%X)
  if(mode == "max"){
    for (tau in thresholds) inits[[paste0("thr", tau)]] <- threshold_init_Y(nrmzXr, tau)
  } else if(mode == "min"){
    for (tau in thresholds) inits[[paste0("thr", tau)]] <- 1 - threshold_init_Y(nrmzXr, tau)
  }
  if(n_random>0){
    etahat0=.get_H(Z)%*%X
    etahat=etahat0+sign_mult*X*1
    muhat=plogis(etahat)
    for (r in 1:n_random) inits[[paste0("rand", r)]] <- rbinom(n, 1, muhat)
    }
  best <- list(r = if (mode == "max") -Inf else Inf, Y = NULL)
  for (nm in 1:length(inits)){ #(nm in names(inits)) {
    if (verbose)
      cat("\n evaluating", names(inits)[nm], "\n")
    res <- greedy_optimize_r(inits[[nm]], Z, X,link=link,
                             sign_mult = sign_mult,
                             max_iter = max_iter, topK = topK,
                             tol = tol, patience = patience,
                             verbose = verbose)
    if ((mode == "max" && res$r > best$r) ||
        (mode == "min" && res$r < best$r)) {
      best <- res
      if (verbose)
        cat("New best", mode, "r =", round(best$r, 6), "from", names(inits)[nm], "\n")
    }
  }
  best
}

# ----------- Brute force functions for small n -----------

bruteforce_r <- function(Z, X, link="logit", mode = "max") {
  n <- nrow(Z)
  sign_mult <- if (mode == "max") 1 else -1
  best_val <- if (mode == "max") -Inf else Inf
  best_Y <- NULL
  for (num in 0:(2^n - 1)) {
    Y <- as.integer(intToBits(num)[1:n])
    info <- fit_null_glm(Y, Z, X,link=link)
    val <- compute_r(Y, info$mu, info$V_diag, as.numeric(info$X_r))
    if ((mode == "max" && val > best_val) || (mode == "min" && val < best_val)) {
      best_val <- val
      best_Y <- Y
    }
  }
  list(Y = best_Y, r = best_val)
}

##############################################################
# exact_threshold_search.R
# Exact optimizer for single-predictor case with intercept-only null.
# Inputs:
#   X : numeric vector (length n)
#   mode: "max_r", "min_r", "max_R2"  (default "max_r")
# Output:
#   list(Y = best binary vector,
#        r = r(Y),
#        R2 = R2(Y),   # for p=1, R2 = r^2
#        k = number of ones)

exact_threshold_search_intercept_only <- function(X) {
  # mode <- match.arg(mode)
  n <- length(X)
  if (n < 2) stop("need n>=2")
  # residualize X wrt intercept: center
  Xr <- as.numeric(X - mean(X))
  normXr <- sqrt(sum(Xr^2))
  if (normXr == 0) stop("X is constant after centering: no information")

  # sort Xr descending and remember original indices
  ord <- order(Xr, decreasing = TRUE)
  Xr_sorted <- Xr[ord]
  # cumulative sums of sorted Xr
  csum <- cumsum(Xr_sorted[-n])  # csum[k] = sum of top-k Xr

  # Prepare formula for r_k for each k=1..n-1
  # We derived (after algebra): ||Y_r|| = sqrt(n) and
  # X_r^T Y_r = (1/sqrt(p(1-p))) * csum[k] where p=k/n
  # So r_k = (X_r^T Y_r) / (||X_r|| * ||Y_r||)
  # = csum[k] / ( sqrt(p(1-p)) * normXr * sqrt(n) )
  ks <- 1:(n-1)
  pks <- ks / n
  denom_k <- sqrt(pks * (1 - pks)) * normXr * sqrt(n)  # >0
  r_k <- csum / denom_k

  # For min, look at negative direction: we can either take bottom-k or simply look for min r_k
  # if (mode == "max_r") { always look for the max
    best_k_idx <- which.max(r_k)
    best_r <- r_k[best_k_idx]
    #best_k <- ks[best_k_idx]
    # construct best Y: ones at top-k positions (in original order)
    #Y_sorted_best <- integer(n)
    #Y_sorted_best[1:best_k] <- 1
    #Y_best <- integer(n)
    #Y_best[ord] <- Y_sorted_best
    #best_R2 <- best_r^2
  # } else if (mode == "min_r") {
  #   # min r could be achieved by taking bottom-k large negative Xr values.
  #   # Instead compute r for reversed order (same formula)
  #   csum_rev <- cumsum(rev(Xr_sorted))      # bottom cumulative sums (smallest first)
  #   r_k_rev <- csum_rev / denom_k        # denom_k same for given k
  #   best_k_idx <- which.min(r_k_rev)
  #   best_r <- r_k_rev[best_k_idx]
  #   best_k <- ks[best_k_idx]
  #   Y_sorted_best <- integer(n)
  #   Y_sorted_best[(n-best_k+1):n] <- 1     # pick bottom-k in original sorted order
  #   Y_best <- integer(n)
  #   Y_best[ord] <- Y_sorted_best
  #   best_R2 <- best_r^2
  # } else { # max_R2: maximize |r|
  #   abs_r_k <- abs(r_k)
  #   best_k_idx <- which.max(abs_r_k)
  #   best_r <- r_k[best_k_idx]
  #   best_k <- ks[best_k_idx]
  #   Y_sorted_best <- integer(n)
  #   # if r positive -> pick top-k, else if r negative -> bottom-k
  #   if (r_k[best_k_idx] >= 0) {
  #     Y_sorted_best[1:best_k] <- 1
  #   } else {
  #     Y_sorted_best[(n-best_k+1):n] <- 1
  #   }
  #   Y_best <- integer(n)
  #   Y_best[ord] <- Y_sorted_best
  #   best_R2 <- best_r^2
  # }

  return(as.numeric(best_r))
  # list(Y = Y_best, r = as.numeric(best_r), R2 = as.numeric(best_R2), k = as.integer(best_k))
}

# Example:
 # set.seed(1)
 # X <- rnorm(10)
 # exact_threshold_search_intercept_only(X)





