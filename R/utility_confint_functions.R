# initialization and bisection algorithms for confint

conf_bound_equitailed <- function(to_be_tested,
                                  mod,
                                  flips = NULL,
                                  alpha = 0.05,
                                  score_type = score_type,
                                  alternative = alternative)
{
  if (alternative == "two.sided")
  {
    work_alpha <- alpha / 2
  } else
  {
    work_alpha <- alpha
  }
  dat <- eval(mod$data)
  if(is.environment(dat) | is.null(dat)) # in case the formula parts come from the environment and not from a "proper" dataset
  {
    frml <- all.vars(mod$terms)
    dat <- sapply(frml, \(.form___) eval(parse(text = .form___)))
    dat <- as.data.frame(dat)
  }
  #___
  # factor/char management
  # build dummy variables for all factors
  dat <- merge(dat, model.matrix(mod, dat))
  #___
  if(substr(mod$family$family, 0, 17) == "Negative Binomial")
  {
    fam <- "negbinom"
  }else
  {
    fam <- as.call(list(
      as.name(mod$family$family),
      link = mod$family$link
    ))
  }

  dat_test_col <- match(names(mod$coefficients[to_be_tested]), colnames(dat))
  bhat <- bb <- unname(mod$coefficients[to_be_tested])

  amp <- max(abs(bhat / 100), 0.2, coef(summary.glm(mod))[to_be_tested, "Std. Error"])

  if (alternative != "greater")
  {
    pp <- 1.1
    alt <- "greater"
    while (pp >= work_alpha)
    {
      bb <- bb - amp
      if (abs(bb - bhat) > (15 * amp))
      {
        warning("A confidence bound was not found in the search interval -- returning Inf")
        bb <- -Inf
        break
      }
      dat$.OFFSET___ <- dat[, dat_test_col] * bb
      args <- list(
        mod,
        data = dat,
        family = fam,
        # family = binomial(link = 'logit'),
        flips = flips,
        score_type = score_type,
        alternative = alt,
        to_be_tested = to_be_tested,
        formula = . ~ . + offset(.OFFSET___)
      )
      pp <- do.call(update, args)$p.values[names(mod$coefficients[to_be_tested])]
    }
    if (bb > -Inf)
    {
      low <- conf_bound_oneside(
        mod,
        to_be_tested,
        bb,
        bb,
        pp,
        score_type,
        flips = flips,
        lower = 1,
        # -1 if upper
        move_direction = 1,
        # -1 if going left
        eps = amp,
        toll = amp * 1e-5,
        alpha = work_alpha
      )
    } else{
      low <- -Inf
      names(low) <- names(mod$coefficients[to_be_tested])
    }
  }else
  {
    low <- bhat
  }

  if (alternative != "less")
  {
    bb <- bhat
    pp <- 1.1
    alt <- "less"
    while (pp >= work_alpha)
    {
      bb <- bb + amp
      if (abs(bb - bhat) > (15 * amp))
      {
        warning("A confidence bound was not found in the search interval -- returning Inf")
        bb <- Inf
        break
      }
      dat$.OFFSET___ <- dat[, dat_test_col] * bb
      args <- list(
        mod,
        data = dat,
        family = fam,
        flips = flips,
        alternative = alt,
        score_type = score_type,
        to_be_tested = to_be_tested,
        formula = . ~ . + offset(.OFFSET___)
      )
      pp <- do.call(update, args)$p.values[names(mod$coefficients[to_be_tested])]
    }
    if (bb < Inf)
    {
      upp <- conf_bound_oneside(
        mod,
        to_be_tested,
        bb,
        bb,
        pp,
        score_type,
        flips = flips,
        lower = -1,
        # -1 if upper
        move_direction = -1,
        # -1 if going left
        eps = amp,
        toll = amp * 1e-5,
        alpha = work_alpha
      )
    } else{
      upp <- Inf
      names(upp) <- names(mod$coefficients[to_be_tested])
    }
  }else
  {
    upp <- bhat
  }
  ci <- rbind(low[1], upp[1])
  attr(ci, "tol") <- if(!is.na(low[2])) low[2] else upp[2]
  # colnames(ci) <- names(dat)[dat_test_col]
  return(ci)
}


# equitailed confidence interval
conf_bound_oneside <- function(mod,
                               to_be_tested,
                               bb,
                               good_bb,
                               good_p,
                               score_type = NULL,
                               flips = NULL,
                               # parms_DV = NULL,
                               lower = 1,
                               # -1 if upper
                               move_direction = 1,
                               # -1 if going left
                               eps = 1,
                               toll = 1e-5,
                               alpha = 0.05)
{
  eps <- eps / 2
  bb <- bb + (move_direction) * eps
  alt <- ifelse(lower > 0, "greater", "less") # change the direction of the test
  fam <- as.call(list(
    as.name(mod$family$family),
    link = mod$family$link
  ))
  # return(list(mod, dat$X * bb))
  dat <- eval(mod$data)
  if(is.environment(dat) | is.null(dat)) # in case the formula parts come from the environment and not from a "proper" dataset
  {
    frml <- all.vars(mod$terms)
    dat <- sapply(frml, \(.form___) eval(parse(text = .form___)))
    dat <- as.data.frame(dat)
  }
  #___
  # factor/char management
  # build dummy variables for all factors
  dat <- merge(dat, model.matrix(mod, dat))
  #___
  if(substr(mod$family$family, 0, 17) == "Negative Binomial")
  {
    fam <- "negbinom"
  }else
  {
    fam <- as.call(list(
      as.name(mod$family$family),
      link = mod$family$link
    ))
  }
  dat_test_col <- match(names(mod$coefficients[to_be_tested]), colnames(dat))
  dat$.OFFSET___ <- dat[, dat_test_col] * bb
  ff <- . ~ . + offset(.OFFSET___)

  args <- list(
    mod,
    formula = ff,
    data = dat,
    alternative = alt,
    flips = flips,
    family = fam,
    to_be_tested = to_be_tested,
    score_type = score_type
  )

  # if (!is.null(parms_DV)) {
  #     args$parms_DV <- parms_DV
  # }
  p_value <- do.call(update, args)$p.values[names(mod$coefficients[to_be_tested])]
  if (p_value < alpha)
  {
    # conservative bound found
    good_bb <- bb
    good_p <- p_value
  }
  if (eps < toll)
  {
    # return the conservative bound and the tolerance
    # attributes(good_bb, "tol") <- eps
    return(c(good_bb, eps))
  } # else
  # if p value smaller than alfa, (lower,lower): move towards beta
  # otherwise, (lower,-lower): move away from beta
  return(
    conf_bound_oneside(
      mod,
      to_be_tested,
      bb,
      good_bb,
      good_p,
      score_type,
      flips,
      lower,
      lower * ifelse(p_value < alpha, 1, -1),
      eps,
      toll,
      alpha
    )
  )
}

# symmetric confidence interval
conf_bound_symmetric <- function(to_be_tested,
                                 mod,
                                 flips = NULL,
                                 alpha = 0.05,
                                 score_type = score_type,
                                 ...)
{
  dat <- eval(mod$data)
  if(is.environment(dat) | is.null(dat)) # in case the formula parts come from the environment and not from a "proper" dataset
  {
    frml <- all.vars(mod$formula)
    dat <- sapply(frml, \(.form___) eval(parse(text = .form___)))
    dat <- as.data.frame(dat)
  }
  #___
  # factor/char management
  # build dummy variables for all factors
  dat <- merge(dat, model.matrix(mod, dat))
  #___
  if(substr(mod$family$family, 0, 17) == "Negative Binomial")
  {
    fam <- "negbinom"
  }else
  {
    fam <- as.call(list(
      as.name(mod$family$family),
      link = mod$family$link
    ))
  }

  dat_test_col <- match(names(mod$coefficients[to_be_tested]), colnames(dat))
  bhat <- bb <- unname(mod$coefficients[to_be_tested])

  amp <- max(abs(bhat / 100), 0.2, coef(summary.glm(mod))[to_be_tested, "Std. Error"])

  bdiff <- 0
  pp <- 1.1
  while (pp >= alpha)
  {
    bdiff <- bdiff + amp
    if (bdiff > (15 * amp))
    {
      warning("A confidence bound was not found in the search interval -- returning Inf")
      bdiff <- Inf
      break
    }
    bl <- bhat - bdiff
    dat$.OFFSET___ <- dat[, dat_test_col] * bl
    args <- list(
      mod,
      data = dat,
      family = fam,
      flips = flips,
      score_type = score_type,
      alternative = "greater",
      to_be_tested = to_be_tested,
      formula = . ~ . + offset(.OFFSET___)
    )
    pl <- do.call(update, args)$p.values[names(mod$coefficients[to_be_tested])]
    bu <- bhat + bdiff
    dat$.OFFSET___ <- dat[, dat_test_col] * bu
    args$alternative <- "less"
    args$data <- dat
    pu <- do.call(update, args)$p.values[names(mod$coefficients[to_be_tested])]
    pp <- pl + pu
  }
  if(bdiff < Inf)
  {
    bdiff <- conf_bound_twoside(
      mod,
      to_be_tested,
      bhat,
      bdiff,
      bdiff,
      pp,
      score_type,
      flips = flips,
      move_direction = -1,
      eps = amp,
      toll = amp * 1e-5,
      alpha = alpha
    )
  }
  low <- bhat - bdiff[1]
  upp <- bhat + bdiff[1]
  names(low) <- names(upp) <- names(mod$coefficients[to_be_tested])
  ci <- rbind(low, upp)
  attr(ci, "tol") <- if(!is.na(bdiff[2])) bdiff[2] else NA
  # colnames(ci) <- names(dat)[dat_test_col]
  return(ci)
}


conf_bound_twoside <- function(mod,
                               to_be_tested,
                               bhat,
                               bdiff,
                               good_diff,
                               good_pp,
                               score_type = NULL,
                               flips = NULL,
                               move_direction = 1,
                               # 1: increase distance, -1:decrease distance
                               eps = 1,
                               toll = 1e-5,
                               alpha = 0.05)
{
  eps <- eps / 2
  bdiff <- bdiff + (move_direction) * eps
  fam <- as.call(list(
    as.name(mod$family$family),
    link = mod$family$link
  ))
  # return(list(mod, dat$X * bb))
  dat <- eval(mod$data)
  if(is.environment(dat) | is.null(dat)) # in case the formula parts come from the environment and not from a "proper" dataset
  {
    frml <- all.vars(mod$formula)
    dat <- sapply(frml, \(.form___) eval(parse(text = .form___)))
    dat <- as.data.frame(dat)
  }
  #___
  # factor/char management
  # build dummy variables for all factors
  dat <- merge(dat, model.matrix(mod, dat))
  #___
  if(substr(mod$family$family, 0, 17) == "Negative Binomial")
  {
    fam <- "negbinom"
  }else
  {
    fam <- as.call(list(
      as.name(mod$family$family),
      link = mod$family$link
    ))
  }
  dat_test_col <- match(names(mod$coefficients[to_be_tested]), colnames(dat))
  # lower
  bl <- bhat - bdiff
  dat$.OFFSET___ <- dat[, dat_test_col] * bl
  ff <- . ~ . + offset(.OFFSET___)
  args <- list(
    mod,
    formula = ff,
    data = dat,
    alternative = "greater",
    flips = flips,
    family = fam,
    to_be_tested = to_be_tested,
    score_type = score_type
  )
  pl <- do.call(update, args)$p.values[names(mod$coefficients[to_be_tested])]

  # upper
  bu <- bhat + bdiff
  dat$.OFFSET___ <- dat[, dat_test_col] * bu
  args$alternative = "less"
  args$data <- dat
  pu <- do.call(update, args)$p.values[names(mod$coefficients[to_be_tested])]

  p_value <- pl + pu
  if (p_value < alpha)
  {
    # conservative interval found
    good_diff <- bdiff
    good_pp <- p_value
  }
  if (eps < toll)
  {
    return(c(good_diff, eps))
  } # else
  # if pvalue < alpha, reduce the distance from beta --> move_dir = -1;
  # otherwise, +1
  return(
    conf_bound_twoside(
      mod,
      to_be_tested,
      bhat,
      bdiff,
      good_diff,
      good_pp,
      score_type,
      flips,
      move_direction = ifelse(p_value < alpha, -1, 1),
      # 1: increase distance, -1:decrease distance
      eps,
      toll,
      alpha
    )
  )
}
