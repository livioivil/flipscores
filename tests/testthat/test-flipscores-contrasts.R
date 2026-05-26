test_that("formula interface reports interaction notes when appropriate", {
  set.seed(11)
  toy <- data.frame(
    trt = factor(rep(c("A", "B"), each = 15)),
    sex = factor(rep(c("F", "M"), length.out = 30))
  )
  toy$y <- rbinom(30, 1, ifelse(toy$trt == "B", 0.65, 0.35))

  fit <- glm(y ~ trt * sex, data = toy, family = binomial, x = TRUE)

  out <- flipscores_contrasts(fit, pairwise ~ trt, n_flips = 20, seed = 1)
  expect_s3_class(out, "flipscores_contrasts")
  expect_equal(nrow(out$table), 1)
  expect_match(out$notes, "trt:sex", fixed = TRUE)

  out_by <- flipscores_contrasts(fit, pairwise ~ trt | sex,
                                 n_flips = 20, seed = 1)
  expect_s3_class(out_by, "flipscores_contrasts")
  expect_equal(nrow(out_by$table), 2)
  expect_length(out_by$notes, 0)
})

test_that("custom coefficient contrasts accept matrix, vector, and list input", {
  set.seed(11)
  toy <- data.frame(trt = factor(rep(c("A", "B"), each = 15)))
  toy$y <- rbinom(30, 1, ifelse(toy$trt == "B", 0.65, 0.35))

  fit <- glm(y ~ trt, data = toy, family = binomial, x = TRUE)
  K <- matrix(c(0, 1), nrow = 1)
  colnames(K) <- names(coef(fit))
  rownames(K) <- "B - A"

  out_matrix <- flipscores_contrasts(fit, linfct = K, n_flips = 20, seed = 1)
  out_vector <- flipscores_contrasts(fit, linfct = c(trtB = 1),
                                     n_flips = 20, seed = 1)
  out_list <- flipscores_contrasts(fit, linfct = list("B - A" = c(trtB = 1)),
                                   n_flips = 20, seed = 1)

  expect_equal(out_matrix$table$estimate, unname(coef(fit)["trtB"]))
  expect_equal(out_vector$table$estimate, out_matrix$table$estimate)
  expect_equal(out_list$table$estimate, out_matrix$table$estimate)
  expect_equal(out_matrix$table$p.value, out_vector$table$p.value)
  expect_equal(out_matrix$table$p.value, out_list$table$p.value)
})

test_that("Dunnett-like formula contrasts compare all levels to one control", {
  set.seed(14)
  toy <- data.frame(
    trt = factor(rep(c("Control", "Low", "High"), each = 12),
                 levels = c("Control", "Low", "High"))
  )
  toy$y <- rbinom(
    nrow(toy), 1,
    ifelse(toy$trt == "High", 0.7, ifelse(toy$trt == "Low", 0.55, 0.35))
  )

  fit <- glm(y ~ trt, data = toy, family = binomial, x = TRUE)

  out <- flipscores_contrasts(fit, dunnett ~ trt,
                              n_flips = 20, seed = 1)
  expect_equal(out$table$contrast, c("Low - Control", "High - Control"))
  expect_equal(nrow(out$table), 2)
  expect_equal(unname(out$linfct[, "trtLow"]), c(1, 0))
  expect_equal(unname(out$linfct[, "trtHigh"]), c(0, 1))

  out_ref_name <- flipscores_contrasts(fit, trt.vs.ctrl ~ trt,
                                       ref = "Low", n_flips = 20, seed = 1)
  expect_equal(out_ref_name$table$contrast, c("Control - Low", "High - Low"))

  out_ref_last <- flipscores_contrasts(fit, trt.vs.ctrlk ~ trt,
                                       n_flips = 20, seed = 1)
  expect_equal(out_ref_last$table$contrast,
               c("Control - High", "Low - High"))
})

test_that("lm and flipscores objects are accepted", {
  set.seed(12)
  lm_dat <- data.frame(
    y = rnorm(30),
    trt = factor(rep(c("A", "B"), each = 15))
  )
  lm_fit <- lm(y ~ trt, data = lm_dat, x = TRUE, y = TRUE)
  lm_out <- flipscores_contrasts(lm_fit, pairwise ~ trt,
                                 n_flips = 20, seed = 1)

  expect_s3_class(lm_out, "flipscores_contrasts")
  expect_equal(nrow(lm_out$table), 1)

  set.seed(13)
  fs_dat <- data.frame(
    x = rnorm(24),
    trt = factor(rep(c("A", "B"), each = 12))
  )
  fs_dat$y <- rpois(24, exp(0.2 + 0.4 * (fs_dat$trt == "B") + 0.3 * fs_dat$x))
  fs_fit <- flipscores(y ~ trt + x, data = fs_dat, family = "poisson",
                       n_flips = 20, seed = 1, x = TRUE)
  fs_out <- flipscores_contrasts(fs_fit, linfct = c(trtB = 1),
                                 n_flips = 20, seed = 1)

  expect_s3_class(fs_out, "flipscores_contrasts")
  expect_equal(nrow(fs_out$table), 1)
  expect_equal(fs_out$table$estimate, unname(coef(fs_fit)["trtB"]))
})
