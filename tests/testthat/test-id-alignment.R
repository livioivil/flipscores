id_models <- function() {
  set.seed(42)
  d1 <- data.frame(ids = c("A", "A", "B", "C", "C", "D", "D", "E"),
                   x = rnorm(8), z = rnorm(8), y = rnorm(8),
                   row.names = paste0("left", seq_len(8)))
  d2 <- data.frame(ids = c("E", "B", "B", "C", "D", "E", "F"),
                   x = rnorm(7), z = rnorm(7), y = rnorm(7),
                   row.names = paste0("right", seq_len(7)))
  list(left = glm(y ~ x + z, data = d1), right = glm(y ~ x + z, data = d2))
}

test_that("model-local IDs align duplicates and retain local scores", {
  models <- id_models()
  out <- flipscores(models, id = ~ids, to_be_tested = "x", n_flips = 25, seed = 3)
  expect_identical(names(out$objects), c("left", "right"))
  expect_identical(vapply(out$objects, function(x) nrow(x$scores), integer(1)),
                   c(left = 8L, right = 7L))
  expect_identical(out$alignment$rows$id,
                   c("A", "A", "B", "B", "C", "C", "D", "D", "E", "E", "F"))
  expect_identical(out$alignment$rows$occurrence, c(1L, 2L, 1L, 2L, 1L, 2L,
                                                 1L, 2L, 1L, 2L, 1L))
  aligned <- .fs_aligned_scores(out)
  expect_equal(dim(aligned), c(11L, 2L))
  expect_identical(unname(which(is.na(aligned[, 1]))), c(4L, 10L, 11L))
  expect_identical(unname(which(is.na(aligned[, 2]))), c(1L, 2L, 6L, 8L))
  expect_equal(unname(aligned[out$alignment$index$left, 1]),
               unname(out$objects$left$scores$x))
  expect_equal(unname(aligned[out$alignment$index$right, 2]),
               unname(out$objects$right$scores$x))
  expect_false(anyNA(out$Tspace))
  string <- flipscores(models, id = "ids", to_be_tested = "x", n_flips = 25, seed = 3)
  expect_equal(string$Tspace, out$Tspace)
  expect_identical(string$alignment, out$alignment)
})

test_that("joint cluster flips produce the same tests as explicit local flips", {
  models <- id_models()
  G <- rbind(c(1, 1, 1, 1, 1, 1), c(-1, 1, -1, 1, 1, -1),
             c(1, -1, -1, 1, -1, 1))
  colnames(G) <- LETTERS[1:6]
  for (type in c("basic", "effective", "standardized")) {
    joint <- flipscores(models, id = "ids", to_be_tested = "x", flips = G,
                         score_type = type)
    for (i in seq_along(models)) {
      ids <- .fs_model_data(models[[i]])$ids
      local <- G[, ids, drop = FALSE]
      colnames(local) <- rownames(model.matrix(models[[i]]))
      expected <- flipscores(models[[i]], id = ids, flips = local,
                              to_be_tested = "x", score_type = type)
      expect_equal(joint$objects[[i]]$Tspace, expected$Tspace)
      expect_equal(joint$objects[[i]]$p.values, expected$p.values)
      expect_identical(attr(joint$objects[[i]]$scores, "id"), ids)
      if (type != "standardized") {
        null <- glm(y ~ z, data = .fs_model_data(models[[i]]))
        residual <- null$y - null$fitted.values
        x <- .fs_model_data(models[[i]])$x
        if (type == "effective") x <- residuals(lm(x ~ null$model$z))
        expected_stat <- as.vector(local %*% (x * residual)) /
          sqrt(summary(null)$dispersion * sum(x^2))
        expect_equal(joint$objects[[i]]$Tspace$x, expected_stat, tolerance = 1e-8)
      }
    }
  }
})

test_that("ID expressions respect subset, omitted rows and model environments", {
  make_model <- function() {
    d <- data.frame(ids = rep(LETTERS[1:5], each = 2),
                    x = seq_len(10), y = c(1, 3, 2, NA, 6, 3, 8, 7, 6, 9),
                    row.names = paste0("r", 1:10))
    glm(y ~ x, data = d, subset = c(9, 2, 4, 7, 5, 10), na.action = na.exclude)
  }
  model <- make_model()
  d <- data.frame(ids = "wrong data")
  out <- flipscores(model, id = ~ids, to_be_tested = "x", n_flips = 10, seed = 1)
  expect_identical(rownames(out$scores), c("r9", "r2", "r7", "r5", "r10"))
  expect_identical(attr(out$scores, "id"), c("E", "A", "D", "C", "E"))
  joint <- flipscores(list(model, model), id = "ids", to_be_tested = "x",
                      n_flips = 10, seed = 1, data = d)
  expect_equal(joint$objects[[1]]$Tspace, out$Tspace)
  expect_identical(attr(joint$objects[[1]]$scores, "id"), attr(out$scores, "id"))
})

test_that("formula interfaces resolve IDs in data and return fitted rows only", {
  set.seed(4)
  d <- data.frame(ids = rep(LETTERS[1:6], each = 2), x = rnorm(12),
                  y = rnorm(12), y2 = rnorm(12))
  d$y[2] <- NA
  d$y2[5] <- NA
  one <- flipscores(y ~ x, data = d, id = "ids", to_be_tested = "x",
                     n_flips = 15, seed = 9)
  by_formula <- flipscores(y ~ x, data = d, id = ~ids, to_be_tested = "x",
                            n_flips = 15, seed = 9)
  expect_equal(one$Tspace, by_formula$Tspace)
  expect_identical(attr(one$scores, "id"), d$ids[-2])
  formulas <- flipscores(list(y ~ x, y2 ~ x), data = d, id = "ids",
                          to_be_tested = "x", n_flips = 15, seed = 9)
  matrix <- flipscores(cbind(y, y2) ~ x, data = d, id = ~ids,
                        to_be_tested = "x", n_flips = 15, seed = 9)
  expect_equal(matrix$Tspace, formulas$Tspace)
  expect_equal(vapply(formulas$objects, function(x) nrow(x$scores), integer(1)),
               c(mod1 = 11L, mod2 = 11L))
  expect_identical(attr(formulas$objects[[2]]$scores, "id"), d$ids[-5])
})

test_that("on-demand flips preserve joint alignment and coefficient synchronization", {
  models <- id_models()
  for (id in list(NULL, "ids")) {
    precomputed <- flipscores(models, id = id, n_flips = 31, seed = 12)
    streamed <- flipscores(models, id = id, n_flips = 31, seed = 12,
                            precompute_flips = FALSE)
    expect_equal(streamed$Tspace, precomputed$Tspace)
    expect_identical(streamed$alignment, precomputed$alignment)
  }
  model <- models[[1]]
  a <- flipscores(model, id = ~ids, n_flips = 31, seed = 12)
  b <- flipscores(model, id = ~ids, n_flips = 31, seed = 12, precompute_flips = FALSE)
  expect_equal(a$Tspace, b$Tspace)
})

test_that("NULL IDs retain row-name matching independently of row ordering", {
  models <- id_models()
  # With no ID, these disjoint row names must remain independent despite overlap in ids.
  out <- flipscores(models, to_be_tested = "x", n_flips = 12, seed = 2)
  expect_identical(out$alignment$by, "rownames")
  expect_equal(nrow(out$alignment$rows), 15L)
  expect_null(attr(out$objects[[1]]$scores, "id"))
  d <- .fs_model_data(models[[1]])
  shuffled <- glm(y ~ x + z, data = d[c(8, 3, 1, 6, 2, 7, 5, 4), ])
  same <- flipscores(list(models[[1]], shuffled), to_be_tested = "x",
                     n_flips = 12, seed = 2)
  expect_equal(nrow(same$alignment$rows), 8L)
  expect_equal(same$objects[[1]]$Tspace, same$objects[[2]]$Tspace)
})

test_that("invalid IDs and inconsistent supplied cluster flips fail clearly", {
  models <- id_models()
  expect_error(flipscores(models[[1]], id = 1:7, n_flips = 5), "one value per score row")
  expect_error(flipscores(models, id = 1:8, n_flips = 5), "one value per score row")
  expect_error(flipscores(models, id = "missing", n_flips = 5), "column missing")
  ids <- rep("external", 8)
  expect_error(flipscores(models, id = ~missing, n_flips = 5), "column.*missing")
  bad <- .fs_model_data(models[[2]])
  bad$ids <- NULL
  expect_error(flipscores(list(models[[1]], glm(y ~ x + z, data = bad)),
                           id = ~ids, n_flips = 5), "column.*missing")
  expect_error(flipscores(models[[1]], id = c(NA, 1:7), n_flips = 5), "missing values")
  expect_error(make_flips(8, 5, id = 1:7), "one value per score row")
  wrong <- make_flips(8, 5, obs_names = rownames(model.matrix(models[[1]])))
  wrong[2, 1:2] <- c(-1, 1)
  expect_error(flipscores(models[[1]], id = ~ids, flips = wrong), "same sign")
  expect_error(flipscores(models, id = ~ids, score_type = "orthogonalized"),
                "not supported")
})

test_that("compute_scores stores resolved IDs without padding score metadata", {
  models <- id_models()
  d <- .fs_model_data(models[[1]])
  null <- glm(y ~ z, data = d, x = TRUE)
  scores <- compute_scores(null, matrix(d$x), score_type = "standardized", id = ~ids,
                            nobservations = 100)
  expect_identical(attr(scores, "id"), d$ids)
  expect_equal(nrow(scores), 8L)
  expect_equal(nrow(attr(scores, "scale_objects")[[1]]$A), 8L)
  expect_error(compute_scores(null, matrix(d$x), score_type = "basic", id = 1:7),
                "one value per score row")
})

test_that("explicit vectors and single-cluster flips work in both generation modes", {
  models <- id_models()
  ids <- c(1, 1, 2, 2, 3, 3, 4, 4)
  vector <- flipscores(list(models[[1]], models[[1]]), id = ids,
                        to_be_tested = "x", n_flips = 11, seed = 2)
  expect_identical(attr(vector$objects[[1]]$scores, "id"), ids)
  expect_equal(vector$objects[[1]]$Tspace, vector$objects[[2]]$Tspace)
  for (precompute in c(TRUE, FALSE)) {
    single <- flipscores(models[[1]], id = rep("one", 8), to_be_tested = "x",
                          n_flips = 1, seed = 2, precompute_flips = precompute)
    expect_equal(nrow(single$Tspace), 1)
    one_cluster <- flipscores(models[[1]], id = rep("one", 8), to_be_tested = "x",
                               n_flips = 11, seed = 2, precompute_flips = precompute)
    expect_equal(abs(one_cluster$Tspace$x), rep(abs(single$Tspace$x), 11))
  }
})

test_that("updates preserve ID expressions and joined calls outlive their temporary frame", {
  models <- id_models()
  single <- flipscores(models[[1]], id = ~ids, to_be_tested = "x", n_flips = 9, seed = 3)
  expect_equal(update(single)$Tspace, single$Tspace)
  joint <- flipscores(models, id = "ids", to_be_tested = "x", n_flips = 9, seed = 3)
  updated <- update(joint$objects[[1]])
  expect_identical(attr(updated$scores, "id"), .fs_model_data(models[[1]])$ids)
  expect_equal(updated$scores, joint$objects[[1]]$scores)
  G <- matrix(1, 3, 6, dimnames = list(NULL, LETTERS[1:6]))
  G[2, ] <- c(-1, 1, -1, 1, -1, 1)
  supplied <- flipscores(models, id = ~ids, to_be_tested = "x", flips = G)
  expect_equal(update(supplied$objects[[2]])$Tspace, supplied$objects[[2]]$Tspace)
})

test_that("ID clustering also works when an fs_lm result is joined", {
  models <- id_models()
  d <- .fs_model_data(models[[1]])
  linear <- flipscores_lm(y ~ x + z, data = d, tested_coeffs = "x", n_flips = 7, seed = 3)
  joint <- flipscores(list(linear = linear, glm = models[[1]]), id = ~ids,
                       to_be_tested = "x", n_flips = 7, seed = 3)
  expect_identical(attr(joint$objects$linear$scores, "id"), d$ids)
  expect_equal(dim(.fs_aligned_scores(joint)), c(8L, 2L))
  expect_identical(joint$objects$linear$info$D, linear$info$D)
})

test_that("cluster IDs are preserved when joining contrast and nested results", {
  models <- id_models()
  contrast <- flipscores_contrasts(models[[1]], linfct = c(x = 1), id = ~ids,
                                    n_flips = 15, seed = 5)
  expect_identical(attr(contrast$scores, "id"), .fs_model_data(models[[1]])$ids)
  joint <- flipscores(list(left = models[[1]], right = models[[2]]),
                       id = "ids", to_be_tested = "x", n_flips = 15, seed = 5)
  rejoined <- flipscores(list(joint = joint, other = models[[1]]), id = ~ids,
                          to_be_tested = "x", n_flips = 15, seed = 5)
  expect_equal(rejoined$objects$joint$Tspace, joint$Tspace)
  expect_equal(dim(.fs_aligned_scores(rejoined)), c(11L, 3L))
  # Identical coefficient and contrast tests must use identical cluster signs.
  mixed <- flipscores(list(coef = models[[1]], contrast = contrast), id = ~ids,
                       to_be_tested = "x", n_flips = 15, seed = 5)
  expect_equal(unname(as.matrix(mixed$objects$coef$Tspace)),
               unname(mixed$objects$contrast$Tspace))
})
