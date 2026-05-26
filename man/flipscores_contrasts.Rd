\name{flipscores_contrasts}
\alias{flipscores_contrasts}
\title{Flip-score tests for custom contrasts}
\usage{
flipscores_contrasts(
  model,
  specs = NULL,
  linfct = NULL,
  ref = NULL,
  score_type = "standardized",
  n_flips = 5000,
  alternative = "two.sided",
  id = NULL,
  seed = NULL,
  flips = NULL,
  precompute_flips = TRUE,
  ...
)
}
\arguments{
\item{model}{a fitted \code{glm}, \code{lm}, or \code{flipscores} object.}

\item{specs}{a contrast formula, coefficient contrast matrix, numeric vector,
named list of numeric vectors, or \code{NULL} when \code{linfct} is supplied.}

\item{linfct}{optional custom coefficient contrast matrix, numeric vector,
or named list of numeric vectors. Rows define linear functions of the model
coefficients, similar to \code{multcomp::glht()}.}

\item{ref}{control level for Dunnett-like formula contrasts. It can be a
level name, cell name, or numeric index. The default is the first level/cell,
except for \code{trt.vs.ctrlk}, where the default is the last level/cell.}

\item{score_type}{type of score, as in \code{\link{flipscores}}.}

\item{n_flips}{number of sign flips.}

\item{alternative}{\code{"two.sided"}, \code{"greater"}, or \code{"less"}.}

\item{id}{optional cluster id.}

\item{seed}{optional random seed.}

\item{flips}{optional precomputed flip matrix.}

\item{precompute_flips}{whether to precompute the flip matrix.}

\item{...}{currently unused.}
}
\value{
An object of class \code{flipscores_contrasts}, containing the resolved
coefficient contrasts, optional reference grid, flip statistic space, p-values,
and printed result table.
}
\description{
Defines factor contrasts from a formula or accepts fully custom coefficient
contrasts, then applies the flip-score test to each contrast.
}
\details{
For formula-based contrasts, an emmeans-style note is printed when the requested
factor is involved in model interactions that are not represented in the
contrast specification. For example, \code{pairwise ~ trt} prints a note when
the fitted model contains \code{trt:sex}; \code{pairwise ~ trt | sex} does not,
because the interaction partner is included as a conditioning factor.

Dunnett-like all-versus-one contrasts can be requested with
\code{dunnett ~ trt}, \code{trt.vs.ctrl ~ trt}, or
\code{trt.vs.ctrlk ~ trt}. Use \code{ref} to choose the control level.
}
\examples{
# Minimal example: one binary factor
set.seed(11)
toy <- data.frame(
  trt = factor(rep(c("A", "B"), each = 15))
)
toy$y <- rbinom(30, 1, ifelse(toy$trt == "B", 0.65, 0.35))

toy_fit <- glm(y ~ trt, data = toy, family = binomial, x = TRUE)

# Same style as emmeans: pairwise comparisons for trt
flipscores_contrasts(toy_fit, pairwise ~ trt,
                     n_flips = 100, seed = 1)

# Dunnett-like all-versus-one comparisons
toy3 <- data.frame(
  trt = factor(rep(c("Control", "Low", "High"), each = 12),
               levels = c("Control", "Low", "High"))
)
toy3$y <- rbinom(nrow(toy3), 1,
                 ifelse(toy3$trt == "High", 0.7,
                        ifelse(toy3$trt == "Low", 0.55, 0.35)))
toy3_fit <- glm(y ~ trt, data = toy3, family = binomial, x = TRUE)

flipscores_contrasts(toy3_fit, dunnett ~ trt,
                     n_flips = 100, seed = 1)
flipscores_contrasts(toy3_fit, trt.vs.ctrl ~ trt, ref = "Low",
                     n_flips = 100, seed = 1)
flipscores_contrasts(toy3_fit, trt.vs.ctrlk ~ trt,
                     n_flips = 100, seed = 1)

# If trt is involved in an interaction, requesting only trt prints a note.
toy$sex <- factor(rep(c("F", "M"), length.out = 30))
toy_int <- glm(y ~ trt * sex, data = toy, family = binomial, x = TRUE)
flipscores_contrasts(toy_int, pairwise ~ trt,
                     n_flips = 100, seed = 1)

# Custom coefficient contrast, similar in spirit to multcomp::glht().
# Here the row tests the trtB coefficient directly.
K <- matrix(c(0, 1), nrow = 1)
colnames(K) <- names(coef(toy_fit))
rownames(K) <- "B - A"
flipscores_contrasts(toy_fit, linfct = K,
                     n_flips = 100, seed = 1)

# The same custom contrast can be supplied as a named vector.
flipscores_contrasts(toy_fit, linfct = c(trtB = 1),
                     n_flips = 100, seed = 1)

set.seed(1)
dat <- data.frame(
  treatment = factor(rep(c("Control", "Drug"), each = 30)),
  sex = factor(rep(c("F", "M"), times = 30)),
  age_group = factor(rep(c("Young", "Old", "Middle"), length.out = 60)),
  severity = runif(60, 0, 1)
)
eta <- with(dat,
  -1 + 0.8 * (treatment == "Drug") + 0.4 * (sex == "M") +
    0.6 * (age_group == "Old") - 0.3 * severity
)
dat$y <- rbinom(nrow(dat), 1, plogis(eta))

fit <- glm(y ~ treatment * sex + age_group + severity,
           data = dat, family = binomial, x = TRUE)

# Test the treatment effect separately for F and M.
flipscores_contrasts(fit, pairwise ~ treatment | sex,
                     n_flips = 100, seed = 1)

# Test selected interaction-cell contrasts. Increase n_flips in real analyses.
flipscores_contrasts(fit, ~ treatment:age_group | sex,
                     n_flips = 100, seed = 1)

# Linear models are accepted too.
lm_dat <- data.frame(
  y = rnorm(30),
  trt = factor(rep(c("A", "B"), each = 15))
)
lm_fit <- lm(y ~ trt, data = lm_dat, x = TRUE, y = TRUE)
flipscores_contrasts(lm_fit, pairwise ~ trt,
                     n_flips = 100, seed = 1)
}
