\name{flipscores_contrasts}
\alias{flipscores_contrasts}
\title{Flip-score tests for custom contrasts}
\usage{
flipscores_contrasts(
  model,
  specs,
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
\item{model}{a fitted \code{glm} or \code{flipscores} object.}

\item{specs}{contrast formula; see \code{\link{custom_contrasts}}.}

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
An object of class \code{flipscores_contrasts}, containing the contrast
definition, reference grid, flip statistic space, p-values, and printed result
table.
}
\description{
Defines custom factor contrasts and applies the flip-score test to each
contrast. This function connects \code{custom_contrasts()} to the flip-score
machinery.
}
\examples{
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
}
