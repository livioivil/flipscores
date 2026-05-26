\name{custom_contrasts}
\alias{custom_contrasts}
\alias{apply_custom_contrasts}
\title{Build custom contrasts for flip-score tests}
\usage{
custom_contrasts(model, specs)

apply_custom_contrasts(contrasts, model = NULL)
}
\arguments{
\item{model}{a fitted \code{glm} or \code{flipscores} object.}

\item{specs}{a one-sided or two-sided formula. Put \code{pairwise} on the
left-hand side to request all pairwise comparisons; otherwise simple treatment
contrasts against the first level combination are generated.}

\item{contrasts}{a \code{custom_contrasts} object.}
}
\value{
\code{custom_contrasts()} returns an object of class
\code{custom_contrasts}. \code{apply_custom_contrasts()} returns the
corresponding contrast matrix.
}
\description{
Creates emmeans-like contrasts for factors in a fitted model. Supported syntax
includes \code{~ factor}, \code{pairwise ~ factor},
\code{~ factor1:factor2}, and
\code{pairwise ~ factor1:factor2 | by_factor}.
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
           data = dat, family = binomial)

# Simple treatment contrast for one factor
custom_contrasts(fit, ~ treatment)

# Pairwise treatment contrasts inside each sex
custom_contrasts(fit, pairwise ~ treatment | sex)

# Interaction-cell contrasts, similar to emmeans pairwise comparisons
cc <- custom_contrasts(fit, pairwise ~ treatment:age_group | sex)
apply_custom_contrasts(cc)
}
