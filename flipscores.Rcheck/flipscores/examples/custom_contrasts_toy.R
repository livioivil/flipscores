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

fit <- glm(
  y ~ treatment * sex + age_group + severity,
  data = dat,
  family = binomial,
  x = TRUE
)

# Build contrasts only.
custom_contrasts(fit, ~ treatment)
custom_contrasts(fit, pairwise ~ treatment | sex)
custom_contrasts(fit, pairwise ~ treatment:age_group | sex)

# Build contrasts and apply flip-score tests.
flipscores_contrasts(
  fit,
  pairwise ~ treatment | sex,
  n_flips = 500,
  seed = 1
)

flipscores_contrasts(
  fit,
  ~ treatment:age_group | sex,
  n_flips = 500,
  seed = 1
)
