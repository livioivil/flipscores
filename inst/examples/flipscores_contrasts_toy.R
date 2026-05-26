set.seed(123)

dat <- data.frame(
  treatment = factor(rep(c("Control", "Drug"), each = 30)),
  sex = factor(rep(c("F", "M"), times = 30)),
  age_group = factor(rep(c("Young", "Old", "Middle"), length.out = 60)),
  severity = runif(60, 0, 1)
)

eta <- with(dat,
  -1 + 1.2 * (treatment == "Drug") * (sex == "M") +
    0.6 * (age_group == "Old") - 0.3 * severity
)
dat$y <- rbinom(nrow(dat), 1, plogis(eta))

fit <- glm(
  y ~ treatment * sex + age_group + severity,
  data = dat,
  family = binomial,
  x = TRUE
)

# Formula interface: build contrasts and apply flip-score tests.
flipscores_contrasts(
  fit,
  pairwise ~ treatment,
  n_flips = 500,
  seed = 1
)

# Formula interface: build contrasts and apply flip-score tests.
flipscores_contrasts(
  fit,
  pairwise ~ treatment | sex,
  n_flips = 500,
  seed = 1
)

# Dunnett-like all-versus-one comparisons.
flipscores_contrasts(
  fit,
  dunnett ~ age_group,
  n_flips = 500,
  seed = 1
)

flipscores_contrasts(
  fit,
  trt.vs.ctrl ~ age_group,
  ref = "Middle",
  n_flips = 500,
  seed = 1
)

flipscores_contrasts(
  fit,
  trt.vs.ctrlk ~ age_group,
  n_flips = 500,
  seed = 1
)

flipscores_contrasts(
  fit,
  ~ treatment:age_group | sex,
  n_flips = 500,
  seed = 1
)

# Custom coefficient contrast matrix, similar in spirit to multcomp::glht().
K <- matrix(0, nrow = 2, ncol = length(coef(fit)))
colnames(K) <- names(coef(fit))
rownames(K) <- c("Drug effect in F", "Drug effect difference M vs F")
K["Drug effect in F", "treatmentDrug"] <- 1
K["Drug effect difference M vs F", "treatmentDrug:sexM"] <- 1

flipscores_contrasts(
  fit,
  linfct = K,
  n_flips = 500,
  seed = 1
)
