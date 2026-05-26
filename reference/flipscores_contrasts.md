# Flip-score tests for custom contrasts

Defines factor contrasts from a formula or accepts fully custom
coefficient contrasts, then applies the flip-score test to each
contrast.

## Usage

``` r
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
```

## Arguments

- model:

  a fitted `glm`, `lm`, or `flipscores` object.

- specs:

  a contrast formula, coefficient contrast matrix, numeric vector, named
  list of numeric vectors, or `NULL` when `linfct` is supplied.

- linfct:

  optional custom coefficient contrast matrix, numeric vector, or named
  list of numeric vectors. Rows define linear functions of the model
  coefficients, similar to `multcomp::glht()`.

- ref:

  control level for Dunnett-like formula contrasts. It can be a level
  name, cell name, or numeric index. The default is the first
  level/cell, except for `trt.vs.ctrlk`, where the default is the last
  level/cell.

- score_type:

  type of score, as in
  [`flipscores`](https://livioivil.github.io/flipscores/reference/flipscores.md).

- n_flips:

  number of sign flips.

- alternative:

  `"two.sided"`, `"greater"`, or `"less"`.

- id:

  optional cluster id.

- seed:

  optional random seed.

- flips:

  optional precomputed flip matrix.

- precompute_flips:

  whether to precompute the flip matrix.

- ...:

  currently unused.

## Value

An object of class `flipscores_contrasts`, containing the resolved
coefficient contrasts, optional reference grid, flip statistic space,
p-values, and printed result table.

## Details

For formula-based contrasts, an emmeans-style note is printed when the
requested factor is involved in model interactions that are not
represented in the contrast specification. For example, `pairwise ~ trt`
prints a note when the fitted model contains `trt:sex`;
`pairwise ~ trt | sex` does not, because the interaction partner is
included as a conditioning factor.

Dunnett-like all-versus-one contrasts can be requested with
`dunnett ~ trt`, `trt.vs.ctrl ~ trt`, or `trt.vs.ctrlk ~ trt`. Use `ref`
to choose the control level.

## Examples

``` r
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
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>  contrast  estimate     Score p.value
#>     A - B -3.332205 -3.409809    0.01

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
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>        contrast     estimate        Score p.value
#>   Low - Control 2.378116e-16 2.266233e-17    0.99
#>  High - Control 6.931472e-01 8.280787e-01    0.44
flipscores_contrasts(toy3_fit, trt.vs.ctrl ~ trt, ref = "Low",
                     n_flips = 100, seed = 1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>       contrast      estimate         Score p.value
#>  Control - Low -2.378116e-16 -2.266233e-17    0.99
#>     High - Low  6.931472e-01  8.280787e-01    0.35
flipscores_contrasts(toy3_fit, trt.vs.ctrlk ~ trt,
                     n_flips = 100, seed = 1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>        contrast   estimate      Score p.value
#>  Control - High -0.6931472 -0.8280787    0.44
#>      Low - High -0.6931472 -0.8280787    0.41

# If trt is involved in an interaction, requesting only trt prints a note.
toy$sex <- factor(rep(c("F", "M"), length.out = 30))
toy_int <- glm(y ~ trt * sex, data = toy, family = binomial, x = TRUE)
flipscores_contrasts(toy_int, pairwise ~ trt,
                     n_flips = 100, seed = 1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>  contrast  estimate     Score p.value
#>     A - B -19.07689 -2.561738    0.02
#> 
#> NOTE: Results may be misleading due to involvement in interactions: trt:sex

# Custom coefficient contrast, similar in spirit to multcomp::glht().
# Here the row tests the trtB coefficient directly.
K <- matrix(c(0, 1), nrow = 1)
colnames(K) <- names(coef(toy_fit))
rownames(K) <- "B - A"
flipscores_contrasts(toy_fit, linfct = K,
                     n_flips = 100, seed = 1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>  contrast estimate    Score p.value
#>     B - A 3.332205 3.409809    0.01

# The same custom contrast can be supplied as a named vector.
flipscores_contrasts(toy_fit, linfct = c(trtB = 1),
                     n_flips = 100, seed = 1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>  contrast estimate    Score p.value
#>      trtB 3.332205 3.409809    0.01

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
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>            contrast  estimate     Score p.value
#>  Control - Drug | F -1.191570 -1.528571    0.13
#>  Control - Drug | M -1.196268 -1.535975    0.17

# Test selected interaction-cell contrasts. Increase n_flips in real analyses.
flipscores_contrasts(fit, ~ treatment:age_group | sex,
                     n_flips = 100, seed = 1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>                            contrast   estimate      Score p.value
#>    Drug:Middle - Control:Middle | F  1.1915705  1.5285708    0.13
#>    Control:Old - Control:Middle | F  0.2302185  0.3368640    0.74
#>       Drug:Old - Control:Middle | F  1.4217889  1.3573139    0.20
#>  Control:Young - Control:Middle | F -0.2226355 -0.3287063    0.72
#>     Drug:Young - Control:Middle | F  0.9689350  0.9355907    0.41
#>    Drug:Middle - Control:Middle | M  1.1962685  1.5359745    0.17
#>    Control:Old - Control:Middle | M  0.2302185  0.3368640    0.74
#>       Drug:Old - Control:Middle | M  1.4264869  1.3547755    0.20
#>  Control:Young - Control:Middle | M -0.2226355 -0.3287063    0.72
#>     Drug:Young - Control:Middle | M  0.9736330  0.9428786    0.37

# Linear models are accepted too.
lm_dat <- data.frame(
  y = rnorm(30),
  trt = factor(rep(c("A", "B"), each = 15))
)
lm_fit <- lm(y ~ trt, data = lm_dat, x = TRUE, y = TRUE)
flipscores_contrasts(lm_fit, pairwise ~ trt,
                     n_flips = 100, seed = 1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 100, alternative = two.sided
#> 
#>  contrast  estimate    Score p.value
#>     A - B 0.8190186 1.761655    0.08
```
