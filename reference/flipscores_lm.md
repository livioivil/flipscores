# Robust score testing for linear models via sign-flipping

Performs robust score tests for one or more coefficients in a linear
model (ordinary least squares) by sign‑flipping score contributions. The
test is valid under heteroscedasticity, overdispersion and, in some
cases, ignored nuisance variables. It is an efficient alternative to
[`flipscores()`](https://livioivil.github.io/flipscores/reference/flipscores.md)
when the model is Gaussian and the link is identity, the `score_type`
implemented is the `"standardized"`.

## Usage

``` r
flipscores_lm(
  formula,
  data = NULL,
  n_flips = 5000,
  alternative = "two.sided",
  seed = NULL,
  tested_coeffs = NULL,
  flips = NULL,
  ...
)
```

## Arguments

- formula:

  A formula describing the linear model (e.g. `y ~ x1 + x2`). The
  response can be a numeric vector or a matrix.

- data:

  An optional data frame containing the variables in the model. If
  `NULL`, variables are taken from the environment of the formula.

- n_flips:

  Number of random sign flips. Defaults to `5000`. When `n_flips >= n^2`
  (with \\n\\ the number of observations) all possible flips are
  performed, giving an exact permutation test.

- alternative:

  Alternative hypothesis. Must be one of `"two.sided"` (default),
  `"greater"` or `"less"`.

- seed:

  Optional integer seed for reproducibility of the random flips.

- tested_coeffs:

  Character vector naming the coefficients to test. If `NULL` (default),
  all coefficients in the model are tested.

- flips:

  Optional matrix of pre‑computed sign flips. It must have `n_flips`
  rows and \\n\\ columns, where each entry is \\+1\\ or \\-1\\. If
  supplied, `n_flips` is ignored.

- ...:

  Additional arguments passed to internal functions (currently not
  used).

## Value

An object of class `c("fs_lm", ...)` with the following components:

- `Tspace`:

  A matrix of dimension `n_flips` by `k`, where
  `k = length(tested_coeffs)`. The first row contains the observed test
  statistics; the remaining rows contain the test statistics obtained
  after random sign flips.

- `summary_table`:

  A data frame with one row per tested coefficient, containing the
  observed score, the two‑sided or one‑sided p‑value, and optionally
  other summary statistics (e.g. standardised score).

- `info`:

  A list with the original formula, the names of the design matrix
  columns, and the names of the response variables.

- `call`:

  The matched call.

## Details

The function works by first projecting the response and the tested
covariate onto the orthogonal complement of the nuisance design matrix.
The resulting score contributions are then randomly sign‑flipped (with
fixed signs across observations) to generate a permutation distribution
of the test statistic.

The observed score statistic for a coefficient is \\X_r^\top Y_r /
\sqrt{\sum X_r^2}\\, where \\X_r\\ and \\Y_r\\ are the residualised
covariate and residualised response, respectively. Under the null
hypothesis the sign of each observation’s contribution is arbitrary, so
repeatedly flipping the signs yields an exact conditional test.

The method is robust to misspecification of the variance structure and
does not require the errors to be normally distributed or homoscedastic.

## References

Hemerik, J., Goeman, J. J., & Finos, L. (2020). Robust testing in
generalized linear models by sign‑flipping score contributions. *Journal
of the Royal Statistical Society Series B: Statistical Methodology*,
82(3), 841‑864.
[doi:10.1111/rssb.12369](https://doi.org/10.1111/rssb.12369)

De Santis, R., Goeman, J. J., Hemerik, J., Davenport, S., & Finos, L.
(2025). Inference in Generalized Linear Models with Robustness to
Misspecified Variances. *Journal of the American Statistical
Association*, 1‑10.
[doi:10.1080/01621459.2025.2491775](https://doi.org/10.1080/01621459.2025.2491775)

## See also

[`flipscores`](https://livioivil.github.io/flipscores/reference/flipscores.md)
for the general GLM version; `summary.fs_lm` for summarising results;
`anova.fs_lm` for joint tests across coefficients.

## Examples

``` r
set.seed(123)
n <- 100
data <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  y  = 2 * rnorm(n, mean = 0, sd = 1)   # pure noise (null holds)
)
data$y2  = data$x1+  rnorm(n, mean = 0, sd = 1)   # pure noise (null holds)

# Test the coefficient of x1
res <- flipscores_lm(cbind(y,y2) ~ x1 + x2, data = data, n_flips = 999, alternative = "two.sided")
flipscores:::summary.fs_lm(res)
#> Flip-score Linear Model
#> Call: flipscores_lm(formula = cbind(y, y2) ~ x1 + x2, data = data, 
#>     n_flips = 999, alternative = "two.sided")
#> score_type = Standardized, n_flips = , alternative = two.sided
#> 
#>  response coefficient   score p.values
#>         y (Intercept)  26.448 0.172172
#>        y2 (Intercept)  -2.651 0.782783
#>         y          x1 -21.917 0.232232
#>        y2          x1  78.354 0.001001
#>         y          x2   4.398 0.805806
#>        y2          x2   4.142 0.666667

# Test only x2
res2 <- flipscores_lm(y ~ x1 + x2, data = data, tested_coeffs = "x2")
flipscores:::summary.fs_lm(res2)
#> Flip-score Linear Model
#> Call: flipscores_lm(formula = y ~ x1 + x2, data = data, tested_coeffs = "x2")
#> score_type = Standardized, n_flips = , alternative = two.sided
#> 
#>  response coefficient score p.values
#>         y          x2 4.398   0.8184
```
