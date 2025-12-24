# Normalized Generalized R-squared for Binomial GLMs

Computes the normalized generalized R-squared coefficient for binomial
GLMs. The normalization scales the R-squared by its maximum possible
value for the specified set of terms.

## Usage

``` r
.socket_compute_gR2_n(terms, full_glm, null_glm = NULL, algorithm, control)
```

## Arguments

- terms:

  Character vector of variable names for which to compute normalized
  R-squared. If \`NULL\` (default), computes for all non-intercept terms
  in the model.

- full_glm:

  A fitted GLM object of class \`glm\` with binomial family.

- null_glm:

  A fitted GLM object of class \`glm\` (the null model). If \`NULL\`
  (default), uses an empty model (intercept-only if intercept is
  present, otherwise no predictors).

- algorithm:

  \`"auto"\` by default. It chooses between \`"brute_force"\` and
  \`"multi_start"\`

- algorithm.control:

  \`list\` of control parameters:

  - \`n_exact\` Integer specifying the sample size threshold for using
    exact methods (brute force). Default is 15.

  - \`thresholds\` Numeric vector of threshold values for multi-start
    initialization.

  - \`n_random\` Integer number of random starts for multi-start
    optimization.

  - \`max_iter\` Integer maximum number of iterations per start.

  - \`topK\` Integer number of top candidates to consider at each
    iteration.

  - \`tol\` Numeric tolerance for convergence.

  - \`patience\` Integer number of iterations without improvement before
    stopping.

## Value

A list with components:

- R2:

  The generalized R-squared coefficient for the set of terms

- R2_n:

  The normalized generalized R-squared coefficient

- algorithm:

  The algorithm used to compute the maximum R-squared

- terms_tested:

  The names of the terms included in the test

## Details

The normalized generalized R-squared is computed as: \$\$ R^2_n =
\frac{R^2}{R^2\_{\max}} \$\$ where \\R^2\_{\max}\\ is the maximum
possible R-squared value for the specified set of terms.

Different algorithms are used based on sample size:

- For small samples (\\n \leq n\_{\text{exact}}\\), brute force search
  finds the exact maximum

- For larger samples, a greedy multi-start algorithm finds approximate
  maximum

## Examples

``` r
set.seed(123)
dt <- data.frame(X = rnorm(20),
                 Z = factor(rep(LETTERS[1:3], length.out = 20)))
dt$Y <- rbinom(n = 20, prob = plogis((dt$Z == "C") * 2), size = 1)
mod <- glm(Y ~ Z + X, data = dt, family = binomial)

# Compute generalized partial correlations for all variables
(results <-  gR2_normalized_binom(mod))
#> Error in gR2_normalized_binom(mod): could not find function "gR2_normalized_binom"
# equivalent to
mod0=glm(Y~1,data=dt,family=binomial)
(results <-  gR2_normalized_binom(mod, mod0))
#> Error in gR2_normalized_binom(mod, mod0): could not find function "gR2_normalized_binom"

# Compute for specific variables only
(results <-  gR2_normalized_binom(mod,terms = c("X","Z")))
#> Error in gR2_normalized_binom(mod, terms = c("X", "Z")): could not find function "gR2_normalized_binom"


```
