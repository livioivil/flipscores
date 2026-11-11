# Flip Scores Test

Extended version of `flipscores()` supporting joint testing when
`formula` is a list of glm objects, a list of formulas, or a formula
with a matrix response.

## Usage

``` r
flipscores(
  formula,
  family = gaussian(),
  data = NULL,
  score_type = "standardized",
  n_flips = 5000,
  alternative = "two.sided",
  id = NULL,
  seed = NULL,
  to_be_tested = NULL,
  tested_coeffs = NULL,
  flips = NULL,
  precompute_flips = TRUE,
  ...
)
```

## Arguments

- formula:

  One of:

  1.  A `formula` object (standard use, returns `flipscores` object)

  2.  A list of fitted `glm` objects (returns `jfs` object)

  3.  A list of `formula` objects (returns `jfs` object)

  4.  A `formula` with matrix response e.g. `cbind(y1,y2) ~ x` (returns
      `jfs` object)

- family:

  Error distribution for `glm`. Used in cases 1, 3, 4.

- data:

  A data frame. Used in cases 1, 3, 4.

- score_type:

  Type of score. One of `"standardized"` (default), `"orthogonalized"`,
  `"effective"`, `"basic"`.

- n_flips:

  Number of random flips. Default `5000`.

- alternative:

  Alternative hypothesis. Default `"two.sided"`.

- id:

  Optional cluster identifiers, also used to align score rows across
  models. Use a vector with one value per score row, a column name, or a
  one-sided formula such as `~ids`. Column names and formulas are
  evaluated in `data`, or separately in each fitted model's original
  data for a list of models. With `NULL`, rows are matched by row names.
  Repeated IDs are aligned by occurrence, with the maximum multiplicity
  across models, and all occurrences share the same random sign. A
  supplied vector is applied unchanged to every model and must match
  each model's number of score rows, after `subset` and missing-value
  removal. A length mismatch is an error; values are never recycled or
  truncated. For models with different fitted row counts, consider
  `id = ~ids` or `id = "ids"`, with an ID column present in every
  model's data.

- seed:

  Random seed for reproducibility.

- to_be_tested:

  Character vector or list of coefficient names to test. If a character
  vector, the same terms are intersected with each model's coefficients.
  If a list, each element applies to the corresponding model. If `NULL`
  all coefficients are tested.

- tested_coeffs:

  Alias for `to_be_tested` for compatibility.

- flips:

  Optional matrix of pre-computed flips. Columns may be named by
  observations or, when `id` is supplied, by cluster identifiers. All
  observations in the same cluster must receive the same sign. For joint
  tests with `id`, observation columns follow the common rows in
  `alignment$rows`; alternatively supply columns named by ID.

- precompute_flips:

  Logical, default `TRUE`. If `FALSE`, GLM flips are generated on demand
  using the same observation or cluster mapping.

- ...:

  Additional arguments passed to
  [`glm()`](https://rdrr.io/r/stats/glm.html) or the internal flipscores
  engine.

## Value

A `flipscores` object (case 1) or a `jfs` object (cases 2, 3, 4). Joint
objects retain local score matrices in `objects`, with resolved
identifiers in `attr(scores, "id")`, and an `alignment` list containing
the common `rows` and local row-to-common-row `index` maps.

## Examples

``` r
set.seed(1)
n <- 50
df <- data.frame(
  y1  = rnorm(n),
  y2  = rnorm(n),
  yb  = rbinom(n, 1, 0.5),
  y1p = rpois(n, lambda = 3),
  y2p = rpois(n, lambda = 5),
  x1  = rnorm(n),
  x2  = rnorm(n)
)

# Case 1: single formula
res=flipscores(y1 ~ x1 + x2, family = gaussian(), data = df)
summary(res)
#> 
#> Call:
#> flipscores(formula = y1 ~ x1 + x2, family = "gaussian(link='identity')", 
#>     data = df)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.10358  5.15889  5.96233  0.86525  0.12489 0.390
#> x1          -0.04491 -2.41629  6.15867 -0.39234 -0.05663 0.709
#> x2           0.02313  1.20293  6.04796  0.19890  0.02871 0.825
#> 
#> (Dispersion parameter for gaussian family taken to be 0.7175968)
#> 
#>     Null deviance: 33.870  on 49  degrees of freedom
#> Residual deviance: 33.727  on 47  degrees of freedom
#> AIC: 130.21
#> 
#> Number of Fisher Scoring iterations: 2
#> 

# Case 2: list of glm objects
models <- list(
  glm(y1 ~ x1 + x2, data = df, family = gaussian()),
  glm(y2 ~ x1 + x2, data = df, family = gaussian()),
  glm(yb ~ x1 + x2, data = df, family = binomial())
)
res=flipscores(models)
summary(res)
#>   model response coefficient estimate   score    se       z     pcor      p
#> 1  mod1       y1 (Intercept)  0.10358   5.159 5.962  0.8652  0.12489 0.3936
#> 2  mod1       y1          x1 -0.04491  -2.416 6.159 -0.3923 -0.05663 0.7040
#> 3  mod1       y1          x2  0.02313   1.203 6.048  0.1989  0.02871 0.8122
#> 4  mod2       y2 (Intercept)  0.13022   6.486 6.650  0.9752  0.14076 0.3394
#> 5  mod2       y2          x1 -0.22619 -12.171 7.066 -1.7225 -0.24863 0.0364
#> 6  mod2       y2          x2 -0.17934  -9.329 6.862 -1.3594 -0.19621 0.1288
#> 7  mod3       yb (Intercept) -0.23693  -2.705 3.385 -0.7992 -0.11326 0.4378
#> 8  mod3       yb          x1 -0.47375  -5.801 3.590 -1.6159 -0.22891 0.0918
#> 9  mod3       yb          x2 -0.38289  -4.558 3.490 -1.3058 -0.18499 0.2236
res=flipscores(models, to_be_tested = "x1")
summary(res)
#>   model response coefficient estimate   score    se       z     pcor      p
#> 1  mod1       y1          x1 -0.04491  -2.416 6.159 -0.3923 -0.05663 0.7000
#> 2  mod2       y2          x1 -0.22619 -12.171 7.066 -1.7225 -0.24863 0.0342
#> 3  mod3       yb          x1 -0.47375  -5.801 3.590 -1.6159 -0.22891 0.0926

# Different dataset sizes: resolve IDs separately in each model's data
df$ids <- rep(seq_len(n / 2), each = 2)
df_short <- df[-c(2, 5, 9), ]
models_by_id <- list(
  full = glm(y1 ~ x1 + x2, data = df),
  shorter = glm(y2 ~ x1 + x2, data = df_short)
)
# A shared vector of length 50 cannot match the shorter model's 47 rows:
tryCatch(
  flipscores(models_by_id, id = df$ids, to_be_tested = "x1", n_flips = 99),
  error = function(e) message(conditionMessage(e))
)
#> The supplied 'id' vector has length 50, but model 'shorter' has 47 score rows. A vector must have one value per score row in every model, after subset and missing-value removal. Consider using id = ~ids or id = "ids" to read an ID column from each model's own data.
# These forms support different row counts and cluster flips by the same IDs:
res_by_id <- flipscores(models_by_id, id = ~ids, to_be_tested = "x1",
                        n_flips = 99, seed = 1)
res_by_name <- flipscores(models_by_id, id = "ids", to_be_tested = "x1",
                          n_flips = 99, seed = 1)

# Case 3: list of formulas, same family
res=flipscores(
  list(y1 ~ x1 + x2, y2 ~ x1 + x2),
  family = gaussian(),
  data   = df
)
summary(res)
#>   model response coefficient estimate   score    se       z     pcor      p
#> 1  mod1       y1 (Intercept)  0.10358   5.159 5.962  0.8652  0.12489 0.3904
#> 2  mod1       y1          x1 -0.04491  -2.416 6.159 -0.3923 -0.05663 0.6980
#> 3  mod1       y1          x2  0.02313   1.203 6.048  0.1989  0.02871 0.8264
#> 4  mod2       y2 (Intercept)  0.13022   6.486 6.650  0.9752  0.14076 0.3394
#> 5  mod2       y2          x1 -0.22619 -12.171 7.066 -1.7225 -0.24863 0.0292
#> 6  mod2       y2          x2 -0.17934  -9.329 6.862 -1.3594 -0.19621 0.1324

# Case 4: matrix response, gaussian
res=flipscores(cbind(y1, y2) ~ x1 + x2, family = gaussian(), data = df)
summary(res)
#>   model response coefficient estimate   score    se       z     pcor      p
#> 1  mod1       y1 (Intercept)  0.10358   5.159 5.962  0.8652  0.12489 0.3884
#> 2  mod1       y1          x1 -0.04491  -2.416 6.159 -0.3923 -0.05663 0.6964
#> 3  mod1       y1          x2  0.02313   1.203 6.048  0.1989  0.02871 0.8140
#> 4  mod2       y2 (Intercept)  0.13022   6.486 6.650  0.9752  0.14076 0.3310
#> 5  mod2       y2          x1 -0.22619 -12.171 7.066 -1.7225 -0.24863 0.0322
#> 6  mod2       y2          x2 -0.17934  -9.329 6.862 -1.3594 -0.19621 0.1350

# Case 4: matrix response, poisson
res=flipscores(cbind(y1p, y2p) ~ x1 + x2, family = poisson(), data = df)
summary(res)
#>   model response coefficient estimate   score     se       z    pcor      p
#> 1  mod1      y1p (Intercept)  0.88061  69.644  6.967  9.9964  0.6718 0.0002
#> 2  mod1      y1p          x1 -0.12451 -16.043 11.393 -1.4082 -0.1930 0.2144
#> 3  mod1      y1p          x2 -0.07236  -9.056 11.132 -0.8135 -0.1136 0.3958
#> 4  mod2      y2p (Intercept)  1.67965 208.330  6.561 31.7534  0.8422 0.0002
#> 5  mod2      y2p          x1  0.08768  25.725 17.094  1.5049  0.1975 0.2462
#> 6  mod2      y2p          x2  0.05057  14.288 16.878  0.8466  0.1140 0.4244

# more examples
n=20
DD=data.frame(X=rnorm(n),Z1=rnorm(n),Z2=rnorm(n))
DD$Y=DD$Z1+DD$X+rnorm(n)
# Run four glms abd combine it in a list
mod1=glm(Y~X+Z1+Z2,data=DD)
mod2=glm(Y~X+poly(Z1,2)+Z2,data=DD)
mod3=glm(Y~X+poly(Z1,2)+poly(Z2,2),data=DD)
mod4=glm(Y~X+Z1+poly(Z2,2),data=DD)
mods=list(mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)
res=flipscores(mods, to_be_tested = "X")
summary(res)
#>   model response coefficient estimate score    se     z   pcor      p
#> 1  mod1        Y           X   1.0095 8.567 3.370 2.542 0.6165 0.0248
#> 2  mod2        Y           X   0.9817 8.209 3.370 2.436 0.6090 0.0452
#> 3  mod3        Y           X   0.9732 8.096 3.447 2.348 0.6063 0.0530
#> 4  mod4        Y           X   1.0040 8.492 3.453 2.459 0.6148 0.0278
```
