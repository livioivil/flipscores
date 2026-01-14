# Normalized Generalized Partial Correlation for Binomial GLMs

Computes the normalized generalized partial correlation coefficient for
binomial GLMs. The normalization scales the correlation by its maximum
possible absolute value.

## Usage

``` r
compute_gcor_normalized_binom(
  model0,
  X,
  algorithm = "auto",
  algorithm.control = list(n_exact = 15, thresholds = c(-0.1, 0, 0.1), n_random = 10,
    max_iter = 1000, topK = 10, tol = 1e-12, patience = 10),
  ...
)
```

## Arguments

- algorithm:

  \`"auto"\` by default. It choose between \`"intercept_only"\`,
  \`"brute_force"\` and \`"multi_start"\`

- algorithm.control:

  \`list\` of coltrol parameters: \`n_exact\` Integer specifying the
  sample size threshold for using exact methods (brute force). Default
  is 15. \`thresholds\` Numeric vector of threshold values for
  multi-start initialization. \`n_random\` Integer number of random
  starts for multi-start optimization. \`max_iter\` Integer maximum
  number of iterations per start. \`topK\` Integer number of top
  candidates to consider at each iteration. \`tol\` Numeric tolerance
  for convergence. \`patience\` Integer number of iterations without
  improvement before stopping.

- full_glm:

  A fitted GLM object of class \`glm\` with binomial family.

- terms:

  Character vector of variable names for which to compute normalized
  correlations. If \`NULL\` (default), computes for all non-intercept
  terms in the model.

- intercept_too:

  Logical indicating whether to include the intercept as a variable.
  Default is FALSE.

## Value

A data frame with five columns:

- terms:

  The variable name

- r:

  The generalized partial correlation coefficient

- r_n:

  The normalized generalized partial correlation coefficient

- null_model:

  The null model used to compute the generalized (partial) correlation

- algorithm:

  The algorithm used to compute the upper/lower bounds of the
  generalized partial correlation coefficient (to compute its normalized
  version)

## Details

The normalized generalized partial correlation is computed as: \$\$ r_n
= \begin{cases} +r / r\_+ & \text{if } r \> 0 \\ -r / r\_- & \text{if }
r \< 0 \end{cases} \$\$ where \\r\_+\\ is the maximum possible
correlation and \\r\_-\\ is the minimum.

When the (full) model has only intercept and only one predictor X, the
generalized (non partial) correlation is computed and the normalization
factor for X is exact. In the more general case with more predictors,
for sample sizes \\n \leq n\_{\text{exact}}\\, brute force search is
used to find the exact extrema. For larger sample sizes, a greedy
multi-start algorithm is employed:

- Multiple starting points are generated using thresholding and random
  sampling

- From each start, coordinates are greedily flipped to improve the
  correlation

- Early stopping is used when no improvements are found for several
  iterations

- The best solution across all starts is returned

This approach provides a good trade-off between computational efficiency
and solution quality for large problems where brute force is infeasible.

## Examples

``` r
set.seed(123)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rbinom(n=20,prob=plogis((dt$Z=="C")*2),size=1)
mod=flipscores(Y~Z+X,data=dt,family="binomial",n_flips=1000)
summary(mod)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "binomial", data = dt, 
#>     n_flips = 1000)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)  
#> (Intercept)  -0.1486  -0.2102     1.1881  -0.1770    -0.067    0.910  
#> ZB          -20.4539  -1.4784     0.7466  -1.9802    -0.530    0.057 .
#> ZC           20.8561   1.8043     0.8180   2.2057     0.615    0.037 *
#> X            -0.4276  -0.3782     0.9574  -0.3951    -0.149    0.735  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 27.5256  on 19  degrees of freedom
#> Residual deviance:  9.4015  on 16  degrees of freedom
#> AIC: 17.401
#> 
#> Number of Fisher Scoring iterations: 19
#> 

(results <- gcor_normalized_binom(mod))
#> Error in gcor_normalized_binom(mod): could not find function "gcor_normalized_binom"
# Compute for specific terms only
gcor_normalized_binom(mod, terms = c("X", "ZC"))
#> Error in gcor_normalized_binom(mod, terms = c("X", "ZC")): could not find function "gcor_normalized_binom"

```
