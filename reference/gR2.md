# Generalized R-squared for GLM Models

Computes the generalized R-squared measure for nested generalized linear
models. The generalized R-squared measures the proportion of "variance"
explained by the additional predictors in the full model compared to the
null model.

## Usage

``` r
gR2(
  full_glm,
  null_glm = NULL,
  terms = NULL,
  normalize = FALSE,
  adjusted = FALSE,
  algorithm = "auto",
  algorithm.control = list(n_exact = 15, thresholds = c(-0.1, 0, 0.1), n_random = max(1,
    13 + log(1/nrow(model.matrix(full_glm)))), max_iter = 1000, topK = max(10, min(100,
    length(nrow(model.matrix(full_glm)))/10)), tol = 1e-12, patience = 10)
)
```

## Arguments

- full_glm:

  A fitted GLM object of class \`glm\` (the full model).

- null_glm:

  A fitted GLM object of class \`glm\` (the null model). If \`NULL\`
  (default), uses an empty model (intercept-only if intercept is
  present, otherwise no predictors).

- terms:

  A character vector of variable names or a formula specifying the
  additional terms in the full model compared to the null. If provided,
  this overrides \`null_glm\` and the null model is refitted excluding
  these terms.

- normalize:

  FALSE by default.

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

- algorith:

  \`"auto"\` by default. It choose between \`"intercept_only"\`,
  \`"brute_force"\` and \`"multi_start"\`

## Value

If normalize==FALSE: A numeric value representing the generalized
R-squared measure. If normalize==TRUE: A list with components:

- R2:

  The generalized R-squared coefficient for the set of terms

- R2_n:

  The normalized generalized R-squared coefficient

- algorithm:

  The algorithm used to compute the maximum R-squared

- terms_tested:

  The names of the terms included in the test

## Details

The generalized R-squared is computed as:

\$\$gR^2 = \frac{Y_r^\top X_r (X_r^\top X_r)^{-1} X_r^\top Y_r}{Y_r^\top
Y_r}\$\$

where:

- \\Y_r = V^{-1/2}(Y - \hat{\mu}\_0)\\ is the standardized residual
  vector from the null model

- \\X_r = (I - H)W^{1/2}X\\ is the residualized additional predictors
  matrix

- \\H\\ is the hat matrix for the null model

- \\W = DV^{-1}D\\ is the weight matrix

This measures the proportion of the standardized residual sum of squares
explained by the additional predictors in the full model.

The normalized generalized R-squared is computed as: \$\$ R^2_n =
\frac{R^2}{R^2\_{\max}} \$\$ where \\R^2\_{\max}\\ is the maximum
possible R-squared value for the specified set of terms.

Different algorithms are used based on sample size:

- For small samples (\\n \leq n\_{\text{exact}}\\), brute force search
  finds the exact maximum

- For larger samples, a greedy multi-start algorithm finds approximate
  maximum

## Author

Livio Finos and Paolo Girardi

## Examples

``` r
set.seed(1)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rpois(n=20,lambda=exp(dt$Z=="C"))
mod=glm(Y~Z+X,data=dt,family="poisson")
summary(mod)
#> 
#> Call:
#> glm(formula = Y ~ Z + X, family = "poisson", data = dt)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) -0.14256    0.40940  -0.348  0.72768   
#> ZB          -0.18558    0.60570  -0.306  0.75931   
#> ZC           1.40981    0.46298   3.045  0.00233 **
#> X           -0.06964    0.20905  -0.333  0.73906   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 28.649  on 19  degrees of freedom
#> Residual deviance: 11.218  on 16  degrees of freedom
#> AIC: 58.102
#> 
#> Number of Fisher Scoring iterations: 5
#> 

# Compute generalized partial correlations for all variables
(results <- gR2(mod))
#>           terms       gR2 null_model
#> 1 ~ ZB + ZC + X 0.6957557      Y ~ 1
# equivalent to
mod0=glm(Y~0,data=dt,family="poisson")
(results <- gR2(mod, mod0))
#>               terms       gR2 null_model
#> 1 ~ 1 + ZB + ZC + X 0.7378818      Y ~ 0
(results <- gR2(mod, mod0,normalize=TRUE))
#> Warning: The Normalized Generalized Partial Correlation (Determination) Coefficient for Count families without interncept in the null model has not implemented, yet. NA will be returned.
#>                         terms       gR2 gR2_n algorithm exact null_model
#> 1 ~ (Intercept) + ZB + ZC + X 0.7378818    NA        NA    NA      Y ~ 0

# Compute for specific variables only
(results <- gR2(mod,terms = c("X","Z")))
#> Error in model.frame.default(formula = Y ~ Z, data = dt, drop.unused.levels = TRUE): 'data' must be a data.frame, environment, or list
(results <- gR2(mod,terms = c("X","Z")))
#> Error in model.frame.default(formula = Y ~ Z, data = dt, drop.unused.levels = TRUE): 'data' must be a data.frame, environment, or list


set.seed(123)
dt <- data.frame(X = rnorm(20),
                 Z = factor(rep(LETTERS[1:3], length.out = 20)))
dt$Y <- rbinom(n = 20, prob = plogis((dt$Z == "C") * 2), size = 1)
mod <- glm(Y ~ Z + X, data = dt, family = binomial)

# Compute generalized partial correlations for all variables
(results <-  gR2(mod,normalize=TRUE))
#>           terms       gR2     gR2_n   algorithm exact null_model
#> 1 ~ ZB + ZC + X 0.6553299 0.6553299 multi_start FALSE      Y ~ 1
# equivalent to
mod0=glm(Y~1,data=dt,family=binomial)
(results <-  gR2(mod, mod0,normalize=TRUE))
#>           terms       gR2     gR2_n   algorithm exact null_model
#> 1 ~ ZB + ZC + X 0.6553299 0.6553299 multi_start FALSE      Y ~ 1

# Compute for specific variables only
(results <-  gR2(mod,terms = c("X","Z"),normalize=TRUE))
#> Error in model.frame.default(formula = Y ~ 1 - 1, data = dt, drop.unused.levels = TRUE): 'data' must be a data.frame, environment, or list


# Compute generalized (non partial!) correlations for all variables
mod <- glm(Y ~ X, data = dt, family = binomial)
(results <-  gR2(mod,normalize=TRUE))
#>   terms          gR2        gR2_n   algorithm exact null_model
#> 1   ~ X 0.0004224238 0.0006806998 multi_start FALSE      Y ~ 1
# note the difference:
(results <-  gR2(mod,normalize=TRUE,algorithm="intercept_only"))
#>   terms          gR2        gR2_n      algorithm exact null_model
#> 1   ~ X 0.0004224238 0.0006806998 intercept_only  TRUE      Y ~ 1
# Despite the result is the same in this case,
# the multi_start algorithm does not ensure exactness (while intercept_only and brute_force do)
(results <-  gR2(mod,normalize=TRUE,algorithm="multi_start"))
#>   terms          gR2        gR2_n   algorithm exact null_model
#> 1   ~ X 0.0004224238 0.0006806998 multi_start FALSE      Y ~ 1
```
