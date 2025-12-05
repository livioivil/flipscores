# Generalized R-squared for GLM Models

Computes the generalized R-squared measure for nested generalized linear
models. The generalized R-squared measures the proportion of "variance"
explained by the additional predictors in the full model compared to the
null model.

## Usage

``` r
gR2(full_glm, null_glm = NULL, terms = NULL)
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

## Value

A numeric value representing the generalized R-squared measure.

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
mod0=glm(Y~1,data=dt,family="poisson")
(results <- gR2(mod, mod0))
#>           terms       gR2 null_model
#> 1 ~ ZB + ZC + X 0.6957557      Y ~ 1

# Compute for specific variables only
(results <- gR2(mod,terms = c("X","Z")))
#> Error in model.frame.default(formula = Y ~ Z, data = dt, drop.unused.levels = TRUE): 'data' must be a data.frame, environment, or list
```
