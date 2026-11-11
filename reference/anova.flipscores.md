# anova.flipscores

This is the `anova` method for `flipscores` object. Remark: it performs
type III deviance decomposition as in `car::Anova`.

## Usage

``` r
# S3 method for class 'flipscores'
anova(object, model1 = NULL, score_type = NULL, n_flips = 5000, id = NULL, ...)
```

## Arguments

- object:

  `glm` (or `flipscores`) object with the model under the null
  hypothesis (i.e. the covariates, the nuisance parameters).

- model1:

  a `glm` (or `flipscores`) or a `matrix` (or `vector`). If it is a
  `glm` object, it has the model under the alternative hypothesis. The
  variables in `model1` are the same variables in `object` plus one or
  more variables to be tested. Alternatively, if `model1` is a `matrix`,
  it contains the tested variables column-wise.

- score_type:

  The type of score that is computed. It can be "orthogonalized",
  "effective" or "basic". Default is "orthogonalized". "effective" and
  "orthogonalized" take into account the nuisance estimation. The
  default is `NULL`, in this case the value is taken from `object`.

- n_flips:

  The number of random flips of the score contributions. When `n_flips`
  is equal or larger than the maximum number of possible flips (i.e.
  n^2), all possible flips are performed. Default is 5000.

- id:

  a `vector` identifying the clustered observations. If `NULL` (default)
  observations are assumed to be independent. NOTE: if `object` is a
  `flipscores` and `model$flip_param_call$id` is not `NULL`, this is
  considered in the inference.

- ...:

  other parameters allowed in
  [`stats::anova`](https://rdrr.io/r/stats/anova.html).

## Examples

``` r
set.seed(1)
dt=data.frame(X=scale(rnorm(50)),
   Z=factor(rep(LETTERS[1:3],length.out=50)))
dt$Y=rpois(n=nrow(dt),lambda=exp(dt$X*(dt$Z=="C")))
mod0=flipscores(Y~Z+X,data=dt,family="poisson")
summary(mod0)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "poisson(link='log')", 
#>     data = dt)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept) -0.35133 -5.42336  4.24036 -1.27899 -0.16912 0.115
#> ZB           0.22762  1.80483  2.81403  0.64137  0.07908 0.449
#> ZC           0.90344 10.33921  3.45268  2.99454  0.34843 0.023
#> X            0.63808 33.11677  8.04700  4.11542  0.47835 0.011
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 82.225  on 49  degrees of freedom
#> Residual deviance: 52.798  on 46  degrees of freedom
#> AIC: 143.12
#> 
#> Number of Fisher Scoring iterations: 5
#> 
anova(mod0)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach (5000 sign flips).
#> 
#> Model: Y ~ Z + X
#>   Ncoeff Mahalanobis      p
#> Z      2      4.5589 0.0944
#> X      1      0.0340 0.0106

mod1=flipscores(Y~Z*X,data=dt,family="poisson")
summary(mod1)
#> 
#> Call:
#> flipscores(formula = Y ~ Z * X, family = "poisson(link='log')", 
#>     data = dt)
#> 
#> Coefficients:
#>             estimate   score      se       z    pcor     p
#> (Intercept)  -0.3549 -4.3865  3.8867 -1.1286 -0.1888 0.095
#> ZB            0.3570  2.3488  2.6093  0.9002  0.1452 0.313
#> ZC            0.6023  3.3952  2.3099  1.4698  0.2347 0.065
#> X             0.6439  8.1773  3.6816  2.2211  0.3310 0.118
#> ZB:X         -0.5979 -3.8933  2.5112 -1.5504 -0.2382 0.346
#> ZC:X          0.4856  2.8440  2.4290  1.1709  0.1836 0.231
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 82.225  on 49  degrees of freedom
#> Residual deviance: 44.750  on 44  degrees of freedom
#> AIC: 139.07
#> 
#> Number of Fisher Scoring iterations: 5
#> 
anova(mod0,model1 = mod1)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach (5000 sign flips).
#> 
#> Model 1: Y ~ Z + X
#> Model 2: Y ~ Z * X
#>                    Ncoeff Mahalanobis      p
#> Model 2 vs Model 1      2      4.6533 0.0634
```
