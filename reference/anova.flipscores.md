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

  (the object) `glm` (or `flipscores`) object with the model under the
  null hypothesis (i.e. the covariates, the nuisance parameters).

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
#> flipscores(formula = Y ~ Z + X, family = "poisson", data = dt)
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value Part. Cor Pr(>|z|)  
#> (Intercept)  -0.3513 -5.4234     4.2404 -1.2790    -0.169   0.1154  
#> ZB            0.2276  1.8048     2.8140  0.6414     0.079   0.4514  
#> ZC            0.9034 10.3392     3.4527  2.9945     0.348   0.0202 *
#> X             0.6381 33.1168     8.0470  4.1154     0.478   0.0124 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
#>   Df  Score Pr(>Score)  
#> Z  2 4.5460     0.0880 .
#> X  1 0.0337     0.0124 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

mod1=flipscores(Y~Z*X,data=dt,family="poisson")
summary(mod1)
#> 
#> Call:
#> flipscores(formula = Y ~ Z * X, family = "poisson", data = dt)
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value Part. Cor Pr(>|z|)  
#> (Intercept)  -0.3549 -4.3865     3.8867 -1.1286    -0.189   0.0982 .
#> ZB            0.3570  2.3488     2.6093  0.9002     0.145   0.3056  
#> ZC            0.6023  3.3952     2.3099  1.4698     0.235   0.0720 .
#> X             0.6439  8.1773     3.6816  2.2211     0.331   0.1294  
#> ZB:X         -0.5979 -3.8933     2.5112 -1.5504    -0.238   0.3418  
#> ZC:X          0.4856  2.8440     2.4290  1.1709     0.184   0.2248  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
#>                    Df  Score Pr(>Score)  
#> Model 2 vs Model 1  2 4.5658     0.0686 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
