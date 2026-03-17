# Confidence intervals for flipscores

This function allows to apply the `confint` method to `flipscores`
objects.

## Usage

``` r
# S3 method for class 'flipscores'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("equitailed", "symmetric"),
  score_type = c("standardized", "effective"),
  alternative = c("two.sided", "greater", "less"),
  flips = NULL,
  n_flips = NULL,
  ...
)
```

## Arguments

- object:

  a `flipscores`) model object.

- parm:

  a specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, they are recovered from the input model (`to_be_tested`).

- level:

  the confidence level required.

- type:

  the type of confidence interval to build, either "equitailed" or
  "symmetric". Default is "equitailed".

- score_type:

  The type of score that is computed. It can be "standardized" or
  "effective". Default is inherited by the input model.

- alternative:

  The direction of the test to invert, and corresponding interval with
  respect to the estimate. It can be "two.sided" (default), "greater",
  or "less".

- flips:

  the fixed flips to use in the interval search. If `NULL` (default), it
  tries to inherit it from the input model.

- n_flips:

  If `flips` are not given or retrieved by the model, the number of
  random flips of the score contributions. Default inherits them from
  the input model.

- ...:

  additional arguments for methods.

## Value

A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter. These will be labelled as (1-level)/2 and 1 -
(1-level)/2 in

## Details

All families allowed in `glm`, as well as Negative Binomial, are
supported, and inherited here by the input model `object`. See the main
function
[`flipscores`](https://livioivil.github.io/flipscores/reference/flipscores.md)
for further information.

`flipscores` is based on randomly sign-flipped score contributions.
Performing the test (and generating the confidence interval) different
times might give slightly different results. Providing a flip matrix
generated with
[`make_flips`](https://livioivil.github.io/flipscores/reference/make_flips.md)
is suggested to uniform the results.

`to_be_tested` works differently for factors and character variables:
the corresponding dummy variables are to be tested separately. So, to
test a numeric variable `X`, one could write `to_be_tested = "X"`, but
to test a factor `Z` with levels `A`, `B` and `C` one should write
`to_be_tested = c("ZB", "ZC")` as they appear in the model summary.

\#' @seealso
[`flipscores`](https://livioivil.github.io/flipscores/reference/flipscores.md),
[`make_flips`](https://livioivil.github.io/flipscores/reference/make_flips.md)

## Examples

``` r
set.seed(1)
dt <- data.frame(X=scale(rnorm(50)),
   Z=factor(rep(LETTERS[1:3],length.out=50)))
dt$Y=rpois(n=nrow(dt),lambda=exp(dt$X*(dt$Z=="C")))
mod0 <- flipscores(Y~Z+X,data=dt,family="poisson")
summary(mod0)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "poisson", data = dt)
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value Part. Cor Pr(>|z|)  
#> (Intercept)  -0.3513 -5.4234     4.2404 -1.2790    -0.169   0.1142  
#> ZB            0.2276  1.8048     2.8140  0.6414     0.079   0.4380  
#> ZC            0.9034 10.3392     3.4527  2.9945     0.348   0.0232 *
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
confint(mod0)
#>                  2.5 %     97.5 %
#> (Intercept) -1.0379192 0.07774207
#> ZB          -0.3950854 0.91736873
#> ZC           0.1714362 1.63156422
#> X            0.1783036 1.14192973

xx <- rnorm(20)
zz <- rnorm(20, 0.2 * xx)
yy <- rnorm(20, 1 + 2 * xx - zz)
flps <- make_flips(20, 5000)
mod1 <- flipscores(yy ~ xx + zz)
summary(mod1)
#> 
#> Call:
#> flipscores(formula = yy ~ xx + zz)
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value Part. Cor Pr(>|z|)    
#> (Intercept)    1.399  15.584      7.568   2.059     0.485    2e-04 ***
#> xx             2.310  16.986     10.918   1.556     0.367    4e-04 ***
#> zz            -1.007 -12.830      5.708  -2.248    -0.530    1e-03 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.940885)
#> 
#>     Null deviance: 116.830  on 19  degrees of freedom
#> Residual deviance:  15.995  on 17  degrees of freedom
#> AIC: 60.288
#> 
#> Number of Fisher Scoring iterations: 2
#> 
confint(mod1, flips = flps)
#> Error in eval(parse(text = .form___)): object 'yy' not found
```
