# flipscores

## The “flipscores” package

The package provides functions for hypothesis tests in GLMs based on
sign-flipping score contributions. The tests are robust against
overdispersion, heteroscedasticity and, in some cases, ignored nuisance
variables.

In this vignette, we provide some use cases and a simulation study to
show possible usage of the package.

The `flipscores` function is the main function of the package. It
supports the standard usage of `R` modeling functions, such as `glm`,
from which it inherits the majority of the arguments.

``` r
nn <- 30
set.seed(123)
dt <- data.frame(X = rnorm(nn), Z = factor(rep(LETTERS[1:3], length.out = nn)))
dt$Y <- rpois(n = nn, lambda = exp(dt$X))
mod <- flipscores(Y ~ Z + X,
                  data = dt,
                  family = "poisson",
                  x = TRUE)
summary(mod)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "poisson", data = dt, 
#>     x = TRUE)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)    
#> (Intercept)  0.02207  0.26383    3.44133  0.07666     0.016   0.9416    
#> ZB           0.02318  0.13383    2.39994  0.05576     0.012   0.9582    
#> ZC           0.15946  1.36568    2.93906  0.46467     0.097   0.6474    
#> X            0.89732 35.59669    6.63014  5.36892     0.726   0.0004 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
```

Functions like `summary` or `anova` can be applied to `flipscores`
objects, obtaining an output in the usual format. The p-values are
obtained by randomly sign-flipping score contributions and comparing the
observed score with the flipped scores.

``` r
# Anova test
anova(mod)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach (5000 sign flips).
#> 
#> Model: Y ~ Z + X
#>   Df    Score Pr(>Score)    
#> Z  2 0.243224     0.9088    
#> X  1 0.042297     0.0004 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# or
mod0 <- flipscores(Y~Z,data=dt,family="poisson",x=TRUE)
anova(mod0,mod)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach (5000 sign flips).
#> 
#> Model 1: Y ~ Z
#> Model 2: Y ~ Z + X
#>                    Df    Score Pr(>Score)    
#> Model 2 vs Model 1  1 0.041756      4e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# and
mod0 <- flipscores(Y~X,data=dt,family="poisson")
anova(mod0,mod)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach (5000 sign flips).
#> 
#> Model 1: Y ~ X
#> Model 2: Y ~ Z + X
#>                    Df   Score Pr(>Score)
#> Model 2 vs Model 1  2 0.27248     0.8886
```

Equivalently, one might apply the `flipscores` function to a precomputed
`glm` model.

``` r
model <- glm(Y~Z+X,data=dt,family="poisson")
mod2 <- flipscores(model)
summary(mod2)
#> 
#> Call:
#> flipscores(formula = model)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)    
#> (Intercept)  0.02207  0.26383    3.44133  0.07666     0.016   0.9364    
#> ZB           0.02318  0.13383    2.39994  0.05576     0.012   0.9600    
#> ZC           0.15946  1.36568    2.93906  0.46467     0.097   0.6506    
#> X            0.89732 35.59669    6.63014  5.36892     0.726   0.0004 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
```

Notice that the p-values are not identical to the ones of the previous
model. In fact, computing all possible combinations of flipped score
contributions would be computationally unfeasible. Only a subset of this
possibilities is computed, regulated by the parameter `n_flips`, with
default at `5000`. Reproducibility can be achieved by setting the
argument `seed` or precomputing flips with the function `make_flips`.

``` r
flps <- make_flips(n_obs = nn, n_flips = 1000)
summary(flipscores(model, flips = flps))
#> 
#> Call:
#> flipscores(formula = model, flips = flps)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)    
#> (Intercept)  0.02207  0.26383    3.44133  0.07666     0.016    0.945    
#> ZB           0.02318  0.13383    2.39994  0.05576     0.012    0.956    
#> ZC           0.15946  1.36568    2.93906  0.46467     0.097    0.653    
#> X            0.89732 35.59669    6.63014  5.36892     0.726    0.001 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
summary(flipscores(model, flips = flps))
#> 
#> Call:
#> flipscores(formula = model, flips = flps)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)    
#> (Intercept)  0.02207  0.26383    3.44133  0.07666     0.016    0.945    
#> ZB           0.02318  0.13383    2.39994  0.05576     0.012    0.956    
#> ZC           0.15946  1.36568    2.93906  0.46467     0.097    0.653    
#> X            0.89732 35.59669    6.63014  5.36892     0.726    0.001 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
```

One-sided testing alternatives are also covered

``` r
summary(flipscores(model, alternative = "less"))
#> 
#> Call:
#> flipscores(formula = model, alternative = "less")
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)
#> (Intercept)  0.02207  0.26383    3.44133  0.07666     0.016    0.533
#> ZB           0.02318  0.13383    2.39994  0.05576     0.012    0.516
#> ZC           0.15946  1.36568    2.93906  0.46467     0.097    0.680
#> X            0.89732 35.59669    6.63014  5.36892     0.726    1.000
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
summary(flipscores(model, alternative = "greater"))
#> 
#> Call:
#> flipscores(formula = model, alternative = "greater")
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)    
#> (Intercept)  0.02207  0.26383    3.44133  0.07666     0.016   0.4654    
#> ZB           0.02318  0.13383    2.39994  0.05576     0.012   0.4816    
#> ZC           0.15946  1.36568    2.93906  0.46467     0.097   0.3228    
#> X            0.89732 35.59669    6.63014  5.36892     0.726   0.0002 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
```

If only a subset of parameters is of interest, the `to_be_tested`
parameter can be set to avoid computing the p-values for all the others

``` r
summary(flipscores(model, to_be_tested = "X", flips = flps))
#> 
#> Call:
#> flipscores(formula = model, to_be_tested = "X", flips = flps)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)    
#> (Intercept)  0.02207       NA         NA       NA        NA       NA    
#> ZB           0.02318       NA         NA       NA        NA       NA    
#> ZC           0.15946       NA         NA       NA        NA       NA    
#> X            0.89732 35.59669    6.63014  5.36892     0.726    0.001 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
summary(flipscores(model, to_be_tested = c(1, 4), flips = flps))
#> 
#> Call:
#> flipscores(formula = model, to_be_tested = c(1, 4), flips = flps)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)    
#> (Intercept)  0.02207  0.26383    3.44133  0.07666     0.016    0.945    
#> ZB           0.02318       NA         NA       NA        NA       NA    
#> ZC           0.15946       NA         NA       NA        NA       NA    
#> X            0.89732 35.59669    6.63014  5.36892     0.726    0.001 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
```

`to_be_tested` accepts vectors with both the name of the variables and
indexing numbers; these indicators refer to the name or the cardinality
of the variable in the summary output, so the Intercept will always
correspond to number 1, while particular attention should be paid if
`factor/character` variables are involved.

``` r
summary(flipscores(model, to_be_tested = "ZB", flips = flps))
#> 
#> Call:
#> flipscores(formula = model, to_be_tested = "ZB", flips = flps)
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value Part. Cor Pr(>|z|)
#> (Intercept)  0.02207      NA         NA      NA        NA       NA
#> ZB           0.02318 0.13383    2.39994 0.05576     0.012    0.956
#> ZC           0.15946      NA         NA      NA        NA       NA
#> X            0.89732      NA         NA      NA        NA       NA
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
summary(flipscores(model, to_be_tested = 2, flips = flps))
#> 
#> Call:
#> flipscores(formula = model, to_be_tested = 2, flips = flps)
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value Part. Cor Pr(>|z|)
#> (Intercept)  0.02207      NA         NA      NA        NA       NA
#> ZB           0.02318 0.13383    2.39994 0.05576     0.012    0.956
#> ZC           0.15946      NA         NA      NA        NA       NA
#> X            0.89732      NA         NA      NA        NA       NA
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 60.819  on 29  degrees of freedom
#> Residual deviance: 26.638  on 26  degrees of freedom
#> AIC: 84.227
#> 
#> Number of Fisher Scoring iterations: 5
```

`flipscores` supports all the available families in `glm`, and also the
Negative Binomial family. This, in particular, should be indicated as
`"negbinom"` between quotes.

``` r
summary(flipscores(Y ~ Z + X, data = dt, family = "negbinom"))
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "negbinom", data = dt)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)    
#> (Intercept)  0.02208  0.26385    3.18354  0.08288     0.016   0.9398    
#> ZB           0.02316  0.13368    2.20995  0.06049     0.012   0.9552    
#> ZC           0.15946  1.36541    2.72310  0.50142     0.096   0.6484    
#> X            0.89733 17.36132    4.68164  3.70839     0.714   0.0004 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(20776.66) family taken to be 0.8801237)
#> 
#>     Null deviance: 60.815  on 29  degrees of freedom
#> Residual deviance: 26.636  on 26  degrees of freedom
#> AIC: 86.227
#> 
#> Number of Fisher Scoring iterations: 1
```

The `confint` method is also applicable to provide confidence intervals
based on the inversion of the sign-flip score tests. These intervals
inherit the properties of robustness to variance misspecification of the
test.

``` r
confint(mod, n_flips = 1000)
#>                  2.5 %    97.5 %
#> (Intercept) -0.8163922 0.4984970
#> ZB          -1.5591010 0.8733270
#> ZC          -0.6914538 0.8176993
#> X            0.6327678 1.2986202
```

As before, various options are allowed, including the possibility to
generate symmetric confidence intervals.

``` r
confint(mod, parm = 4, flips = flps) # parm == to_be_tested!
#>       2.5 %   97.5 %
#> X 0.6467617 1.279425
confint(mod, parm = 4, level = 0.99, flips = flps)
#>       0.5 %   99.5 %
#> X 0.5443523 1.824792
confint(mod, parm = 4, alternative = "greater", flips = flps)
#>        Est.     95 %
#> X 0.8973217 1.183059
confint(mod, parm = 4, type = "symmetric", flips = flps)
#>       2.5 %   97.5 %
#> X 0.5751064 1.219537
```
