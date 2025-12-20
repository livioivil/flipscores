# flipscores: Robust Score Testing in GLMs, by Sign-Flip Contributions

Provides robust tests for testing in GLMs, by sign-flipping score
contributions. The tests are robust against overdispersion,
heteroscedasticity and, in some cases, ignored nuisance variables. See
Hemerik, Goeman and Finos (2020)
[doi:10.1111/rssb.12369](https://doi.org/10.1111/rssb.12369) .

It provides robust tests for testing in GLMs, by sign-flipping score
contributions. The tests are often robust against overdispersion,
heteroscedasticity and, in some cases, ignored nuisance variables.

## See also

Useful links:

- <https://livioivil.github.io/flipscores/>

Useful links:

- <https://livioivil.github.io/flipscores/>

## Author

**Maintainer**: Livio Finos <livio.finos@unipd.it>
([ORCID](https://orcid.org/0000-0003-3181-8078))

Other contributors:

- Jelle J. Goeman \[contributor\]

- Jesse Hemerik \[contributor\]

- Riccardo De Santis \[contributor\]

Livio Finos, Jelle Goeman and Jesse Hemerik, with contribution of
Riccardo De Santis.

## Examples

``` r
set.seed(1)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rpois(n=20,lambda=exp(dt$X))
mod=flipscores(Y~Z+X,data=dt,family="poisson",x=TRUE)
summary(mod)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "poisson", data = dt, 
#>     x = TRUE)
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value Part. Cor Pr(>|z|)  
#> (Intercept)  -0.1026 -0.7229     2.7127 -0.2665    -0.088    0.752  
#> ZB           -0.1501 -0.7125     2.1789 -0.3270    -0.104    0.658  
#> ZC            0.1633  0.8106     2.2232  0.3646     0.117    0.697  
#> X             0.9439 16.2058     4.7272  3.4282     0.671    0.011 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 27.135  on 19  degrees of freedom
#> Residual deviance: 12.888  on 16  degrees of freedom
#> AIC: 57.459
#> 
#> Number of Fisher Scoring iterations: 5
#> 

# Anova test
anova(mod)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach (5000 sign flips).
#> 
#> Model: Y ~ Z + X
#>   Df   Score Pr(>Score)  
#> Z  2 0.73321     0.7136  
#> X  1 0.02957     0.0110 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# or
mod0=flipscores(Y~Z,data=dt,family="poisson",x=TRUE)
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
#> Model 2 vs Model 1  1 0.029251     0.0098 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# and
mod0=flipscores(Y~X,data=dt,family="poisson")
anova(mod0,mod)
#> Error in model.frame.default(formula = Y ~ X, data = dt, drop.unused.levels = TRUE): 'data' must be a data.frame, environment, or list
```
