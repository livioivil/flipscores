# flipscores: Robust Score Testing in GLMs, by Sign-Flip Contributions

Provides robust tests for testing in GLMs, by sign-flipping score
contributions. The tests are robust against overdispersion,
heteroscedasticity and, in some cases, ignored nuisance variables. See
Hemerik, Goeman and Finos (2020)
[doi:10.1111/rssb.12369](https://doi.org/10.1111/rssb.12369) , De
Santis, Goeman, Hemerik, Davenport, Finos (2025)
[doi:10.1080/01621459.2025.2491775](https://doi.org/10.1080/01621459.2025.2491775)
and De Santis, Goeman, Davenport, Hemerik, Finos (2025)
[doi:10.1214/25-EJS2409](https://doi.org/10.1214/25-EJS2409) .

It provides robust tests for testing in GLMs, by sign-flipping score
contributions. The tests are often robust against overdispersion,
heteroscedasticity and, in some cases, ignored nuisance variables.

## See also

Useful links:

- <https://livioivil.github.io/flipscores/>

## Author

**Maintainer**: Livio Finos <livio.finos@unipd.it>
([ORCID](https://orcid.org/0000-0003-3181-8078))

Authors:

- Livio Finos <livio.finos@unipd.it>
  ([ORCID](https://orcid.org/0000-0003-3181-8078))

- Angela Andreella

- Anna Vesely

Other contributors:

- Jelle J. Goeman \[contributor\]

- Jesse Hemerik \[contributor\]

- Riccardo De Santis \[contributor\]

- Andrea Panarotto \[contributor\]

- Filippo Gambarota \[contributor\]

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
#> flipscores(formula = Y ~ Z + X, family = "poisson(link='log')", 
#>     data = dt, x = TRUE)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept) -0.10259 -0.72287  2.71274 -0.26647 -0.08824 0.767
#> ZB          -0.15010 -0.71249  2.17890 -0.32699 -0.10415 0.660
#> ZC           0.16332  0.81062  2.22317  0.36462  0.11735 0.700
#> X            0.94385 16.20584  4.72718  3.42823  0.67119 0.012
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
#>   Ncoeff Mahalanobis      p
#> Z      2     0.75809 0.7048
#> X      1     0.02940 0.0120
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
#>                    Ncoeff Mahalanobis      p
#> Model 2 vs Model 1      1    0.029845 0.0084
# and
mod0=flipscores(Y~X,data=dt,family="poisson")
anova(mod0,mod)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach (5000 sign flips).
#> 
#> Model 1: Y ~ X
#> Model 2: Y ~ Z + X
#>                    Ncoeff Mahalanobis     p
#> Model 2 vs Model 1      2      1.4527 0.516
```
