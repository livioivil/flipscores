<!-- README.md is generated from README.Rmd. Please edit that file -->
Welcome to the dev-version of the
=================================

[library flip on CRAN](http://cran.r-project.org/web/packages/flip/index.html)
==============================================================================

------------------------------------------------------------------------

Set up
------

To **install** this github version type (in R):

    #if devtools is not installed yet: 
    # install.packages("devtools") 
    library(devtools)
    install_github("livioivil/flip")

------------------------------------------------------------------------

Some examples
-------------

``` r
library(flipscores)

set.seed(1)
x=(rep(0:1,20))
D=data.frame(y=rbinom(40,1,.25+x*.5),x=x,
             z=rnorm(40),id=rep(1:20,each=2))

mod_par=glm(y~x*z,data=D,family = binomial)
summary(mod_par)
#> 
#> Call:
#> glm(formula = y ~ x * z, family = binomial, data = D)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.9510  -0.7475   0.5615   0.7202   2.1981  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)  -0.7001     0.5195  -1.348  0.17778   
#> x             2.1557     0.7862   2.742  0.00611 **
#> z            -1.1320     0.8147  -1.389  0.16470   
#> x:z           1.5070     1.0329   1.459  0.14456   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 55.051  on 39  degrees of freedom
#> Residual deviance: 41.830  on 36  degrees of freedom
#> AIC: 49.83
#> 
#> Number of Fisher Scoring iterations: 4
mod_par
#> 
#> Call:  glm(formula = y ~ x * z, family = binomial, data = D)
#> 
#> Coefficients:
#> (Intercept)            x            z          x:z  
#>     -0.7001       2.1557      -1.1320       1.5070  
#> 
#> Degrees of Freedom: 39 Total (i.e. Null);  36 Residual
#> Null Deviance:       55.05 
#> Residual Deviance: 41.83     AIC: 49.83

mod=glm_flipscores(y~x*z,data=D,family = binomial,score_type = "ortho")
summary(mod)
#> 
#> Flip Score Test: 
#>   score_type = orthogonalized 
#>   n_flips= 1000
#> 
#> Call:
#> glm(formula = y ~ x * z, family = binomial, data = D, x = TRUE)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.9510  -0.7475   0.5615   0.7202   2.1981  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error z value Pr(>|z|)   
#> (Intercept) -0.70008 -0.06760    0.05194  -1.301    0.202   
#> x            2.15567  0.10335    0.04041   2.557    0.008 **
#> z           -1.13197 -0.05254    0.04660  -1.127    0.302   
#> x:z          1.50703  0.04013    0.03397   1.181    0.258   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 55.051  on 39  degrees of freedom
#> Residual deviance: 41.830  on 36  degrees of freedom
#> AIC: 49.83
#> 
#> Number of Fisher Scoring iterations: 4

print(mod)
#> 
#> Flip Score Test: 
#>   score_type = orthogonalized 
#>   n_flips= 1000 
#> Call: glm(formula = y ~ x * z, family = binomial, data = D, x = TRUE)
#> 
#> Coefficients:
#> (Intercept)           x           z         x:z 
#>  -0.7000829   2.1556675  -1.1319694   1.5070293


mod=glm_flipscores(y~x*z,data=D,family = binomial,score_type = "ortho",
                   id=D$id)
summary(mod)
#> 
#> Flip Score Test: 
#>   score_type = orthogonalized 
#>   n_flips= 1000
#> 
#> Call:
#> glm(formula = y ~ x * z, family = binomial, data = D, x = TRUE)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.9510  -0.7475   0.5615   0.7202   2.1981  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error z value Pr(>|z|)  
#> (Intercept) -0.70008 -0.13520    0.10374  -1.303    0.210  
#> x            2.15567  0.20670    0.08413   2.457    0.020 *
#> z           -1.13197 -0.10507    0.08889  -1.182    0.250  
#> x:z          1.50703  0.08026    0.06643   1.208    0.226  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 55.051  on 39  degrees of freedom
#> Residual deviance: 41.830  on 36  degrees of freedom
#> AIC: 49.83
#> 
#> Number of Fisher Scoring iterations: 4
```

------------------------------------------------------------------------

References
----------

Bug reports
-----------

If you encounter a bug, please file a [reprex](https://github.com/tidyverse/reprex) (minimal reproducible example) on [github](https://github.com/livioivil/flipscores/issues).
