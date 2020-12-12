<!-- README.md is generated from README.Rmd. Please edit that file -->

Welcome to the dev-version of the
=================================

\#[library flipscores on
CRAN](http://cran.r-project.org/web/packages/flipscores/index.html)

------------------------------------------------------------------------

Set up
------

To **install** this github version type (in R):

    ##if devtools is not installed yet: 
    ## install.packages("devtools") 
    library(devtools)
    install_github("livioivil/flipscores")

Some examples
-------------

Binomial
--------

``` r
library(flipscores)
set.seed(1)
x=(rep(0:1,20))
D=data.frame(y=rbinom(40,1,.25+x*.5), x=x,z=rnorm(40),id=rep(1:20,each=2))
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
mod=flipscores(y~x*z,data=D,family = binomial, score_type = "effective")
summary(mod)
#> 
#> Call:
#> flipscores(formula = y ~ x * z, family = binomial, data = D, 
#>     score_type = "effective")
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.9510  -0.7475   0.5615   0.7202   2.1981  
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value Pr(>|z|)   
#> (Intercept)  -0.7001 -2.7040     1.9441  -1.391    0.164   
#> x             2.1557  4.1340     1.3697   3.018    0.004 **
#> z            -1.1320 -2.1015     1.6562  -1.269    0.238   
#> x:z           1.5070  1.6052     1.2372   1.297    0.238   
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
#> Flip Score Test: 
#>          score_type = effective , n_flips = 5000 
#> Call: flipscores(formula = y ~ x * z, family = binomial, data = D, 
#>     score_type = "effective")
#> 
#> Coefficients:
#> (Intercept)           x           z         x:z 
#>  -0.7000829   2.1556675  -1.1319694   1.5070293

# with clustered observations:
mod=flipscores(y~x*z,data=D,family = binomial, id=D$id, score_type = "ort")
#> Warning in print(warning("WARNING: Use of id is allowed only with score_type=='effective', yet. \n Nothoing done.")): WARNING: Use of id is allowed only with score_type=='effective', yet. 
#>  Nothoing done.
#> [1] "WARNING: Use of id is allowed only with score_type=='effective', yet. \n Nothoing done."
summary(mod)
#> Length  Class   Mode 
#>      0   NULL   NULL
```

Poisson
-------

``` r
set.seed(1)
y=(sample(1:200,20))
x=(rep(0:1,20))
D=data.frame(y=y, x=rbinom(40,1,.25+x*.5), z=rnorm(40),id=rep(1:20,each=2))
mod_par=glm(y~x*z,data=D, family = poisson)
summary(mod_par)
#> 
#> Call:
#> glm(formula = y ~ x * z, family = poisson, data = D)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -12.1384   -5.2048   -0.7249    3.4922   11.3918  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  4.44058    0.02393 185.571  < 2e-16 ***
#> x            0.12091    0.03407   3.549 0.000386 ***
#> z           -0.31402    0.02715 -11.568  < 2e-16 ***
#> x:z          0.45883    0.03595  12.762  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 1468.7  on 39  degrees of freedom
#> Residual deviance: 1289.4  on 36  degrees of freedom
#> AIC: 1541.4
#> 
#> Number of Fisher Scoring iterations: 5
mod=flipscores(y~x*z, data=D, family = poisson, score_type = "effective")
summary(mod)
#> 
#> Call:
#> flipscores(formula = y ~ x * z, family = poisson, data = D, score_type = "effective")
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -12.1384   -5.2048   -0.7249    3.4922   11.3918  
#> 
#> Coefficients:
#>              Estimate     Score Std. Error z value Pr(>|z|)   
#> (Intercept)    4.4406 1532.4218   375.5455   4.081    0.002 **
#> x              0.1209  104.0738   172.8491   0.602    0.538   
#> z             -0.3140 -440.9885   255.9402  -1.723    0.092 . 
#> x:z            0.4588  357.5901   159.8972   2.236    0.048 * 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 1468.7  on 39  degrees of freedom
#> Residual deviance: 1289.4  on 36  degrees of freedom
#> AIC: 1541.4
#> 
#> Number of Fisher Scoring iterations: 5
```

Linear model
------------

``` r
set.seed(1)
D$y=rnorm(40)
mod=flipscores(y~x*z, data=D, family = gaussian, score_type = "effective")
summary(mod)
#> 
#> Call:
#> flipscores(formula = y ~ x * z, family = gaussian, data = D, 
#>     score_type = "effective")
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -1.82822  -0.45485   0.05528   0.48941   1.44802  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error z value Pr(>|t|)   
#> (Intercept)  0.03272  0.71313    3.82746   0.186    0.874   
#> x            0.06963  0.68635    2.18839   0.314    0.750   
#> z           -0.50382 -8.06234    3.30916  -2.436    0.010 **
#> x:z          1.11710  9.46884    2.80226   3.379    0.002 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.5501491)
#> 
#>     Null deviance: 30.661  on 39  degrees of freedom
#> Residual deviance: 19.805  on 36  degrees of freedom
#> AIC: 95.398
#> 
#> Number of Fisher Scoring iterations: 2
```

Negative Binomial
-----------------

``` r
set.seed(1)
D=data.frame(x=rnorm(40), z=rnorm(40))
D$y=rnbinom(40,mu=exp(-D$x),size=3)

library(MASS)
#> Warning: package 'MASS' was built under R version 3.6.3
mod_par=glm.nb(y~x+z,data=D, link="log")
#> Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
#> control$trace > : iteration limit reached
#> Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
#> control$trace > : iteration limit reached
summary(mod_par)
#> 
#> Call:
#> glm.nb(formula = y ~ x + z, data = D, link = "log", init.theta = 5040.070663)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.6713  -1.1157  -0.7917   0.5558   2.5341  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.3371     0.1969  -1.712 0.086934 .  
#> x            -0.6175     0.1831  -3.373 0.000743 ***
#> z             0.1428     0.1920   0.744 0.456954    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(5040.071) family taken to be 1)
#> 
#>     Null deviance: 57.275  on 39  degrees of freedom
#> Residual deviance: 46.578  on 37  degrees of freedom
#> AIC: 98.234
#> 
#> Number of Fisher Scoring iterations: 1
#> 
#> 
#>               Theta:  5040 
#>           Std. Err.:  161311 
#> Warning while fitting theta: iteration limit reached 
#> 
#>  2 x log-likelihood:  -90.234
mod=flipscores(y~x+z, data=D, family = "negbinom",score_type = "effective") 
#> Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
#> control$trace > : iteration limit reached

#> Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
#> control$trace > : iteration limit reached
summary(mod)
#> 
#> Call:
#> flipscores(formula = y ~ x + z, family = "negbinom", data = D, 
#>     score_type = "effective")
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.6713  -1.1157  -0.7917   0.5558   2.5341  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error z value Pr(>|t|)  
#> (Intercept)  -0.3371  -8.9318     5.1061  -1.749    0.066 .
#> x            -0.6175 -11.0407     4.8196  -2.291    0.024 *
#> z             0.1428   3.4697     5.5809   0.622    0.556  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(5040.071) family taken to be 1.301545)
#> 
#>     Null deviance: 57.275  on 39  degrees of freedom
#> Residual deviance: 46.578  on 37  degrees of freedom
#> AIC: 98.234
#> 
#> Number of Fisher Scoring iterations: 1
```

References
----------

J Hemerik, JJ Goeman and L Finos (2019) Robust testing in generalized
linear models by sign-flipping score contributions, Submitted, arXiv:
1909.03796

<a href="https://arxiv.org/abs/1909.03796" class="uri">https://arxiv.org/abs/1909.03796</a>

Bug reports
-----------

If you encounter a bug, please file a
[reprex](https://github.com/tidyverse/reprex) (minimal reproducible
example) on [github](https://github.com/livioivil/flipscores/issues).

Bug reports
-----------

If you encounter a bug, please file a
[reprex](https://github.com/tidyverse/reprex) (minimal reproducible
example) on [github](https://github.com/livioivil/flipscores/issues).
