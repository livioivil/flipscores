<!-- README.md is generated from README.Rmd. Please edit that file -->
Welcome to the dev-version of the
=================================

[library flipscores on CRAN](http://cran.r-project.org/web/packages/flipscores/index.html)
==========================================================================================

------------------------------------------------------------------------

Set up
------

To **install** this github version type (in R):

    ##if devtools is not installed yet: 
    ## install.packages("devtools") 
    library(devtools)
    install_github("livioivil/flipscores")

------------------------------------------------------------------------

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
mod=flipscores_glm(y~x*z,data=D,family = binomial,score_type = "ortho")
summary(mod)
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000
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
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000 
#> Call: glm(formula = y ~ x * z, family = binomial, data = D, x = TRUE)
#> 
#> Coefficients:
#> (Intercept)           x           z         x:z 
#>  -0.7000829   2.1556675  -1.1319694   1.5070293
mod=flipscores_glm(y~x*z,data=D,family = binomial,score_type = "ortho", id=D$id)
summary(mod)
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000
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

Poisson
-------

``` r
set.seed(2)
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
#> -11.4824   -6.6399    0.2325    4.7232    8.0694  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  4.78149    0.02190 218.297  < 2e-16 ***
#> x           -0.15845    0.03048  -5.198 2.01e-07 ***
#> z           -0.07039    0.01991  -3.535 0.000407 ***
#> x:z          0.09269    0.02679   3.460 0.000539 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 1546.2  on 39  degrees of freedom
#> Residual deviance: 1499.0  on 36  degrees of freedom
#> AIC: 1759.5
#> 
#> Number of Fisher Scoring iterations: 5
mod=flipscores_glm(y~x*z, data=D, family = poisson, score_type = "ortho")
summary(mod)
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000
#> 
#> Call:
#> glm(formula = y ~ x * z, family = poisson, data = D, x = TRUE)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -11.4824   -6.6399    0.2325    4.7232    8.0694  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error z value Pr(>|z|)   
#> (Intercept)  4.78149 47.22525   12.97956   3.638    0.002 **
#> x           -0.15845 -4.23970    4.50268  -0.942    0.442   
#> z           -0.07039 -4.46547    7.61931  -0.586    0.540   
#> x:z          0.09269  3.22026    1.80400   1.785    0.088 . 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 1546.2  on 39  degrees of freedom
#> Residual deviance: 1499.0  on 36  degrees of freedom
#> AIC: 1759.5
#> 
#> Number of Fisher Scoring iterations: 5
mod=flipscores_glm(y~x*z, data=D, family = poisson, score_type = "ortho", id=D$id)
summary(mod)
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000
#> 
#> Call:
#> glm(formula = y ~ x * z, family = poisson, data = D, x = TRUE)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -11.4824   -6.6399    0.2325    4.7232    8.0694  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error z value Pr(>|z|)   
#> (Intercept)  4.78149 94.45051   34.55977   2.733    0.006 **
#> x           -0.15845 -8.47940    8.77391  -0.966    0.398   
#> z           -0.07039 -8.93094   12.80239  -0.698    0.536   
#> x:z          0.09269  6.44052    3.43586   1.875    0.058 . 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 1546.2  on 39  degrees of freedom
#> Residual deviance: 1499.0  on 36  degrees of freedom
#> AIC: 1759.5
#> 
#> Number of Fisher Scoring iterations: 5
```

Linear model
------------

``` r
set.seed(3)
y=rnorm(40)
x=(rep(0:1,20))
D=data.frame(y=y, x=rbinom(40,1,.25+x*.5), z=rnorm(40),id=rep(1:20,each=2))

mod_par=glm(y~x*z,data=D, family = gaussian)
summary(mod_par)
#> 
#> Call:
#> glm(formula = y ~ x * z, family = gaussian, data = D)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -1.80932  -0.55385   0.01882   0.59155   1.44258  
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)  0.06305    0.20034   0.315    0.755
#> x           -0.10949    0.27058  -0.405    0.688
#> z           -0.40916    0.27361  -1.495    0.144
#> x:z          0.11649    0.36169   0.322    0.749
#> 
#> (Dispersion parameter for gaussian family taken to be 0.6747523)
#> 
#>     Null deviance: 27.063  on 39  degrees of freedom
#> Residual deviance: 24.291  on 36  degrees of freedom
#> AIC: 103.56
#> 
#> Number of Fisher Scoring iterations: 2
mod=flipscores_glm(y~x*z, data=D, family = gaussian, score_type = "ortho")
summary(mod)
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000
#> 
#> Call:
#> glm(formula = y ~ x * z, family = gaussian, data = D, x = TRUE)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -1.80932  -0.55385   0.01882   0.59155   1.44258  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error z value Pr(>|t|)  
#> (Intercept)  0.06305  0.02650    0.10146   0.261    0.798  
#> x           -0.10949 -0.02523    0.06916  -0.365    0.714  
#> z           -0.40916 -0.09220    0.04859  -1.898    0.080 .
#> x:z          0.11649  0.01502    0.03557   0.422    0.694  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.6747523)
#> 
#>     Null deviance: 27.063  on 39  degrees of freedom
#> Residual deviance: 24.291  on 36  degrees of freedom
#> AIC: 103.56
#> 
#> Number of Fisher Scoring iterations: 2
mod=flipscores_glm(y~x*z, data=D, family = gaussian, score_type = "ortho", id=D$id)
summary(mod)
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000
#> 
#> Call:
#> glm(formula = y ~ x * z, family = gaussian, data = D, x = TRUE)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -1.80932  -0.55385   0.01882   0.59155   1.44258  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error z value Pr(>|t|)  
#> (Intercept)  0.06305  0.05300    0.19313   0.274    0.778  
#> x           -0.10949 -0.05045    0.12576  -0.401    0.728  
#> z           -0.40916 -0.18440    0.10435  -1.767    0.094 .
#> x:z          0.11649  0.03004    0.06858   0.438    0.684  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.6747523)
#> 
#>     Null deviance: 27.063  on 39  degrees of freedom
#> Residual deviance: 24.291  on 36  degrees of freedom
#> AIC: 103.56
#> 
#> Number of Fisher Scoring iterations: 2
```

Gamma
-----

``` r
set.seed(4)
y=rgamma(40, shape=9, scale = 0.5)
x=(rep(0:1,20))
D=data.frame(y=y, x=rbinom(40,1,.25+x*.5), z=rnorm(40), x2=rpois(40, lambda = 4) ,id=rep(1:20,each=2))

mod_par=glm(y~x*z+x2,data=D, family = Gamma)
summary(mod_par)
#> 
#> Call:
#> glm(formula = y ~ x * z + x2, family = Gamma, data = D)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -0.54437  -0.17918  -0.04148   0.10429   0.61267  
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.203829   0.028924   7.047 3.32e-08 ***
#> x           -0.001971   0.021768  -0.091   0.9284    
#> z           -0.031979   0.016694  -1.916   0.0636 .  
#> x2           0.002689   0.005229   0.514   0.6103    
#> x:z          0.032547   0.023328   1.395   0.1717    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.09268343)
#> 
#>     Null deviance: 3.4693  on 39  degrees of freedom
#> Residual deviance: 3.1289  on 35  degrees of freedom
#> AIC: 144.21
#> 
#> Number of Fisher Scoring iterations: 5
mod=flipscores_glm(y~x*z+x2, data=D, family = Gamma, score_type = "ortho") #> senza +x2 non va?
summary(mod)
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000
#> 
#> Call:
#> glm(formula = y ~ x * z + x2, family = Gamma, data = D, x = TRUE)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -0.54437  -0.17918  -0.04148   0.10429   0.61267  
#> 
#> Coefficients:
#>              Estimate     Score Std. Error z value Pr(>|t|)   
#> (Intercept)  0.203829  1.021417   0.538284   1.898    0.002 **
#> x           -0.001971 -0.009645   0.051459  -0.187    0.838   
#> z           -0.031979 -0.260645   0.137029  -1.902    0.066 . 
#> x2           0.002689  0.233317   0.550261   0.424    0.718   
#> x:z          0.032547  0.142528   0.140503   1.014    0.392   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.09268343)
#> 
#>     Null deviance: 3.4693  on 39  degrees of freedom
#> Residual deviance: 3.1289  on 35  degrees of freedom
#> AIC: 144.21
#> 
#> Number of Fisher Scoring iterations: 5
mod=flipscores_glm(y~x*z+x2, data=D, family = Gamma, score_type = "ortho", id=D$id)
summary(mod)
#> Flip Score Test: 
#>          score_type = orthogonalized , n_flips= 1000
#> 
#> Call:
#> glm(formula = y ~ x * z + x2, family = Gamma, data = D, x = TRUE)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -0.54437  -0.17918  -0.04148   0.10429   0.61267  
#> 
#> Coefficients:
#>              Estimate     Score Std. Error z value Pr(>|t|)   
#> (Intercept)  0.203829  2.042835   1.239281   1.648    0.002 **
#> x           -0.001971 -0.019290   0.103982  -0.186    0.860   
#> z           -0.031979 -0.521289   0.332218  -1.569    0.176   
#> x2           0.002689  0.466634   0.991815   0.470    0.636   
#> x:z          0.032547  0.285055   0.278291   1.024    0.424   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.09268343)
#> 
#>     Null deviance: 3.4693  on 39  degrees of freedom
#> Residual deviance: 3.1289  on 35  degrees of freedom
#> AIC: 144.21
#> 
#> Number of Fisher Scoring iterations: 5
```

Negative Binomial
-----------------

``` r
set.seed(2)
y=(sample(1:200,20))
x=(rep(0:1,20))
D=data.frame(y=y, x=rbinom(40,1,.25+x*.5), z=rnorm(40),id=rep(1:20,each=2))

mod_par=MASS::glm.nb(y~x*z,data=D, link="log")
summary(mod_par)
#> 
#> Call:
#> MASS::glm.nb(formula = y ~ x * z, data = D, link = "log", init.theta = 2.450862935)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -2.26290  -1.14663   0.02337   0.64921   1.07088  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  4.78302    0.15342  31.176   <2e-16 ***
#> x           -0.16008    0.20662  -0.775    0.439    
#> z           -0.06200    0.14045  -0.441    0.659    
#> x:z          0.08336    0.18310   0.455    0.649    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(2.4509) family taken to be 1)
#> 
#>     Null deviance: 43.679  on 39  degrees of freedom
#> Residual deviance: 42.706  on 36  degrees of freedom
#> AIC: 452.34
#> 
#> Number of Fisher Scoring iterations: 1
#> 
#> 
#>               Theta:  2.451 
#>           Std. Err.:  0.533 
#> 
#>  2 x log-likelihood:  -442.339
#>mod=flipscores::flipscores_glm(y~x*z, data=D, family = "negbinom", score_type = "ortho") 
#>summary(mod)
#>mod=flipscores::flipscores_glm(y~x*z, data=D, family = "negbinom", score_type = "ortho", id=D$id)
#>summary(mod)
```

------------------------------------------------------------------------

References
----------

Bug reports
-----------

If you encounter a bug, please file a [reprex](https://github.com/tidyverse/reprex) (minimal reproducible example) on [github](https://github.com/livioivil/flipscores/issues).

