<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to the dev-version of the

\#[library flipscores on
CRAN](http://cran.r-project.org/web/packages/flipscores/index.html)

------------------------------------------------------------------------

## Set up

To **install** this github version type (in R):

    ##if devtools is not installed yet: 
    ## install.packages("devtools") 
    library(devtools)
    install_github("livioivil/flipscores")

## Some examples

``` r
library(flipscores)
set.seed(1)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rpois(n=20,lambda=exp(dt$X))
mod=flipscores(Y~Z+X,data=dt,family="poisson",x=TRUE,)
summary(mod)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "poisson", data = dt, 
#>     x = TRUE)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.6910  -0.5792   0.1012   0.4900   1.0440  
#> 
#> Coefficients:
#>             Estimate   Score Std. Error z value eff_size Pr(>|z|)   
#> (Intercept)  -0.1026 -0.1616     0.4751 -0.3402   -0.088   0.7588   
#> ZB           -0.1501 -0.1593     0.3108 -0.5127   -0.104   0.6544   
#> ZC            0.1633  0.1813     0.3989  0.4544    0.117   0.6974   
#> X             0.9439  3.6238     1.4409  2.5149    0.671   0.0092 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 27.135  on 19  degrees of freedom
#> Residual deviance: 12.888  on 16  degrees of freedom
#> AIC: 57.459
#> 
#> Number of Fisher Scoring iterations: 5

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
#> Z  2 0.80075     0.6904   
#> X  1 0.02974     0.0092 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
#> Model 2 vs Model 1  1 0.029084      0.011 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
#>                    Df  Score Pr(>Score)
#> Model 2 vs Model 1  2 1.4679      0.506
```

### Negative Binomial

``` r
set.seed(1)
D=data.frame(x=(1:40)/20, z=rnorm(40))
D$y=rnbinom(40,mu=exp(D$x),size=3)

library(MASS)
mod_par=glm.nb(y~x+z,data=D, link="log")
summary(mod_par)
#> 
#> Call:
#> glm.nb(formula = y ~ x + z, data = D, link = "log", init.theta = 7.972747099)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -2.0746  -0.7748  -0.1086   0.4617   2.0435  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.15365    0.29358  -0.523    0.601    
#> x            0.92089    0.21481   4.287 1.81e-05 ***
#> z           -0.01282    0.13606  -0.094    0.925    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(7.9727) family taken to be 1)
#> 
#>     Null deviance: 61.217  on 39  degrees of freedom
#> Residual deviance: 41.312  on 37  degrees of freedom
#> AIC: 154.29
#> 
#> Number of Fisher Scoring iterations: 1
#> 
#> 
#>               Theta:  7.97 
#>           Std. Err.:  6.95 
#> 
#>  2 x log-likelihood:  -146.286
mod=flipscores(y~x+z, data=D, family = "negbinom") 
summary(mod)
#> 
#> Call:
#> flipscores(formula = y ~ x + z, family = "negbinom", data = D)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -2.0746  -0.7748  -0.1086   0.4617   2.0435  
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value eff_size Pr(>|t|)    
#> (Intercept) -0.15365 -0.29119    0.50440 -0.57729   -0.087   0.5712    
#> x            0.92089  2.36084    0.78171  3.02011    0.547   0.0008 ***
#> z           -0.01282 -0.11457    1.31945 -0.08683   -0.016   0.9588    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(7.9727) family taken to be 0.9960228)
#> 
#>     Null deviance: 61.217  on 39  degrees of freedom
#> Residual deviance: 41.312  on 37  degrees of freedom
#> AIC: 154.29
#> 
#> Number of Fisher Scoring iterations: 1
```

## References

J Hemerik, JJ Goeman and L Finos (2019) Robust testing in generalized
linear models by sign-flipping score contributions, Submitted, arXiv:
1909.03796

<https://arxiv.org/abs/1909.03796>

## Bug reports

If you encounter a bug, please file a
[reprex](https://github.com/tidyverse/reprex) (minimal reproducible
example) on [github](https://github.com/livioivil/flipscores/issues).
