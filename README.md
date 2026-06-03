<!-- README.md is generated from README.Rmd. Please edit that file -->

# The flipscores package <a href="https://livioivil.github.io/flipscores/"><img src="man/figures/logo.png" align="right" height="139" alt="flipscores website" /></a>

[On CRAN
version](http://cran.r-project.org/web/packages/flipscores/index.html)

------------------------------------------------------------------------

## Set up

To **install** this github version type (in R):

    ##if devtools is not installed yet: 
    ## install.packages("devtools") 
    library(devtools)
    install_github("livioivil/flipscores")

## References

J Hemerik, JJ Goeman and L Finos (2019) Robust testing in generalized
linear models by sign-flipping score contributions. Journal of the Royal
Statistical Society Series B: Statistical Methodology, Volume 82, Issue
3, July 2020, Pages 841–864.  
<https://doi.org/10.1111/rssb.12369>

De Santis, R., Goeman, J. J., Hemerik, J., Davenport, S., & Finos, L.
(2025). Inference in Generalized Linear Models with Robustness to
Misspecified Variances. Journal of the American Statistical Association,
1–10. <https://doi.org/10.1080/01621459.2025.2491775>

De Santis, R., Goeman, J. J., Davenport, S., Hemerik, J., & Finos, L.
(2025) Permutation-based multiple testing when fitting many generalized
linear models Electronic Journal of Statistics 19 (2), 3317-3332
<https://doi.org/10.1214/25-EJS2409>

## Some examples

``` r
library(flipscores)
set.seed(1)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rpois(n=20,lambda=exp(dt$X))
mod=flipscores(Y~Z+X,data=dt,family="poisson",x=TRUE)
summary(mod)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "poisson(link='log')", 
#>     data = data, score_type = "standardized", n_flips = 5000, 
#>     alternative = "two.sided", to_be_tested = c("(Intercept)", 
#>     "ZB", "ZC", "X"), precompute_flips = TRUE, x = TRUE)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept) -0.10259 -0.72287  2.71274 -0.26647 -0.08824 0.755
#> ZB          -0.15010 -0.71249  2.17890 -0.32699 -0.10415 0.667
#> ZC           0.16332  0.81062  2.22317  0.36462  0.11735 0.678
#> X            0.94385 16.20584  4.72718  3.42823  0.67119 0.011
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
#>   Ncoeff Mahalanobis      p
#> Z      2     0.75546 0.6994
#> X      1     0.02941 0.0114
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
#> Model 2 vs Model 1      1     0.02898 0.0134
# and
mod0=flipscores(Y~X,data=dt,family="poisson")
anova(mod0,mod)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach ( sign flips).
#> 
#> Model 1: Y ~ X
#> Model 2: Y ~ Z + X
#>                    Ncoeff Mahalanobis     p
#> Model 2 vs Model 1      2      1.4212 0.522
```

### Custom contrasts

The `flipscores_contrasts()` helper defines factor contrasts with a
syntax similar to `emmeans`, or accepts a custom coefficient contrast
matrix in the spirit of `multcomp::glht()`, then applies the flip-score
test.

``` r
set.seed(11)
toy=data.frame(
  trt=factor(rep(c("A","B"), each=15))
)
toy$y=rbinom(30, 1, ifelse(toy$trt=="B", .65, .35))

toy_fit=glm(y~trt, data=toy, family=binomial, x=TRUE)

# Formula interface
flipscores_contrasts(toy_fit, pairwise ~ trt,
                     n_flips=500, seed=1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 500, alternative = two.sided
#> 
#>  contrast  estimate     Score p.value
#>     A - B -3.332205 -3.409809   0.002

# Dunnett-like all-versus-one comparisons
toy3=data.frame(
  trt=factor(rep(c("Control","Low","High"), each=12),
             levels=c("Control","Low","High"))
)
toy3$y=rbinom(nrow(toy3), 1,
              ifelse(toy3$trt=="High", .7,
                     ifelse(toy3$trt=="Low", .55, .35)))
toy3_fit=glm(y~trt, data=toy3, family=binomial, x=TRUE)
flipscores_contrasts(toy3_fit, dunnett ~ trt,
                     n_flips=500, seed=1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 500, alternative = two.sided
#> 
#>        contrast estimate    Score p.value
#>   Low - Control 1.609438 1.732051   0.092
#>  High - Control 1.945910 2.108185   0.042
flipscores_contrasts(toy3_fit, trt.vs.ctrl ~ trt, ref="Low",
                     n_flips=500, seed=1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 500, alternative = two.sided
#> 
#>       contrast   estimate      Score p.value
#>  Control - Low -1.6094379 -1.7320508   0.092
#>     High - Low  0.3364722  0.4096732   0.714
flipscores_contrasts(toy3_fit, trt.vs.ctrlk ~ trt,
                     n_flips=500, seed=1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 500, alternative = two.sided
#> 
#>        contrast   estimate      Score p.value
#>  Control - High -1.9459101 -2.1081851   0.042
#>      Low - High -0.3364722 -0.4096732   0.670

# If trt is involved in an interaction, the output includes a note unless
# the interaction partner is included in the contrast specification.
toy$sex=factor(rep(c("F","M"), length.out=nrow(toy)))
toy_int=glm(y~trt*sex, data=toy, family=binomial, x=TRUE)
flipscores_contrasts(toy_int, pairwise ~ trt,
                     n_flips=500, seed=1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 500, alternative = two.sided
#> 
#>  contrast  estimate     Score p.value
#>     A - B -19.07689 -2.561738   0.012
#> 
#> NOTE: Results may be misleading due to involvement in interactions: trt:sex
flipscores_contrasts(toy_int, pairwise ~ trt | sex,
                     n_flips=500, seed=1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 500, alternative = two.sided
#> 
#>   contrast   estimate     Score p.value
#>  A - B | F -19.076894 -2.561738   0.012
#>  A - B | M  -2.862201 -2.324174   0.034

# Custom coefficient contrast matrix
K=matrix(c(0,1), nrow=1)
colnames(K)=names(coef(toy_fit))
rownames(K)="B - A"
flipscores_contrasts(toy_fit, linfct=K,
                     n_flips=500, seed=1)
#> Flip-score custom contrasts
#> score_type = standardized, n_flips = 500, alternative = two.sided
#> 
#>  contrast estimate    Score p.value
#>     B - A 3.332205 3.409809   0.002
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
#> flipscores(formula = y ~ x + z, family = "Negative Binomial(7.9727)(link='log')", 
#>     data = data, score_type = "standardized", n_flips = 5000, 
#>     alternative = "two.sided", to_be_tested = c("(Intercept)", 
#>     "x", "z"), precompute_flips = TRUE)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept) -0.15365 -1.84162  3.42399 -0.53786 -0.08725 0.570
#> x            0.92089 14.93128  4.42610  3.37346  0.54725 0.001
#> z           -0.01282 -0.72457  7.24491 -0.10001 -0.01622 0.952
#> 
#> (Dispersion parameter for Negative Binomial(7.9727) family taken to be 0.9960228)
#> 
#>     Null deviance: 61.217  on 39  degrees of freedom
#> Residual deviance: 41.312  on 37  degrees of freedom
#> AIC: 154.29
#> 
#> Number of Fisher Scoring iterations: 1
```

## Bug reports

If you encounter a bug, please file a
[reprex](https://github.com/tidyverse/reprex) (minimal reproducible
example) on [github](https://github.com/livioivil/flipscores/issues).
