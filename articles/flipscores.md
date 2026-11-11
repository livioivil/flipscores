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
#> flipscores(formula = Y ~ Z + X, family = "poisson(link='log')", 
#>     data = dt, x = TRUE)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02207  0.26383  3.44133  0.07666  0.01595 0.937
#> ZB           0.02318  0.13383  2.39994  0.05576  0.01165 0.957
#> ZC           0.15946  1.36568  2.93906  0.46467  0.09651 0.645
#> X            0.89732 35.59669  6.63014  5.36892  0.72553 0.000
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
#>   Ncoeff Mahalanobis      p
#> Z      2    0.247415 0.8968
#> X      1    0.042301 0.0004
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
#>                    Ncoeff Mahalanobis     p
#> Model 2 vs Model 1      1    0.041869 2e-04
# and
mod0 <- flipscores(Y~X,data=dt,family="poisson",x=TRUE)
anova(mod0,mod)
#> Analysis of Deviance Table (Type III test)
#> 
#> Model: poisson, link: log
#> 
#> Inference is provided by FlipScores approach (5000 sign flips).
#> 
#> Model 1: Y ~ X
#> Model 2: Y ~ Z + X
#>                    Ncoeff Mahalanobis      p
#> Model 2 vs Model 1      2     0.28261 0.8828
```

Equivalently, one might apply the `flipscores` function to a precomputed
`glm` model.

``` r

model <- glm(Y~Z+X,data=dt,family="poisson")
mod2 <- flipscores(model)
summary(mod2)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, data = data.frame_30x3, family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02207  0.26383  3.44133  0.07666  0.01595 0.936
#> ZB           0.02318  0.13383  2.39994  0.05576  0.01165 0.956
#> ZC           0.15946  1.36568  2.93906  0.46467  0.09651 0.642
#> X            0.89732 35.59669  6.63014  5.36892  0.72553 0.000
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
#> flipscores(formula = Y ~ Z + X, flips = flps, data = data.frame_30x3, 
#>     family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02207  0.26383  3.44133  0.07666  0.01595 0.938
#> ZB           0.02318  0.13383  2.39994  0.05576  0.01165 0.955
#> ZC           0.15946  1.36568  2.93906  0.46467  0.09651 0.631
#> X            0.89732 35.59669  6.63014  5.36892  0.72553 0.001
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
#> flipscores(formula = Y ~ Z + X, flips = flps, data = data.frame_30x3, 
#>     family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02207  0.26383  3.44133  0.07666  0.01595 0.938
#> ZB           0.02318  0.13383  2.39994  0.05576  0.01165 0.955
#> ZC           0.15946  1.36568  2.93906  0.46467  0.09651 0.631
#> X            0.89732 35.59669  6.63014  5.36892  0.72553 0.001
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
#> flipscores(formula = Y ~ Z + X, alternative = "less", data = data.frame_30x3, 
#>     family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02207  0.26383  3.44133  0.07666  0.01595 0.930
#> ZB           0.02318  0.13383  2.39994  0.05576  0.01165 0.959
#> ZC           0.15946  1.36568  2.93906  0.46467  0.09651 0.641
#> X            0.89732 35.59669  6.63014  5.36892  0.72553 0.000
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
#> flipscores(formula = Y ~ Z + X, alternative = "greater", data = data.frame_30x3, 
#>     family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02207  0.26383  3.44133  0.07666  0.01595 0.936
#> ZB           0.02318  0.13383  2.39994  0.05576  0.01165 0.959
#> ZC           0.15946  1.36568  2.93906  0.46467  0.09651 0.665
#> X            0.89732 35.59669  6.63014  5.36892  0.72553 0.000
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
#> flipscores(formula = Y ~ Z + X, to_be_tested = "X", flips = flps, 
#>     data = data.frame_30x3, family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02207       NA       NA       NA       NA    NA
#> ZB           0.02318       NA       NA       NA       NA    NA
#> ZC           0.15946       NA       NA       NA       NA    NA
#> X            0.89732 35.59669  6.63014  5.36892  0.72553 0.001
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
#> flipscores(formula = Y ~ Z + X, to_be_tested = c(1, 4), flips = flps, 
#>     data = data.frame_30x3, family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02207  0.26383  3.44133  0.07666  0.01595 0.938
#> ZB           0.02318       NA       NA       NA       NA    NA
#> ZC           0.15946       NA       NA       NA       NA    NA
#> X            0.89732 35.59669  6.63014  5.36892  0.72553 0.001
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
#> flipscores(formula = Y ~ Z + X, to_be_tested = "ZB", flips = flps, 
#>     data = data.frame_30x3, family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate   score      se       z    pcor     p
#> (Intercept)  0.02207      NA      NA      NA      NA    NA
#> ZB           0.02318 0.13383 2.39994 0.05576 0.01165 0.955
#> ZC           0.15946      NA      NA      NA      NA    NA
#> X            0.89732      NA      NA      NA      NA    NA
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
#> flipscores(formula = Y ~ Z + X, to_be_tested = 2, flips = flps, 
#>     data = data.frame_30x3, family = "poisson(link='log')")
#> 
#> Coefficients:
#>             estimate   score      se       z    pcor     p
#> (Intercept)  0.02207      NA      NA      NA      NA    NA
#> ZB           0.02318 0.13383 2.39994 0.05576 0.01165 0.955
#> ZC           0.15946      NA      NA      NA      NA    NA
#> X            0.89732      NA      NA      NA      NA    NA
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
#> flipscores(formula = Y ~ Z + X, family = "Negative Binomial(20776.65)(link='log')", 
#>     data = dt)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.02208  0.26385  3.18354  0.08288  0.01595 0.935
#> ZB           0.02316  0.13368  2.20995  0.06049  0.01164 0.959
#> ZC           0.15946  1.36541  2.72310  0.50142  0.09650 0.645
#> X            0.89733 17.36132  4.68164  3.70839  0.71368 0.000
#> 
#> (Dispersion parameter for Negative Binomial(20776.65) family taken to be 0.8801237)
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
#> Warning in FUN(X[[i]], ...): A confidence bound was not found in the search
#> interval -- returning Inf
#>                  2.5 %    97.5 %
#> (Intercept) -1.0782344 0.6344857
#> ZB                -Inf 0.9538402
#> ZC          -0.8221180 1.0590446
#> X            0.5996654 1.4233455
```

As before, various options are allowed, including the possibility to
generate symmetric confidence intervals.

``` r

confint(mod, parm = 4, flips = flps) # parm == to_be_tested!
#>       2.5 %   97.5 %
#> X 0.5780193 1.446592
confint(mod, parm = 4, level = 0.99, flips = flps)
#>       0.5 %   99.5 %
#> X 0.4604931 1.766795
confint(mod, parm = 4, alternative = "greater", flips = flps)
#>        Est.    95 %
#> X 0.8973217 1.31611
confint(mod, parm = 4, type = "symmetric", flips = flps)
#>       2.5 %   97.5 %
#> X 0.4586757 1.335968
```

## Other examples

``` r

set.seed(1)

dat <- data.frame(
  treatment = factor(rep(c("Control", "Drug"), each = 30)),
  sex = factor(rep(c("F", "M"), times = 30)),
  age_group = factor(rep(c("Young", "Old", "Middle"), length.out = 60)),
  severity = runif(60, 0, 1)
)

eta <- with(dat,
  -1 + 1.2 * (treatment == "Drug") * (sex == "M") +
    0.6 * (age_group == "Old") - 0.3 * severity
)
dat$y <- rbinom(nrow(dat), 1, plogis(eta))


flips=make_flips(60,5000)
fit1 <- glm(
  y ~ treatment * sex + age_group + severity,
  data = dat,
  family = binomial
)

# Formula interface: build contrasts and apply flip-score tests.
cntr1=flipscores_contrasts(
  fit1,
  pairwise ~ treatment,
  flips = flips,
  seed = 1
)

fit2 <- flipscores(
  y ~ treatment + sex * age_group + severity,
  data = dat,
  family = binomial,
  flips = flips
)
summary(fit2)
#> 
#> Call:
#> flipscores(formula = y ~ treatment + sex * age_group + severity, 
#>     family = "binomial(link='logit')", data = dat, flips = flips)
#> 
#> Coefficients:
#>                     estimate    score       se        z     pcor     p
#> (Intercept)         -1.47850 -1.67124  1.11276 -1.50189 -0.19460 0.145
#> treatmentDrug        0.60981  2.01506  1.82720  1.10281  0.14221 0.299
#> sexM                 1.39795  1.47267  1.05204  1.39982  0.18139 0.164
#> age_groupOld         0.55042  0.50052  0.95861  0.52212  0.06747 0.624
#> age_groupYoung       1.00299  1.00111  1.01392  0.98736  0.12736 0.366
#> severity            -0.45260 -0.41470  0.95782 -0.43297 -0.05617 0.669
#> sexM:age_groupOld   -0.14754 -0.07530  0.71503 -0.10531 -0.01368 0.916
#> sexM:age_groupYoung -1.42334 -0.77349  0.74401 -1.03962 -0.13456 0.332
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 80.761  on 59  degrees of freedom
#> Residual deviance: 75.082  on 52  degrees of freedom
#> AIC: 91.082
#> 
#> Number of Fisher Scoring iterations: 4
# Formula interface: build contrasts and apply flip-score tests.
cntr2=flipscores_contrasts(
  fit2,
  pairwise ~ treatment,
  flips = flips,
  seed = 1
)

cmbn=combine_contrasts(fit2)
summary(cmbn) #perchè non coimbina tra livelli di group??
#>      model         coefficient        stat ntests        S      p
#> 1 model1.0         (Intercept) Mahalanobis      1 0.020974 0.1446
#> 2 model8.1       treatmentDrug Mahalanobis      1 0.014919 0.2992
#> 3 model5.2                sexM Mahalanobis      1 0.020069 0.1642
#> 4 model2.3        age_groupOld Mahalanobis      1 0.007012 0.6244
#> 5 model3.3      age_groupYoung Mahalanobis      1 0.012951 0.3662
#> 6 model4.4            severity Mahalanobis      1 0.006092 0.6694
#> 7 model6.5   sexM:age_groupOld Mahalanobis      1 0.001423 0.9164
#> 8 model7.5 sexM:age_groupYoung Mahalanobis      1 0.013818 0.3322


# jfs è la vecchia classe join_flipscores

flipscores(list(cntr1,cntr2,fit2))
#> == Joining n = 3 objects
#> 
#>   model .assign response    coefficient   estimate      score      p se  z pcor
#> 1  mod1      NA        y Control - Drug  0.3173374  0.4928128 0.6998 NA NA   NA
#> 2  mod2      NA        y Control - Drug -0.6098113 -2.0150610 0.3006 NA NA   NA
#> ... 6 rows ...
#>    model .assign response         coefficient   estimate       score      p
#> 9   mod3       5        y   sexM:age_groupOld -0.1475364 -0.07529954 0.9202
#> 10  mod3       5        y sexM:age_groupYoung -1.4233363 -0.77348703 0.3320
#>           se          z        pcor
#> 9  0.7150335 -0.1053091 -0.01368444
#> 10 0.7440070 -1.0396233 -0.13456011

cres=flipscores:::c.jfs(cntr1,cntr2,fit2,cmbn)
#deve funzionare anche per flipscores, fs_contrasts, fs_cobine,jfs
#forse meglio chiamarla join_fs() o fualcosa del genere?
#rivedere il nome degli objects: mode1 (evitare duplicati, cambia nome in caso) e credo anche .assign della summary_table

#aggiungere pcor nel summary dei contrasti?

#Le classi fs_contrast (e clip, flip2sss etc) devono avere i seguenti methods:
#model.matrix per gli obs_names
#update indispensabile per refit con stessa matrice flip
#as.jfs  convertire ciascuno così poi il c() è più facile. c'è già un default che dovrebbe coprire quasi tutti i casi: basta che nell'oggetto ci sia una $call, una $summary_table e un $Tspace, il resto viene cacciato in $objects)


#fs_contrast (e clip, flip2sss etc devono avere gli argomenti n_flips (forse questo non serve in realtà) e flips in modo che poi l'update(flips=FLIP) funzioni


# fs_combine è in pratica un offetto jfs (che non serva nemmeno farne un metodo a parte?)


#IMPORTANTE: c() o come lo chiamiamo non controlla la coerenza dei flip (altrimenti fs_combine non può essere agganciato). il controllo lo fa flipscores() - usando quindi model.matrix e update. 



combine_contrasts(fit2) #non combine per livelli di age_group perchè non è un oggetto joint_flipscores
#> == Joining n = 8 objects
#> 
#>      model   coefficient        stat ntests          S      p
#> 1 model1.0   (Intercept) Mahalanobis      1 0.02097440 0.1446
#> 2 model8.1 treatmentDrug Mahalanobis      1 0.01491949 0.2992
#> ... 4 rows ...
#>      model         coefficient        stat ntests           S      p
#> 7 model6.5   sexM:age_groupOld Mahalanobis      1 0.001423184 0.9164
#> 8 model7.5 sexM:age_groupYoung Mahalanobis      1 0.013817997 0.3322

combine_contrasts(cres) #invece funziona
#> == Joining n = 15 objects
#> 
#>        model    coefficient        stat ntests         S      p
#> 1 object1.NA Control - Drug Mahalanobis      2 3.0410563 0.2204
#> 2 object1. 0    (Intercept) Mahalanobis      1 0.0209744 0.1446
#> ... 11 rows ...
#>          model         coefficient        stat ntests           S      p
#> 14 model6.5.NA   sexM:age_groupOld Mahalanobis      1 0.001423184 0.9164
#> 15 model7.5.NA sexM:age_groupYoung Mahalanobis      1 0.013817997 0.3322

combine_tests(cres)
#> == Joining n = 1 objects
#> 
#>     model coefficient stat ntests        S     p
#> 1 overall        many maxT     18 1.501888 0.666
combine_tests(cres,by="coefficient",comb_funct = "Fisher")
#> == Joining n = 9 objects
#> 
#>            model    coefficient   stat ntests        S      p
#> 1 Control - Drug Control - Drug Fisher      2 1.558757 0.4882
#> 2    (Intercept)    (Intercept) Fisher      2 3.867568 0.1446
#> ... 5 rows ...
#>                 model         coefficient   stat ntests         S      p
#> 8   sexM:age_groupOld   sexM:age_groupOld Fisher      2 0.1746047 0.9164
#> 9 sexM:age_groupYoung sexM:age_groupYoung Fisher      2 2.2040362 0.3322

#ora construire la selezione [[]] usa .assign
length(cres$objects)
#> [1] 4

###############
set.seed(1)
n <- 50
df <- data.frame(
  y1  = rnorm(n),
  y2  = rnorm(n),
  yb  = rbinom(n, 1, 0.5),
  y1p = rpois(n, lambda = 3),
  y2p = rpois(n, lambda = 5),
  x1  = rnorm(n),
  x2  = rnorm(n)
)
#'
# Case 1: single formula
res=flipscores(y1 ~ x1 + x2, family = gaussian(), data = df)
summary(res)
#> 
#> Call:
#> flipscores(formula = y1 ~ x1 + x2, family = "gaussian(link='identity')", 
#>     data = df)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept)  0.10358  5.15889  5.96233  0.86525  0.12489 0.393
#> x1          -0.04491 -2.41629  6.15867 -0.39234 -0.05663 0.699
#> x2           0.02313  1.20293  6.04796  0.19890  0.02871 0.819
#> 
#> (Dispersion parameter for gaussian family taken to be 0.7175968)
#> 
#>     Null deviance: 33.870  on 49  degrees of freedom
#> Residual deviance: 33.727  on 47  degrees of freedom
#> AIC: 130.21
#> 
#> Number of Fisher Scoring iterations: 2
#'
# Case 2: list of glm objects

res_lm=flipscores_lm(y1 ~ x1 + x2, data = df, tested_coeffs = "x2")
summary(res_lm)
#> Flip-score Linear Model
#> Call: flipscores_lm(formula = y1 ~ x1 + x2, data = df, tested_coeffs = "x2")
#> score_type = Standardized, n_flips = , alternative = two.sided
#> 
#>  response coefficient score p.values
#>        y1          x2 1.203   0.8272

models <- list(
  glm(y1 ~ x1 + x2, data = df, family = gaussian()),
  cntr1,
  cntr2,
  fit2,
  res_lm
)

res=flipscores(models, family = gaussian(), data = df)

summary(res)
#>    model response         coefficient estimate   score     se       z     pcor
#> 1   mod1       y1         (Intercept)  0.10358  5.1589 5.9623  0.8652  0.12489
#> 2   mod1       y1                  x1 -0.04491 -2.4163 6.1587 -0.3923 -0.05663
#> 3   mod1       y1                  x2  0.02313  1.2029 6.0480  0.1989  0.02871
#> 4   mod2        y      Control - Drug  0.31734  0.4928     NA      NA       NA
#> 5   mod3        y      Control - Drug -0.60981 -2.0151     NA      NA       NA
#> 6   mod4        y         (Intercept) -1.47850 -1.6712 1.1128 -1.5019 -0.19460
#> 7   mod4        y       treatmentDrug  0.60981  2.0151 1.8272  1.1028  0.14221
#> 8   mod4        y                sexM  1.39795  1.4727 1.0520  1.3998  0.18139
#> 9   mod4        y        age_groupOld  0.55042  0.5005 0.9586  0.5221  0.06747
#> 10  mod4        y      age_groupYoung  1.00299  1.0011 1.0139  0.9874  0.12736
#> 11  mod4        y            severity -0.45260 -0.4147 0.9578 -0.4330 -0.05617
#> 12  mod4        y   sexM:age_groupOld -0.14754 -0.0753 0.7150 -0.1053 -0.01368
#> 13  mod4        y sexM:age_groupYoung -1.42334 -0.7735 0.7440 -1.0396 -0.13456
#> 14  mod5       y1                <NA>       NA  1.2029     NA      NA       NA
#>         p p.values
#> 1  0.3852       NA
#> 2  0.7062       NA
#> 3  0.8228       NA
#> 4  0.7170       NA
#> 5  0.3020       NA
#> 6  0.1332       NA
#> 7  0.3020       NA
#> 8  0.1550       NA
#> 9  0.6206       NA
#> 10 0.3524       NA
#> 11 0.6696       NA
#> 12 0.9202       NA
#> 13 0.3274       NA
#> 14     NA   0.8228


flipscores_lm(y1 ~ x1 + x2, data = df)
#> two.sided
```
