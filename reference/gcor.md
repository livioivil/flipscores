# Compute Generalized Partial Correlations for GLM terms

This function computes the generalized partial correlation coefficient
\\r\\ for each specified variable in a generalized linear model. For
each variable, it refits the null model excluding that variable and
computes the cosine similarity between the residualized predictor and
standardized residuals.

## Usage

``` r
gcor(full_glm, terms = NULL, normalize = FALSE, intercept_too = FALSE)
```

## Arguments

- full_glm:

  A fitted GLM object of class \`glm\`.

- terms:

  Character vector of variable names (referred to the model.matrix, that
  is pay attention for factors) for which to compute generalized partial
  correlations. If \`NULL\` (default), computes for all non-intercept
  terms in the model.

- normalize:

  FALSE by default.

- intercept_too:

  Logical indicating whether to include the intercept as a variable.
  Default is FALSE.

## Value

if `normalize` is `FALSE`, a data frame with five columns:

- variable:

  The variable name

- r:

  The generalized partial correlation coefficient

while if `normalize` is `TRUE`

- terms:

  The variable name

- r:

  The generalized partial correlation coefficient

- r_n:

  The normalized generalized partial correlation coefficient

- null_model:

  The null model used to compute the generalized (partial) correlation

- algorithm:

  The algorithm used to compute the upper/lower bounds of the
  generalized partial correlation coefficient (to compute its normalized
  version)

- exact:

  logical

## Details

The generalized partial correlation \\\r\\ measures the association
between a predictor and response after adjusting for all other terms in
the model. It is defined as the cosine similarity between the
residualized predictor \\X_r\\ and standardized residuals \\Y_r\\:

\$\$\r = \frac{X_r^\top Y_r}{\\X_r\\ \\Y_r\\}\$\$

where:

- \\X_r = (I - H)W^{1/2}X\\ is the residualized predictor

- \\Y_r = V^{-1/2}(Y - \hat{\mu})\\ is the standardized residual vector

- \\H\\ is the hat matrix for the nuisance covariates

- \\W = DV^{-1}D\\ is the weight matrix

- \\V\\ is the variance matrix and \\D\\ is the derivative matrix

The function uses \`flipscores:::get_par_expo_fam()\` to compute \\V\\
and \\D\\ consistently with the flipscores package methodology.

The normalized generalized partial correlation is computed as: \$\$ r_n
= \begin{cases} +r / r\_+ & \text{if } r \> 0 \\ -r / r\_- & \text{if }
r \< 0 \end{cases} \$\$ where \\r\_+\\ is the maximum possible
correlation and \\r\_-\\ is the minimum.

When the (full) model has only intercept and only one predictor X, the
generalized (non partial) correlation is computed and the normalization
factor for X is exact. In the more general case with more predictors,
for sample sizes \\n \leq n\_{\text{exact}}\\, brute force search is
used to find the exact extrema. For larger sample sizes, a greedy
multi-start algorithm is employed:

- Multiple starting points are generated using thresholding and random
  sampling

- From each start, coordinates are greedily flipped to improve the
  correlation

- Early stopping is used when no improvements are found for several
  iterations

- The best solution across all starts is returned

This approach provides a good trade-off between computational efficiency
and solution quality for large problems where brute force is infeasible.

## Author

Livio Finos and Paolo Girardi

## Examples

``` r
set.seed(1)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rpois(n=20,lambda=exp(dt$Z=="C"))
mod=flipscores(Y~Z+X,data=dt,family="poisson",n_flips=1000)
summary(mod)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "poisson(link='log')", 
#>     data = dt, n_flips = 1000)
#> 
#> Coefficients:
#>             estimate    score       se        z     pcor     p
#> (Intercept) -0.14256 -0.91360  2.62144 -0.34851 -0.12658 0.734
#> ZB          -0.18558 -0.50868  1.65785 -0.30683 -0.10810 0.672
#> ZC           1.40981  8.55380  2.58950  3.30326  0.76519 0.001
#> X           -0.06964 -1.56935  4.70999 -0.33320 -0.11745 0.655
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 28.649  on 19  degrees of freedom
#> Residual deviance: 11.218  on 16  degrees of freedom
#> AIC: 58.102
#> 
#> Number of Fisher Scoring iterations: 5
#> 

# Compute generalized partial correlations for all terms
(results <- gcor(mod))
#>   terms       gcor null_model
#> 1   ~ZB -0.1080979    ~1+ZC+X
#> 2   ~ZC  0.7651942    ~1+ZB+X
#> 3    ~X -0.1174534   ~1+ZB+ZC

# Compute for specific terms only
gcor(mod, terms = c("X", "ZC"))
#>   terms       gcor null_model
#> 1    ~X -0.1174534   ~1+ZB+ZC
#> 2   ~ZC  0.7651942    ~1+ZB+X

gcor(mod, terms = c("X", "ZC"),normalize=TRUE)
#>    terms          r        r_n null_model
#> X     ~X -0.1174534 -0.2755351   ~1+ZB+ZC
#> ZC   ~ZC  0.7651942  1.0442471    ~1+ZB+X


gcor(mod, intercept_too=TRUE, normalize=TRUE)
#>                    terms          r        r_n null_model
#> (Intercept) ~(Intercept) -0.1265802 -0.2077374 ~1+ZB+ZC+X
#> ZB                   ~ZB -0.1080979 -0.1532418    ~1+ZC+X
#> ZC                   ~ZC  0.7651942  1.0442471    ~1+ZB+X
#> X                     ~X -0.1174534 -0.2755351   ~1+ZB+ZC
set.seed(123)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rbinom(n=20,prob=plogis((dt$Z=="C")*2),size=1)
mod=flipscores(Y~Z+X,data=dt,family="binomial",n_flips=1000)
summary(mod)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "binomial(link='logit')", 
#>     data = dt, n_flips = 1000)
#> 
#> Coefficients:
#>              estimate     score        se         z      pcor     p
#> (Intercept)  -0.14858  -0.21025   1.18815  -0.17696  -0.06718 0.870
#> ZB          -20.45393  -1.47837   0.74659  -1.98016  -0.52957 0.059
#> ZC           20.85614   1.80427   0.81801   2.20568   0.61463 0.033
#> X            -0.42757  -0.37825   0.95742  -0.39507  -0.14932 0.740
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 27.5256  on 19  degrees of freedom
#> Residual deviance:  9.4015  on 16  degrees of freedom
#> AIC: 17.401
#> 
#> Number of Fisher Scoring iterations: 19
#> 

(results <- gcor(mod,normalize=TRUE))
#>    terms          r        r_n null_model
#> ZB   ~ZB -0.5295718        NaN    ~1+ZC+X
#> ZC   ~ZC  0.6146304        NaN    ~1+ZB+X
#> X     ~X -0.1493213 -0.3669599   ~1+ZB+ZC
# Compute for specific terms only
gcor(mod, terms = c("X", "ZC"),normalize=TRUE)
#>    terms          r        r_n null_model
#> X     ~X -0.1493213 -0.3669599   ~1+ZB+ZC
#> ZC   ~ZC  0.6146304        NaN    ~1+ZB+X


############################
set.seed(1)
Z=rnorm(20)
X=Z+rnorm(20)-10
Y=rpois(n=20,lambda=exp(Z+.25*X))
model1=glm(Y~Z+X,family="poisson")
gcor(model1)
#>   terms       gcor null_model
#> 1    ~Z -0.2379952       ~1+X
#> 2    ~X  0.3139891       ~1+Z
gcor(model1, normalize=TRUE)
#>   terms          r        r_n null_model
#> Z    ~Z -0.2379952 -0.4708683       ~1+X
#> X    ~X  0.3139891  0.4983668       ~1+Z
#normalized_gcor_poisson(model1, "X", verbose = TRUE)[1:3]
model0=glm(Y~Z,family="poisson")
if (FALSE) sapply(c('v','kl','sse','lr','n'),function(type) rsq.partial(model1,model0,type=type)$partial.rsq) # \dontrun{}

set.seed(1)
Z=rnorm(20)
X=Z+rnorm(20)
Y=rbinom(n=20,1,prob=plogis((-.5*Z+1*X)))
model1=glm(Y~Z+X,family="binomial")
gcor(model1)
#>   terms       gcor null_model
#> 1    ~Z -0.1852201       ~1+X
#> 2    ~X  0.3229866       ~1+Z
gcor(model1, normalize=TRUE)
#>   terms          r        r_n null_model
#> Z    ~Z -0.1852201 -0.2467218       ~1+X
#> X    ~X  0.3229866  0.4340634       ~1+Z
#normalized_gcor_binomial(model1, "X", verbose = TRUE)[1:3]
model0=glm(Y~Z,family="binomial")
if (FALSE) sapply(c('v','kl','sse','lr','n'),function(type) rsq.partial(model1,model0,type=type)$partial.rsq) # \dontrun{}#'

set.seed(1)
Z=rnorm(20)
X=10+Z+rnorm(20)
Y=rnorm(n=20,mean=(Z+2*X))
model1=glm(Y~Z+X)
gcor(model1)
#>   terms      gcor null_model
#> 1    ~Z 0.7986106       ~1+X
#> 2    ~X 0.8963620       ~1+Z
gcor(model1, normalize=TRUE)
#>   terms         r       r_n null_model
#> Z    ~Z 0.7986106 0.7986106       ~1+X
#> X    ~X 0.8963620 0.8963620       ~1+Z

#library(rsq)
model0x=glm(X~Z)
model0y=glm(Y~Z)
#t(.nrmz(resid(model0x)))%*% .nrmz(resid(model0y))
#library(ppcor)
#pcor(cbind(X,Z,Y))
if (FALSE) sapply(c('v','kl','sse','lr','n'),function(type) rsq.partial(model1,model0,type=type)$partial.rsq) # \dontrun{}
```
