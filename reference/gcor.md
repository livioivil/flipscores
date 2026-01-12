# Compute Generalized Partial Correlations for GLM terms

This function computes the generalized partial correlation coefficient
\\r\\ for each specified variable in a generalized linear model. For
each variable, it refits the null model excluding that variable and
computes the cosine similarity between the residualized predictor and
standardized residuals.

## Usage

``` r
gcor(
  full_glm,
  terms = NULL,
  normalize = FALSE,
  intercept_too = FALSE,
  algorithm = "auto",
  algorithm.control = list(n_exact = 15, thresholds = c(-0.1, 0, 0.1), n_random = max(1,
    13 + log(1/nrow(model.matrix(full_glm)))), max_iter = 1000, topK = max(10, min(100,
    length(nrow(model.matrix(full_glm)))/10)), tol = 1e-12, patience = 10)
)
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

- algorithm.control:

  Only used if `normalize` is `TRUE`. \`list\` of control parameters:
  \`n_exact\` Integer specifying the sample size threshold for using
  exact methods (brute force). Default is 15. \`thresholds\` Numeric
  vector of threshold values for multi-start initialization.
  \`n_random\` Integer number of random starts for multi-start
  optimization. \`max_iter\` Integer maximum number of iterations per
  start. \`topK\` Integer number of top candidates to consider at each
  iteration. \`tol\` Numeric tolerance for convergence. \`patience\`
  Integer number of iterations without improvement before stopping.

- algorith:

  Only used if `normalize` is `TRUE`. \`"auto"\` by default. It choose
  between \`"intercept_only"\`, \`"brute_force"\` and \`"multi_start"\`

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
#> flipscores(formula = Y ~ Z + X, family = "poisson", data = dt, 
#>     n_flips = 1000)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)   
#> (Intercept) -0.14256 -0.91360    2.62144 -0.34851    -0.127    0.722   
#> ZB          -0.18558 -0.50868    1.65785 -0.30683    -0.108    0.665   
#> ZC           1.40981  8.55380    2.58950  3.30326     0.765    0.006 **
#> X           -0.06964 -1.56935    4.70999 -0.33320    -0.117    0.682   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
#>    terms          r        r_n null_model   algorithm is.exact
#> X     ~X -0.1174534 -0.1174534   ~1+ZB+ZC from theory     TRUE
#> ZC   ~ZC  0.7651942  0.7651942    ~1+ZB+X from theory     TRUE


gcor(mod, intercept_too=TRUE, normalize=TRUE)
#> Warning: The Normalized Generalized Partial Correlation (Determination) Coefficient for Count families without interncept in the null model has not implemented, yet. NA will be returned.
#>                    terms          r        r_n null_model   algorithm is.exact
#> (Intercept) ~(Intercept) -0.1265802         NA ~1+ZB+ZC+X        <NA>       NA
#> ZB                   ~ZB -0.1080979 -0.1080979    ~1+ZC+X from theory     TRUE
#> ZC                   ~ZC  0.7651942  0.7651942    ~1+ZB+X from theory     TRUE
#> X                     ~X -0.1174534 -0.1174534   ~1+ZB+ZC from theory     TRUE
set.seed(123)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rbinom(n=20,prob=plogis((dt$Z=="C")*2),size=1)
mod=flipscores(Y~Z+X,data=dt,family="binomial",n_flips=1000)
summary(mod)
#> 
#> Call:
#> flipscores(formula = Y ~ Z + X, family = "binomial", data = dt, 
#>     n_flips = 1000)
#> 
#> Coefficients:
#>             Estimate    Score Std. Error  z value Part. Cor Pr(>|z|)  
#> (Intercept)  -0.1486  -0.2102     1.1881  -0.1770    -0.067    0.911  
#> ZB          -20.4539  -1.4784     0.7466  -1.9802    -0.530    0.065 .
#> ZC           20.8561   1.8043     0.8180   2.2057     0.615    0.039 *
#> X            -0.4276  -0.3782     0.9574  -0.3951    -0.149    0.754  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
#>    terms          r        r_n null_model   algorithm is.exact
#> ZB   ~ZB -0.5295718 -0.5304792    ~1+ZC+X multi_start    FALSE
#> ZC   ~ZC  0.6146304  0.6148542    ~1+ZB+X multi_start    FALSE
#> X     ~X -0.1493213 -0.1891066   ~1+ZB+ZC multi_start    FALSE
# Compute for specific terms only
gcor(mod, terms = c("X", "ZC"),normalize=TRUE)
#>    terms          r        r_n null_model   algorithm is.exact
#> X     ~X -0.1493213 -0.1891066   ~1+ZB+ZC multi_start    FALSE
#> ZC   ~ZC  0.6146304  0.6148542    ~1+ZB+X multi_start    FALSE

```
