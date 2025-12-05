# Compute Generalized Partial Correlations for GLM terms

This function computes the generalized partial correlation coefficient
\\r\\ for each specified variable in a generalized linear model. For
each variable, it refits the null model excluding that variable and
computes the cosine similarity between the residualized predictor and
standardized residuals.

## Usage

``` r
gcor(full_glm, terms = NULL, intercept_too = FALSE)
```

## Arguments

- full_glm:

  A fitted GLM object of class \`glm\`.

- terms:

  Character vector of variable names (referred to the model.matrix, that
  is pay attention for factors) for which to compute generalized partial
  correlations. If \`NULL\` (default), computes for all non-intercept
  terms in the model.

## Value

A data frame with two columns:

- variable:

  The variable name

- r:

  The generalized partial correlation coefficient

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
#> (Intercept) -0.14256 -0.91360    2.62144 -0.34851    -0.127    0.737   
#> ZB          -0.18558 -0.50868    1.65785 -0.30683    -0.108    0.666   
#> ZC           1.40981  8.55380    2.58950  3.30326     0.765    0.003 **
#> X           -0.06964 -1.56935    4.70999 -0.33320    -0.117    0.654   
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
```
