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
#> (Intercept) -0.14256 -0.91360    2.62144 -0.34851    -0.127    0.726   
#> ZB          -0.18558 -0.50868    1.65785 -0.30683    -0.108    0.642   
#> ZC           1.40981  8.55380    2.58950  3.30326     0.765    0.002 **
#> X           -0.06964 -1.56935    4.70999 -0.33320    -0.117    0.669   
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
# equivalent to
mod0=glm(Y~1,data=dt,family="poisson")
(results <- gcor(mod, mod0))
#> Error in gcor(mod, mod0): terms not found in model: c(`(Intercept)` = 0.470003629247622), c(`1` = 0.249999999997642, `2` = -0.375000000001179, `3` = 1.49999999999528, `4` = -0.375000000001179, `5` = -0.375000000001179, `6` = 1.49999999999528, `7` = -1, `8` = -0.375000000001179, `9` = 1.49999999999528, `10` = -0.375000000001179, `11` = -0.375000000001179, `12` = 2.1249999999941, `13` = -0.375000000001179, `14` = -1, `15` = -0.375000000001179, `16` = -1, `17` = -1, `18` = 0.874999999996462, `19` = -0.375000000001179, `20` = -0.375000000001179), c(`1` = 1.60000000000302, `2` = 1.60000000000302, `3` = 1.60000000000302, `4` = 1.60000000000302, `5` = 1.60000000000302, `6` = 1.60000000000302, `7` = 1.60000000000302, `8` = 1.60000000000302, `9` = 1.60000000000302, `10` = 1.60000000000302, `11` = 1.60000000000302, `12` = 1.60000000000302, `13` = 1.60000000000302, `14` = 1.60000000000302, `15` = 1.60000000000302, `16` = 1.60000000000302, `17` = 1.60000000000302, `18` = 1.60000000000302, `19` = 1.60000000000302, `20` = 1.60000000000302), c(`(Intercept)` = -2.6587446098628, -0.53212985828345, 1.8395760831688, -0.53212985828345, -0.53212985828345, 1.8395760831688, -1.3226985054342, -0.53212985828345, 1.8395760831688, -0.53212985828345, -0.53212985828345, 2.63014473031955, -0.53212985828345, -1.3226985054342, -0.53212985828345, -1.3226985054342, -1.3226985054342, 1.04900743601805, -0.53212985828345, -0.53212985828345), -5.6568597440809, 1, list(qr = c(-5.6568597440809, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979), rank = 1, qraux = 1.22360679774998, pivot = 1, tol = 1e-11), list(family = "poisson", link = "log", linkfun = function (mu) 
#> log(mu), linkinv = function (eta) 
#> pmax(exp(eta), .Machine$double.eps), variance = function (mu) 
#> mu, dev.resids = function (y, mu, wt) 
#> {
#>     r <- mu * wt
#>     p <- which(y > 0)
#>     r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
#>     2 * r
#> }, aic = function (y, n, mu, wt, dev) 
#> -2 * sum(dpois(y, mu, log = TRUE) * wt), mu.eta = function (eta) 
#> pmax(exp(eta), .Machine$double.eps), initialize = expression({
#>     if (any(y < 0)) 
#>         stop("negative values not allowed for the 'Poisson' family")
#>     n <- rep.int(1, nobs)
#>     mustart <- y + 0.1
#> }), validmu = function (mu) 
#> all(is.finite(mu)) && all(mu > 0), valideta = function (eta) 
#> TRUE, simulate = function (object, nsim) 
#> {
#>     wts <- object$prior.weights
#>     if (any(wts != 1)) 
#>         warning("ignoring prior weights")
#>     ftd <- fitted(object)
#>     rpois(nsim * length(ftd), ftd)
#> }, dispersion = 1), c(`1` = 0.470003629247622, `2` = 0.470003629247622, `3` = 0.470003629247622, `4` = 0.470003629247622, `5` = 0.470003629247622, `6` = 0.470003629247622, `7` = 0.470003629247622, `8` = 0.470003629247622, `9` = 0.470003629247622, `10` = 0.470003629247622, `11` = 0.470003629247622, `12` = 0.470003629247622, `13` = 0.470003629247622, `14` = 0.470003629247622, `15` = 0.470003629247622, `16` = 0.470003629247622, `17` = 0.470003629247622, `18` = 0.470003629247622, `19` = 0.470003629247622, `20` = 0.470003629247622
#> ), 28.6494739737397, 69.5328874955007, 5, c(`1` = 1.60000310821015, `2` = 1.60000310821015, `3` = 1.60000310821015, `4` = 1.60000310821015, `5` = 1.60000310821015, `6` = 1.60000310821015, `7` = 1.60000310821015, `8` = 1.60000310821015, `9` = 1.60000310821015, `10` = 1.60000310821015, `11` = 1.60000310821015, `12` = 1.60000310821015, `13` = 1.60000310821015, `14` = 1.60000310821015, `15` = 1.60000310821015, `16` = 1.60000310821015, `17` = 1.60000310821015, `18` = 1.60000310821015, `19` = 1.60000310821015, `20` = 1.60000310821015), c(`1` = 1, `2` = 1, `3` = 1, `4` = 1, `5` = 1, `6` = 1, `7` = 1, `8` = 1, `9` = 1, `10` = 1, `11` = 1, `12` = 1, `13` = 1, `14` = 1, `15` = 1, `16` = 1, `17` = 1, `18` = 1, `19` = 1, `20` = 1), 19, c(`1` = 2, `2` = 1, `3` = 4, `4` = 1, `5` = 1, `6` = 4, `7` = 0, `8` = 1, `9` = 4, `10` = 1, `11` = 1, `12` = 5, `13` = 1, `14` = 0, `15` = 1, `16` = 0, `17` = 0, `18` = 3, `19` = 1, `20` = 1), TRUE, FALSE, list(Y = c(2, 1, 4, 1, 1, 4, 0, 1, 4, 1, 1, 5, 1, 0, 1, 0, 0, 3, 1, 1)), glm(formula = Y ~ 1, family = "poisson", data = dt), Y ~ 1, Y ~ 1, list(X = c(-0.626453810742332, 0.183643324222082, -0.835628612410047, 1.59528080213779, 0.329507771815361, -0.820468384118015, 0.487429052428485, 0.738324705129217, 0.575781351653492, -0.305388387156356, 1.51178116845085, 0.389843236411431, -0.621240580541804, -2.2146998871775, 1.12493091814311, -0.0449336090152309, -0.0161902630989461, 0.943836210685299, 0.821221195098089, 0.593901321217509), Z = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2), Y = c(2, 1, 4, 1, 1, 4, 0, 1, 4, 1, 1, 
#> 5, 1, 0, 1, 0, 0, 3, 1, 1)), NULL, list(epsilon = 1e-08, maxit = 25, trace = FALSE), glm.fit
# Compute for specific terms only
gcor(mod, terms = c("X", "ZC"))
#>   terms       gcor null_model
#> 1    ~X -0.1174534   ~1+ZB+ZC
#> 2   ~ZC  0.7651942    ~1+ZB+X
```
