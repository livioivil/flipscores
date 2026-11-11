# Nonparametric combination of `jfs` objects

Methods for combining `jfs` objects.

`combine_tests` combines the tests derived from multiverse models.

`combine_contrasts` combines the tests derived from the contrasts of a
factor variable to get a global test for the factor (i.e. categorical
predictor). It has strong analogies with ANOVA test.

## Usage

``` r
combine_tests(mods, comb_funct = "maxT", by = NULL, by_list=NULL, tail = 0)

combine_contrasts(mods, comb_funct = "Mahalanobis", tail = 0)
```

## Arguments

- mods:

  a `jfs` object.

- comb_funct:

  combining function to be used. Several functions are implemented:
  "mean", "median", "Fisher", "Liptak", (equal to) "Stoufer", "Tippet",
  (equal to) "minp", "maxT", "Mahalanobis". Alternatively it can be a
  custom function that has a Tspace matrix as input. For `combine_tests`
  the default is `comb_funct="maxT"`, while for `combine_contrasts` the
  default is `comb_funct="Mahalanobis"`.

- by:

  if `NULL` (default), it combines all test statistics. If a characters,
  it refers to the column's name of `summary_table` (and printed by
  something like `summary(mods)`). The elements with the same value will
  be combined. If `by` is a vector, the values are defined by row-wise
  concatenation of the values of the columns in `by`. The argument is
  inactive if `by_list` is not `NULL`.

- by_list:

  NULL (default) or a list of vectors. For each vector of the list it
  combines test statistics with position given by the element of the
  vector. If the vectors in the list are characters, these refer to
  names(mods\$Tspace).

- tail:

  direction of the alternative hypothesis. It can be "two.sided" (or 0,
  the default), "less" (or -1) or "greater" (or +1).

## Value

The function returns a `jfs`-object.

## Examples

``` r
#First example
library(flipscores)
set.seed(123)

#Simulate data
n=20
D=data.frame(X=rnorm(n),Z1=rnorm(n),Z2=rnorm(n))
D$Y=D$Z1+D$X+rnorm(n)

# Run four glms abd combine it in a list
mod1=glm(Y~X+Z1+Z2,data=D)
mod2=glm(Y~X+poly(Z1,2)+Z2,data=D)
mod3=glm(Y~X+poly(Z1,2)+poly(Z2,2),data=D)
mod4=glm(Y~X+Z1+poly(Z2,2),data=D)
mods=list(mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)

# Let us analyze the tests related to coefficient "X" and combine them
res=flipscores(mods,n_flips = 5000, seed = 1, tested_coeffs = "X")
summary(combine_tests(res))
#>     model coefficient stat ntests     S      p
#> 1 overall           X maxT      4 3.136 0.0232
# Second (continued) example
# flipscores jointly on all models and all coefficients
mod1=glm(Y~X+Z1+Z2,data=D)
mod2=glm(Y~X+poly(Z1,2)+Z2,data=D)
mod3=glm(Y~X+poly(Z1,2)+poly(Z2,2),data=D)
mod4=glm(Y~X+Z1+poly(Z2,2),data=D)
mods=list(mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)
res=flipscores(mods,n_flips = 2000)
summary(combine_tests(res))
#>     model coefficient stat ntests     S      p
#> 1 overall        many maxT     20 3.247 0.0105
summary(combine_tests(res, by="model"))
#>   model coefficient stat ntests     S      p
#> 1  mod1        many maxT      4 3.136 0.0145
#> 2  mod2        many maxT      5 3.247 0.0060
#> 3  mod3        many maxT      6 3.091 0.0095
#> 4  mod4        many maxT      5 3.001 0.0125
summary(combine_tests(res, by="coefficient"))
#>          model  coefficient stat ntests      S      p
#> 1  (Intercept)  (Intercept) maxT      4 0.8948 0.5755
#> 2            X            X maxT      4 3.1356 0.0205
#> 3           Z1           Z1 maxT      2 3.1222 0.0010
#> 4           Z2           Z2 maxT      2 0.6681 0.5875
#> 5 poly(Z1, 2)1 poly(Z1, 2)1 maxT      2 3.2475 0.0010
#> 6 poly(Z1, 2)2 poly(Z1, 2)2 maxT      2 2.3372 0.0585
#> 7 poly(Z2, 2)1 poly(Z2, 2)1 maxT      2 0.6175 0.6505
#> 8 poly(Z2, 2)2 poly(Z2, 2)2 maxT      2 0.7280 0.5260
res2=combine_contrasts(res)
summary(res2)
#>     model coefficient        stat ntests        S      p
#> 1  mod1.0 (Intercept) Mahalanobis      1 0.009994 0.6785
#> 2  mod1.1           X Mahalanobis      1 0.044775 0.0140
#> 3  mod1.2          Z1 Mahalanobis      1 0.057614 0.0010
#> 4  mod1.3          Z2 Mahalanobis      1 0.019639 0.3950
#> 5  mod2.0 (Intercept) Mahalanobis      1 0.019358 0.4015
#> 6  mod2.1           X Mahalanobis      1 0.045576 0.0135
#> 7  mod2.2 poly(Z1, 2) Mahalanobis      2 6.697510 0.0045
#> 8  mod2.3          Z2 Mahalanobis      1 0.003766 0.8775
#> 9  mod3.0 (Intercept) Mahalanobis      1 0.020063 0.3870
#> 10 mod3.1           X Mahalanobis      1 0.047143 0.0165
#> 11 mod3.2 poly(Z1, 2) Mahalanobis      2 6.187591 0.0140
#> 12 mod3.3 poly(Z2, 2) Mahalanobis      2 0.955023 0.6485
#> 13 mod4.0 (Intercept) Mahalanobis      1 0.012247 0.6115
#> 14 mod4.1           X Mahalanobis      1 0.047215 0.0120
#> 15 mod4.2          Z1 Mahalanobis      1 0.059404 0.0010
#> 16 mod4.3 poly(Z2, 2) Mahalanobis      2 1.079869 0.6140
#custom combinations:
coeffs=c("(Intercept)","X","Z1","Z2")
coeffs_ids=lapply(coeffs,grep,res2$summary_table$coeff)
names(coeffs_ids)=coeffs
summary(combine_tests(res2,by_list =   coeffs_ids))
#>         model coefficient stat ntests       S      p
#> 1 (Intercept) (Intercept) maxT      4 0.02006 0.5625
#> 2           X           X maxT      4 0.04721 0.0255
#> 3          Z1        many maxT      4 6.69751 0.0100
#> 4          Z2        many maxT      4 1.07987 0.7595
```
