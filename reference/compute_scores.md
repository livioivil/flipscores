# compute_scores

Same usage as `anova.glm`. The parameter `id` is used too, if present in
`model0` (with priority) or in `model1`.

## Usage

``` r
compute_scores(model0, model1, score_type = "standardized", ...)
```

## Arguments

- model0:

  a `glm` object with the model under the null hypothesis (i.e. the
  covariates, the nuisance parameters).

- model1:

  a `glm` or a `matrix` (or `vector`). If it is a `glm` object, it has
  the model under the alternative hypothesis. The variables in `model1`
  are the same variables in `model0` plus one or more variables to be
  tested. Alternatively, if `model1` is a `matrix`, it contains the
  tested variables column-wise.

- score_type:

  The type of score that is computed. It is "orthogonalized",
  "effective" or "basic". "effective" and "orthogonalized" take into
  account the nuisance estimation.

- ...:

  other arguments.

## Author

Jesse Hemerik, Riccardo De Santis, Vittorio Giatti, Jelle Goeman and
Livio Finos

## Examples

``` r
set.seed(1)
Z=rnorm(20)
X=Z+rnorm(20)
Y=rpois(n=20,lambda=exp(Z+X))
mod0=glm(Y~Z,family="poisson")
X=data.frame(X=X)
scr0=compute_scores(model0 = mod0, model1 = X)
head(scr0)
#>             X
#> 1 -0.09157255
#> 2  1.29487548
#> 3  0.14949028
#> 4 46.19093971
#> 5  1.30167541
#> 6  0.18155826
```
