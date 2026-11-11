# Flip-score tests for custom contrasts

Defines factor contrasts from a formula or accepts fully custom
coefficient contrasts, then applies the flip-score test to each
contrast.

## Usage

``` r
flipscores_contrasts(
  model,
  specs = NULL,
  linfct = NULL,
  ref = NULL,
  score_type = "standardized",
  n_flips = 5000,
  alternative = "two.sided",
  id = NULL,
  seed = NULL,
  flips = NULL,
  ...
)
```

## Arguments

- model:

  a fitted `glm`, `lm`, or `flipscores` object.

- specs:

  a contrast formula, coefficient contrast matrix, numeric vector, named
  list of numeric vectors, or `NULL` when `linfct` is supplied.

- linfct:

  optional custom coefficient contrast matrix, numeric vector, or named
  list of numeric vectors. This follows the same idea as
  `multcomp::glht()`: rows define linear functions of the model
  coefficients.

- ref:

  control level for Dunnett-like formula contrasts. It can be a level
  name, cell name, or numeric index. The default is the first
  level/cell, except for `trt.vs.ctrlk`, where the default is the last
  level/cell.

- score_type:

  type of score, as in
  [`flipscores()`](https://livioivil.github.io/flipscores/reference/flipscores.md).

- n_flips:

  number of sign flips.

- alternative:

  `"two.sided"`, `"greater"`, or `"less"`.

- id:

  optional cluster id.

- seed:

  optional random seed.

- flips:

  optional precomputed flip matrix.

- ...:

  currently unused.

## Value

an object of class `fs_contrasts`.
