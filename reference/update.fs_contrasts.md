# Convert an fs_contrasts object to a jfs object

This method creates a \`jfs\` object from a fitted model and a set of
contrasts. The resulting object can be used with \`flipscores\` for
joint permutation testing across multiple contrasts or models.

## Usage

``` r
# S3 method for class 'fs_contrasts'
update(object, formula., ..., evaluate = TRUE)
```

## Arguments

- object:

  An object of class \`fs_contrasts\`.

- formula.:

  Changes to the formula – see \[update.formula()\] for details.

- ...:

  Additional arguments passed to \`update()\` for the underlying model,
  or to the contrast‑generating function.

- data:

  New data frame.

- family:

  New family for the GLM.

- linfct:

  New contrast matrix (or formula to generate it).

## Value

An object of class \`jfs\` (joint flipscores) that contains the original
model, the contrast matrix, and other relevant information.

An updated object of class \`fs_contrasts\`.
