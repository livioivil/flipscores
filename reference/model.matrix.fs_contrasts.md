# Model matrix for fs_contrasts objects

Extract the model matrix from the underlying fitted model stored in an
\`fs_contrasts\` object.

## Usage

``` r
# S3 method for class 'fs_contrasts'
model.matrix(object, ...)
```

## Arguments

- object:

  An object of class \`fs_contrasts\`.

- ...:

  Additional arguments passed to \`model.matrix\` (e.g., \`data\`).

## Value

The model matrix (design matrix) of the original model.
