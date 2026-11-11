# Update method for flipscores

Ensures that [`update()`](https://rdrr.io/r/stats/update.html) on a
`flipscores` object returns a proper `flipscores` object rather than a
plain `glm`.

## Usage

``` r
# S3 method for class 'flipscores'
update(object, ...)
```

## Arguments

- object:

  A `flipscores` object.

- ...:

  Additional arguments to update in the call, such as `formula`, `data`,
  `family`, `score_type`, etc. Supports the `.` shorthand in formula
  updates, e.g. `formula = . ~ . + offset(.OFFSET___)`.

## Value

A `flipscores` object.
