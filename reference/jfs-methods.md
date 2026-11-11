# Methods for jfs objects

Collection of methods for objects of class `jfs`.

`p.adjust` method for class "`jfs`" and "`flipscores`". Add adjusted
p-values into the `jfs`and `flipscores` object.

## Usage

``` r
p.adjust(object, method = "maxT", tail = 0, ...)
```

## Arguments

- object:

  an object of class `jfs` and `flipscores`.

- method:

  any method implemented in `flip.adjust` or a custom function. In the
  last case it must be a function that uses a matrix as input and
  returns a vector of adjusted p-values equal to the number of columns
  of the inputed matrix.

- tail:

  argument: expresses the tail direction of the alternative hypothesis.
  It can be "two.sided" (or 0, the default), "less" (or -1) or "greater"
  (or +1).

- ...:

  additional arguments to be passed
