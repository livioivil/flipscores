# It creates a `n_flips`x`n_obs` matrix of random +1 and -1. The first row is made by ones (i.e. the observed test statistic is computed)

It creates a `n_flips`x`n_obs` matrix of random +1 and -1. The first row
is made by ones (i.e. the observed test statistic is computed)

## Usage

``` r
make_flips(n_obs, n_flips)
```

## Arguments

- n_obs:

  number of observations

- n_flips:

  number of flips

## Examples

``` r
# example code
make_flips(n_obs=10,n_flips=20)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,]    1    1    1    1    1    1    1    1    1     1
#>  [2,]    1    1   -1    1    1    1    1    1    1     1
#>  [3,]    1   -1   -1    1   -1    1    1    1    1    -1
#>  [4,]   -1    1    1   -1   -1   -1   -1   -1    1    -1
#>  [5,]   -1   -1   -1   -1    1    1    1    1   -1    -1
#>  [6,]   -1    1    1   -1   -1    1    1   -1    1    -1
#>  [7,]    1   -1   -1    1    1    1   -1    1    1    -1
#>  [8,]    1    1    1    1   -1    1    1    1   -1    -1
#>  [9,]    1   -1   -1   -1    1    1   -1   -1    1     1
#> [10,]   -1    1    1    1   -1   -1    1    1   -1    -1
#> [11,]   -1    1    1   -1    1    1    1    1    1    -1
#> [12,]    1   -1   -1    1   -1   -1    1    1   -1    -1
#> [13,]   -1   -1   -1    1    1   -1    1   -1    1    -1
#> [14,]    1    1    1    1   -1   -1   -1   -1   -1    -1
#> [15,]   -1    1    1   -1   -1    1   -1   -1    1    -1
#> [16,]    1   -1    1   -1   -1    1    1   -1    1     1
#> [17,]    1    1    1   -1   -1    1   -1    1   -1     1
#> [18,]    1    1    1   -1   -1    1    1   -1   -1    -1
#> [19,]   -1    1    1   -1   -1    1    1   -1   -1    -1
#> [20,]    1    1   -1    1   -1    1    1   -1   -1    -1
```
