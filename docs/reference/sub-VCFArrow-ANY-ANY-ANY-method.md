# Subset method for VCFArrow

Subset a VCFArrow object by variants (rows) and samples (columns)

## Usage

``` r
# S4 method for class 'VCFArrow,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A VCFArrow object

- i:

  Row indices (numeric or logical)

- j:

  Column indices: numeric, logical, or sample name character vector

- drop:

  Ignored; kept for S4 compatibility

## Value

A new VCFArrow object containing the selected variants and samples

## Details

This function is a method of the VCFArrow S4 class Method to subset by
row and column of GT

## Author

Tomas Hrbek April 2026
