# Garbage-collect VCFArrow temp directories

Triggers R's garbage collector (three full passes to handle the
finalizer → pending-queue → unlink chain), then flushes any directories
whose reference count has already reached zero.

## Usage

``` r
vcf_gc(force = FALSE, verbose = TRUE)
```

## Arguments

- force:

  Logical. If TRUE, also force-deregister and delete directories for
  objects that are still nominally live in the registry. Use this when
  you have called [`rm()`](https://rdrr.io/r/base/rm.html) on all
  VCFArrow objects but the temp directories have not been cleaned up
  (common in RStudio, which can hold display references that delay GC).

- verbose:

  Logical. Print a status message.

## Typical workflow


      rm(vcf1, vcf2, vcf3)
      vcf_gc()              # usually sufficient
      vcf_gc(force = TRUE)  # if directories are still present after rm()
