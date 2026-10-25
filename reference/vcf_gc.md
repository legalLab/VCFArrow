# Garbage-collect VCFArrow temp directories

Triggers R's garbage collector (three full passes to handle the
finalizer -\> pending-queue -\> unlink chain), then flushes any
directories whose reference count has already reached zero.

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

## Value

Invisibly returns `NULL`; called for its side effect of deleting
temporary directories.

## Typical workflow


      rm(vcf1, vcf2, vcf3)
      vcf_gc()              # usually sufficient
      vcf_gc(force = TRUE)  # if directories remain after rm()

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f)
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
rm(vcf)
vcf_gc()
#> VCFArrow: 2 temp directories are still live (VCFArrow objects not yet GC'd):
#>   /tmp/Rtmpmnning/arrow_vcf_3670398f295f
#>   /tmp/Rtmpmnning/arrow_vcf_36704dc0020a
#>   → rm() all VCFArrow objects, then call vcf_gc() again.
#>   → Or call vcf_gc(force = TRUE) to delete regardless.
```
