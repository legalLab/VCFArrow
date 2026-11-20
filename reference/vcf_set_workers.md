# vcf_set_workers

Set the number of parallel worker processes used by VCFArrow

## Usage

``` r
vcf_set_workers(workers = 1L)
```

## Arguments

- workers:

  -\> number of worker processes, default 1 (no parallelism) (integer)

## Value

Invisibly, the previous number of workers.

## Details

Functions that process the on-disk genotype chunks independently
distribute the chunks over `workers` background R processes: read_vcf(),
the genotype-based filters, sample and group extraction, vcf_stats() and
the assess functions, vcf_bind(), write_vcf(), the vcf2\*() exporters
and vcf2gt_long() with CSV output (with feather or Parquet output it
runs serially). Results do not depend on the number of workers. The
workers are started once and reused until `vcf_set_workers(1)` is called
or the session ends. Each worker holds one chunk in memory at a time, so
memory use grows with the number of workers (see vcf_memory_estimate()).
Useful values are up to the number of physical CPU cores; disk speed
limits the gain beyond that.

## Author

Tomas Hrbek September 2026

## Examples

``` r
# \donttest{
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
# 2 workers (CRAN examples may use at most 2 cores)
vcf_set_workers(2)
#> ℹ VCFArrow will use 2 worker processes.
vcf <- read_vcf(f)
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_set_workers(1)
# }
```
