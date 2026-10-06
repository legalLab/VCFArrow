# Estimate RAM required for a VCFArrow export operation

Estimate RAM required for a VCFArrow export operation

## Usage

``` r
vcf_memory_estimate(
  vcf_arrow,
  keep_groups = NULL,
  format = c("individual", "pop", "chunk"),
  chunk_size = 100000L,
  lowmem = FALSE
)
```

## Arguments

- vcf_arrow:

  A VCFArrow object.

- keep_groups:

  Groups to export (NULL = all).

- format:

  One of "individual" (Structure, Arlequin, FASTA, ...) or "pop"
  (BayesScan, Treemix, Migrate-N C, ...) or "chunk" (SmartSNP,
  fineRADstructure, sNMF, EIGENSTRAT, ...).

- chunk_size:

  Feather chunk size used at read_vcf() time.

- lowmem:

  If TRUE, estimate uses raw-byte matrices (vcf2\*() lowmem variants);
  otherwise integer matrices.

## Value

Invisibly, a list with elements:

- `n_var`: number of variants

- `n_samples`: number of samples retained

- `n_pops`: number of populations

- `chunk_arrow_bytes`: Arrow memory per chunk read

- `matrix_bytes`: size of the accumulation matrices

- `peak_bytes`: estimated peak memory

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_memory_estimate(vcf, format = "individual")
#> 
#> ── VCFArrow memory estimate ──
#> 
#> • Variants (filtered): 9313
#> • Samples retained: 18
#> • Populations: 4
#> • Chunk size: 100,000 variants
#> • Low-memory mode: FALSE
#> 
#> ── Per-operation footprint 
#> • Arrow pool per chunk read: 34.3 MiB
#> • Accumulation matrices (a1+a2): 1.3 MiB [4 B/cell]
#> • Estimated peak RAM: 35.6 MiB
```
