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

  One of "individual" (Structure, Arlequin, FASTA, …) or "pop"
  (BayesScan, Treemix, Migrate-N C, …) or "chunk" (SmartSNP,
  fineRADstructure, sNMF, EIGENSTRAT, …).

- chunk_size:

  Feather chunk size used at read_vcf() time.

- lowmem:

  If TRUE, estimate uses raw-byte matrices (vcf2\*() lowmem variants);
  otherwise integer matrices.
