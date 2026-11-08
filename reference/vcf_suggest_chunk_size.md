# Suggest a chunk_size for read_vcf() given available RAM

Suggest a chunk_size for read_vcf() given available RAM

## Usage

``` r
vcf_suggest_chunk_size(
  available_gb,
  n_samples,
  fmt_chars = 25,
  sample_chars = 10,
  workers = .vcf_workers(),
  n_columns = NULL
)
```

## Arguments

- available_gb:

  RAM available in gigabytes.

- n_samples:

  Number of samples.

- fmt_chars:

  Average length in characters of a sample's FORMAT text in the VCF
  (e.g. 24 for "0/1:12,8:20:99:200,0,300"), default 25. Longer FORMAT
  fields (e.g. PL with many alleles) need more memory.

- sample_chars:

  Average length of the sample names, default 10.

- workers:

  Number of worker processes read_vcf() will use (see
  vcf_set_workers()), default the number currently set.

- n_columns:

  Deprecated and ignored: the size of a genotype row is now computed
  from its columns.

## Value

Invisibly, the suggested chunk size (integer).

## Details

Half of `available_gb` is left for the operating system, R and later
operations; the other half is shared by the chunks read_vcf() holds at
once (one per worker). A chunk holds one row per variant and sample, of
about 44 + sample_chars + fmt_chars bytes in memory: .row_id, a1 and a2
(4 bytes each), the sample name and the sample's FORMAT text (each 4
bytes plus its characters), phased (1 bit), and DP, GQ and ADR (8 bytes
each). The FORMAT text is kept so write_vcf() can write the sample
fields unchanged; it is often the largest column. With workers, each
chunk's raw VCF lines (about fmt_chars + 1 bytes per genotype) are also
held, by the main process and by the worker parsing them. The result is
limited to between 1,000 and 500,000 variants.

## Examples

``` r
vcf_suggest_chunk_size(available_gb = 8, n_samples = 100)
#> ℹ For 8 GB RAM, 100 samples and 1 worker, recommended chunk_size = 500,000
#> (use: read_vcf(vcf_file, chunk_size = 500000))
```
