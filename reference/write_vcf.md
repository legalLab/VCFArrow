# write_vcf

Write a VCFArrow object to an external VCF file

## Usage

``` r
write_vcf(
  vcf_arrow,
  out_file,
  gzip = grepl("\\.gz$", out_file),
  update_info = TRUE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of the VCF file to be written to, no default (character)

- gzip:

  -\> a flag to GZIP VCF when writing, default TRUE if out_file ends in
  ".gz" (Boolean)

- update_info:

  -\> recompute the INFO fields AC, AN, AF and NS from the genotypes
  written, where present, default TRUE (Boolean)

## Value

Invisibly returns the path of the written file.

## Details

This function writes a VCFArrow object to an external VCF file. Writing
occurs in chunks whose size is determined by the read_vcf() function.
Larger chunks result in faster writing speeds. It writes both
uncompressed and gz compressed files. Compressing increases writing time
be about 50%. For large files, it is recommended to output an
uncompressed VCF file, and then compress with GZIP or PIGZ (or BGZIP,
which tabix needs to index the file). The sample fields are written as
stored; genotypes changed by VCFArrow (e.g. by vcf_filter_adr() or
vcf_filter_coverage()) have their GT updated. With update_info = TRUE
(the default), the INFO fields AC (count of each ALT allele), AN (number
of called alleles), AF (AC / AN) and NS (number of samples with a called
genotype) are recomputed from the samples and genotypes written,
wherever a variant's INFO has them, so they stay correct after samples
or genotypes were removed or changed; other INFO fields are written
unchanged. Haploid genotypes, which read_vcf() stores as homozygous,
count as two alleles. Variants without a valid position (POS) are not
written, with a warning. A warning is also given if the variants are not
sorted by position (e.g. after vcf_bind(mode = "union")), since indexing
tools such as tabix need sorted files; `bcftools sort` sorts them. A
'##source=VCFArrow' line with the package version is added to the
header.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
write_vcf(vcf, out_file = tempfile(fileext = ".vcf"))
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ VCF is being written in 1 chunk
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ✔ VCFArrow object successfully written to /tmp/RtmpVdrUAs/file3153bb76e28.vcf
write_vcf(vcf, out_file = tempfile(fileext = ".vcf.gz"))
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ VCF is being written in 1 chunk
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ✔ VCFArrow object successfully written to /tmp/RtmpVdrUAs/file3153dc8f12c.vcf.gz
```
