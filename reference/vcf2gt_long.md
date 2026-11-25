# vcf2gt_long

Converts a VCFArrow object to tidy long format infile

## Usage

``` r
vcf2gt_long(
  vcf_arrow,
  out_file,
  keep_groups = NULL,
  format = c("feather", "parquet", "csv"),
  col_select = NULL
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- format:

  -\> output format, "feather" (Arrow IPC), "parquet" or "csv", default
  "feather" (character)

- col_select:

  -\> optional selection of columns to save, default ALL

## Value

Invisibly returns the path of the written file.

## Details

This function converts a VCFArrow object to an external tidy long-format
genotype file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Only
biallelic SNVs are exported (multiallelic sites and indels are left
out), and every sample must be assigned to a group (see
set_vcf_groups()), also when keep_groups is NULL. The tidy data can be
saved in either Arrow, Parquet or CSV formats. File extension is added
automatically if missing. Optionally, specific columns can be saved, by
default all columns are saved. The gt long slot contains pre-calculated
metrics in addition to just genotypes.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2gt_long(vcf, out_file = tempfile(), format = "csv")
#> ℹ Writing gt_long: 9313 variants x 18 samples
#> ✔ gt_long table written to /tmp/Rtmpq4egqB/file1db26c968565.csv
```
