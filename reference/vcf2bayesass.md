# vcf2bayesass

Converts a VCFArrow object to BayesAss format infile

## Usage

``` r
vcf2bayesass(vcf_arrow, out_file, keep_groups = NULL)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external BayesAss
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. The
file is in BA3's native format: one line per individual and locus,
'indivID popID locID allele1 allele2', separated by spaces, with the
alleles labelled 1 (REF) and 2 (ALT) and a missing genotype as '0 0'.
Sample names, group labels and locus names must not contain spaces. BA3
accepts at most 100,000 loci, 100 populations and 20,000 individuals; a
warning is given above these limits.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
vcf2bayesass(vcf, out_file = tempfile(fileext = ".immanc"))
#> ℹ Accumulating BayesAss: 9313 variants x 18 samples (0 MiB raw storage)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Writing BayesAss file...
#> ✔ BayesAss file written to /tmp/RtmpVdrUAs/file315335cce7c9.immanc
```
