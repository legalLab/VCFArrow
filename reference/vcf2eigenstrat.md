# vcf2eigenstrat

Converts a VCFArrow object to Eigenstrat format infile

## Usage

``` r
vcf2eigenstrat(
  vcf_arrow,
  out_file,
  keep_groups = NULL,
  sex = NULL,
  chrom_code = c("auto", "index", "keep")
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> output file prefix; extensions are added automatically, no default
  (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- sex:

  -\> sex of the individual, default = U (undefined) (character)

- chrom_code:

  -\> how CHROM is written in the .snp file: "auto" (keep CHROM if it is
  numeric, otherwise number the contigs 1, 2, 3, ...), "index" (number
  the contigs 1, 2, 3, ... in order of first appearance) or "keep"
  (verbatim), default "auto" (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external EIGENSTRAT
formatted fileset (.geno, .ind and .snp). Writing occurs in chunks whose
size is determined by the read_vcf() function. Larger chunks result in
faster writing speeds. If no groups are defined, the default behavior is
to use all groups. In the .geno file, each digit is the number of copies
of the reference allele, as EIGENSTRAT defines it: the REF allele, given
in column 5 of the .snp file (2 = REF homozygote, 1 = heterozygote, 0 =
ALT homozygote, 9 = missing). EIGENSOFT needs integer chromosome codes
and removes SNPs on chromosome 0, so contig names (e.g. from de novo SNP
calling) are numbered 1, 2, 3, ... in order of first appearance
(chrom_code = "auto" or "index"), and the correspondence is written to
\<out_file\>.chrommap. EIGENSOFT reads chromosomes above its numchrom
parameter (default 22) as X, Y, ..., so with more than 22 codes set
numchrom to at least the largest code in the parameter file; a warning
gives the value.

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
vcf2eigenstrat(vcf, out_file = file.path(tempdir(), "eigenstrat_in"))
#> ℹ Building EIGENSTRAT: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing EIGENSTRAT files...
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Chromosome names were recoded as integers; mapping written to /tmp/RtmpVdrUAs/eigenstrat_in.chrommap
#> Warning: Chromosome codes go up to 9309: EIGENSOFT reads chromosomes above its numchrom
#> parameter (default 22) as X, Y, ...
#> ℹ Set `numchrom: 9309` (or higher) in the parameter file.
#> ✔ EIGENSTRAT fileset written to /tmp/RtmpVdrUAs/eigenstrat_in.geno, /tmp/RtmpVdrUAs/eigenstrat_in.ind, /tmp/RtmpVdrUAs/eigenstrat_in.snp
```
