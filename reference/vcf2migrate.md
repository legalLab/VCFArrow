# vcf2migrate

Converts a VCFArrow object to a MIGRATE-N format infile

## Usage

``` r
vcf2migrate(
  vcf_arrow,
  out_file,
  keep_groups = NULL,
  block_size = 100L,
  method = c("S", "N", "C")
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- block_size:

  -\> number of SNPs (method S) or base pairs (method N) per linked
  block, default 100 (integer)

- method:

  -\> "C" (allele counts, HapMap format), "S" (sequences, loci of
  block_size SNPs, "(sK)") or "N" (SNPs, loci of block_size bp per
  CHROM, "(nK)"), default "S" (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external MIGRATE-N
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. The
data file follows the Migrate-N manual (section 4, "Data file
specification",
<https://peterbeerli.com/downloads/migrate/docs/migratedoc.pdf>). Method
"S" writes sequence data: each sample as two haplotypes (\_1, \_2) of
the SNPs' nucleotides ('?' if missing), cut into loci of block_size
consecutive SNPs, labelled "(sK)" (a sequence locus of K sites). Method
"N" writes the same sequences cut into loci of the SNPs on one CHROM
within the same interval of block_size bp, labelled "(nK)" (K linked
SNPs), with the data type N (SNPs) on the first line. Method "C" writes
the allele counts of each population in Migrate-N's HapMap format (data
type H): for each SNP its position, REF and ALT nucleotides with their
counts, and the total; the manual describes this format as preliminary.
SNPs without called alleles in some population are left out (with a
warning), as Migrate-N cannot compute their allele frequencies there.
The first line gives the number of populations and of loci. Individual
names are written in a fixed field of 10 characters, Migrate-N's
default; if longer names are needed, a warning gives the nmlength to set
in the parmfile. Population names must not start with a number. Whether
loci of SNPs should be treated as linked depends on whether the SNPs
were mapped (so that linkage is known) or called de novo.

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
vcf2migrate(vcf, out_file = tempfile(fileext = ".txt"))
#> ℹ Accumulating Migrate-N (S): 9313 variants x 18 samples (0 MiB raw storage)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Writing Migrate-N (S) file: 94 loci...
#> ✔ Migrate-N file written to /tmp/RtmpVdrUAs/file31535f02c0cf.txt
```
