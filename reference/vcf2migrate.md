# vcf2migrate

Converts a VCFArrow object to a MIGRATE-N format infile

## Usage

``` r
vcf2migrate(
  vcf_arrow,
  out_file,
  keep_groups = NULL,
  block_size = 100L,
  method = "S"
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

  -\> "C" (allele counts), "S" (sequences, fixed blocks) or "N"
  (sequences, chromosome intervals), default "S" (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external MIGRATE-N
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. The
function implements two different SNP formats, new (un)linked SNPs (S) -
fixed blocks, and new (un)linked SNPs (N) - chromosome intervals.
Alleles are sequence data when SNPs are mapped to a reference. The S
option generates blocks of a specific number of SNPs and treats them as
linked; this is appropriate if SNPs are extracted without a reference.
The N format extracts all SNPs within a chromosome, within the block
size and treats them as linked; this is appropriate if SNPs are mapped
against a reference. The size of the linked block is determined by the
block_size parameter. See
<https://peterbeerli.com/programs/migrate/distribution_4.x/migratedoc4.x.pdf>
for format detail.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2migrate(vcf, out_file = tempfile(fileext = ".txt"))
#> ℹ Accumulating Migrate-N (S): 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing Migrate-N (S) file...
#> ✔ Migrate-N file written to /tmp/RtmpFwXb6F/file1f6b6dacf909.txt
```
