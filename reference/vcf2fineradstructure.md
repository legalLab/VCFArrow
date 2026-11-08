# vcf2fineradstructure

Converts a VCFArrow object to FineRadStructure format infile

## Usage

``` r
vcf2fineradstructure(
  vcf_arrow,
  out_file,
  keep_groups = NULL,
  block_size = 10000
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

  -\> size in bp of the blocks that define RAD loci, as in
  vcf_filter_multiSNV(), default 10000 (integer)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external fineRADstructure
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. The
first line holds the sample names; each further line is one RAD locus,
defined as in vcf_filter_multiSNV(): the SNVs on one CHROM within the
same block of block_size bases (counted from the CHROM's first SNV).
Each cell holds a sample's two haplotypes over the locus's SNVs in
RADpainter's syntax, e.g. 'TGAT/CGAT', or 'TGAT' if both are the same; a
sample missing at any SNV of the locus has an empty cell. The genotypes
are taken as phased: the first allele of each SNV (a1) forms one
haplotype and the second (a2) the other. This holds for phased VCFs, and
for loci reconstructed as phased genotypes such as the output of
vcf_filter_multiSNV(); for unphased heterozygotes at several SNVs of a
locus the haplotypes are arbitrary. Use the same block_size as in
vcf_filter_multiSNV(). All genotypes are held in memory while writing (2
bytes per sample and SNV, see vcf_memory_estimate(format =
"individual")).

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
vcf_linked <- vcf_filter_multiSNV(vcf)
#> ℹ Applying linked SNV filter
#> ℹ Retained 8 / 10000 variants (linked SNVs)
vcf2fineradstructure(vcf_linked, out_file = tempfile(fileext = ".txt"))
#> ℹ Accumulating fineRADstructure: 8 variants x 18 samples (0 MiB raw storage)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Writing fineRADstructure file: 4 loci from 8 SNVs...
#> ✔ fineRADstructure file written to /tmp/RtmpVdrUAs/file31532eb153b2.txt
```
