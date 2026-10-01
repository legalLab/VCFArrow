# vcf2genlight

Converts a VCFArrow object to Genlight format infile

## Usage

``` r
vcf2genlight(
  vcf_arrow,
  out_file = NULL,
  keep_groups = NULL,
  ploidy = 2L,
  save = FALSE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of .rds file to write, required when save = TRUE, default
  NULL (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- ploidy:

  -\> ploidy level, default = 2 (integer)

- save:

  -\> save as R data object, default = FALSE (Boolean)

## Value

An adegenet `genlight` object.

## Details

This function converts a VCFArrow object to an adegenet Genlight object.
Genotypes are read in chunks whose size is determined by the read_vcf()
function. If no groups are defined, the default behavior is to use all
groups. Genlight objects can encode polyploid genomes, by default
diploid genomes are assumed. Genlight objects are in-memory S4 objects
and are always returned; if save = TRUE, the object is also written to
out_file with saveRDS(). Unlike the other exporters, out_file is
optional because nothing is written unless save = TRUE.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
gl <- vcf2genlight(vcf)
#> ℹ Accumulating Genlight: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Building Genlight object...
gl
#>  /// GENLIGHT OBJECT /////////
#> 
#>  // 18 genotypes,  9,313 binary SNPs, size: 1.9 Mb
#>  79069 (47.17 %) missing data
#> 
#>  // Basic content
#>    @gen: list of 18 SNPbin
#>    @ploidy: ploidy of each individual  (range: 2-2)
#> 
#>  // Optional content
#>    @ind.names:  18 individual labels
#>    @loc.names:  9313 locus labels
#>    @chromosome: factor storing chromosomes of the SNPs
#>    @position: integer storing positions of the SNPs
#>    @pop: population of each individual (group size range: 3-5)
#>    @strata: a data frame with 1 columns ( pop )
#>    @other: a list containing: elements without names 
#> 
gl <- vcf2genlight(vcf, out_file = tempfile(fileext = ".rds"),
                   save = TRUE)
#> ℹ Accumulating Genlight: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Building Genlight object...
#> ℹ Writing Genlight object...
#> ✔ Genlight object written to /tmp/RtmpodaNXU/file1ffe19ecdaf4.rds
```
