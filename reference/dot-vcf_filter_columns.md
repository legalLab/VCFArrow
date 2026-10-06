# vcf_filter_columns

Unified API for filtering VCFArrow objects by row IDs.

## Usage

``` r
.vcf_filter_columns(vcf_arrow, keep, f_invar = TRUE, verbose = TRUE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep:

  -\> rows to keep (numeric or logical)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

- verbose:

  -\> report removed samples and final % missing data, default TRUE
  (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes all rows in the 'keep' parameter from a VCFArrow
object, returning a new VCFArrow object. Optionally will remove any loci
that may have become invariant as the result of the removal of samples.
Optionally will report removed samples, final % missing data, and number
of retained samples after sample filtering.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
.vcf_filter_columns(vcf, keep = vcf@groups != "OG")
#> ℹ Compacting GT: 18 -> 15 samples across 1 chunk
#> ℹ Applying invariant filter
#> ℹ Removed 1705 invariant variants; 8295 retained.
#> ℹ Removed samples: Pb2Jp, Pb2Scx, and Pb1Rd
#> ℹ Variants retained: 8295 | Samples retained: 15
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 8295 
#>   Samples:  15 
#> 
#> Quick stats:
#>   Non-missing variants: 8295 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpFwXb6F/arrow_vcf_samp_1f6b4fdfea38 
#> 
#> Genotype storage (Arrow):
#> FileSystemDataset with 1 Feather file
#> 9 columns
#> .row_id: int32
#> sample: string
#> a1: int32
#> a2: int32
#> phased: bool
#> fmt: string
#> DP: double
#> GQ: double
#> ADR: double
#> 
#> See $metadata for additional Schema metadata
#> 
#> Variants (first 5 rows):
#>                     CHROM POS        ID REF ALT QUAL FILTER Rk RPT n_alt
#> 2 SNP_higher_path_9984432 105   9984432   C   T    .      .  1  NA     1
#> 5 SNP_higher_path_9803974  33   9803974   A   G    .      .  1  NA     1
#> 6 SNP_higher_path_9791272  31   9791272   T   G    .      .  1  NA     1
#> 7  SNP_higher_path_976062  43    976062   A   G    .      .  1  NA     1
#> 8 SNP_higher_path_9736643  72 9736643_3   G   A    .      .  1  NA     1
#>   is_biallelic is_indel .row_id
#> 2         TRUE    FALSE       2
#> 5         TRUE    FALSE       5
#> 6         TRUE    FALSE       6
#> 7         TRUE    FALSE       7
#> 8         TRUE    FALSE       8
#>   ... 8290 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=3;UR=1;CL=75;CR=7;Genome=.;Sd=.;Cluster=17296599;ClSize=1" 
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"  
#> [3] "Ty=SNP;Rk=1.0;UL=1;UR=7;CL=1;CR=42;Genome=.;Sd=.;Cluster=14628692;ClSize=2" 
#> [4] "Ty=SNP;Rk=1.0;UL=13;UR=8;CL=13;CR=8;Genome=.;Sd=.;Cluster=17727234;ClSize=4"
#> [5] "Ty=SNP;Rk=1.0;UL=16;UR=0;CL=16;CR=0;Genome=.;Sd=.;Cluster=18305039;ClSize=1"
#> 
#> FORMAT (first 5):
#>           FORMAT .row_id
#> 1 GT:DP:PL:AD:HQ       1
#> 2 GT:DP:PL:AD:HQ       2
#> 3 GT:DP:PL:AD:HQ       3
#> 4 GT:DP:PL:AD:HQ       4
#> 5 GT:DP:PL:AD:HQ       5
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 10 more
```
