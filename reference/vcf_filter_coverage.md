# vcf_filter_coverage

Remove genotypes below read DP threshold from a VCFArrow object

## Usage

``` r
vcf_filter_coverage(vcf_arrow, threshold = 10)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> DP threshold, default 10 (integer)

## Value

subsetted VCFArrow object

## Details

This function removes genotypes below a DP threshold from a VCFArrow
object, returning a new VCFArrow object.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_filter_coverage(vcf, threshold = 10)
#> ℹ Applying read coverage filter
#> ℹ Retained 5817 / 10000 variants (polymorphic with DP >= 10)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 5817 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 5817 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpOwogNL/arrow_vcf_1eff567167b1 
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
#>                      CHROM POS        ID REF ALT QUAL FILTER Rk RPT n_alt
#> 4   SNP_higher_path_993510  88  993510_4   A   G    .      .  1  NA     1
#> 5  SNP_higher_path_9803974  33   9803974   A   G    .      .  1  NA     1
#> 6  SNP_higher_path_9791272  31   9791272   T   G    .      .  1  NA     1
#> 18 SNP_higher_path_9393899  61 9393899_3   T   A    .      .  1  NA     1
#> 21  SNP_higher_path_927319  45    927319   A   C    .      .  1  NA     1
#>    is_biallelic is_indel .row_id
#> 4          TRUE    FALSE       4
#> 5          TRUE    FALSE       5
#> 6          TRUE    FALSE       6
#> 18         TRUE    FALSE      18
#> 21         TRUE    FALSE      21
#>   ... 5812 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=1;UR=0;CL=1;CR=0;Genome=.;Sd=.;Cluster=16400006;ClSize=2"     
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"     
#> [3] "Ty=SNP;Rk=1.0;UL=1;UR=7;CL=1;CR=42;Genome=.;Sd=.;Cluster=14628692;ClSize=2"    
#> [4] "Ty=SNP;Rk=1.0;UL=7;UR=0;CL=7;CR=0;Genome=.;Sd=.;Cluster=18375190;ClSize=2"     
#> [5] "Ty=SNP;Rk=1.0;UL=12;UR=27;CL=15;CR=27;Genome=.;Sd=.;Cluster=10461812;ClSize=23"
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
#>   ... 13 more
```
