# vcf_extract_groups

Extract groups of samples from a VCFArrow object

## Usage

``` r
vcf_extract_groups(
  vcf_arrow,
  groups,
  keep = TRUE,
  f_invar = TRUE,
  verbose = TRUE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- groups:

  -\> groups to retain/drop (character)

- keep:

  -\> retain (TRUE) or drop (FALSE) individuals flag, default TRUE
  (Boolean)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

- verbose:

  -\> report filtering stats, default TRUE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes groups of samples from a VCFArrow object,
returning a new VCFArrow object. It uses the 'keep' flag to either keep
or drop the groups in the list. By default will remove any loci that may
have become invariant as the result of the removal of samples. By
default will report removed samples, final % missing data, and number of
retained samples after sample filtering.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_extract_groups(vcf, groups = c("GS", "BS"))
#> ℹ Compacting GT: 18 -> 10 samples across 1 chunk
#> ℹ Applying invariant filter
#> ℹ Removed 4217 invariant variants; 5783 retained.
#> ℹ Removed samples: Pv120, Pv2, Pv31, Pv56, Pv93, Pb2Jp, Pb2Scx, and Pb1Rd
#> ℹ Variants retained: 5783 | Samples retained: 10
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 5783 
#>   Samples:  10 
#> 
#> Quick stats:
#>   Non-missing variants: 5783 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpSpcr9J/arrow_vcf_samp_1ff82b5a0dea 
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
#> 2  SNP_higher_path_9984432 105   9984432   C   T    .      .  1  NA     1
#> 5  SNP_higher_path_9803974  33   9803974   A   G    .      .  1  NA     1
#> 13 SNP_higher_path_9568904  69   9568904   C   T    .      .  1  NA     1
#> 14 SNP_higher_path_9516492  34 9516492_2   T   G    .      .  1  NA     1
#> 17 SNP_higher_path_9462243 105   9462243   A   G    .      .  1  NA     1
#>    is_biallelic is_indel .row_id
#> 2          TRUE    FALSE       2
#> 5          TRUE    FALSE       5
#> 13         TRUE    FALSE      13
#> 14         TRUE    FALSE      14
#> 17         TRUE    FALSE      17
#>   ... 5778 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=3;UR=1;CL=75;CR=7;Genome=.;Sd=.;Cluster=17296599;ClSize=1"   
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"    
#> [3] "Ty=SNP;Rk=1.0;UL=39;UR=0;CL=39;CR=0;Genome=.;Sd=.;Cluster=18460267;ClSize=2"  
#> [4] "Ty=SNP;Rk=1.0;UL=2;UR=0;CL=2;CR=0;Genome=.;Sd=.;Cluster=17232656;ClSize=1"    
#> [5] "Ty=SNP;Rk=1.0;UL=12;UR=35;CL=75;CR=35;Genome=.;Sd=.;Cluster=18316722;ClSize=7"
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
#> [1] "Pv126" "Pv13"  "Pv14"  "Pv27"  "Pv28" 
#>   ... 5 more
vcf_extract_groups(vcf, groups = "OG", keep = FALSE)
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
#>   Path: /tmp/RtmpSpcr9J/arrow_vcf_samp_1ff81be38557 
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
