# vcf_filter_paralogs

Remove loci above REPEAT and below RANK threshold from a VCFArrow object

## Usage

``` r
vcf_filter_paralogs(
  vcf_arrow,
  threshold_rk = 0.2,
  threshold_rpt = 0.5,
  keep_na = FALSE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold_rk:

  -\> decimal RANK threshold, loci below it are removed, default 0.2
  (numeric)

- threshold_rpt:

  -\> decimal REPEAT threshold, loci above it are removed, default 0.5
  (numeric)

- keep_na:

  -\> retain loci without a RANK or REPEAT value, default FALSE
  (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes loci below the RANK threshold and above the REPEAT
threshold from a VCFArrow object, returning a new VCFArrow object. It
applies vcf_filter_rank() and then vcf_filter_rpt(). RANK and REPEAT
metrics focus on different paralog signals, and used jointly maximize
true positive and minimize false positive detections. Simulations
indicate that removing loci with Rk \< 0.2 or RPT \> 0.5 gives the best
results, with Rk \< 0.4 or RPT \> 0.5 being slightly more conservative.

## Author

Tomas Hrbek September 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_filter_paralogs(vcf, threshold_rk = 0.2, threshold_rpt = 0.5)
#> ℹ Applying RANK filter
#> ℹ Retained 10000 / 10000 variants (Rank >= 0.2)
#> ! RPT field not present; returning unfiltered VCFArrow
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 10000 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 10000 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpSpcr9J/arrow_vcf_1ff8cb1bdeb 
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
#> 1 SNP_higher_path_9994239  41 9994239_1   C   G    .      .  1  NA     1
#> 2 SNP_higher_path_9984432 105   9984432   C   T    .      .  1  NA     1
#> 3 SNP_higher_path_9967574  50   9967574   A   C    .      .  1  NA     1
#> 4  SNP_higher_path_993510  88  993510_4   A   G    .      .  1  NA     1
#> 5 SNP_higher_path_9803974  33   9803974   A   G    .      .  1  NA     1
#>   is_biallelic is_indel .row_id
#> 1         TRUE    FALSE       1
#> 2         TRUE    FALSE       2
#> 3         TRUE    FALSE       3
#> 4         TRUE    FALSE       4
#> 5         TRUE    FALSE       5
#>   ... 9995 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=11;UR=8;CL=11;CR=27;Genome=.;Sd=.;Cluster=20022615;ClSize=3"
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=1;CL=75;CR=7;Genome=.;Sd=.;Cluster=17296599;ClSize=1"  
#> [3] "Ty=SNP;Rk=1.0;UL=20;UR=4;CL=20;CR=4;Genome=.;Sd=.;Cluster=18266172;ClSize=3" 
#> [4] "Ty=SNP;Rk=1.0;UL=1;UR=0;CL=1;CR=0;Genome=.;Sd=.;Cluster=16400006;ClSize=2"   
#> [5] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"   
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
