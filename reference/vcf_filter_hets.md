# vcf_filter_hets

Remove loci above a heterozygosity threshold from a VCFArrow object

## Usage

``` r
vcf_filter_hets(vcf_arrow, threshold = 0.5)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> heterozygosity threshold; loci with a proportion of heterozygous
  called genotypes above it are removed, those at it are kept, default
  0.5 (numeric)

## Value

subsetted VCFArrow object

## Details

This function removes loci above a heterozygosity threshold from a
VCFArrow object (loci at the threshold are kept), as well as loci
without called genotypes, returning a new VCFArrow object. High
heterozygosities are indicative of potential paralogs.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_filter_hets(vcf, threshold = 0.5)
#> ℹ Applying heterozygosity filter
#> ℹ Retained 8984 / 10000 variants (heterozygosity rate <= 0.5)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 8984 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 8436 
#>   Multiallelic:   0 
#>   Indels:         548 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> 
#> Genotypes (first chunk with variants): 161712 
#>   Missing: 48.39% 
#>   Phased:  39474 
#> 
#> Storage:
#>   Path: /tmp/RtmppaKiJH/arrow_vcf_1c7d7497dabc 
#>   Genotype chunk files: 1 
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
#>   ... 8979 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=11;UR=8;CL=11;CR=27;Genome=.;Sd=.;Cluster=20022615;ClSize=3"
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=1;CL=75;CR=7;Genome=.;Sd=.;Cluster=17296599;ClSize=1"  
#> [3] "Ty=SNP;Rk=1.0;UL=20;UR=4;CL=20;CR=4;Genome=.;Sd=.;Cluster=18266172;ClSize=3" 
#> [4] "Ty=SNP;Rk=1.0;UL=1;UR=0;CL=1;CR=0;Genome=.;Sd=.;Cluster=16400006;ClSize=2"   
#> [5] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"   
#> 
#> FORMAT (first 5):
#>  .row_id         FORMAT
#>        1 GT:DP:PL:AD:HQ
#>        2 GT:DP:PL:AD:HQ
#>        3 GT:DP:PL:AD:HQ
#>        4 GT:DP:PL:AD:HQ
#>        5 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 13 more
```
