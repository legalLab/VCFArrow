# vcf_sub_loci

Randomly subsets SNVs from a VCFArrow object

## Usage

``` r
vcf_sub_SNVs(vcf_arrow, n_SNVs = 10000, seed = NULL)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- n_SNVs:

  -\> number of SNVs to subset, default 10000 (integer)

- seed:

  -\> random number generator seed, default NULL (integer)

## Value

VCFArrow object

## Details

This function subsets a VCFArrow object to specific number of SNVs,
returning new VCFArrow object. The seed for random number generator is
automatically generated unless otherwise specified.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_sub_SNVs(vcf, n_SNVs = 1000, seed = 42)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 1000 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 1000 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpFwXb6F/arrow_vcf_1f6b313934c7 
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
#> 16 SNP_higher_path_9468236  33   9468236   A   T    .      .  1  NA     1
#> 19 SNP_higher_path_9359469  56 9359469_1   A   C    .      .  1  NA     1
#> 27 SNP_higher_path_9160465 158   9160465   C   T    .      .  1  NA     1
#> 35 SNP_higher_path_8822666  30   8822666   C   G    .      .  1  NA     1
#>    is_biallelic is_indel .row_id
#> 2          TRUE    FALSE       2
#> 16         TRUE    FALSE      16
#> 19         TRUE    FALSE      19
#> 27         TRUE    FALSE      27
#> 35         TRUE    FALSE      35
#>   ... 995 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=3;UR=1;CL=75;CR=7;Genome=.;Sd=.;Cluster=17296599;ClSize=1"   
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=22;CL=3;CR=91;Genome=.;Sd=.;Cluster=20063772;ClSize=5"  
#> [3] "Ty=SNP;Rk=1.0;UL=26;UR=12;CL=26;CR=12;Genome=.;Sd=.;Cluster=17801571;ClSize=1"
#> [4] "Ty=SNP;Rk=1.0;UL=98;UR=0;CL=128;CR=0;Genome=.;Sd=.;Cluster=19875537;ClSize=5" 
#> [5] "Ty=SNP;Rk=1.0;UL=0;UR=1;CL=0;CR=1;Genome=.;Sd=.;Cluster=17143733;ClSize=2"    
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
