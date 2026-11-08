# vcf_sub_SNVs_stratified

Randomly subsets SNVs from a VCFArrow object

## Usage

``` r
vcf_sub_SNVs_stratified(vcf_arrow, n_SNVs = 1000, seed = NULL)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- n_SNVs:

  -\> number of SNVs to subset, default 1000 (integer)

- seed:

  -\> random number generator seed, default NULL (integer)

## Value

VCFArrow object

## Details

This function subsets a VCFArrow object to specific number of SNVs,
returning new VCFArrow object. The subsampling is stratified, i.e. the
same proportion of SNVs per CHROM. Each CHROM contributes its
proportional share of n_SNVs, rounded down; the SNVs still needed to
reach n_SNVs come one each from the CHROMs with the largest remainders
(ties broken at random), so exactly n_SNVs SNVs are returned. When most
CHROMs hold a single variant (e.g. DiscoSnp-RAD output), this amounts to
simple random sampling. The seed for random number generator is
automatically generated unless specified.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
vcf_sub_SNVs_stratified(vcf, n_SNVs = 1000, seed = 42)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 1000 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 930 
#>   Multiallelic:   0 
#>   Indels:         70 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 18000 
#>   Missing: 47.37% 
#>   Phased:  4788 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_31535b005909 
#>   Genotype chunk files: 1 
#> 
#> Variants (first 5 rows):
#>                      CHROM POS        ID REF ALT QUAL FILTER Rk RPT n_alt
#> 25 SNP_higher_path_9205279  41   9205279   T   G    .      .  1  NA     1
#> 35 SNP_higher_path_8822666  30   8822666   C   G    .      .  1  NA     1
#> 37  SNP_higher_path_869715  36    869715   A   G    .      .  1  NA     1
#> 43 SNP_higher_path_8540263  44 8540263_2   A   G    .      .  1  NA     1
#> 55 SNP_higher_path_8280643  41   8280643   A   G    .      .  1  NA     1
#>    is_biallelic is_indel .row_id
#> 25         TRUE    FALSE      25
#> 35         TRUE    FALSE      35
#> 37         TRUE    FALSE      37
#> 43         TRUE    FALSE      43
#> 55         TRUE    FALSE      55
#>   ... 995 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=11;UR=0;CL=11;CR=0;Genome=.;Sd=.;Cluster=12523492;ClSize=1"
#> [2] "Ty=SNP;Rk=1.0;UL=0;UR=1;CL=0;CR=1;Genome=.;Sd=.;Cluster=17143733;ClSize=2"  
#> [3] "Ty=SNP;Rk=1.0;UL=6;UR=4;CL=6;CR=4;Genome=.;Sd=.;Cluster=15447647;ClSize=3"  
#> [4] "Ty=SNP;Rk=1.0;UL=8;UR=0;CL=8;CR=0;Genome=.;Sd=.;Cluster=18769025;ClSize=1"  
#> [5] "Ty=SNP;Rk=1.0;UL=11;UR=0;CL=11;CR=0;Genome=.;Sd=.;Cluster=18237634;ClSize=1"
#> 
#> FORMAT (first 5):
#>  .row_id         FORMAT
#>       25 GT:DP:PL:AD:HQ
#>       35 GT:DP:PL:AD:HQ
#>       37 GT:DP:PL:AD:HQ
#>       43 GT:DP:PL:AD:HQ
#>       55 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 13 more
```
