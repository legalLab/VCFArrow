# vcf_filter_missingness

Remove loci above missingness threshold from a VCFArrow object

## Usage

``` r
vcf_filter_missingness(vcf_arrow, threshold = 0.1, verbose = TRUE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> decimal missingness threshold, default 0.1 (numeric)

- verbose:

  -\> flag to report total % missing data after filtering, default TRUE
  (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes loci above missingness threshold from a VCFArrow
object, returning a new VCFArrow object. Missingness is a locus focused
metric, i.e. missing data per locus. With verbose = TRUE, the number of
retained variants and the final % missing data (missing genotypes over
the retained variants and samples) are reported.

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
vcf_filter_missingness(vcf, threshold = 0.2)
#> ℹ Applying locus missingness filter
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Retained 2317 / 10000 variants (per-variant missingness <= 0.2) | Missing data: 7.68%
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 2317 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 2139 
#>   Multiallelic:   0 
#>   Indels:         178 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 41706 
#>   Missing: 7.68% 
#>   Phased:  9756 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_315316fafbb4 
#>   Genotype chunk files: 1 
#> 
#> Variants (first 5 rows):
#>                       CHROM POS        ID REF ALT QUAL FILTER Rk RPT n_alt
#> 4    SNP_higher_path_993510  88  993510_4   A   G    .      .  1  NA     1
#> 77  SNP_higher_path_7688542  44 7688542_4   T   C    .      .  1  NA     1
#> 136 SNP_higher_path_6271717  30   6271717   A   C    .      .  1  NA     1
#> 193 SNP_higher_path_4759738  42   4759738   A   G    .      .  1  NA     1
#> 195 SNP_higher_path_4663868  32   4663868   C   T    .      .  1  NA     1
#>     is_biallelic is_indel .row_id
#> 4           TRUE    FALSE       4
#> 77          TRUE    FALSE      77
#> 136         TRUE    FALSE     136
#> 193         TRUE    FALSE     193
#> 195         TRUE    FALSE     195
#>   ... 2312 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=1;UR=0;CL=1;CR=0;Genome=.;Sd=.;Cluster=16400006;ClSize=2"  
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=0;CL=3;CR=0;Genome=.;Sd=.;Cluster=13843457;ClSize=1"  
#> [3] "Ty=SNP;Rk=1.0;UL=0;UR=36;CL=0;CR=36;Genome=.;Sd=.;Cluster=17102126;ClSize=3"
#> [4] "Ty=SNP;Rk=1.0;UL=12;UR=0;CL=12;CR=0;Genome=.;Sd=.;Cluster=19240024;ClSize=1"
#> [5] "Ty=SNP;Rk=1.0;UL=2;UR=5;CL=2;CR=5;Genome=.;Sd=.;Cluster=18749679;ClSize=2"  
#> 
#> FORMAT (first 5):
#>  .row_id         FORMAT
#>        4 GT:DP:PL:AD:HQ
#>       77 GT:DP:PL:AD:HQ
#>      136 GT:DP:PL:AD:HQ
#>      193 GT:DP:PL:AD:HQ
#>      195 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 13 more
```
