# vcf_filter_oneSNV

Subset a VCFArrow object keeping only 1 SNV per locus

## Usage

``` r
vcf_filter_oneSNV(vcf_arrow, block_size = 10000, by_quality = FALSE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- block_size:

  -\> size of linked SNV blocks, default 10000 bp (integer)

- by_quality:

  -\> keep the SNV with the highest QUAL of each block instead of the
  first, default FALSE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function subsets a VCFArrow object keeping only 1 SNV per locus,
returning a new VCFArrow object. Locus is defined as a different
chromosome or a block of the 'block_size' parameter value within a
chromosome. By default the first SNV (lowest POS) of each locus is kept,
whatever its quality. With by_quality = TRUE the SNV with the highest
QUAL is kept; SNVs without QUAL ('.') rank below all others, and ties
(or a locus without any QUAL) go to the first SNV. VCFs without QUAL
values, such as those of DiscoSnpRad, therefore give the same result as
the default, with a warning.

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
vcf_filter_oneSNV(vcf)
#> ℹ Applying unlinked SNV filter
#> ℹ Retained 9996 / 10000 variants (unlinked SNVs)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 9996 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 9309 
#>   Multiallelic:   0 
#>   Indels:         687 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 179928 
#>   Missing: 47.08% 
#>   Phased:  45540 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_31533497538a 
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
#>   ... 9991 more
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
# keep the SNV with the highest QUAL of each locus
vcf_filter_oneSNV(vcf, by_quality = TRUE)
#> ℹ Applying unlinked SNV filter
#> ! No QUAL values: keeping the first SNV of each block.
#> ℹ Retained 9996 / 10000 variants (unlinked SNVs, highest QUAL)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 9996 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 9309 
#>   Multiallelic:   0 
#>   Indels:         687 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 179928 
#>   Missing: 47.08% 
#>   Phased:  45540 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_31533497538a 
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
#>   ... 9991 more
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
