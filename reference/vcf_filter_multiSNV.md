# vcf_filter_multiSNV

Subset a VCFArrow object keeping only loci with 2+ SNVs per locus

## Usage

``` r
vcf_filter_multiSNV(vcf_arrow, block_size = 10000, minSNV = 2, maxSNV = 5)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- block_size:

  -\> size of linked SNV blocks, default 10000 bp (integer)

- minSNV:

  -\> minimum linked block size, default 2 (integer)

- maxSNV:

  -\> maximum number of selected linked SNVs per block, default 5
  (integer)

## Value

subsetted VCFArrow object

## Details

This function subsets a VCFArrow object keeping only loci with between
min and max \# of SNVs per locus, returning a new VCFArrow object.
Default min = 2 and max = 5 SNVs per locus (recommended as input for
fineRADstructure analyses). Locus is defined as a different chromosome
or a block of the 'block_size' parameter value within a chromosome.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_filter_multiSNV(vcf, block_size = 10000, minSNV = 2, maxSNV = 5)
#> ℹ Applying linked SNV filter
#> ℹ Retained 8 / 10000 variants (linked SNVs)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 8 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 8 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpSpcr9J/arrow_vcf_1ff841ab4118 
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
#>                         CHROM POS         ID REF ALT QUAL FILTER      Rk RPT
#> 1128 SNP_higher_path_22880229  30 22880229_1   T   G    .      . 1.00000  NA
#> 1129 SNP_higher_path_22880229  34 22880229_3   T   C    .      . 1.00000  NA
#> 3055 SNP_higher_path_35991569  33 35991569_2   T   C    .      . 0.93928  NA
#> 3056 SNP_higher_path_35991569  47 35991569_3   G   T    .      . 0.93928  NA
#> 5500  SNP_higher_path_4000943  65  4000943_4   T   G    .      . 0.77253  NA
#>      n_alt is_biallelic is_indel .row_id
#> 1128     1         TRUE    FALSE    1128
#> 1129     1         TRUE    FALSE    1129
#> 3055     1         TRUE    FALSE    3055
#> 3056     1         TRUE    FALSE    3056
#> 5500     1         TRUE    FALSE    5500
#>   ... 3 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=0;UR=0;CL=0;CR=0;Genome=.;Sd=.;Cluster=14533803;ClSize=1"     
#> [2] "Ty=SNP;Rk=1.0;UL=0;UR=0;CL=0;CR=0;Genome=.;Sd=.;Cluster=14533803;ClSize=1"     
#> [3] "Ty=SNP;Rk=0.93928;UL=2;UR=3;CL=2;CR=3;Genome=.;Sd=.;Cluster=3995598;ClSize=5"  
#> [4] "Ty=SNP;Rk=0.93928;UL=2;UR=3;CL=2;CR=3;Genome=.;Sd=.;Cluster=3995598;ClSize=5"  
#> [5] "Ty=SNP;Rk=0.77253;UL=25;UR=0;CL=25;CR=0;Genome=.;Sd=.;Cluster=8839140;ClSize=1"
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
