# vcf_filter_multiSNV

Subset a VCFArrow object keeping only loci with several linked SNVs

## Usage

``` r
vcf_filter_multiSNV(
  vcf_arrow,
  block_size = 10000,
  minSNV = 2,
  maxSNV = 5,
  over_max = c("remove", "truncate")
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- block_size:

  -\> size of linked SNV blocks, default 10000 bp (integer)

- minSNV:

  -\> minimum number of SNVs per block; blocks with fewer are removed,
  default 2 (integer)

- maxSNV:

  -\> maximum number of SNVs per block, default 5 (integer); Inf keeps
  every SNV of each block

- over_max:

  -\> what to do with blocks with more than maxSNV SNVs: "remove" the
  whole block or "truncate" it to its first maxSNV SNVs (by position),
  default "remove" (character)

## Value

subsetted VCFArrow object

## Details

This function subsets a VCFArrow object keeping only linked SNVs, i.e.
blocks (loci) with at least minSNV SNVs, returning a new VCFArrow
object. A block is a CHROM, or within a CHROM a window of block_size bp
counted from its first SNV; with de novo RAD data, where each locus is
its own CHROM, a block is a locus. Blocks with fewer than minSNV SNVs
are removed. Blocks with more than maxSNV SNVs are removed (over_max =
"remove") or truncated to their first maxSNV SNVs (over_max =
"truncate"); blocks with minSNV to maxSNV SNVs are kept whole. With RAD
loci and de novo SNV discovery, an unusually high number of SNVs in a
short locus often indicates misassembled paralogs, so such loci are best
removed (the default; fineRADstructure's Stacks2fineRAD.py script
likewise removes loci with more than a maximum number of SNPs). With
SNVs from whole genome data mapped to a reference, a block is an
arbitrary window and many SNVs indicate a variable region rather than an
artefact; over_max = "truncate" keeps the first, most tightly linked,
maxSNV SNVs of such blocks, which suits unphased data, where the
haplotypes of loci with several heterozygous SNVs are arbitrary, and
maxSNV = Inf keeps all SNVs, e.g. for phased data.

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
#> ℹ 0 blocks with more than 5 SNVs removed
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 8 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 8 
#>   Multiallelic:   0 
#>   Indels:         0 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> 
#> Genotypes (first chunk with variants): 144 
#>   Missing: 66.67% 
#>   Phased:  126 
#> 
#> Storage:
#>   Path: /tmp/Rtmpq4egqB/arrow_vcf_1db272a520c0 
#>   Genotype chunk files: 1 
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
#>  .row_id         FORMAT
#>     1128 GT:DP:PL:AD:HQ
#>     1129 GT:DP:PL:AD:HQ
#>     3055 GT:DP:PL:AD:HQ
#>     3056 GT:DP:PL:AD:HQ
#>     5500 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 13 more
# whole genome data: truncate blocks with many SNVs instead of removing them
vcf_filter_multiSNV(vcf, minSNV = 2, maxSNV = 5, over_max = "truncate")
#> ℹ Applying linked SNV filter
#> ℹ Retained 8 / 10000 variants (linked SNVs)
#> ℹ 0 blocks with more than 5 SNVs truncated
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 8 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 8 
#>   Multiallelic:   0 
#>   Indels:         0 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> 
#> Genotypes (first chunk with variants): 144 
#>   Missing: 66.67% 
#>   Phased:  126 
#> 
#> Storage:
#>   Path: /tmp/Rtmpq4egqB/arrow_vcf_1db272a520c0 
#>   Genotype chunk files: 1 
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
#>  .row_id         FORMAT
#>     1128 GT:DP:PL:AD:HQ
#>     1129 GT:DP:PL:AD:HQ
#>     3055 GT:DP:PL:AD:HQ
#>     3056 GT:DP:PL:AD:HQ
#>     5500 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 13 more
```
