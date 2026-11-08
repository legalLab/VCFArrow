# vcf_filter_coverage

Set genotypes below a read depth (DP) threshold to missing, or remove
variants without enough genotypes at or above it, in a VCFArrow object

## Usage

``` r
vcf_filter_coverage(
  vcf_arrow,
  threshold = 10,
  strict = TRUE,
  strict_poly = FALSE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> DP threshold, default 10 (integer)

- strict:

  -\> set called genotypes with DP below the threshold to missing (TRUE)
  or remove variants with fewer than two polymorphic called genotypes
  with DP \>= threshold (FALSE), default TRUE (Boolean)

- strict_poly:

  -\> with strict = FALSE, polymorphism at the genotype level (TRUE) or
  the allele level (FALSE), default FALSE (Boolean)

## Value

VCFArrow object with genotypes set to missing (strict = TRUE) or
subsetted VCFArrow object (strict = FALSE)

## Details

With strict = TRUE (the default), called genotypes with DP below the
threshold, or without DP, are set to missing (as vcftools –minDP does),
so every called genotype of the result has DP \>= threshold; the number
of genotypes set to missing is reported. No variants are removed:
variants can become invariant or have much missing data, so follow with
e.g. vcf_filter_invariant() and vcf_filter_missingness(). With strict =
FALSE the genotypes are not changed, and only the variants that have at
least two called genotypes with DP \>= threshold and are polymorphic
among them are kept; the kept variants can still contain called
genotypes with DP below the threshold. With strict_poly = FALSE (the
default), polymorphic means that not all of these genotypes' alleles are
the same, e.g. both REF and ALT alleles occur. With strict_poly = TRUE,
it means that not all of these genotypes are the same, so a variant
whose genotypes with DP \>= threshold are all heterozygous is removed
although both alleles occur. strict_poly is ignored with strict = TRUE.

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
# set genotypes with DP < 10 to missing
vcf_filter_coverage(vcf, threshold = 10)
#> ℹ Applying read coverage filter
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> ℹ Set 34935 / 95245 called genotypes with DP < 10 to missing
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 10000 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 9313 
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
#> Genotypes (first chunk with variants): 180000 
#>   Missing: 66.49% 
#>   Phased:  45612 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_coverage_31533e4a3457 
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
# keep the genotypes, remove variants with fewer than two polymorphic
# genotypes with DP >= 10
vcf_filter_coverage(vcf, threshold = 10, strict = FALSE)
#> ℹ Applying read coverage filter
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Retained 5983 / 10000 variants (polymorphic at the allele level with DP >= 10)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 5983 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 5534 
#>   Multiallelic:   0 
#>   Indels:         449 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 107694 
#>   Missing: 34.71% 
#>   Phased:  26622 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_3153643841ae 
#>   Genotype chunk files: 1 
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
#>   ... 5978 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=1;UR=0;CL=1;CR=0;Genome=.;Sd=.;Cluster=16400006;ClSize=2"     
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"     
#> [3] "Ty=SNP;Rk=1.0;UL=1;UR=7;CL=1;CR=42;Genome=.;Sd=.;Cluster=14628692;ClSize=2"    
#> [4] "Ty=SNP;Rk=1.0;UL=7;UR=0;CL=7;CR=0;Genome=.;Sd=.;Cluster=18375190;ClSize=2"     
#> [5] "Ty=SNP;Rk=1.0;UL=12;UR=27;CL=15;CR=27;Genome=.;Sd=.;Cluster=10461812;ClSize=23"
#> 
#> FORMAT (first 5):
#>  .row_id         FORMAT
#>        4 GT:DP:PL:AD:HQ
#>        5 GT:DP:PL:AD:HQ
#>        6 GT:DP:PL:AD:HQ
#>       18 GT:DP:PL:AD:HQ
#>       21 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 13 more
# also remove variants whose well-covered genotypes are all heterozygous
vcf_filter_coverage(vcf, threshold = 10, strict = FALSE, strict_poly = TRUE)
#> ℹ Applying read coverage filter
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Retained 5817 / 10000 variants (polymorphic at the genotype level with DP >= 10)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 5817 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 5388 
#>   Multiallelic:   0 
#>   Indels:         429 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 104706 
#>   Missing: 34.06% 
#>   Phased:  25506 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_3153643841ae 
#>   Genotype chunk files: 1 
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
#>  .row_id         FORMAT
#>        4 GT:DP:PL:AD:HQ
#>        5 GT:DP:PL:AD:HQ
#>        6 GT:DP:PL:AD:HQ
#>       18 GT:DP:PL:AD:HQ
#>       21 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 13 more
```
