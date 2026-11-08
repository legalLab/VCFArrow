# vcf_extract_samples

Extract samples from a VCFArrow object

## Usage

``` r
vcf_extract_samples(
  vcf_arrow,
  samples,
  keep = TRUE,
  f_invar = TRUE,
  verbose = TRUE,
  strict_poly = FALSE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- samples:

  -\> individuals to retain/drop (character)

- keep:

  -\> retain (TRUE) or drop (FALSE) individuals flag, default TRUE
  (Boolean)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

- verbose:

  -\> report filtering stats, default TRUE (Boolean)

- strict_poly:

  -\> when removing invariant loci, polymorphism at the genotype level
  (TRUE) or the allele level (FALSE), see vcf_filter_invariant(),
  default FALSE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes samples from a VCFArrow object, returning a new
VCFArrow object. It uses the 'keep' flag to either keep or drop the
samples in the list. By default will remove any loci that may have
become invariant as the result of the removal of samples. By default
will report removed samples, final % missing data, and number of
retained samples after sample filtering.

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
vcf_extract_samples(vcf, samples = c("Pv14", "Pv27", "Pv28"))
#> ℹ Compacting GT: 18 -> 3 samples across 1 chunk
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> ℹ Applying invariant filter
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Retained 2516 / 10000 variants (polymorphic at the allele level)
#> ℹ Removed samples: Pv120, Pv126, Pv13, Pv2, Pv31, Pv56, Pv62, Pv68, Pv73, Pv78, Pv79, Pv93, Pb2Jp, Pb2Scx, and Pb1Rd
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Variants retained: 2516 | Samples retained: 3 | Missing data: 14.97%
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 2516 
#>   Samples:  3 
#> 
#> Variant types:
#>   Biallelic SNVs: 2290 
#>   Multiallelic:   0 
#>   Indels:         226 
#> 
#> Groups: 1 
#> 
#> GS 
#>  3 
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 7548 
#>   Missing: 14.97% 
#>   Phased:  1971 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_samp_315331e55608 
#>   Genotype chunk files: 1 
#> 
#> Variants (first 5 rows):
#>                       CHROM POS      ID REF ALT QUAL FILTER Rk RPT n_alt
#> 5   SNP_higher_path_9803974  33 9803974   A   G    .      .  1  NA     1
#> 72  SNP_higher_path_7805510  39 7805510   A   G    .      .  1  NA     1
#> 73  SNP_higher_path_7790901  42 7790901   A   G    .      .  1  NA     1
#> 97  SNP_higher_path_7171957  38 7171957   T   G    .      .  1  NA     1
#> 121 SNP_higher_path_6557883  95 6557883   C   T    .      .  1  NA     1
#>     is_biallelic is_indel .row_id
#> 5           TRUE    FALSE       5
#> 72          TRUE    FALSE      72
#> 73          TRUE    FALSE      73
#> 97          TRUE    FALSE      97
#> 121         TRUE    FALSE     121
#>   ... 2511 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"    
#> [2] "Ty=SNP;Rk=1.0;UL=9;UR=35;CL=9;CR=35;Genome=.;Sd=.;Cluster=13557893;ClSize=5"  
#> [3] "Ty=SNP;Rk=1.0;UL=12;UR=15;CL=12;CR=65;Genome=.;Sd=.;Cluster=16814114;ClSize=3"
#> [4] "Ty=SNP;Rk=1.0;UL=8;UR=5;CL=8;CR=5;Genome=.;Sd=.;Cluster=18980179;ClSize=3"    
#> [5] "Ty=SNP;Rk=1.0;UL=16;UR=4;CL=65;CR=4;Genome=.;Sd=.;Cluster=20793194;ClSize=3"  
#> 
#> FORMAT (first 5):
#>  .row_id         FORMAT
#>        5 GT:DP:PL:AD:HQ
#>       72 GT:DP:PL:AD:HQ
#>       73 GT:DP:PL:AD:HQ
#>       97 GT:DP:PL:AD:HQ
#>      121 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv14" "Pv27" "Pv28"
vcf_extract_samples(vcf, samples = c("Pb2Jp", "Pb2Scx", "Pb1Rd"),
                    keep = FALSE)
#> ℹ Compacting GT: 18 -> 15 samples across 1 chunk
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> ℹ Applying invariant filter
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Retained 8421 / 10000 variants (polymorphic at the allele level)
#> ℹ Removed samples: Pb2Jp, Pb2Scx, and Pb1Rd
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Variants retained: 8421 | Samples retained: 15 | Missing data: 40.00%
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 8421 
#>   Samples:  15 
#> 
#> Variant types:
#>   Biallelic SNVs: 7816 
#>   Multiallelic:   0 
#>   Indels:         605 
#> 
#> Groups: 3 
#> 
#> BS GS WA 
#>  5  5  5 
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 126315 
#>   Missing: 40.00% 
#>   Phased:  29610 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_samp_3153667986f5 
#>   Genotype chunk files: 1 
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
#>   ... 8416 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=3;UR=1;CL=75;CR=7;Genome=.;Sd=.;Cluster=17296599;ClSize=1" 
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"  
#> [3] "Ty=SNP;Rk=1.0;UL=1;UR=7;CL=1;CR=42;Genome=.;Sd=.;Cluster=14628692;ClSize=2" 
#> [4] "Ty=SNP;Rk=1.0;UL=13;UR=8;CL=13;CR=8;Genome=.;Sd=.;Cluster=17727234;ClSize=4"
#> [5] "Ty=SNP;Rk=1.0;UL=16;UR=0;CL=16;CR=0;Genome=.;Sd=.;Cluster=18305039;ClSize=1"
#> 
#> FORMAT (first 5):
#>  .row_id         FORMAT
#>        2 GT:DP:PL:AD:HQ
#>        5 GT:DP:PL:AD:HQ
#>        6 GT:DP:PL:AD:HQ
#>        7 GT:DP:PL:AD:HQ
#>        8 GT:DP:PL:AD:HQ
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 10 more
```
