# vcf_extract_samples

Extract samples from a VCFArrow object

## Usage

``` r
vcf_extract_samples(
  vcf_arrow,
  samples,
  keep = TRUE,
  f_invar = TRUE,
  verbose = TRUE
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
#> ✔ VCF successfully read into a VCFArrow object
vcf_extract_samples(vcf, samples = c("Pv14", "Pv27", "Pv28"))
#> ℹ Compacting GT: 18 -> 3 samples across 1 chunk
#> ℹ Applying invariant filter
#> ℹ Removed 8275 invariant variants; 1725 retained.
#> ℹ Removed samples: Pv120, Pv126, Pv13, Pv2, Pv31, Pv56, Pv62, Pv68, Pv73, Pv78, Pv79, Pv93, Pb2Jp, Pb2Scx, and Pb1Rd
#> ℹ Variants retained: 1725 | Samples retained: 3
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 1725 
#>   Samples:  3 
#> 
#> Quick stats:
#>   Non-missing variants: 1725 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpodaNXU/arrow_vcf_samp_1ffe24b47990 
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
#>   ... 1720 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"    
#> [2] "Ty=SNP;Rk=1.0;UL=9;UR=35;CL=9;CR=35;Genome=.;Sd=.;Cluster=13557893;ClSize=5"  
#> [3] "Ty=SNP;Rk=1.0;UL=12;UR=15;CL=12;CR=65;Genome=.;Sd=.;Cluster=16814114;ClSize=3"
#> [4] "Ty=SNP;Rk=1.0;UL=8;UR=5;CL=8;CR=5;Genome=.;Sd=.;Cluster=18980179;ClSize=3"    
#> [5] "Ty=SNP;Rk=1.0;UL=16;UR=4;CL=65;CR=4;Genome=.;Sd=.;Cluster=20793194;ClSize=3"  
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
#> [1] "Pv14" "Pv27" "Pv28"
vcf_extract_samples(vcf, samples = c("Pb2Jp", "Pb2Scx", "Pb1Rd"),
                    keep = FALSE)
#> ℹ Compacting GT: 18 -> 15 samples across 1 chunk
#> ℹ Applying invariant filter
#> ℹ Removed 1705 invariant variants; 8295 retained.
#> ℹ Removed samples: Pb2Jp, Pb2Scx, and Pb1Rd
#> ℹ Variants retained: 8295 | Samples retained: 15
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 8295 
#>   Samples:  15 
#> 
#> Quick stats:
#>   Non-missing variants: 8295 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpodaNXU/arrow_vcf_samp_1ffe43561978 
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
#>   ... 8290 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=3;UR=1;CL=75;CR=7;Genome=.;Sd=.;Cluster=17296599;ClSize=1" 
#> [2] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"  
#> [3] "Ty=SNP;Rk=1.0;UL=1;UR=7;CL=1;CR=42;Genome=.;Sd=.;Cluster=14628692;ClSize=2" 
#> [4] "Ty=SNP;Rk=1.0;UL=13;UR=8;CL=13;CR=8;Genome=.;Sd=.;Cluster=17727234;ClSize=4"
#> [5] "Ty=SNP;Rk=1.0;UL=16;UR=0;CL=16;CR=0;Genome=.;Sd=.;Cluster=18305039;ClSize=1"
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
#>   ... 10 more
```
