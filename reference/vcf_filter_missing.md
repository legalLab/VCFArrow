# vcf_filter_missing

Remove samples with \> % missing data from a VCFArrow object

## Usage

``` r
vcf_filter_missing(vcf_arrow, threshold = 0.5, f_invar = TRUE, verbose = TRUE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> decimal missing threshold, default 0.5 (numeric)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

- verbose:

  -\> report filtering stats, default TRUE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes samples from a VCFArrow object if they have above
threshold missing loci, returning a new VCFArrow object. By default will
remove any loci that may have become invariant as the result of the
removal of samples. By default will report removed samples, final %
missing data, and number of retained samples after sample filtering.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf |>
  vcf_filter_missingness(threshold = 0.2) |>
  vcf_filter_missing(threshold = 0.3)
#> ℹ Applying locus missingness filter
#> ℹ Retained 2317 / 10000 variants (per-variant missingness <= 0.2)
#> ℹ Applying sample missingness filter
#> ℹ Variants retained: 2317 | Samples retained: 18
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 2317 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 2317 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpOwogNL/arrow_vcf_1eff299c224c 
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
