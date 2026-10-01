# vcf_filter_adr

Correct or remove genotypes with \> % ALT/REF ratio from a VCFArrow
object

## Usage

``` r
vcf_filter_adr(
  vcf_arrow,
  mode = c("correct", "remove"),
  threshold = 0.1,
  f_invar = TRUE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- mode:

  -\> switch between 'correct' or 'remove' mode (character)

- threshold:

  -\> decimal missing threshold, default 0.1 (numeric)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function either changes genotypes or makes genotypes missing in a
VCFArrow object if they have above/below threshold normalized ADR ratio,
returning a new VCFArrow object. The ADR is calculated as ADR = ALT /
(REF + ALT). If ADR \> threshold, the genotype becomes ALT homozygous.
If ADR \< threshold, the genotype becomes REF homozygous. By default
will remove any loci that may have become invariant as the result of the
change/removal of genotypes

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_filter_adr(vcf, mode = "correct")
#> ℹ Applying ADR filter
#> ℹ Applying ADR correct (threshold = 0.1) across 1 chunk...
#> ℹ Applying invariant filter
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 10000 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 10000 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpodaNXU/arrow_vcf_adr_1ffe323db5d4 
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
