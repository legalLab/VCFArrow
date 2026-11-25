# vcf_filter_maf

Remove loci below MAF threshold from a VCFArrow object

## Usage

``` r
vcf_filter_maf(vcf_arrow, threshold = 0.01)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> decimal MAF threshold, default 0.01 (numeric)

## Value

subsetted VCFArrow object

## Details

This function removes loci below MAF threshold from a VCFArrow object,
returning a new VCFArrow object. Loci at the threshold are kept.

The default threshold of 0.01 follows the conventional definition of a
polymorphism as a variant whose minor allele has a frequency of at least
1% in a population (Karki et al. 2015). Karki et al. note that this
cutoff is a convention rather than a biological boundary, so set
`threshold` to suit the data and the analysis.

## References

Karki R, Pandya D, Elston RC, Ferlini C (2015). Defining "mutation" and
"polymorphism" in the era of personal genomics. *BMC Medical Genomics*
8:37.
[doi:10.1186/s12920-015-0115-z](https://doi.org/10.1186/s12920-015-0115-z)

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_filter_maf(vcf, threshold = 0.05)
#> ℹ Applying MAF filter
#> ℹ Retained 8955 / 10000 variants (MAF >= 0.05)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 8955 
#>   Samples:  18 
#> 
#> Variant types:
#>   Biallelic SNVs: 8331 
#>   Multiallelic:   0 
#>   Indels:         624 
#> 
#> Groups: 4 
#> 
#> BS GS OG WA 
#>  5  5  3  5 
#> 
#> Genotypes (first chunk with variants): 161190 
#>   Missing: 50.35% 
#>   Phased:  41688 
#> 
#> Storage:
#>   Path: /tmp/Rtmpq4egqB/arrow_vcf_1db226eed1e5 
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
#>   ... 8950 more
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
