# read_vcf

Read VCF file and store its content within a VCFArrow object

## Usage

``` r
read_vcf(vcf_file, chunk_size = 50000)
```

## Arguments

- vcf_file:

  -\> VCF file

- chunk_size:

  -\> number of variants to read in at a time, default 50000 (integer)

## Value

VCFArrow object

## Details

This function read a VCF file into an S4 class object in chunks,
returning a VCFArrow object. It accepts both uncompressed and gz
compressed files. The GT field is stored as an Apache Arrow in Long
format. Genotypes are read from the GT field, wherever it is in FORMAT,
as two allele numbers (0 = REF, 1 = first ALT, ...; alleles of any
number of digits). Haploid genotypes (e.g. "1") are stored as homozygous
("1/1"). Genotypes with one missing allele (e.g. "./1"), with more than
two alleles, or that cannot be read are stored as missing, as are all
genotypes of lines without a GT field; a warning reports how many there
were, and how many lines had fewer sample columns than samples. The
allele depth ratio ADR = ALT / (REF + ALT) is computed from AD at
biallelic sites only (NA elsewhere). Various metrics are precalculated
for fast and easy filtering. GT slot content stored in a TEMP directory
for lazy loading.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f)
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
vcf
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
#> Groups: not set (use set_vcf_groups())
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> 
#> Genotypes (first chunk with variants): 180000 
#>   Missing: 47.09% 
#>   Phased:  45612 
#> 
#> Storage:
#>   Path: /tmp/RtmpVdrUAs/arrow_vcf_31533f3682c8 
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
```
