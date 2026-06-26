
<!-- README.md is generated from vignettes/VCFArrow.Rmd. Please edit that file -->

# VCFArrow <img src="man/VCFArrow_logo.png" align="right" height="138"  />

<!-- badges: start -->

# Installation

This package needs to be installed from GitHub.
`devtools::install_github("legalLab/VCFArrow")`

# Introduction

This package has six main functionalities: (1) reading in and writing
out VCF data files, including precalculation of diverse metrics used in
the filtering functions; (2) data filtering functions; (3) merging VCFs;
(4) subsetting VCFs; (5) generation of summary statistics and graphs;
and (6) data transformation and writing out of other data formats (
population genetic and phylogenetic formats) for downstream analyses.
This package is a reimagination of the `vcf2others` package. Same as
`vcf2others`, it uses the Apache Arrow library functions to manipulate
SNP data via dplyr verb equivalents. However, unlike `vcf2others`,
`VCFArrow` now reads in the VCF data file into an S4 R object
constructed around an Apache Arrow data structure. It precalculate a
number of metrics and stores them in the S4 object. The VCF genotype
matrix is stored as a long tidy table, and the data in the S4 object are
lazy loaded. Many common tasks have been delegated to shared, C++ backed
APIs and functions. These modifications and optimizations now permit
extremely fast data filtering, manipulation and transformation. It is
possible to process VCFs comprising tens of millions of SNVs and
hundreds of individuals in matter minutes with minimal RAM overhead.

# How to use the functions of this package

Following are examples of the usage of the functions of this package.

## Load VCF and associated files of individuals and groups

First we need to read in the VCF. If one has information on the grouping
of individuals, this information can also be put into the VCFArrow
object. This grouping information is stored in a ‘strata’ file (a
tab-separated text file with two columns with an ‘id’ and ‘pop’
headers). The grouping of individuals can be used in filtering steps,
and it is necessary for transformation of the VCFArrow object to other
population genetic formats, many of which require it.

``` r
library(VCFArrow)
library(dplyr)

# set path to example files and project name
data_path <- paste0(system.file("extdata", package="VCFArrow"), "/")
project <- "crocs_"
postfix <- "discosnp_sub"

# load vcf and assign individuals to groups based on 'strata'
vcf <- read_vcf(paste0(data_path, project, postfix, ".vcf.gz")) |>
  set_vcf_groups(data_path)
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object

# visualize the VCFArrow object
vcf
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 10000 
#>   Samples:  17 
#> 
#> Quick stats:
#>   Non-missing variants: 10000 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/Rtmp486OmM/arrow_vcf_2e60338573d9d 
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
#>                    CHROM POS       ID REF ALT QUAL FILTER Rk n_alt is_biallelic is_indel
#> 1 SNP_higher_path_999984  79 999984_1   A   G    .      .  1     1         TRUE    FALSE
#> 2 SNP_higher_path_999984  90 999984_2   T   C    .      .  1     1         TRUE    FALSE
#> 3 SNP_higher_path_999984  92 999984_3   C   T    .      .  1     1         TRUE    FALSE
#> 4 SNP_higher_path_999976  47   999976   C   T    .      .  1     1         TRUE    FALSE
#> 5 SNP_higher_path_999973 132   999973   A   T    .      .  1     1         TRUE    FALSE
#>   .row_id
#> 1       1
#> 2       2
#> 3       3
#> 4       4
#> 5       5
#>   ... 9995 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=49;UR=142;CL=49;CR=142;Genome=.;Sd=.;Cluster=3792743;ClSize=1"
#> [2] "Ty=SNP;Rk=1.0;UL=49;UR=142;CL=49;CR=142;Genome=.;Sd=.;Cluster=3792743;ClSize=1"
#> [3] "Ty=SNP;Rk=1.0;UL=49;UR=142;CL=49;CR=142;Genome=.;Sd=.;Cluster=3792743;ClSize=1"
#> [4] "Ty=SNP;Rk=1.0;UL=17;UR=5;CL=17;CR=88;Genome=.;Sd=.;Cluster=3790187;ClSize=1"   
#> [5] "Ty=SNP;Rk=1.0;UL=42;UR=6;CL=102;CR=19;Genome=.;Sd=.;Cluster=3876718;ClSize=4"  
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
#> [1] "CTGA_H4639" "CTGA_H4640" "CTGA_H4641" "CTGA_H4644" "CTGA_H4646"
#>   ... 12 more

# read individuals to include
# if indivs is blank, the default is to use all individuals
indivs <- read.table(paste0(data_path, "indivs_b"), header = TRUE)$id
  
# check if all samples are in VCF sample names
if (any(!(indivs %in% vcf@samples))) stop(paste("Some individuals in list not in VCF"))
```

## Assessment of Missing Data

The function `assess_vcf_missing_data()` is used to generate a table and
graph of missing data for each sample in a VCFArrow object. It is useful
for visualizing data before and after filtering to evaluate the effect
of filtering. The function accepts parameters which form part of the
name of the output file. The idea is that the output file name contains
information on the name of the project (project), how the VCF was
extracted (postfix), and how it was filtered (fltr). The “postfix” and
“fltr” can be left blank if the VCF name does not contain this
information. Samples are ordered and colored by group assignment. The
‘details’ flag is species whether or not filtering information is
presented in the figure or if only the species name is reported. The
default is to report species name only.

``` r
# define results path
res_path <- data_path
# species for plot title
species <- "Paleosuchus/Caiman"
# filter - filter parameters in file name
fltr <- ""

assess_vcf_missing_data(vcf, res_path, species, paste0(project, postfix, fltr), details = TRUE)
```

## Basic stats for individuals in VCF

The function `vcf_stats()` is used to generate a table with basic stats
for each sample in the VCFArrow object. These stats include per sample
average read depth, heterozygosity, number of heterozygotes,
homozygotes, REF homozygotes and ALT homogygotes, percent missing SNVs,
total missing SNVs, total non-missing SNVs and total SNVs, Watterson’s
Theta and Pi for the entire dataset, and per group Watterson’s Theta and
Pi, plus the number of individuals in each group.

``` r
# get a table of basic sample stats

vcf_stats(vcf, res_path, paste0(project, postfix, fltr))
#> ℹ Computing per-sample stats: 10000 variants x 17 samples, reading 1 chunk(s) directly
#> ℹ Accumulating theta/pi: 10000 variants x 3 pops (0 MiB raw storage, vs 0 MiB with integer matrices)
```

## Filtering, subsetting, merging and otherwise wrangling VCF files

Most but not all functions can be used by the user directly on the
`VCFArrow` object. Some functions, such as `vcf_filter_rows()` and
`vcf_filter_columns()` are common APIs used by other functions to
perform filtering. Other functions, such as `vcf_filter_invariant()`
will remove invariant SNPs from a VCFArrow object; however, VCF by
definition should not have invariant SNPs. So this function is primarily
called by other functions to remove loci that may have become invariant
following filtering.

### Subsetting samples and variants

There are two types of subsetting functions. The first type is used to
subset the VCFArrow object by specific samples `vcf_extract_samples()`
or groups (all samples of a group) of samples `vcf_extract_groups()`.
The ‘keep’ flag in both functions determines whether the indicated
samples should be extracted from the VCF (‘keep = TRUE’) or excluded
from the VCF (‘keep = FALSE’). Once subsetted, the invariant SNVs are
removed by default, however, this behavior is controlled by the
‘f_invar’ flag. The second type subsets a random number of SNVs from the
VCF. The subsetting is done either across the entire VCF
`vcf_sub_SNVs()` or proportionately within each chromosome/linkage group
`vcf_sub_SNVs_stratified()`. For repeatability, both functions accept a
‘seed’; if ‘seed’ is not specified, it is randomly generated.

``` r
# subset VCF by individuals (default keep = TRUE, remove invariants f_invar = TRUE)
indivs1 <- c("CTGA_H4635", "CTGA_H4667")
vcf1 <- vcf_extract_samples(vcf, indivs1)
#> ℹ Applying invariant filter
#> ℹ Removed 9205 invariant variants; 795 retained.
#> ℹ Removed samples: CTGA_H4639, CTGA_H4640, CTGA_H4641, CTGA_H4644, CTGA_H4646, CTGA_H4661, CTGA_H4662, CTGA_H4663, CTGA_H4666, CTGA_H4668, CTGA_H4669, CTGA_H4637, CTGA_H4643, CTGA_H4645, and CTGA_H4647
#> ℹ Variants retained: 795 | Samples retained: 2

# subset VCF by groups (default keep = TRUE, remove invariants f_invar = TRUE)
groups1 <- c("paleosuchus_1", "paleosuchus_2")
vcf1 <- vcf_extract_groups(vcf, groups1)
#> ℹ Applying invariant filter
#> ℹ Removed 9627 invariant variants; 373 retained.
#> ℹ Removed samples: CTGA_H4635, CTGA_H4637, CTGA_H4643, CTGA_H4645, and CTGA_H4647
#> ℹ Variants retained: 373 | Samples retained: 12

# subset VCF by a random  number of variants (default n_SNVs = 1000)
vcf1 <- vcf_sub_SNVs(vcf, n_SNVs = 5000)

# subset VCF by a random  number of variants proportional to chromosome length (default n_SNVs = 1000)
vcf1 <- vcf_sub_SNVs_stratified(vcf, n_SNVs = 5000)
```

### Filtering

Filtering consists of removing samples and variants based on specific
properties. The function `vcf_filter_missing()` is the only function
that filters samples based on % missing SNV threshold. The other
functions filter on properties of individuals SNVs across samples The
function `vcf_filter_missingness()` removes SNVs above a % missing
threshold, `vcf_filter_quality()` removes SNVs below a quality
threshold, `vcf_filter_pass()` removes SNVs that did not PASS filters,
`vcf_filter_maf()` removes SNVs with minor allele frequency below a %
threshold, `vcf_filter_hets()` removes SNVs with heterozigosity above a
% threshold, `vcf_filter_rank()` removes SNPs with rank below a %
threshold (rank is specific to VCFs generated by DiscoSnpRad and
measures the likelihood of observed SNP variation being due to
paralogs), `vcf_filter_adr()` either removes (make genotype missing) or
corrects (make genotype homozygous) genotypes with allele ratio
imbalance, `vcf_filter_coverage()` removes SNVs below a specific
coverage (make genotypes missing), and `vcf_filter_biallelic()` removes
all variants that are not biallelic. Finally the
`vcf_filter_invariant()` removes invariant SNVs (SNVs may become
invariant after subsetting and filtering). The function
`vcf_filter_indels()` removes indels from the VCF. The functions
`vcf_filter_oneSNV()` and `vcf_filter_multiSNV()` filter the VCF to
remove linked SNVs (retain only unliked SNPs) or retain only linked
SNVs, respectively. In general, all analyses assume that SNVs are
unlinked, however, some analyses such as FineRadStructure work with
linked SNVs.

``` r
# filter VCF for analyses (unlinked SNPs)
vcf_oneSNP <- vcf_extract_samples(vcf, indivs) %>%
  vcf_filter_missing(.8) |>
  vcf_filter_maf(.03) |>
  vcf_filter_coverage(6) |>
  vcf_filter_oneSNV() |>
  vcf_filter_missingness(.2) |>
  vcf_filter_missing(.3)
#> ℹ No samples to keep or remove - keeping all samples
#> ℹ Applying sample missingness filter
#> ℹ Applying invariant filter
#> ℹ Removed 73 invariant variants; 9927 retained.
#> ℹ Removed samples: CTGA_H4663 and CTGA_H4668
#> ℹ Variants retained: 9927 | Samples retained: 15
#> ℹ Applying MAF filter
#> ℹ Retained 9927 / 9927 variants (MAF >= 0.03)
#> ℹ Applying read coverage filter
#> ℹ Retained 6548 / 9927 variants (polymorphic with DP >= 6)
#> ℹ Applying unlinked SNV filter
#> ℹ Applying locus missingness filter
#> ℹ Retained 605 / 3833 variants (per-variant missingness <= 0.2)
#> ℹ Applying sample missingness filter
#> ℹ Applying invariant filter
#> ℹ Removed 2 invariant variants; 603 retained.
#> ℹ Removed samples: CTGA_H4667 and CTGA_H4647
#> ℹ Variants retained: 603 | Samples retained: 13

# filter VCF for analyses (linked SNPs)
vcf_multiSNP <- vcf_extract_samples(vcf, indivs) |>
  vcf_filter_missing(.8) |>
  vcf_filter_maf(.03) |>
  vcf_filter_coverage(6) |>
  vcf_filter_multiSNV() |>
  vcf_filter_missingness(.2) |>
  vcf_filter_missing(.3)
#> ℹ No samples to keep or remove - keeping all samples
#> ℹ Applying sample missingness filter
#> ℹ Applying invariant filter
#> ℹ Removed 73 invariant variants; 9927 retained.
#> ℹ Removed samples: CTGA_H4663 and CTGA_H4668
#> ℹ Variants retained: 9927 | Samples retained: 15
#> ℹ Applying MAF filter
#> ℹ Retained 9927 / 9927 variants (MAF >= 0.03)
#> ℹ Applying read coverage filter
#> ℹ Retained 6548 / 9927 variants (polymorphic with DP >= 6)
#> ℹ Applying linked SNV filter
#> ℹ Applying locus missingness filter
#> ℹ Retained 642 / 4359 variants (per-variant missingness <= 0.2)
#> ℹ Applying sample missingness filter
#> ℹ Applying invariant filter
#> ℹ Removed samples: CTGA_H4667 and CTGA_H4669
#> ℹ Variants retained: 642 | Samples retained: 13
```

### Extracting, merging and adding

The functions `vcf__merge_bind()`, `vcf_bind_sparse()`,
`vcf_extract_samples()` and `vcf_extract_groups()` are used to merge two
or more VCFArrow objects, and to extract samples/groups of samples from
a VCFArrow objects storing them in a separate VCFArrow object.

``` r
# extract individuals from VCFArrow object (keep all loci)
indivs1 <- c("CTGA_H4635", "CTGA_H4637", "CTGA_H4643", "CTGA_H4645", "CTGA_H4647")
vcf_outgrp <- vcf_extract_samples(vcf, indivs1, f_invar = FALSE)
#> ℹ Removed samples: CTGA_H4639, CTGA_H4640, CTGA_H4641, CTGA_H4644, CTGA_H4646, CTGA_H4661, CTGA_H4662, CTGA_H4663, CTGA_H4666, CTGA_H4667, CTGA_H4668, and CTGA_H4669
#> ℹ Variants retained: 10000 | Samples retained: 5
vcf_ingrp <- vcf_extract_samples(vcf, indivs1, keep = FALSE, f_invar = FALSE)
#> ℹ Removed samples: CTGA_H4635, CTGA_H4637, CTGA_H4643, CTGA_H4645, and CTGA_H4647
#> ℹ Variants retained: 10000 | Samples retained: 12

# extract groups of individuals from VCFArrow object (keep all loci)
groups1 <- "caiman"
vcf_outgrp <- vcf_extract_groups(vcf, groups1, f_invar = FALSE)
#> ℹ Removed samples: CTGA_H4639, CTGA_H4640, CTGA_H4641, CTGA_H4644, CTGA_H4646, CTGA_H4661, CTGA_H4662, CTGA_H4663, CTGA_H4666, CTGA_H4667, CTGA_H4668, and CTGA_H4669
#> ℹ Variants retained: 10000 | Samples retained: 5
vcf_ingrp <- vcf_extract_groups(vcf, groups1, keep = FALSE, f_invar = FALSE)
#> ℹ Removed samples: CTGA_H4635, CTGA_H4637, CTGA_H4643, CTGA_H4645, and CTGA_H4647
#> ℹ Variants retained: 10000 | Samples retained: 12

# bind vcf_outgrp and vcf_ingrp
vcf1 <- vcf_bind(vcf_ingrp, vcf_outgrp)
```

### Filtering with outgroup

Outgroups, which tend to have many fewer individuals than ingroups, also
have more missing data than ingroups due to phylogenetic effects.
Filtering the entire dataset often results in the removal of the
outgroup taxa when `vcf_filter_missing()` is used during filtering to
remove individuals with % missing data above some threshold.The solution
is to extract the outgroup taxa with `vcf_extract_samples()` or
`vcf_extract_groups()`, filter the ingroup VCFArrow object, and add the
outgroup taxa with `vcf_bind()`. It is important to NOT filter invariant
loci during the filtering of the ingroup, and only filter invariants
from the final datasets after the outgroups have been bound.

``` r
# extract outgroup from VCFArrow (keep all loci)
groups1 <- "caiman"
vcf_outgrp <- vcf_extract_groups(vcf, groups1, f_invar = FALSE)
#> ℹ Removed samples: CTGA_H4639, CTGA_H4640, CTGA_H4641, CTGA_H4644, CTGA_H4646, CTGA_H4661, CTGA_H4662, CTGA_H4663, CTGA_H4666, CTGA_H4667, CTGA_H4668, and CTGA_H4669
#> ℹ Variants retained: 10000 | Samples retained: 5

# filter ingroup VCFArrow for filtering then bind outgroups
vcf1 <- vcf_extract_groups(vcf, groups1, keep = FALSE, f_invar = FALSE) |>
  vcf_filter_missing(.8, f_invar = FALSE) |>
  vcf_filter_maf(.03) |>
  vcf_filter_coverage(6) |>
  vcf_filter_oneSNV() |>
  vcf_filter_missingness(.2) |>
  vcf_filter_missing(.3, f_invar = FALSE) |>
  vcf_bind_sparse(vcf_outgrp) |>
  vcf_filter_invariant()
#> ℹ Removed samples: CTGA_H4635, CTGA_H4637, CTGA_H4643, CTGA_H4645, and CTGA_H4647
#> ℹ Variants retained: 10000 | Samples retained: 12
#> ℹ Applying sample missingness filter
#> ℹ Removed samples: CTGA_H4663 and CTGA_H4668
#> ℹ Variants retained: 10000 | Samples retained: 10
#> ℹ Applying MAF filter
#> ℹ Retained 328 / 10000 variants (MAF >= 0.03)
#> ℹ Applying read coverage filter
#> ℹ Retained 94 / 328 variants (polymorphic with DP >= 6)
#> ℹ Applying unlinked SNV filter
#> ℹ Applying locus missingness filter
#> ℹ Retained 8 / 86 variants (per-variant missingness <= 0.2)
#> ℹ Applying sample missingness filter
#> ℹ Removed samples: CTGA_H4667
#> ℹ Variants retained: 8 | Samples retained: 9
#> ℹ Binding 2 VCFArrow objects: 8 common variants, 14 total samples.
#> ℹ Applying invariant filter
#> Warning: Invalid metadata$r
```

## Converting a VCF file to other population genetic and phylogenetic formats

All the following functions will take a VCFArrow object, and convert it
other populations genetic and phylogenetic formats, exporting/writing a
file of this format. Group/population information is extracted from the
VCFArrow object for those formats that require it. The function
`vcf2genlight()` automatically returns a genlight object and optionally
can also export/write a file of this format to the working directory.
Generally the `vcf2genlight()` function is called within a script using
functions of the `adegenet` and `poppr` packages rather than importing
the genlight object.

``` r
res_path <- data_path
project <- "trigonatus_"
postfix <- "discosnp_sub"
fltr <- ""

##########
# datasets for analyses with unlinked SNPs
vcf <- vcf_oneSNP

##########
# export data formats
# migrate-n https://peterbeerli.com/migrate-html5/
vcf2migrate(vcf, out_file = paste0(res_path, project, postfix, fltr, '_migrate.txt'))
#> ℹ Accumulating Migrate-N (S): 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing Migrate-N (S): 603 variants x 13 samples, 7 blocks...
#> ✔ Migrate-N file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub_migrate.txt'
# arlequin http://cmpg.unibe.ch/software/arlequin35/
vcf2arlequin(vcf, out_file = paste0(res_path, project, postfix, fltr, '.arp'))
#> ℹ Accumulating Arlequin: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing Arlequin file...
#> ✔ Arlequin file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.arp'
# structure https://web.stanford.edu/group/pritchardlab/structure.html
vcf2structure(vcf, out_file = paste0(res_path, project, postfix, fltr, '.str'))
#> ℹ Accumulating Structure: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing Structure file...
#> ✔ Structure file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.str'
# faststucture http://rajanil.github.io/fastStructure/
vcf2structure(vcf, out_file = paste0(res_path, project, postfix, fltr, '.fstr'), method = "F")
#> ℹ Accumulating Structure: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing Structure file...
#> ✔ Structure file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.fstr'
# genepop https://gitlab.mbb.univ-montp2.fr/francois/genepop
vcf2genepop(vcf, out_file = paste0(res_path, project, postfix, fltr, '.gen'))
#> ℹ Accumulating Genepop: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing Genepop file...
#> ✔ Genepop file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.gen'
# smartsnp https://github.com/ChristianHuber/smartsnp
vcf2smartsnp(vcf, out_file = paste0(res_path, project, postfix, fltr, '.smartsnp'))
#> ℹ Writing SmartSNP: 603 variants x 13 samples
#> ✔ SmartSNP file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.smartsnp'
# eigenstrat https://github.com/DReichLab/EIG/tree/master
vcf2eigenstrat(vcf, out_file = paste0(res_path, project, postfix, fltr, '_eigenstrat'))
#> ℹ Writing EIGENSTRAT .geno: 603 variants x 13 samples...
#> ✔ EIGENSTRAT files written: '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub_eigenstrat.geno', '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub_eigenstrat.ind', '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub_eigenstrat.snp'
# bayescan https://github.com/mfoll/BayeScan
vcf2bayescan(vcf, out_file = paste0(res_path, project, postfix, fltr, '.bayescan'))
#> ℹ Accumulating BayesScan: 603 variants x 3 pops (0 MiB raw storage, vs 0 MiB with integer matrices)
#> ℹ Writing BayesScan file...
#> ✔ BayesScan file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.bayescan'
# bayesass https://github.com/brannala/BA3
vcf2bayesass(vcf, out_file = paste0(res_path, project, postfix, fltr, '.bayesass'))
#> ℹ Accumulating BayesAss: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing BayesAss file...
#> ✔ BayesAss file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.bayesass'
# treemix https://bitbucket.org/nygcresearch/treemix/wiki/Home
vcf2treemix(vcf, out_file = paste0(res_path, project, postfix, fltr, '.treemix'))
#> ℹ Writing Treemix: 603 variants x 3 pops
#> ✔ Treemix file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.treemix'
# apparent https://github.com/halelab/apparent/tree/master
vcf2apparent(vcf, out_file = paste0(res_path, project, postfix, fltr, '.apparent'))
#> ℹ Accumulating Apparent: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing Apparent file...
#> ✔ Apparent file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.apparent'
# related https://github.com/timothyfrasier/related
vcf2related(vcf, out_file = paste0(res_path, project, postfix, fltr, '.related'))
#> ℹ Accumulating Related: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing Related file...
#> ✔ Related file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.related'
# long tidy dataframe of genotypes
vcf2gt_long(vcf, out_file = paste0(res_path, project, postfix, fltr, '.csv'), format = 'csv')
#> ℹ Exporting gt_long: 603 variants x 13 samples -> '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.csv' (csv)
#> ℹ Combining and writing...
#> ✔ gt_long table written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.csv'
# snapp https://www.beast2.org/snapp/
vcf2snapp(vcf, out_file = paste0(res_path, project, postfix, fltr, '_snapp.nex'))
#> ℹ Accumulating SNAPP: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing SNAPP file...
#> ✔ SNAPP file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub_snapp.nex'
# nexus - only SNPs, meant for SDVq analyses https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/
vcf2nexus(vcf, out_file = paste0(res_path, project, postfix, fltr, '_sdvq.nex'))
#> ℹ Accumulating Nexus: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing Nexus file...
#> ✔ Nexus file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub_sdvq.nex'
# fasta https://www.ncbi.nlm.nih.gov/genbank/fastaformat/
vcf2fasta(vcf, out_file = paste0(res_path, project, postfix, fltr, '.fna'))
#> ℹ Accumulating FASTA: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Writing FASTA file...
#> ✔ FASTA file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.fna'

# genlight object https://www.rdocumentation.org/packages/adegenet/versions/2.0.0/topics/genlight-class
genlight <- vcf2genlight(vcf)
#> ℹ Accumulating Genlight: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Building genlight object...
# genlight object with an optional save
genlight <- vcf2genlight(vcf, out_file = paste0(res_path, project, postfix, fltr, '_genlight.rds'), save = TRUE)
#> ℹ Accumulating Genlight: 603 variants x 13 samples (0 MiB raw storage)
#> ℹ Building genlight object...
#> ✔ Genlight object saved to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub_genlight.rds'

##########
# datasets for analyses with linked SNPs
vcf <- vcf_multiSNP

##########
# export data formats
# fineRadStructure - expects VCF of linked SNPs
vcf2fineradstructure(vcf, out_file = paste0(res_path, project, postfix, fltr, '.finerad'))
#> ℹ Writing fineRADstructure: 642 variants x 13 samples
#> ✔ fineRADstructure file written to '/home/tomas/git/legal_public/packages/VCFArrow/inst/extdata/trigonatus_discosnp_sub.finerad'
```
