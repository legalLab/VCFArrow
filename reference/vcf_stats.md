# vcf_stats

Calculates basic stats of each samples from VCFArrow format data.
Includes average read depth per individual, missing data per individual,
the inbreeding coefficient F, Watterson's theta and pi.

## Usage

``` r
vcf_stats(
  vcf_arrow,
  res_path,
  project,
  theta = FALSE,
  f_ref = c("all", "group")
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- res_path:

  -\> directory where to write results

- project:

  -\> base name of the project file

- theta:

  -\> flag to perform theta and pi calculation, default FALSE (Boolean)

- f_ref:

  -\> allele frequencies for the inbreeding coefficient F: "all" (all
  samples, default, as vcftools and PLINK) or "group" (the sample's own
  group) (character)

## Value

table of statistics

## Details

This function calculates average read depth, heterozygosity number of
heterozygotes, number of reference and alternative homozygotes, missing
data and total number SNPs of each sample in an VCFArrow object.
Optionally calls vcf_theta() to get total and group Watterson's theta
and pi.

The inbreeding coefficient F of each sample is the method-of-moments
estimate of vcftools –het and plink –het, F = (O - E) / (N - E). It uses
the biallelic SNVs at which the sample is called, at least two samples
of the reference set are called and the reference set is polymorphic (N
sites; monomorphic sites would not change F); O is the number of
homozygous genotypes of the sample at these sites, and E the expected
number, the sum of 1 - n / (n - 1) \* 2p(1 - p) over the sites, where p
is the ALT allele frequency and n the number of called alleles in the
reference set. With f_ref = "all" the reference set is all samples; with
f_ref = "group" it is the sample's group, which avoids the excess
homozygosity that differences in allele frequency between groups cause
when they are pooled (Wahlund effect). F is NA for samples without such
sites (e.g. a group of one sample with f_ref = "group", or a sample
without a group). The table gives F with N, O and E (F_sites, F_hom_obs,
F_hom_exp), and for each sample the mean F of its group (F_group_mean;
over the group's samples with an F).

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_stats(vcf, res_path = tempdir(), project = "vaillantii",
          theta = TRUE)
#> ℹ Computing per-sample stats: 10000 variants x 18 samples, reading 1 chunk directly
```
