# vcf_bind

Bind two or more VCFArrow objects into a new VCFArrow object

## Usage

``` r
vcf_bind(
  ...,
  mode = c("intersect", "union"),
  absent_as = c("missing", "hom_ref"),
  outgroup = FALSE,
  recover_loci = FALSE,
  .check = NULL
)
```

## Arguments

- ...:

  -\> a collection of VCFArrow objects

- mode:

  -\> "intersect" keeps variants present in all objects, "union" keeps
  variants present in any object, default "intersect" (character)

- absent_as:

  -\> for mode = "union" or outgroup = TRUE, how genotypes absent from a
  source object are coded, "missing" or "hom_ref", default "missing"
  (character)

- outgroup:

  -\> add the other objects (e.g. outgroups) to the first object (e.g. a
  filtered ingroup), keeping exactly the variants of the first object;
  mode and recover_loci are then not used, default FALSE (Boolean)

- recover_loci:

  -\> for mode = "union", restore the genotypes of variants an object
  removed by filtering, where its genotype files still hold them,
  instead of coding them according to absent_as, default FALSE (Boolean)

- .check:

  -\> deprecated and ignored; objects no longer need to share the same
  variants

## Value

VCFArrow object

## Details

This function binds two or more VCFArrow objects with different sets of
individuals, returning a new VCFArrow object. Sample names must be
unique across objects. The objects need not have the same variants: mode
= "intersect" keeps the variants present in all objects, mode = "union"
keeps the variants present in any object and codes the genotypes an
object lacks according to absent_as. For objects with the same variants
both modes give the same result. In mode = "union", variants an object
removed stay removed for it: a variant that one object removed by
filtering but that another object contains is included, with the
genotypes of the object that removed it coded according to absent_as.
Binding in steps and binding all objects in one call then give the same
result. With recover_loci = TRUE, the real genotypes of such variants
are restored instead, read from the genotype files of the object that
removed them (which keep every variant the object had when its files
were written). This suits joining separately filtered datasets; binding
in steps can then differ from binding in one call.

outgroup = TRUE is meant for adding outgroups to an ingroup that was
filtered on its own (so that properties of the outgroup, such as its
missingness, cannot remove ingroup variants): the result keeps exactly
the variants of the first object, and variants removed from it stay
removed. The other objects contribute their genotypes for these
variants, read from their genotype files wherever present (including
variants they removed, for example as invariant); genotypes they do not
have are coded according to absent_as.

Variants are matched by CHROM, POS, REF and ALT, so the objects may come
from different VCF files (e.g. separate sequencing runs); REF and ALT
must be written the same way in all of them (including the order of
alternative alleles). This is assured when all datasets were called
against the same reference genome assembly, which gives the same CHROM,
POS and REF for the same locus. It is generally not possible for de novo
(reference-free) SNP calling, for example with DiscoSnp, where locus
names and positions are generated anew in each calling run, so the same
locus will usually not match between separately called datasets; such
datasets should be called together in a single run instead.
Multi-allelic sites can list their ALT alleles in a different order in
different files, so it is best to normalize the VCFs first, for example
by splitting multi-allelic sites with `bcftools norm -m -any`. Variants
keep the .row_id of the first object where its source VCF has them,
other variants get new ids; the variant metadata of the first object
and, in mode = "intersect", its record of removed invariant variants are
kept, so bound objects can be bound again.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
ingroup <- vcf_extract_groups(vcf, c("GS", "BS", "WA"))
#> ℹ Compacting GT: 18 -> 15 samples across 1 chunk
#> ℹ Applying invariant filter
#> ℹ Removed 1705 invariant variants; 8295 retained.
#> ℹ Removed samples: Pb2Jp, Pb2Scx, and Pb1Rd
#> ℹ Variants retained: 8295 | Samples retained: 15
outgroup <- vcf_extract_groups(vcf, "OG")
#> ℹ Compacting GT: 18 -> 3 samples across 1 chunk
#> ℹ Applying invariant filter
#> ℹ Removed 8904 invariant variants; 1096 retained.
#> ℹ Removed samples: Pv120, Pv126, Pv13, Pv14, Pv27, Pv28, Pv2, Pv31, Pv56, Pv62, Pv68, Pv73, Pv78, Pv79, and Pv93
#> ℹ Variants retained: 1096 | Samples retained: 3
vcf_bind(ingroup, outgroup, mode = "intersect")
#> ℹ Binding 2 VCFArrow objects (intersect): 421 variants, 18 total samples.
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 421 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 421 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpOwogNL/arrow_vcf_bind_1eff6b2739ce 
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
#> Variants (first 5 rows):
#>                         CHROM POS       ID REF ALT QUAL FILTER      Rk RPT
#> 1707 SNP_higher_path_17804743  30 17804743   A   T    .      . 0.99936  NA
#> 1837 SNP_higher_path_12811050  47 12811050   A   G    .      . 0.99124  NA
#> 2088 SNP_higher_path_24014087  41 24014087   A   G    .      . 0.98306  NA
#> 2315  SNP_higher_path_6283988  47  6283988   A   G    .      . 0.97461  NA
#> 2400 SNP_higher_path_32152907 102 32152907   A   G    .      . 0.97066  NA
#>      n_alt is_biallelic is_indel .row_id
#> 1707     1         TRUE    FALSE    1707
#> 1837     1         TRUE    FALSE    1837
#> 2088     1         TRUE    FALSE    2088
#> 2315     1         TRUE    FALSE    2315
#> 2400     1         TRUE    FALSE    2400
#>   ... 416 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=0.99936;UL=0;UR=0;CL=0;CR=0;Genome=.;Sd=.;Cluster=11121508;ClSize=9"     
#> [2] "Ty=SNP;Rk=0.99124;UL=17;UR=14;CL=17;CR=14;Genome=.;Sd=.;Cluster=10708221;ClSize=8" 
#> [3] "Ty=SNP;Rk=0.98306;UL=11;UR=12;CL=11;CR=12;Genome=.;Sd=.;Cluster=18551936;ClSize=3" 
#> [4] "Ty=SNP;Rk=0.97461;UL=17;UR=18;CL=17;CR=52;Genome=.;Sd=.;Cluster=11896588;ClSize=12"
#> [5] "Ty=SNP;Rk=0.97066;UL=12;UR=0;CL=72;CR=0;Genome=.;Sd=.;Cluster=12236649;ClSize=22"  
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
vcf_bind(ingroup, outgroup, mode = "union", absent_as = "missing")
#> ℹ Binding 2 VCFArrow objects (union): 10000 variants, 18 total samples.
#> ℹ Absent genotypes will be filled as: missing
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
#>   Path: /tmp/RtmpOwogNL/arrow_vcf_bind_1eff5dc5aa11 
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
vcf_bind(vcf_filter_missingness(ingroup, 0.2), outgroup, outgroup = TRUE,
         absent_as = "missing")
#> ℹ Applying locus missingness filter
#> ℹ Retained 3005 / 8295 variants (per-variant missingness <= 0.2)
#> ℹ Binding 2 VCFArrow objects (outgroup): 3005 variants, 18 total samples.
#> ℹ Absent genotypes will be filled as: missing
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 3005 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 3005 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpOwogNL/arrow_vcf_bind_1eff56c38283 
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
#> Variants (first 5 rows):
#>                      CHROM POS      ID REF ALT QUAL FILTER Rk RPT n_alt
#> 5  SNP_higher_path_9803974  33 9803974   A   G    .      .  1  NA     1
#> 27 SNP_higher_path_9160465 158 9160465   C   T    .      .  1  NA     1
#> 44 SNP_higher_path_8507041  30 8507041   C   T    .      .  1  NA     1
#> 52  SNP_higher_path_834734  31  834734   A   T    .      .  1  NA     1
#> 54 SNP_higher_path_8329240  42 8329240   T   G    .      .  1  NA     1
#>    is_biallelic is_indel .row_id
#> 5          TRUE    FALSE       5
#> 27         TRUE    FALSE      27
#> 44         TRUE    FALSE      44
#> 52         TRUE    FALSE      52
#> 54         TRUE    FALSE      54
#>   ... 3000 more
#> 
#> INFO (first 5):
#> [1] "Ty=SNP;Rk=1.0;UL=3;UR=7;CL=3;CR=7;Genome=.;Sd=.;Cluster=10880986;ClSize=5"   
#> [2] "Ty=SNP;Rk=1.0;UL=98;UR=0;CL=128;CR=0;Genome=.;Sd=.;Cluster=19875537;ClSize=5"
#> [3] "Ty=SNP;Rk=1.0;UL=0;UR=7;CL=0;CR=7;Genome=.;Sd=.;Cluster=20284570;ClSize=3"   
#> [4] "Ty=SNP;Rk=1.0;UL=1;UR=3;CL=1;CR=3;Genome=.;Sd=.;Cluster=16426942;ClSize=4"   
#> [5] "Ty=SNP;Rk=1.0;UL=12;UR=3;CL=12;CR=3;Genome=.;Sd=.;Cluster=10901796;ClSize=10"
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
