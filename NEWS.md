# VCFArrow 0.10.0

* New function `vcf2fstructure()`, an alias of `vcf2structure()` that
  writes the fastSTRUCTURE format by default (`method = "F"`).
* New function `rename_vcf_samples()` renames samples, from a named vector
  (`c(old = "new")`), a data frame or a file with columns `id` and `new_id`;
  by default each sample keeps its group, and `update_groups = TRUE` takes
  the groups from a `pop` column. The genotype files are rewritten, as
  sample names are stored with every genotype.
* `vcf_filter_adr()` now changes only called heterozygotes: in `"correct"`
  mode an imbalanced heterozygote becomes the homozygote of its majority
  allele, and in `"remove"` mode it becomes missing. Previously any genotype
  with an ADR outside the threshold was affected, so `"correct"` filled in
  missing genotypes that had allele depths (AD), and `"remove"` set
  homozygotes to missing. It also reports the number of genotypes corrected
  or removed.
* `vcf_stats()` (`read_depth`) and `assess_vcf_coverage()` use the read depth
  (DP) of called genotypes only. DP is sometimes written for missing genotypes
  which lowered the mean depth of samples with much missing data and added 
  low depths to the coverage plot.
* `vcf_sub_SNVs_stratified()` now returns exactly `n_SNVs` SNVs. Each CHROM
  contributes its proportional share rounded down, and the remaining SNVs
  come from the CHROMs with the largest remainders. Previously each share was
  only rounded down, so CHROMs with few variants contributed none: with
  DiscoSnp-RAD data, where nearly every variant has its own CHROM, no SNVs
  were returned. Subsamples drawn with a given `seed` differ from earlier
  versions.
* `vcf_filter_coverage()` gains `strict` (default `TRUE`): called genotypes
  with DP below `threshold`, or without DP, are set to missing, as vcftools
  `--minDP` does, and their number is reported; no variants are removed, so
  follow with e.g. `vcf_filter_invariant()` and `vcf_filter_missingness()`.
  `strict = FALSE` gives the previous behavior: the genotypes are not
  changed, and only the variants with at least two called genotypes with
  DP >= `threshold`, polymorphic among them, are kept.
* `vcf_filter_oneSNV()` gains `by_quality` (default `FALSE`): with `TRUE` it
  keeps the SNV with the highest QUAL of each block instead of the first.
  SNVs without QUAL rank last and ties go to the first SNV, so VCFs without
  QUAL (such as DiscoSnp-RAD output) give the default result, with a warning.
* `vcf2plink_bed()` (and so `vcf2admixture()`) and `vcf2plink_ped()` with
  `allele_code = "12"` now write ALT as allele 1 and REF as allele 2, as
  PLINK does when it converts a VCF (`plink --vcf`): ALT is in column 5 of
  the .bim file and is coded 1 in the .ped file. PLINK's A1 counts and
  frequencies therefore refer to ALT. Previously allele 1 was REF.
* `assess_vcf_coverage()`: each violin now ends at the sample's lowest and
  highest plotted depth (`trim = TRUE`) instead of extending beyond them.
* `vcf_suggest_chunk_size()` now computes the memory of a genotype row from
  its columns, including the samples' FORMAT text (new argument
  `fmt_chars`) and sample names (`sample_chars`), and shares the memory
  between workers (`workers`). It previously assumed 20 bytes per genotype,
  about a quarter of the real size, and so suggested chunks about four times
  too large. `n_columns` is deprecated and ignored.
* `vcf_memory_estimate()`: the default `chunk_size` is now 50,000, as in
  `read_vcf()` (was 100,000).
* `vcf_filter_invariant()` and `vcf_filter_coverage()` (with `strict = FALSE`) 
  have a new argument `strict_poly` that sets what counts as polymorphic. 
  With `strict_poly = FALSE` (the default) a variant is polymorphic if not 
  all of its counted alleles are the same, e.g. both REF and ALT occur. 
  With `strict_poly = TRUE` it is polymorphic if not all of its counted 
  genotypes are the same, the previous behavior, which also removes variants 
  where every sample is heterozygous. The functions that remove invariant 
  variants after changing samples or genotypes (`vcf_filter_missing()`, 
  `vcf_extract_samples()`, `vcf_extract_groups()`, `vcf_filter_adr()`) take 
  the same argument, also with the default `FALSE`, so by default they now 
  keep such all-heterozygous variants.
* `show()` (printing a VCFArrow object) now reports the numbers of
  biallelic SNVs, multiallelic variants and indels, the samples per group,
  the number of chunk files, and the percentage of missing genotypes and the
  number of phased genotypes in the first chunk file that holds current
  variants. It previously always reported "Phased genotypes: FALSE", and
  after filtering its FORMAT preview showed the first variants of the source
  VCF instead of the object's own. The Arrow schema is no longer printed.
* `set_vcf_groups()` is more robust:
  * samples not in the strata file get the group `NA` with a warning, and the
    other samples are still assigned; previously a single missing sample left
    every group `NA`;
  * an id assigned to more than one group, a missing `id` or `pop` column, or
    a missing file now stops with a clear error; previously a duplicated id or
    a misnamed column silently gave `@groups` the wrong length;
  * the file is read as tab-separated when its header line has a tab (so
    fields may contain spaces or be empty), otherwise as space-separated, and
    all columns are read as text, so numeric ids keep their leading zeros.
* `vcf_filter_missing()` and `vcf_filter_missingness()` report the final
  percentage of missing genotypes over the retained variants and samples
  (with `verbose = TRUE`), as their documentation described; so do
  `vcf_extract_samples()` and `vcf_extract_groups()`. For the sample filters
  this reads the genotypes once more, which `verbose = FALSE` avoids.
* `read_vcf()` reads genotypes correctly in cases it misread before:
  * GT is found by name in FORMAT, and all its alleles are read, including
    alleles of several digits (`0/11`, `10/2`); previously the first three
    characters of the sample field were taken as the genotype, so haploid
    calls and lines without GT gave alleles taken from the next field (e.g.
    `1:5:...` was read as `1/5`) and `10/2` as `1/-1`;
  * haploid genotypes are stored as homozygous (`1` as `1/1`);
  * genotypes with one missing allele (`./1`) or more than two alleles, and
    unreadable ones, are stored as missing, and a warning reports how many
    there were, together with lines without GT or with fewer sample columns
    than samples;
  * ADR is computed at biallelic sites only, so `vcf_filter_adr()` no longer
    changes genotypes at multiallelic sites, where the first two AD values
    need not be the genotype's alleles;
  * DP, GQ and AD values in other notations (e.g. `1e2`) are read correctly,
    and negative values are treated as missing;
  * `chunk_size` must be a positive whole number.
* `write_vcf()`:
  * the INFO fields AC, AN, AF and NS are recomputed from the samples and
    genotypes written, where present (new argument `update_info = TRUE`), so
    they stay correct after samples or genotypes were removed or changed;
  * variants without a valid position are not written (with a warning)
    instead of being written with the position -2147483648;
  * warns if the variants are not sorted by position (e.g. after
    `vcf_bind(mode = "union")`);
  * compresses by default when `out_file` ends in `.gz`;
  * adds a `##source=VCFArrow` line to the header.
* POS is stored as a double, so positions beyond 2,147,483,647 (chromosomes
  longer than 2.1 Gb) are kept instead of becoming `NA`. This adds 4 bytes
  per variant to the variant tables. `vcf2genlight()` stores such positions
  as `NA` (genlight positions are integers), with a warning.
* `vcf2eigenstrat()` writes positions in the `.snp` file without scientific
  notation (e.g. `0.00004100`, not `4.1e-05`).
* `vcf_bind()` (`mode = "union"` or `outgroup = TRUE`, `absent_as =
  "missing"`): genotypes absent from a source object are written by
  `write_vcf()` as `./.` instead of `.`.
* Exporters now write the formats their programs document:
  * `vcf2bayescan()`: each line starts with the locus number, as BayeScan
    requires (it was missing, shifting every column);
  * `vcf2bayesass()`: writes BA3's native format, one line per individual
    and locus (`indivID popID locID allele1 allele2`, alleles 1 = REF,
    2 = ALT, missing `0 0`), and warns above BA3's limits;
  * `vcf2related()`: each line starts with the sample name, and missing
    genotypes are `0 0` (were `NA`);
  * `vcf2apparent()`: each line starts with the sample name, and missing
    genotypes are `-/-` (were `?/?`);
  * `vcf2fineradstructure()`: writes one line per RAD locus, defined as in
    `vcf_filter_multiSNV()` (new argument `block_size`), with each sample's
    two haplotypes over the locus's SNVs in RADpainter's syntax (`TGAT/CGAT`,
    `TGAT` if both are the same, empty if missing), the genotypes taken as
    phased; it wrote one line per SNV with allele pairs such as `A,C` and
    `NA,NA`;
  * `vcf2eigenstrat()`: the `.geno` digits count copies of the reference
    allele (REF, column 5 of `.snp`), as EIGENSTRAT defines them; they
    counted ALT copies;
  * `vcf2eigenstrat()`: new argument `chrom_code`; by default non-numeric
    CHROM (e.g. contig names) is numbered 1, 2, 3, ... in the `.snp` file,
    as EIGENSOFT needs integer chromosomes, with the mapping written to
    `<out_file>.chrommap`, and a warning gives the `numchrom` to set when
    there are more than 22 codes;
  * `vcf2migrate()` follows the Migrate-N manual's data file format:
    - the first line gives the number of loci (blocks), not of SNVs, and for
      method `"N"` starts with the data type `N` (SNPs);
    - the loci are labelled `(sK)` for method `"S"` and `(nK)` for method
      `"N"` throughout (method `"S"` labelled its last block `(nK)`, method
      `"N"` all blocks `(sK)`), and for method `"N"` the block lengths follow
      the order of the sequences (they were listed in sorted CHROM order);
    - individual names are written in a fixed field of 10 characters, with a
      warning giving the `nmlength` to set for longer names (they were
      followed by a tab, which Migrate-N would read as part of the name);
    - population names starting with a number are an error;
    - method `"C"` writes allele counts in Migrate-N's HapMap format (data
      type `H`: position, REF and ALT with their counts, and the total); its
      previous format is not one Migrate-N documents;
    - method `"C"` leaves out SNPs without called alleles in some
      population (with a warning): Migrate-N crashes on them;
  * `vcf2genepop()`: each population starts with `Pop` alone on its line, as
    GENEPOP documents it (`adegenet::read.genepop()` failed on `pop WA`);
    the new argument `pop_names = TRUE` writes `pop <group>` as before;
  * `vcf2smartsnp()`: the genotype file has no header line, as smartsnp
    requires (`smart_pca()` failed on it); the samples and their groups, in
    column order, are written to `<out_file without extension>_samples.smartsnp`
    for smartsnp's `sample_group`;
  * `vcf2plink_ped()`: new argument `allele_code`; `"12"` writes alleles as
    1 (REF) and 2 (ALT), the only `.ped` coding ADMIXTURE reads;
  * `vcf2treemix()`: compresses the file with gzip, as TreeMix expects, when
    `out_file` ends in `.gz` (new argument `gzip`);
  * `vcf2structure()`: the documentation describes the file layout and notes
    that fastStructure adds `.str` to its `--input` name;
  * locus names (GENEPOP, BA3, EIGENSTRAT, Migrate-N, genlight) are
    `CHROM_POS` wherever the ID is missing, empty or `.`; previously only
    when every ID was missing, so a VCF with `.` IDs gave `.` for all loci;
  * `vcf2genlight()`: `ploidy` other than 2 gives a warning and is set to 2,
    as the genotypes are diploid.
    * `vcf2snapp()` declares the data as `DATATYPE=INTEGERDATA SYMBOLS="012"`
    instead of `DATATYPE=DNA`: BEAST 2 reads DNA data as nucleotides, whose
    alphabet has no 0, 1 or 2.

# VCFArrow 0.9.1

* `vcf_theta()` (and `vcf_stats(theta = TRUE)`): per-site nucleotide diversity
  now uses the unbiased estimator n/(n - 1) * 2p(1 - p) for n called alleles,
  as vcftools (`--site-pi`) and tidypopgen (`loci_pi()`) do. Values of pi are
  therefore higher than in 0.9.0 by the factor n/(n - 1), most for small
  groups; Watterson's theta is unchanged.
* `vcf_stats()` reports each sample's inbreeding coefficient F, the
  method-of-moments estimate of vcftools and PLINK (`--het`), with the numbers
  of sites and of observed and expected homozygotes it is based on, and the
  mean F of each sample group (columns `F`, `F_sites`, `F_hom_obs`,
  `F_hom_exp` and `F_group_mean`). The new argument `f_ref` takes the allele
  frequencies from all samples (`"all"`, the default) or from each sample's
  group (`"group"`).

# VCFArrow 0.9.0

* Initial CRAN submission.

## Reading and storage

* `read_vcf()` reads plain or gzip-compressed VCF files (including large
  DiscoSnp-RAD output) into a `VCFArrow` S4 object. Variant information is
  kept in memory; genotypes are written to disk as an 'Apache Arrow' dataset
  in chunks and read lazily, so datasets larger than memory can be processed.
* `set_vcf_groups()` assigns samples to groups (populations); `vcf_copy()`
  makes an independent copy; `vcf_gc()` deletes the on-disk data of objects
  that are no longer referenced.
* `vcf_memory_estimate()` estimates the memory an export needs, and
  `vcf_suggest_chunk_size()` suggests a chunk size for `read_vcf()` for the
  available memory.

## Filtering, subsetting and merging

* Variant filters: `vcf_filter_quality()`, `vcf_filter_pass()`,
  `vcf_filter_biallelic()`, `vcf_filter_indels()`, `vcf_filter_coverage()`,
  `vcf_filter_maf()`, `vcf_filter_hets()`, `vcf_filter_missingness()`,
  `vcf_filter_invariant()`, `vcf_filter_oneSNV()` and `vcf_filter_multiSNV()`
  (SNV blocks), and, for DiscoSnp-RAD data, `vcf_filter_rank()`,
  `vcf_filter_rpt()` and `vcf_filter_paralogs()`.
* Sample filters and subsets: `vcf_filter_missing()`, `vcf_extract_samples()`,
  `vcf_extract_groups()`, `vcf_sub_SNVs()` and `vcf_sub_SNVs_stratified()`.
* `vcf_filter_adr()` corrects or removes genotypes with an extreme allele
  depth ratio (ALT / (REF + ALT)): they are changed to homozygous or set to
  missing.
* `vcf_bind()` merges objects with the same or different samples, matching
  variants by CHROM, POS, REF and ALT, in `"intersect"` or `"union"` mode,
  with an ingroup/outgroup mode (`outgroup = TRUE`) and optional recovery of
  previously filtered loci (`recover_loci = TRUE`). `vcf_bind_sparse()` is a
  deprecated alias of `vcf_bind(..., recover_loci = TRUE)`.

## Statistics and plots

* `vcf_stats()` computes per-sample statistics (heterozygosity, missing data,
  depth) and, optionally, Watterson's theta and nucleotide diversity per group
  (`vcf_theta()`), accumulated chunk by chunk with little memory.
* `assess_vcf_missing_data()` writes a table and plots of missing data per
  sample, and `assess_vcf_coverage()` plots of read depth per sample (PDF,
  PNG and SVG).

## Writing and conversion

* `write_vcf()` writes plain or gzip-compressed VCF files.
* Converters to the input formats of population-genetic and phylogenetic
  software: `vcf2admixture()`, `vcf2apparent()`, `vcf2arlequin()`,
  `vcf2bayesass()`, `vcf2bayescan()`, `vcf2eigenstrat()`, `vcf2fasta()`,
  `vcf2fineradstructure()`, `vcf2genepop()`, `vcf2migrate()`, `vcf2nexus()`,
  `vcf2plink_bed()`, `vcf2plink_ped()`, `vcf2related()`, `vcf2smartsnp()`,
  `vcf2snapp()`, `vcf2snmf()`, `vcf2structure()` and `vcf2treemix()`;
  `vcf2gt_long()` writes a long genotype table (CSV, feather or Parquet) and
  `vcf2genlight()` returns an 'adegenet' genlight object.

## Performance

* `vcf_set_workers()` runs the per-chunk work of reading, filtering,
  statistics, merging and exporting on background R processes. Results do
  not depend on the number of workers.
* Performance-critical steps (VCF parsing, per-sample counts, theta and pi,
  writers of the exported formats) are implemented in C++.

## Dependencies

* 'adegenet' (for `vcf2genlight()`) and 'svglite' (for SVG plots) are
  suggested rather than imported; without 'svglite', SVG plots are written
  with R's Cairo device.
