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
