# Changelog

## VCFArrow 0.9.0

- Initial CRAN submission.

### Reading and storage

- [`read_vcf()`](https://legallab.github.io/VCFArrow/reference/read_vcf.md)
  reads plain or gzip-compressed VCF files (including large DiscoSnp-RAD
  output) into a `VCFArrow` S4 object. Variant information is kept in
  memory; genotypes are written to disk as an ‘Apache Arrow’ dataset in
  chunks and read lazily, so datasets larger than memory can be
  processed.
- [`set_vcf_groups()`](https://legallab.github.io/VCFArrow/reference/set_vcf_groups.md)
  assigns samples to groups (populations);
  [`vcf_copy()`](https://legallab.github.io/VCFArrow/reference/vcf_copy.md)
  makes an independent copy;
  [`vcf_gc()`](https://legallab.github.io/VCFArrow/reference/vcf_gc.md)
  deletes the on-disk data of objects that are no longer referenced.
- [`vcf_memory_estimate()`](https://legallab.github.io/VCFArrow/reference/vcf_memory_estimate.md)
  estimates the memory an export needs, and
  [`vcf_suggest_chunk_size()`](https://legallab.github.io/VCFArrow/reference/vcf_suggest_chunk_size.md)
  suggests a chunk size for
  [`read_vcf()`](https://legallab.github.io/VCFArrow/reference/read_vcf.md)
  for the available memory.

### Filtering, subsetting and merging

- Variant filters:
  [`vcf_filter_quality()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_quality.md),
  [`vcf_filter_pass()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_pass.md),
  [`vcf_filter_biallelic()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_biallelic.md),
  [`vcf_filter_indels()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_indels.md),
  [`vcf_filter_coverage()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_coverage.md),
  [`vcf_filter_maf()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_maf.md),
  [`vcf_filter_hets()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_hets.md),
  [`vcf_filter_missingness()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_missingness.md),
  [`vcf_filter_invariant()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_invariant.md),
  [`vcf_filter_oneSNV()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_oneSNV.md)
  and
  [`vcf_filter_multiSNV()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_multiSNV.md)
  (SNV blocks), and, for DiscoSnp-RAD data,
  [`vcf_filter_rank()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_rank.md),
  [`vcf_filter_rpt()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_rpt.md)
  and
  [`vcf_filter_paralogs()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_paralogs.md).
- Sample filters and subsets:
  [`vcf_filter_missing()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_missing.md),
  [`vcf_extract_samples()`](https://legallab.github.io/VCFArrow/reference/vcf_extract_samples.md),
  [`vcf_extract_groups()`](https://legallab.github.io/VCFArrow/reference/vcf_extract_groups.md),
  [`vcf_sub_SNVs()`](https://legallab.github.io/VCFArrow/reference/vcf_sub_SNVs.md)
  and
  [`vcf_sub_SNVs_stratified()`](https://legallab.github.io/VCFArrow/reference/vcf_sub_SNVs_stratified.md).
- [`vcf_filter_adr()`](https://legallab.github.io/VCFArrow/reference/vcf_filter_adr.md)
  corrects or removes genotypes with an extreme allele depth ratio (ALT
  / (REF + ALT)): they are changed to homozygous or set to missing.
- [`vcf_bind()`](https://legallab.github.io/VCFArrow/reference/vcf_bind.md)
  merges objects with the same or different samples, matching variants
  by CHROM, POS, REF and ALT, in `"intersect"` or `"union"` mode, with
  an ingroup/outgroup mode (`outgroup = TRUE`) and optional recovery of
  previously filtered loci (`recover_loci = TRUE`).
  [`vcf_bind_sparse()`](https://legallab.github.io/VCFArrow/reference/vcf_bind_sparse.md)
  is a deprecated alias of `vcf_bind(..., recover_loci = TRUE)`.

### Statistics and plots

- [`vcf_stats()`](https://legallab.github.io/VCFArrow/reference/vcf_stats.md)
  computes per-sample statistics (heterozygosity, missing data, depth)
  and, optionally, Watterson’s theta and nucleotide diversity per group
  ([`vcf_theta()`](https://legallab.github.io/VCFArrow/reference/vcf_theta.md)),
  accumulated chunk by chunk with little memory.
- [`assess_vcf_missing_data()`](https://legallab.github.io/VCFArrow/reference/assess_vcf_missing_data.md)
  writes a table and plots of missing data per sample, and
  [`assess_vcf_coverage()`](https://legallab.github.io/VCFArrow/reference/assess_vcf_coverage.md)
  plots of read depth per sample (PDF, PNG and SVG).

### Writing and conversion

- [`write_vcf()`](https://legallab.github.io/VCFArrow/reference/write_vcf.md)
  writes plain or gzip-compressed VCF files.
- Converters to the input formats of population-genetic and phylogenetic
  software:
  [`vcf2admixture()`](https://legallab.github.io/VCFArrow/reference/vcf2admixture.md),
  [`vcf2apparent()`](https://legallab.github.io/VCFArrow/reference/vcf2apparent.md),
  [`vcf2arlequin()`](https://legallab.github.io/VCFArrow/reference/vcf2arlequin.md),
  [`vcf2bayesass()`](https://legallab.github.io/VCFArrow/reference/vcf2bayesass.md),
  [`vcf2bayescan()`](https://legallab.github.io/VCFArrow/reference/vcf2bayescan.md),
  [`vcf2eigenstrat()`](https://legallab.github.io/VCFArrow/reference/vcf2eigenstrat.md),
  [`vcf2fasta()`](https://legallab.github.io/VCFArrow/reference/vcf2fasta.md),
  [`vcf2fineradstructure()`](https://legallab.github.io/VCFArrow/reference/vcf2fineradstructure.md),
  [`vcf2genepop()`](https://legallab.github.io/VCFArrow/reference/vcf2genepop.md),
  [`vcf2migrate()`](https://legallab.github.io/VCFArrow/reference/vcf2migrate.md),
  [`vcf2nexus()`](https://legallab.github.io/VCFArrow/reference/vcf2nexus.md),
  [`vcf2plink_bed()`](https://legallab.github.io/VCFArrow/reference/vcf2plink_bed.md),
  [`vcf2plink_ped()`](https://legallab.github.io/VCFArrow/reference/vcf2plink_ped.md),
  [`vcf2related()`](https://legallab.github.io/VCFArrow/reference/vcf2related.md),
  [`vcf2smartsnp()`](https://legallab.github.io/VCFArrow/reference/vcf2smartsnp.md),
  [`vcf2snapp()`](https://legallab.github.io/VCFArrow/reference/vcf2snapp.md),
  [`vcf2snmf()`](https://legallab.github.io/VCFArrow/reference/vcf2snmf.md),
  [`vcf2structure()`](https://legallab.github.io/VCFArrow/reference/vcf2structure.md)
  and
  [`vcf2treemix()`](https://legallab.github.io/VCFArrow/reference/vcf2treemix.md);
  [`vcf2gt_long()`](https://legallab.github.io/VCFArrow/reference/vcf2gt_long.md)
  writes a long genotype table (CSV, feather or Parquet) and
  [`vcf2genlight()`](https://legallab.github.io/VCFArrow/reference/vcf2genlight.md)
  returns an ‘adegenet’ genlight object.

### Performance

- [`vcf_set_workers()`](https://legallab.github.io/VCFArrow/reference/vcf_set_workers.md)
  runs the per-chunk work of reading, filtering, statistics, merging and
  exporting on background R processes. Results do not depend on the
  number of workers.
- Performance-critical steps (VCF parsing, per-sample counts, theta and
  pi, writers of the exported formats) are implemented in C++.

### Dependencies

- ‘adegenet’ (for
  [`vcf2genlight()`](https://legallab.github.io/VCFArrow/reference/vcf2genlight.md))
  and ‘svglite’ (for SVG plots) are suggested rather than imported;
  without ‘svglite’, SVG plots are written with R’s Cairo device.
