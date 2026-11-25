## Submission

This is an update of VCFArrow 0.9.0 (on CRAN since 16 October 2026) to
version 1.0.0.

* It is compatible with 'arrow' 26.0.0, which deprecated `read_feather()`,
  `write_feather()` and `format = "feather"` in `open_dataset()`. VCFArrow
  now uses `read_ipc_file()`, `write_ipc_file()` and `format = "ipc"`, so it
  no longer gives deprecation warnings with 'arrow' 26, and requires
  'arrow' (>= 10.0.0).
* It adds the functions `rename_vcf_samples()` and `vcf2fstructure()`.
* It fixes errors in `read_vcf()` (genotypes that were misread), in the
  exporters (file formats now as their programs document them), in
  `vcf_filter_adr()`, `vcf_filter_coverage()`, `vcf_filter_invariant()`
  (multiallelic sites), `vcf_filter_multiSNV()`, `vcf_sub_SNVs_stratified()`
  and `vcf_bind()`, makes `vcf_filter_missing()` and `vcf_filter_hets()`
  keep values equal to their thresholds, and adds options to several filters.
* The default `threshold` of `vcf_filter_maf()` is now 0.01 (was 0.05), the
  conventional minor allele frequency of a polymorphism (Karki et al. 2015,
  cited on its help page).

NEWS.md lists all changes.

## Test environments

* Local: Ubuntu 22.04, R 4.6.1, with 'arrow' 26.0.1.1, 26.0.1 and 25.0.1
  (`R CMD check --as-cran --run-donttest`, including the remote CRAN
  incoming checks).
* GitHub Actions (r-lib/actions): ubuntu-latest (R devel, release and
  oldrel-1), macos-latest (R release), windows-latest (R release).

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

## Notes for the reviewers

* Examples and the vignette write files only to `tempdir()`.
* Examples use at most 2 cores: the `vcf_set_workers()` example starts 2
  background R processes and stops them again (in `\donttest{}`, as it starts
  processes).
* The package sets the environment variable `ARROW_DEFAULT_MEMORY_POOL` to
  `"system"` when it is loaded, only if the user has not set it, so that
  memory freed by 'Apache Arrow' is returned to the operating system. This
  must happen before 'arrow' initializes its memory pool.
* Temporary genotype data of `VCFArrow` objects is kept in subdirectories of
  `tempdir()` and removed when the objects are garbage-collected or the
  session ends.
