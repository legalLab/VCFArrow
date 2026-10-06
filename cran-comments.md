## Submission

This is a new submission.

VCFArrow reads Variant Call Format (VCF) files into an S4 object whose
genotypes are stored on disk as an 'Apache Arrow' dataset, and provides
filters, statistics and converters to the input formats of population-genetic
software. It is designed for datasets that do not fit in memory.

## Test environments

* Local: Ubuntu 22.04, R 4.6.1 (`R CMD check --as-cran`, including the CRAN
  incoming checks and `--run-donttest`).
* GitHub Actions (r-lib/actions): ubuntu-latest (R devel, release and
  oldrel-1), macos-latest (R release), windows-latest (R release).

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

## Notes for the reviewers

* Examples and the vignette write files only to `tempdir()`.
* Examples use at most 2 cores: the `vcf_set_workers()` example starts 2
  background R processes and stops them again (in `\donttest{}`, as it starts
  processes).
* The package sets the environment variable `ARROW_DEFAULT_MEMORY_POOL` to
  `"system"` when it is loaded, only if the user has not set it, so that
  memory freed by 'Apache Arrow' is returned to the operating system. This
  must happen before 'arrow' initialises its memory pool.
* Temporary genotype data of `VCFArrow` objects is kept in subdirectories of
  `tempdir()` and removed when the objects are garbage-collected or the
  session ends.
