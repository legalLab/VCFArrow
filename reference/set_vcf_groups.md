# set_vcf_groups

Associate samples with group info based on "strt" file

## Usage

``` r
set_vcf_groups(vcf_arrow, data_path, strt = "strata")
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- data_path:

  -\> path to data (where strata are located)

- strt:

  -\> file with 2+ columns, one id column and one pop column (tab- or
  space-separated, with a header line), default 'strata'

## Value

vcf_arrow object with groups slot filled

## Details

This function assigns each sample of a VCFArrow object to a group, based
on a strata file. The strt is a file with a header line and 2+ columns,
two of them named 'id' (sample name) and 'pop' (group); other columns
are ignored. The file is read as tab-separated when its header line
contains a tab, so fields may contain spaces and may be empty; otherwise
its columns are separated by spaces (any whitespace). All columns are
read as text, so ids such as '007' keep their leading zeros. Ids in the
file that are not samples of the object are ignored. An id listed more
than once with different groups is an error. Samples not in the file, or
with an empty group, get the group NA and a warning lists them;
functions that need groups (e.g. the exporters) stop until every sample
has one. By default strt is named 'strata' and found in the datapath
where the VCF file is found.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- set_vcf_groups(read_vcf(f), data_path = dirname(f),
                      strt = "strata")
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
table(vcf@groups)
#> 
#> BS GS OG WA 
#>  5  5  3  5 
```
