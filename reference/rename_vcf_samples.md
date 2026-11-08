# rename_vcf_samples

Rename samples of a VCFArrow object

## Usage

``` r
rename_vcf_samples(vcf_arrow, new_names, update_groups = FALSE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- new_names:

  -\> the new sample names: a named character vector (c(old = "new",
  ...)), a data frame, or the path of a tab- or space-separated file
  with a header line; the data frame or file has the columns id (current
  name) and new_id (new name), and optionally pop (group) (character,
  data.frame or file path)

- update_groups:

  -\> also update the group assignments, from the pop column of
  new_names, default FALSE (keep the current group of each sample)
  (Boolean)

## Value

VCFArrow object with the samples renamed

## Details

This function renames samples of a VCFArrow object, returning a new
VCFArrow object. Samples not listed in new_names keep their names; names
in new_names that are not samples of the object are ignored, with a
warning. The new names must be unique (also with respect to the samples
not renamed) and must not be empty or contain tabs or line breaks; a
sample listed twice with different new names is an error. Names can also
be exchanged (e.g. c(A = "B", B = "A")). By default (update_groups =
FALSE) each sample keeps its group. With update_groups = TRUE, the
samples listed get the group in the pop column of new_names (a listed
sample with an empty pop keeps its group); this needs a data frame or
file with a pop column. A file is read as in set_vcf_groups():
tab-separated if its header line has a tab (fields may then contain
spaces), otherwise space-separated, all columns as text. Sample names
are stored with every genotype, so the genotype files are rewritten (in
parallel with workers, see vcf_set_workers()); the genotypes themselves
are unchanged.

## Author

Tomas Hrbek November 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
vcf2 <- rename_vcf_samples(vcf, c(Pv120 = "WA_Pv120", Pv13 = "WA_Pv13"))
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> ℹ Renamed 2 samples; group assignments kept.
head(vcf2@samples)
#> [1] "WA_Pv120" "Pv126"    "WA_Pv13"  "Pv14"     "Pv27"     "Pv28"    
# rename and move a sample to another group
new <- data.frame(id = "Pv120", new_id = "OG_Pv120", pop = "OG")
vcf3 <- rename_vcf_samples(vcf, new, update_groups = TRUE)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> ℹ Renamed 1 sample; 1 group assignment(s) updated.
```
