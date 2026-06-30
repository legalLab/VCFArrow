# write_vcf

Write a VCFArrow object to an external VCF file

## Usage

``` r
write_vcf(vcf_arrow, out_file = "output.vcf", gzip = FALSE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of the VCF file to be written to, default 'output.vcf'
  (character)

- gzip:

  -\> a flag to GZIP VCF when writing, default FALSE (Boolean)

## Details

This function writes a VCFArrow object to an external VCF file. Writing
occurs in chunks whose size is determined by the read_vcf() function.
Larger chunks result in faster writing speeds. It writes both
uncompressed and gz compressed files. Compressing increases writing time
be about 50%. For large files, it is recommended to output an
uncompressed VCF file, and then compress with GZIP or PIGZ.

## Author

Tomas Hrbek April 2026

## Examples

``` r
write_vcf(vcf = my_vcf, out_file = "output.vcf", gzip = FALSE)
#> Error: object 'my_vcf' not found
write_vcf(my_vcf, "output.vcf", FALSE)
#> Error: object 'my_vcf' not found
write_vcf(my_vcf)
#> Error: object 'my_vcf' not found
```
