# VCFArrow: Fast and Memory-Efficient Manipulation and Conversion of VCF Files

Reads Variant Call Format (VCF) files into an S4 object whose genotypes
are stored on disk as an 'Apache Arrow' dataset and loaded lazily,
allowing large datasets to be processed with little memory. Provides
functions to filter variants and samples, subset and merge datasets,
compute summary statistics and diagnostic plots, write VCF files, and
convert data to input formats of population genetic and phylogenetic
software such as 'STRUCTURE', 'ADMIXTURE', 'PLINK', 'EIGENSTRAT',
'Arlequin', 'GENEPOP', 'BayeScan', 'TreeMix', 'SNAPP' and 'adegenet'.

## See also

Useful links:

- <https://legallab.github.io/VCFArrow/>

- <https://github.com/legalLab/VCFArrow>

- Report bugs at <https://github.com/legalLab/VCFArrow/issues>

## Author

**Maintainer**: Tomas Hrbek <hrbek@evoamazon.net>
([ORCID](https://orcid.org/0000-0003-3239-7068)) \[copyright holder\]

Authors:

- Tomas Hrbek <hrbek@evoamazon.net>
  ([ORCID](https://orcid.org/0000-0003-3239-7068)) \[copyright holder\]
