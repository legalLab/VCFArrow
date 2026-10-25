<!-- Converted from VCFArrow_tool_comparison_online.odt with
     pandoc -t gfm --wrap=none (empty anchor tags and bold markers removed). -->

# VCFArrow, vcftools, vcfR and tidypopgen compared

*Complete* P. vaillantii *dataset (18.6 M variants × 56 samples, 26 GB VCF); 50 GB memory and 60-minute limits. Tables follow.*

**Where the genotypes live** explains most of the differences in run time and memory:

  - ***VCFArrow*** (R) keeps the variant table in R memory and the genotypes on disk as chunked Arrow files. Every filter writes a new, complete object chunk by chunk, optionally on parallel workers. It ran all 45 functions: 1.8 h with 4 workers, peaking at 42–47 GB, mostly for the 18.6 M-row variant table.
  - ***vcftools*** (C++, command line) streams the VCF and needs little memory (4.3 GB peak), but each command re-parses the entire file (1.1–2.3 min per filter); filters only list the kept sites. No merging or converters beyond 012 and PLINK ped. 19 functions in 61 min.
  - ***tidypopgen*** (R) first converts the VCF to a file-backed 0/1/2 matrix (30 min), which drops read depth and other FORMAT fields. After that, its filters are lazy (they only record what is kept; 20–25 s). Writing a VCF or PLINK ped did not finish in 60 min, and 4 threads were as fast as 15. 15 of 17 functions in 42 min.
  - ***vcfR*** (R) holds the whole VCF in memory as strings. Subsetting samples is instant (3 s; *VCFArrow* 4 min). Every genotype operation, however, first parses a full sample × locus string matrix, so filters took 6–11 min, and two functions (the heterozygosity filter and the FASTA export) exceeded 50 GB. 15 of 19 functions in 94 min (*VCFArrow*: 54 min).


# What each tool does for the same function

Functions with the same name in the tables do different amounts of work, which explains much of the difference in run time.

Extract samples and groups.

  - **VCFArrow** (*vcf\_extract\_samples()*, *vcf\_extract\_groups()*) writes a new set of genotype chunk files with only the kept samples, then by default (*f\_invar = TRUE*) removes the loci that have become invariant among them, in a second pass over the genotypes (*vcf\_filter\_invariant()*), and reports the remaining missing data. The result is a complete, independent object.
  - **vcftools** (*--keep … --recode*) reads the whole VCF and writes a new VCF with only the kept samples. All sites are kept, including those that have become invariant, and without *--recode-INFO-all* the INFO field is dropped.
  - **tidypopgen** (*gt\[1:10, \]*, *dplyr::filter(gt, population %in% …)*) subsets the rows of the gen\_tibble. The genotypes stay in the shared file-backed matrix, so only the indices of the kept samples are stored, and loci are not filtered.
  - **vcfR** (*vcf\[, c(TRUE, keep)\]*) copies the kept columns of the genotype string matrix in memory. The fixed fields, including INFO, are unchanged and invariant loci are kept.

Genotype table.

  - **VCFArrow** (*vcf2gt\_long()*) writes a long table with one row per sample and variant (on the complete dataset, 1.04 billion rows). Each row holds the variant’s row id, the sample, both allele codes, the phasing, the raw FORMAT string, DP, GQ and the allele-depth ratio. It is written chunk by chunk, here as feather.
  - **vcftools** (*--012*) writes a wide matrix with one row per sample and one column per site. Each cell is the number of non-reference alleles (0, 1 or 2; −1 = missing). The samples and positions are listed in separate files, and the table holds genotypes only.
  - **tidypopgen** has no equivalent.
  - **vcfR** (*vcfR2tidy(single\_frame = TRUE, format\_fields = "GT")*) builds a long table in memory, with one row per sample and variant. Every row repeats the variant’s fixed fields and all parsed INFO fields, and adds the genotype as codes and as alleles. The table is written as CSV, but only after all of it has been built.

θ / pi.

  - **VCFArrow** (*vcf\_theta()*) uses biallelic SNVs only. A site is valid in a group if at least two samples in the group are called. At each valid site:
      - pi = n/(n − 1) · 2p(1 − p), where p is the ALT allele frequency;
      - the Watterson contribution is 1/a<sub>1</sub> if the site is segregating, where a<sub>1</sub> = Σ<sub>i=1..n−1</sub> 1/i and n is the number of called alleles.
    Watterson’s θ and pi are the means of these values over the valid sites. They are computed for each group and for all samples together, giving two numbers per group.
  - **vcftools** (*--site-pi*) computes pi for each site over all samples, as the proportion of differing pairs among all pairs of called alleles: n/(n − 1) · 2p(1 − p) at a biallelic site. It writes one line per site (18.6 M lines on the complete dataset), with no Watterson’s θ, no groups and no average.
  - **tidypopgen** (*loci\_pi()*) computes pi for each locus over all samples, as n/(n − 1) · 2p(1 − p) at a biallelic site. It returns a vector in memory, with no Watterson’s θ (a grouped gen\_tibble gives pi per group).
  - **vcfR** has no equivalent.

Sample statistics.

  - **VCFArrow** (*vcf\_stats()*) computes per sample, over all variants:
      - the mean read depth;
      - the observed heterozygosity (heterozygotes / called genotypes);
      - the numbers of heterozygotes, of homozygotes (REF and ALT) and of missing and called genotypes;
      - the proportion missing and the number of loci.
    It adds each sample’s group and writes a CSV file. θ and π can be added (*theta = TRUE*).
  - **vcftools** (*--het*) gives per sample the observed and expected numbers of homozygous genotypes, the number of sites, and the inbreeding coefficient F. F is a method-of-moments estimate, with the expected homozygosity taken from the allele frequencies of all samples. It reads the whole VCF to do so.
  - **tidypopgen** (*qc\_report\_indiv()*) gives per sample only the observed heterozygosity and the missingness, from the 0/1/2 matrix, in memory. There is no read depth, which the matrix does not hold.
  - **vcfR** has no equivalent.

Merge two VCFs.

  - **VCFArrow** (*vcf\_bind()*) matches variants by CHROM/POS/REF/ALT. Can bind VCFs with non-overlapping or partially overlapping variants as union or intersection.
  - **vcftools** has no equivalent.
  - **tidypopgen** (*rbind()*) row binds/stacks matrices, no matching of loci.
  - **vcfR** has no equivalent.


# Run time by group of functions

Median run time and range of run times per tool and group: “Read” (read VCF),  “Extract” (extract samples, extract groups, subsample loci), “Filters”, and “Writers” (writers/transformers: “→” exporters, and writing VCF). "n": completed functions. "Not complete": functions running out of RAM (\> 50 GB) or time (\> 60 min).

| Tool            | Group   | n  | Median   | Range             | Not complete                                              |
| --------------- | ------- | -- | -------- | ----------------- | --------------------------------------------------------- |
| VCFArrow (1)    | Read    | 1  | 7.0 min  |                   |                                                           |
| VCFArrow (1)    | Extract | 3  | 3.7 min  | 0.5 s–4.1 min     |                                                           |
| VCFArrow (1)    | Filters | 15 | 18.0 s   | 0.4 s–9.0 min     |                                                           |
| VCFArrow (1)    | Writers | 19 | 4.4 min  | 3.1 min–26.3 min  |                                                           |
| VCFArrow (4)    | Read    | 1  | 6.9 min  |                   |                                                           |
| VCFArrow (4)    | Extract | 3  | 55.6 s   | 0.5 s–1.0 min     |                                                           |
| VCFArrow (4)    | Filters | 15 | 17.8 s   | 0.5 s–3.6 min     |                                                           |
| VCFArrow (4)    | Writers | 19 | 2.3 min  | 1.0 min–24.6 min  |                                                           |
| VCFArrow (8)    | Read    | 1  | 7.1 min  |                   |                                                           |
| VCFArrow (8)    | Extract | 3  | 48.9 s   | 0.6 s–56.2 s      |                                                           |
| VCFArrow (8)    | Filters | 15 | 14.6 s   | 0.6 s–2.7 min     |                                                           |
| VCFArrow (8)    | Writers | 19 | 2.7 min  | 1.1 min–25.8 min  |                                                           |
| vcftools        | Read    | 1  | 1.1 min  |                   |                                                           |
| vcftools        | Extract | 2  | 5.6 min  | 5.2 min–5.9 min   |                                                           |
| vcftools        | Filters | 9  | 1.4 min  | 1.1 min–2.3 min   |                                                           |
| vcftools        | Writers | 2  | 3.9 min  | 3.0 min–4.8 min   |                                                           |
| tidypopgen (4)  | Read    | 1  | 30.0 min |                   |                                                           |
| tidypopgen (4)  | Extract | 3  | 0.0 s    | 0.0 s–6.1 min     |                                                           |
| tidypopgen (4)  | Filters | 4  | 21.0 s   | 4.7 s–25.1 s      |                                                           |
| tidypopgen (4)  | Writers | 2  | 1.7 min  | 4.2 s–3.3 min     | → PLINK ped (\> 60 min), → VCF (\> 60 min)                |
| tidypopgen (15) | Read    | 1  | 30.0 min |                   |                                                           |
| tidypopgen (15) | Extract | 3  | 0.0 s    | 0.0 s–6.1 min     |                                                           |
| tidypopgen (15) | Filters | 4  | 20.5 s   | 4.7 s–24.1 s      |                                                           |
| tidypopgen (15) | Writers | 2  | 1.7 min  | 4.1 s–3.3 min     | → PLINK ped (\> 60 min), → VCF (\> 60 min)                |
| vcfR            | Read    | 1  | 12.9 min |                   |                                                           |
| vcfR            | Extract | 3  | 2.7 s    | 0.0 s–3.0 s       |                                                           |
| vcfR            | Filters | 8  | 4.0 min  | 11.3 s–11.1 min   | heterozygosity (out of memory), MAF (\> 60 min)           |
| vcfR            | Writers | 2  | 18.3 min | 15.0 min–21.7 min | → FASTA (out of memory), → genotype table (out of memory) |

  - **Extract and subsample:** VCFArrow rewrites the genotype files (3.7 min with 1 worker, 49–56 s with 4 workers); vcfR (2.7 s) and tidypopgen (0.0 s) only index existing data; vcftools re-reads the file (5.6 min).
  - **Filters:** VCFArrow 15–18 s (variant-table filters 0.4–8 s, genotype filters 1–2 min with 1–8 workers), tidypopgen 20–25 s (lazy), vcftools 1.4 min, vcfR 4.0 min.
  - **Writers/transformers:** VCFArrow 4.4 min with 1 worker and 2.3–2.7 min with 4 workers (19 converters; high outlier: the feather genotype table, 25–26 min); tidypopgen 1.7 min (2 of 3 completed), vcftools 3.9 min (2), vcfR 18.3 min (2 of 4 completed).


# Run time

“VCFArrow (*n*)”: *n* worker processes. “tidypopgen (*n*)”: *n* threads. “—”: no equivalent function. “\> 60 min”: nor within 60 min when run again. “out of memory”: the R session was killed at the 50 GB limit. Each function was run once, from the unfiltered data.

| Function                    | VCFArrow (1) | VCFArrow (4) | VCFArrow (8) | vcftools | tidypopgen (4) | tidypopgen (15) | vcfR          |
| --------------------------- | ------------ | ------------ | ------------ | -------- | -------------- | --------------- | ------------- |
| read VCF                    | 7.0 min      | 6.9 min      | 7.1 min      | 1.1 min  | 30.0 min       | 30.1 min        | 12.9 min      |
| set groups                  | 0.0 s        | 0.0 s        | 0.0 s        | —        | 0.0 s          | 0.0 s           | —             |
| extract samples             | 4.1 min      | 55.6 s       | 56.2 s       | 5.9 min  | 0.0 s          | 0.0 s           | 3.0 s         |
| extract groups              | 3.7 min      | 1.0 min      | 48.9 s       | 5.2 min  | 0.0 s          | 0.0 s           | 2.7 s         |
| subsample loci (1000)       | 0.5 s        | 0.5 s        | 0.6 s        | —        | 6.1 min        | 6.1 min         | 0.0 s         |
| filter coverage (DP 10)     | 2.0 min      | 41.0 s       | 48.1 s       | 2.3 min  | —              | —               | 7.2 min       |
| filter heterozygosity (0.4) | 1.5 min      | 29.6 s       | 28.3 s       | —        | —              | —               | out of memory |
| filter MAF (0.05)           | 1.7 min      | 33.7 s       | 28.9 s       | 2.2 min  | 24.1 s         | 25.1 s          | \> 60 min     |
| filter locus missing (0.2)  | 56.4 s       | 19.0 s       | 14.6 s       | 2.1 min  | 21.0 s         | 21.1 s          | 5.7 min       |
| filter sample missing (0.8) | 9.0 min      | 3.6 min      | 2.7 min      | —        | 4.7 s          | 4.7 s           | 5.9 min       |
| filter quality (30)         | 5.0 s        | 5.9 s        | 8.0 s        | 1.1 min  | —              | —               | 11.3 s        |
| filter rank (0.4)           | 3.9 s        | 4.7 s        | 6.8 s        | —        | —              | —               | —             |
| filter repeat (0.5)         | 0.4 s        | 0.5 s        | 0.6 s        | —        | —              | —               | —             |
| filter paralogs (0.2, 0.5)  | 5.0 s        | 6.0 s        | 8.2 s        | —        | —              | —               | —             |
| filter PASS                 | 3.9 s        | 4.8 s        | 6.9 s        | 1.4 min  | —              | —               | 11.5 s        |
| filter biallelic            | 3.8 s        | 4.4 s        | 4.6 s        | 1.3 min  | —              | —               | 19.9 s        |
| filter indels               | 3.7 s        | 4.3 s        | 4.3 s        | 1.2 min  | —              | —               | 2.3 min       |
| filter invariant            | 2.0 min      | 41.8 s       | 56.3 s       | 2.3 min  | 20.1 s         | 21.0 s          | 11.1 min      |
| filter one SNV per block    | 19.2 s       | 18.5 s       | 18.9 s       | 1.2 min  | —              | —               | —             |
| filter multi SNV blocks     | 18.0 s       | 17.8 s       | 17.4 s       | —        | —              | —               | —             |
| merge two VCFs              | 20.7 min     | 10.3 min     | 9.3 min      | —        | 1.4 min        | 1.5 min         | —             |
| sample statistics           | 1.3 min      | 1.0 min      | 1.4 min      | 2.2 min  | 4.7 s          | 4.8 s           | —             |
| theta / pi                  | 3.7 min      | 1.1 min      | 38.7 s       | 2.5 min  | 0.6 s          | 1.1 s           | —             |
| missing data plot + table   | 49.2 s       | 46.0 s       | 11.9 s       | 1.9 min  | 4.7 s          | 4.7 s           | —             |
| coverage plot               | 2.3 min      | 41.7 s       | 36.2 s       | 2.1 min  | —              | —               | —             |
| write VCF                   | 12.2 min     | 8.5 min      | 10.0 min     | 17.0 min | \> 60 min      | \> 60 min       | 11.6 min      |
| → apparent                  | 5.7 min      | 3.9 min      | 3.5 min      | —        | —              | —               | —             |
| → Arlequin                  | 4.3 min      | 2.8 min      | 2.5 min      | —        | —              | —               | —             |
| → BayesAss                  | 5.4 min      | 3.1 min      | 4.9 min      | —        | —              | —               | —             |
| → BayeScan                  | 6.5 min      | 3.1 min      | 3.4 min      | —        | —              | —               | —             |
| → EIGENSTRAT                | 3.5 min      | 1.2 min      | 1.1 min      | —        | —              | —               | —             |
| → FASTA                     | 3.8 min      | 2.3 min      | 2.2 min      | —        | —              | —               | out of memory |
| → fineRADstructure          | 3.1 min      | 1.0 min      | 1.2 min      | —        | —              | —               | —             |
| → GENEPOP                   | 5.7 min      | 2.8 min      | 3.7 min      | —        | —              | —               | —             |
| → genlight                  | 8.9 min      | 6.9 min      | 6.9 min      | —        | 3.3 min        | 3.3 min         | 15.0 min      |
| → genotype table            | 26.3 min     | 24.6 min     | 25.8 min     | 3.0 min  | —              | —               | out of memory |
| → Migrate-N                 | 4.0 min      | 2.1 min      | 2.4 min      | —        | —              | —               | 21.7 min      |
| → NEXUS                     | 3.9 min      | 2.1 min      | 2.8 min      | —        | —              | —               | —             |
| → PLINK ped                 | 5.7 min      | 2.9 min      | 3.8 min      | 4.8 min  | \> 60 min      | \> 60 min       | —             |
| → PLINK bed                 | 3.7 min      | 1.7 min      | 1.9 min      | —        | 4.1 s          | 4.2 s           | —             |
| → related                   | 4.4 min      | 2.4 min      | 2.7 min      | —        | —              | —               | —             |
| → SmartSNP                  | 3.8 min      | 1.4 min      | 1.2 min      | —        | —              | —               | —             |
| → SNAPP                     | 3.6 min      | 2.1 min      | 2.0 min      | —        | —              | —               | —             |
| → STRUCTURE                 | 4.6 min      | 2.1 min      | 3.1 min      | —        | —              | —               | —             |
| → TreeMix                   | 5.1 min      | 2.3 min      | 2.2 min      | —        | —              | —               | —             |


# Memory

Median/peak memory in GB (PSS of all processes of the run, including the data already loaded) while each function ran.

| Function                    | VCFArrow (1) | VCFArrow (4) | VCFArrow (8) | vcftools  | tidypopgen (4) | tidypopgen (15) | vcfR        |
| --------------------------- | ------------ | ------------ | ------------ | --------- | -------------- | --------------- | ----------- |
| read VCF                    | 13.13/22.42  | 8.94/14.10   | 12.88/17.91  | 0.00/0.00 | 23.38/32.29    | 23.52/32.29     | 10.71/23.84 |
| set groups                  | 20.61/20.62  | 13.06/13.06  | 16.18/16.25  | —         | 30.96/30.96    | 30.96/30.96     | —           |
| extract samples             | 21.56/22.45  | 16.65/19.05  | 23.65/26.32  | 0.00/0.00 | 30.96/30.96    | 30.96/30.96     | 16.37/23.84 |
| extract groups              | 21.65/24.61  | 19.20/21.61  | 25.79/29.03  | 0.00/0.00 | 30.96/30.96    | 30.96/30.96     | 18.72/19.90 |
| subsample loci (1000)       | 24.61/24.61  | 21.06/21.11  | 27.74/27.74  | —         | 30.96/30.96    | 30.96/30.96     | 17.13/19.90 |
| filter coverage (DP 10)     | 23.56/24.61  | 16.58/20.70  | 23.03/27.35  | 0.00/0.00 | —              | —               | 25.24/44.83 |
| filter heterozygosity (0.4) | 23.77/25.47  | 17.82/22.94  | 24.91/29.61  | —         | —              | —               | 48.78/49.54 |
| filter MAF (0.05)           | 24.76/25.53  | 18.47/23.73  | 25.20/29.33  | 0.00/0.00 | 30.96/38.25    | 30.96/38.24     | 33.83/34.21 |
| filter locus missing (0.2)  | 24.66/25.51  | 18.27/23.94  | 25.22/27.08  | 0.00/0.00 | 30.96/38.25    | 30.96/38.24     | 25.73/32.05 |
| filter sample missing (0.8) | 24.21/27.37  | 20.44/24.04  | 28.31/30.39  | —         | 30.96/30.96    | 30.96/38.24     | 26.80/35.38 |
| filter quality (30)         | 27.37/28.06  | 21.26/25.13  | 29.81/32.31  | 0.00/0.00 | —              | —               | 34.87/36.11 |
| filter rank (0.4)           | 27.99/28.06  | 24.72/25.13  | 29.13/32.36  | —         | —              | —               | —           |
| filter repeat (0.5)         | 27.66/27.66  | 24.86/24.86  | 28.54/32.36  | —         | —              | —               | —           |
| filter paralogs (0.2, 0.5)  | 25.77/27.66  | 19.80/23.94  | 27.50/32.41  | —         | —              | —               | —           |
| filter PASS                 | 27.07/27.66  | 24.52/24.80  | 28.18/32.41  | 0.00/0.00 | —              | —               | 35.63/36.88 |
| filter biallelic            | 27.66/27.66  | 24.52/24.80  | 28.75/32.36  | 0.00/0.00 | —              | —               | 29.66/37.29 |
| filter indels               | 27.66/27.66  | 24.24/24.52  | 32.14/32.23  | 0.00/0.00 | —              | —               | 27.97/37.10 |
| filter invariant            | 25.68/28.72  | 19.63/24.21  | 26.90/32.05  | 0.00/0.00 | 30.96/38.25    | 30.96/38.24     | 34.74/45.85 |
| filter one SNV per block    | 26.48/28.72  | 21.46/25.54  | 28.59/31.40  | 0.00/0.00 | —              | —               | —           |
| filter multi SNV blocks     | 25.13/26.55  | 20.62/25.17  | 28.20/31.31  | —         | —              | —               | —           |
| merge two VCFs              | 25.93/28.37  | 31.79/32.53  | 38.92/41.52  | —         | 32.65/38.25    | 32.71/38.24     | —           |
| sample statistics           | 28.83/29.85  | 23.35/32.27  | 38.06/41.21  | 0.00/0.00 | 30.96/33.82    | 30.96/31.87     | —           |
| theta/pi                    | 22.50/23.12  | 24.21/25.10  | 27.65/29.38  | 0.00/0.00 | 30.96/30.96    | 30.96/30.96     | —           |
| missing data plot + table   | 23.92/24.15  | 13.64/13.64  | 16.99/17.65  | 0.00/0.00 | 30.96/30.96    | 30.96/30.96     | —           |
| coverage plot               | 23.65/24.15  | 13.75/14.55  | 17.53/18.35  | 0.00/0.00 | —              | —               | —           |
| write VCF                   | 23.80/24.05  | 15.87/16.64  | 19.21/20.95  | 0.00/0.00 | 30.96/30.96    | 30.96/30.96     | 29.58/29.58 |
| → apparent                  | 26.79/42.05  | 22.65/34.29  | 26.66/39.95  | —         | —              | —               | —           |
| → Arlequin                  | 30.51/42.30  | 31.30/35.71  | 35.06/39.55  | —         | —              | —               | —           |
| → BayesAss                  | 30.64/44.49  | 32.13/35.76  | 36.19/39.93  | —         | —              | —               | —           |
| → BayeScan                  | 30.13/42.34  | 30.94/35.76  | 34.99/41.08  | —         | —              | —               | —           |
| → EIGENSTRAT                | 26.63/40.04  | 23.43/32.99  | 27.95/36.25  | —         | —              | —               | —           |
| → FASTA                     | 26.93/42.20  | 26.28/35.16  | 31.37/40.22  | —         | —              | —               | 21.04/48.90 |
| → fineRADstructure          | 31.37/42.20  | 24.23/35.16  | 29.58/40.22  | —         | —              | —               | —           |
| → GENEPOP                   | 27.09/42.25  | 27.03/35.28  | 39.59/40.27  | —         | —              | —               | —           |
| → genlight                  | 39.44/45.13  | 36.18/42.33  | 39.64/46.85  | —         | 38.25/39.01    | 38.25/39.00     | 35.17/47.49 |
| → genotype table            | 25.59/40.48  | 25.51/37.41  | 29.82/40.87  | 0.00/0.06 | —              | —               | 26.93/49.35 |
| → Migrate-N                 | 30.69/42.34  | 30.72/40.29  | 33.80/44.40  | —         | —              | —               | 22.82/43.50 |
| → NEXUS                     | 30.69/45.31  | 31.59/40.29  | 34.07/44.51  | —         | —              | —               | —           |
| → PLINK ped                 | 30.44/45.51  | 29.10/41.02  | 35.99/45.57  | 2.51/4.34 | 32.53/38.64    | 32.79/38.51     | —           |
| → PLINK bed                 | 33.70/42.81  | 27.28/41.02  | 31.71/45.57  | —         | 30.84/31.10    | 30.84/31.10     | —           |
| → related                   | 27.10/42.36  | 31.60/40.23  | 39.06/44.74  | —         | —              | —               | —           |
| → SmartSNP                  | 31.48/42.36  | 23.96/40.23  | 28.59/44.74  | —         | —              | —               | —           |
| → SNAPP                     | 27.25/42.40  | 27.62/39.24  | 32.20/43.70  | —         | —              | —               | —           |
| → STRUCTURE                 | 30.81/42.41  | 31.50/40.30  | 39.07/44.60  | —         | —              | —               | —           |
| → TreeMix                   | 30.88/42.41  | 23.84/40.30  | 28.35/44.60  | —         | —              | —               | —           |


# Entire run

| Run                    | Total run time | Completed | Time out | Memory limit | Peak memory |
| ---------------------- | -------------- | --------- | -------- | ------------ | ----------- |
| VCFArrow, 1 worker     | 3.10 h         | 45 of 45  | 0        | 0            | 45.51       |
| VCFArrow, 4 workers    | 1.83 h         | 45 of 45  | 0        | 0            | 42.33       |
| VCFArrow, 8 workers    | 1.92 h         | 45 of 45  | 0        | 0            | 46.85       |
| vcftools               | 1.02 h         | 19 of 19  | 0        | 0            | 4.34        |
| tidypopgen, 4 threads  | 42.4 min       | 15 of 17  | 2        | 0            | 39.00       |
| tidypopgen, 15 threads | 42.2 min       | 15 of 17  | 2        | 0            | 39.01       |
| vcfR                   | 1.57 h         | 15 of 19  | 1        | 3            | 49.54       |
