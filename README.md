# WGSregionalPlot
**Generate regional views of GWAS results from whole genome sequence data.**

Maintainer: Tim Majarian

Version: 0.1

## Description:
------------------------------------------
This package provides functions for plotting summary statistics from GWAS with a particular focus on whole genome sequence and functional data. P-values, genomic context, and multiple bed file tracks can be plotted over small or large regions of the genome. LD information can also be displayed.

### What data can this package be used with?
**Summary statistics**

The minimum data required to generate a plot are per-variant GWAS summary statistics. Four values must be defined for each variant:
- Variant id: a unique identifier for each variant
- Chromosome: the chromosome location of each variant
- Position: the genomic position of each variant
- P-value: a p-value for each variant  (or some other data to plot on the Y-axis)


| MarkerName | chr | pos | pvalue |
-----------|--------|----|---------
| rs12345 |	1	| 234030 |	0.803 |
| rs54321 |	10	| 1268346 |	0.001 |
| rs98765 |	22	| 68263057 |	0.283 |

**LD**

LD information can also be incorperated into the plot. To do so, a file with LD-values per variant needs to be input. This package is designed to be used in conjuction with and handles outputs from our [LDGds pipeline](https://github.com/AnalysisCommons/LDGds). LD data should be stored in a plain, delimited text file as either a matrix or a row vector. Either format must have variant identifiers as row and column names. In conjunction with the LD data, a reference variant must be specified. See the package vignette and accompanying data for an example of input LD data.

| | rs12345 | rs54321 | rs98765 |
---|-------|----------|-------------
| rs12345 |	1	| 0.0001 |	0.1001 |

## Installation
------------------------------------------
To install this package, ensure that [*devtools*](https://cran.r-project.org/web/packages/devtools/readme/README.html) is installed. Then run:

```R
devtools::install_github("manning-lab/WGSregionalPlot")
library(WGSregionalPlot)
```

## Main Functions
------------------------------------------
**Note:** See further documentation for advanced options.

### make_regional_plot
Generates a regional plot.

Inputs:
- chr: [string] chromosome of desired plotting interval
- start: [int] start position of desired plotting interval
- end: [int] end position of desired plotting interval
- variant_data: [data frame] variant-level summary statistics organized with one row per variant, including columns as follows.
- variant_chr_column: [string] column name within *variant_data* of chromosome values for each variant
- variant_pos_column: [string] column name within *variant_data* of position values for each variant
- variant_y_column: [string] column name within *variant_data* of p-values values for each variant

## Tutorial
------------------------------------------
See [vignettes/regional_plot_vignette.html](https://github.com/manning-lab/WGSregionalPlot/blob/master/vignettes/regional_plot_vingette.html)
