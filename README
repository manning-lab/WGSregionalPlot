# WGSregionalPlot
## Generate regional views of GWAS results from whole genome sequence data.
Maintainer: Tim Majarian
Version: 0.1

## Description:
This package provides functions for plotting summary statistics from GWAS with a particular focus on whole genome sequence and functional data. P-values, genomic context, and multiple bed file tracks can be plotted over small or large regions of the genome. LD information can also be displayed.

#### What data can this package be used with?
The minimum data required to generate a plot are per-variant GWAS summary statistics. Four values must be defined for each variant:
- Variant id: a unique identifier for each variant
- Chromosome: the chromosome location of each variant
- Position: the genomic position of each variant
- P-value: a p-value for each variant  (or some other data to plot on the Y-axis)

## Main Functions
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
See vignettes/regional_plot_vignette.html
