#' Mock summary statistics from 1kg and fake phenotypes
#'
#' Variant level summary statistics as output from the Analysis Commons Genesis_wdl workflow
#' diamonds.
#'
#' @format A data frame with 426 rows and 13 variables:
#' \describe{
#'   \item{MarkerName}{unique variant identifier, shared with LD data}
#'   \item{chr}{chromosome of the variant}
#'   \item{pos}{genomic position of the variant}
#'   \item{ref}{reference allele}
#'   \item{alt}{alternate allele}
#'   \item{pvalue}{pvalue from GWAS}
#'   ...
#' }
"variant_data"

#' LD r^2 values for a small number of variants from 1kg phase 3
#'
#' @format A data frame with 1 row and 115 variables:
#'
"variant_ld_data"

#' Chromatin state predictions from Varshney, et al 2017
#'
#' A list of data frames from bed files with genomic regions and chromatin state predictions
#'
#' @format A list with one data frame with 53 rows and 5 variables:
#' \describe{
#'   \item{V1}{Chromosome}
#'   \item{V2}{genomic start coordinate}
#'   \item{V3}{genomic end coordinate}
#'   \item{V4}{chromatin state label}
#'   \item{V5}{chromatin state color for plotting}
#'   ...
#' }
"bed_data"
