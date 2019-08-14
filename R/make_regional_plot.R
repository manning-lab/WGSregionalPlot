#' Generate a regional plot from summary statistics for a desired interval
#' @param chr interval chromosome
#' @param start interval start
#' @param end interval end
#' @param variant_data variant-level summary statistics
#' @param variant_chr_column chromosome column in \code{variant_data}
#' @param variant_pos_column variant position column in \code{variant_data}
#' @param variant_y_column p-value column in \code{variant_data}
#' @param variant_marker_column unique variant identifier column in \code{variant_data}
#' @param variant_ld_data variant LD values, the output of LDGds workflow
#' @param variant_ld_ref reference variant unique identifier for LD
#' @param genome_build genome build of variants
#' @param variant_should_log should the pvalues in \code{variant_data} be -log10'ed?
#' @param variant_horizontal_value Y-intercept for horizontal line in Manhattan plot
#' @param variant_horizontal_color Color for horizontal line in Manhattan plot
#' @param variant_horizontal_style Style for horizontal line in Manhattan plot
#' @param variant_background_color Manhattan plot title box color
#' @param variant_background_frame Should a frame be placed around the Manhattan plot?
#' @param variant_point_color Color of points if no LD information is provided
#' @param variant_title Manhattan plot title
#' @param gene_highlight Should a gene be highlighted? hgnc symbol
#' @param gene_title Title for gene plot
#' @param gene_background_color Gene plot title box color
#' @param gene_frame Should a frame be placed around the Gene plot?
#' @param bed_data List of dataframes in bed format: chr start end label color
#' @param bed_titles List of titles for bed plots
#' @param bed_background_colors List of background colors for bed plot title bars
#' @param bed_frame Should frames be added around each bed plot?

#' @importFrom Gviz GenomeAxisTrack IdeogramTrack plotTracks
#' @export

# Make a regional plot
####

# Notes:
# Function for fixing chromosome name convention is needed (chr1 vs 1)
# Check that the number of titles matches the number of bed files
# Assign default background colors to bed tracks
# Move from biomart query to tabix/htslib for gene regions

# Required packages
# library(GenomicRanges)
# library(Gviz)
# library(biomaRt)
# # library(RColorBrewer)
# library(data.table)

make_regional_plot <- function(chr, start, end, variant_data, variant_chr_column, variant_pos_column, variant_y_column,
                     variant_marker_column = NULL,
                     variant_ld_data = NULL,
                     variant_ld_ref = NULL,
                     genome_build = "hg19",
                     variant_should_log = T,
                     variant_horizontal_value = 5e-8,
                     variant_horizontal_color = "red",
                     variant_horizontal_style = "dashed",
                     variant_background_color = "#e6550d",
                     variant_background_frame = T,
                     variant_point_color = "#252525",
                     variant_title = " ",
                     gene_highlight = NULL,
                     gene_title = NULL,
                     gene_background_color = "#31a354",
                     gene_frame = T,
                     bed_data = NULL,
                     bed_titles = NULL,
                     bed_background_colors = NULL,
                     bed_frame = T){

  # make the data track for variant statistics
  variant.track <- make_variant_track(variant_data,
                                    variant_chr_column,
                                    variant_pos_column,
                                    variant_y_column,
                                    variant_marker_column,
                                    variant_ld_data,
                                    variant_ld_ref,
                                    genome_build,
                                    variant_should_log,
                                    variant_horizontal_value,
                                    variant_horizontal_color,
                                    variant_horizontal_style,
                                    variant_background_color,
                                    variant_background_frame,
                                    variant_point_color,
                                    variant_title)

  # make a mart object for gene coordinates
  mart <- make_mart(genome_build)

  # make gene track
  gene.track <- make_gene_track(chr,
                              start,
                              end,
                              mart,
                              genome_build,
                              gene_highlight,
                              gene_title)

  # make bed tracks
  if (!is.null(bed_data)){
    bed.tracks <- list()
    for (bed_index in 1:length(bed_data)){
      bed.tracks[[bed_index]] <- make_bed_track(bed_data[[bed_index]],
                                              bed_titles[bed_index],
                                              genome_build,
                                              bed_background_colors[bed_index],
                                              bed_frame)
    }
  }

  # make scale track
  scale.track <- GenomeAxisTrack(genome = genome_build,
                                 chromosome = chr)

  # make idiogram track
  idiogram.track <- IdeogramTrack(genome = genome_build,
                                  chromosome = chr)

  # generate the plot
  if (!is.null(bed_data)){
    plotTracks(
      c(list(idiogram.track, scale.track, variant.track, gene.track), bed.tracks),
      showTitle = TRUE,
      sizes=c(1, 2, 5, 4, rep(2, length(bed.tracks))),
      from = start,
      to = end
    )
  } else {
    plotTracks(
      list(idiogram.track, scale.track, variant.track, gene.track),
      showTitle = TRUE,
      sizes=c(1, 2, 5, 4),
      from = start,
      to = end
    )
  }
}
