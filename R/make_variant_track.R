#' Create a manhattan plot
#' @param chr interval chromosome
#' @param start interval start
#' @param end interval end
#' @param variant_data variant-level summary statistics
#' @param chr_column chromosome column in \code{variant_data}
#' @param pos_column variant position column in \code{variant_data}
#' @param y_column p-value column in \code{variant_data}
#' @param marker_column unique variant identifier column in \code{variant_data}
#' @param ld_data variant LD values, the output of LDGds workflow
#' @param ld_ref reference variant unique identifier for LD
#' @param genome_build genome build of variants
#' @param should_log should the pvalues in \code{variant_data} be -log10'ed?
#' @param horizontal_value Y-intercept for horizontal line in Manhattan plot
#' @param horizontal_color Color for horizontal line in Manhattan plot
#' @param horizontal_style Style for horizontal line in Manhattan plot
#' @param background_color Manhattan plot title box color
#' @param background_frame Should a frame be placed around the Manhattan plot?
#' @param point_color Color of points if no LD information is provided
#' @param title Manhattan plot title
#' @return \code{variant.track} gviz track

#' @importFrom Gviz DataTrack OverlayTrack


make_variant_track <- function(variant_data, chr_column, pos_column, y_column,
                             marker_column = NULL,
                             ld_data = NULL,
                             ld_ref = NULL,
                             genome_build = "hg19",
                             should_log = T,
                             horizontal_value = 5e-8,
                             horizontal_color = "red",
                             horizontal_style = "dashed",
                             background_color = "#e6550d",
                             background_frame = T,
                             point_color = "#252525",
                             title = " "){

  # log the data if we need to
  if (should_log){
    variant_data[,y_column] <- -log10(variant_data[,y_column])
    horizontal_value <- -log10(horizontal_value)
  }

  # merge with LD information if we need to
  if (!is.null(marker_column) & !is.null(ld_data) & !is.null(ld_ref)){
    # first make sure the markername column is correct in the ld data file
    names(ld_data)[1] <- marker_column

    # make sure that we have the same sep character in the marker names
    ld_data[,marker_column] <- gsub("[[:punct:]]", "_", ld_data[,marker_column])
    ld_ref <- gsub("[[:punct:]]", "_", ld_ref)
    variant_data[,marker_column] <- gsub("[[:punct:]]", "_", variant_data[,marker_column])

    # then merge
    variant_data <- merge(variant_data,
                          ld_data,
                          by = marker_column,
                          all.x = T)

    # column number with ld information
    ld_col <- ncol(variant_data)

    # generate the colors
    variant_data$color <- "#B8B8B8"
    v.na <- variant_data[is.na(variant_data[,ld_col]),]
    v.notna <- variant_data[!is.na(variant_data[,ld_col]),]
    v.notna[v.notna[,ld_col] >= 0, "color"] <- "#357EBD"
    v.notna[v.notna[,ld_col] >= 0.2, "color"] <- "#46B8DA"
    v.notna[v.notna[,ld_col] >= 0.4, "color"] <- "#5CB85C"
    v.notna[v.notna[,ld_col] >= 0.6, "color"] <- "#EEA236"
    v.notna[v.notna[,ld_col] >= 0.8, "color"] <- "#D43F3A"

    variant_data <- rbind(v.na, v.notna)
    variant_data[variant_data[,marker_column] == ld_ref, "color"] <- "#9632B8"

    colors <- c("#B8B8B8", "#357EBD", "#46B8DA", "#5CB85C", "#EEA236", "#D43F3A", "#9632B8")

    # get y bounds for plot
    ylims <- c(0, max(horizontal_value, max(variant_data$pvalue)))

    # generate the variant track list, to combine later
    variant.tracks <- list()
    for (v in 1:length(colors)){
      val <- colors[v]
      cur_data <- variant_data[variant_data$color == val,]
      if (nrow(cur_data) == 0) next
      variant.tracks[[val]] <- DataTrack(
        data = cur_data[,y_column],
        start = cur_data[,pos_column],
        end = cur_data[,pos_column],
        chromosome = cur_data[,chr_column],
        genome = genome_build,
        name = title,
        background.title = background_color,
        col.frame = background_color,
        frame = background_frame,
        baseline = horizontal_value,
        col.baseline = horizontal_color,
        lty.baseline = horizontal_style,
        col = val,
        ylim = ylims)
    }

    variant.track <- OverlayTrack(variant.tracks,
                                  background.title = background_color,
                                  col.frame = background_color,
                                  frame = background_frame,
    )

  } else {
    print("No LD information found.")
    # generate the variant track
    variant.track <- DataTrack(
      data = variant_data[,y_column],
      start = variant_data[,pos_column],
      end = variant_data[,pos_column],
      chromosome = variant_data[,chr_column],
      genome = genome_build,
      name = title,
      background.title = background_color,
      col.frame = background_color,
      frame = background_frame,
      baseline = horizontal_value,
      col.baseline = horizontal_color,
      lty.baseline = horizontal_style,
      col = point_color)
  }
  return(variant.track)
}
