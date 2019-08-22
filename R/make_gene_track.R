#' Make a mart object for plotting
#' @param chr interval chromosome
#' @param start interval start
#' @param end interval end
#' @param mart mart object
#' @param genome_build either hg19 or hg38
#' @param highlight what gene to higlight
#' @param title title for sidebar
#' @param background_color color of sidebar
#' @param background_frame whether to add frame
#' @return \code{gene.track} gviz track

#' @importFrom Gviz BiomartGeneRegionTrack
#' @importFrom biomaRt getBM

make_gene_track <- function(chr, start, end, mart,
                          genome_build = "hg19",
                          highlight = NULL,
                          title = "Genomic Context",
                          background_color = "#31a354",
                          background_frame = T){

  depricated <- T
  if (!depricated){
    # query biomart for the gene positions and symbols
    gene.data <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),filters = c("chromosome_name"), values=list(chr), mart=mart)

    # subset to genes in the region
    gene.data <- gene.data[(gene.data$end_position >= start) & (gene.data$start_position <= end),]
    # gene.data[gene.data$hgnc_symbol == "", "hgnc_symbol"] <- gene.data[gene.data$hgnc_symbol == "", "ensembl_gene_id"]
    gene.data <- gene.data[gene.data$hgnc_symbol != "",]

    # define a color for non highlighted genes (tan)
    gene.data$color <- "#ffd470"

    # change the color of highlighted genes
    if (!is.null(highlight)){
      highlight <- c(highlight)
      gene.data[gene.data$hgnc_symbol %in% highlight, "color"] <- "#de2d26"
    }

    gene.track <- AnnotationTrack(
      name = title,
      start = gene.data$start_position,
      end = gene.data$end_position,
      id = gene.data$hgnc_symbol,
      feature = gene.data$color,
      chromosome = gene.data$chromosome_name[1],
      genome = genome_build,
      col="transparent",
      background.title = background_color,
      col.frame = background_color,
      frame = background_frame,
      groupAnnotation="id",
      "#ffd470" = "#ffd470",
      "#0fa3ff" = "#0fa3ff",
      "#de2d26" = "#de2d26"
    )
  } else {
    gene.track <- BiomartGeneRegionTrack(
      name = title,
      start = start,
      end = end,
      chromosome = chr,
      genome = genome_build,
      background.title = background_color,
      col.frame = background_color,
      frame = background_frame,
      collapseTranscripts="meta",
      transcriptAnnotation="symbol",
      biomart = mart,
      cex = 1.2
      # shape="arrow",
    )
  }
  return(gene.track)
}
