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

#' @importFrom Gviz GeneRegionTrack
#' @importFrom EnsDb.Hsapiens.v75 EnsDb.Hsapiens.v75
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom ensembldb getGeneRegionTrackForGviz


make_gene_track <- function(chr, start, end, mart,
                          genome_build = "hg19",
                          highlight = NULL,
                          title = "Ensembl",
                          background_color = "#E69F00",
                          background_frame = T){

  # query to get gene coordinates
  if (genome_build %in% c("hg38", "grch38")){
    gene.db <- EnsDb.Hsapiens.v86
  } else {
    gene.db <- EnsDb.Hsapiens.v75
  }

  grt <- getGeneRegionTrackForGviz(
    gene.db,
    start = start,
    end = end,
    chromosome = chr
  )

  gene.track <- GeneRegionTrack(
    grt,
    name = title,
    genome = genome_build,
    stacking = "squish",
    collapseTranscripts = "meta",
    transcriptAnnotation = "symbol",
    just.group = "above",
    background.title = background_color,
    frame = background_frame,
    col.frame = background_color,
    stackHeight = 0.5,
    cex.title = 1.5,
    cex.group = 1.5
  )

  return(gene.track)
}
