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
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene


make_gene_track <- function(chr, start, end, mart,
                          genome_build = "hg19",
                          highlight = NULL,
                          title = "Ensembl",
                          background_color = "#E69F00",
                          background_frame = T){

  # query to get gene coordinates
  if (genome_build %in% c("hg38", "grch38")){
    gene.db <- TxDb.Hsapiens.UCSC.hg38.knownGene
  } else {
    gene.db <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }

  gene.track <- GeneRegionTrack(
    gene.db,
    name = title,
    start = start,
    end = end,
    chromosome = chr,
    genome = genome_build,
    stacking = "squish",
    collapseTranscripts="meta",
    transcriptAnnotation="symbol",
    just.group = "above",
    background.title = background_color,
    col.frame = background_color,
    stackHeight=0.4,
    cex.title = 1.5
  )

  return(gene.track)
}
