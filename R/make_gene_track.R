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

#' @importFrom Gviz BiomartGeneRegionTrack GeneRegionTrack
#' @importFrom biomaRt getBM
#' @importFrom AnnotationDbi loadDb

make_gene_track <- function(chr, start, end, mart,
                          genome_build = "hg19",
                          highlight = NULL,
                          title = "Ensembl",
                          background_color = "#E69F00",
                          background_frame = T){

  # query to get gene coordinates
  if (genome_build %in% c("hg19", "grch37")){
    sql.file <- paste0(system.file('sql', package='WGSregionalPlot'), "/ensembl_grch37.sqlite")
  } else {
    sql.file <- paste0(system.file('sql', package='WGSregionalPlot'), "/ensembl_grch38.sqlite")
  }

  gene.db  <- loadDb(sql.file)

  gene.track <- GeneRegionTrack(
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
