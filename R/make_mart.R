#' Make a mart object for plotting
#' @param genome_build either hg19 or hg38
#' @return \code{mart} mart object

#' @importFrom biomaRt useMart

make_mart <- function(genome_build = "hg19"){
  genome_build <- tolower(genome_build)
  if (genome_build == "hg19" | genome_build == "grch37"){
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  } else if (genome_build == "hg38" | genome_build == "grch38") {
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  } else {
    warning("Unrecognized genome build, defaulting to hg19.")
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  }
  return(mart)
}
