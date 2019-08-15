#' Fix chromosome encoding
#' @param chr1 vector of reference chromosome names
#' @param chr2 chromosome encodings to change

fix_chr <- function(chr1, chr2){
  if (startsWith(chr1[1], "chr") & !startsWith(chr2[1], "chr")){
    chr2 <- paste0("chr", chr2)
  } else if (!startsWith(chr1[1], "chr") & startsWith(chr2[1], "chr")){
    chr2 <- sub("^chr", "", chr2)
  }
  return(chr2)
}
