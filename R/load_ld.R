#' Load and parse an LD matrix in the format returned by LDGds workflow
#' @param file filepath to LD data, stored as a plain, delimited text file
#' @param df dataframe of preloaded ld data
#' @param ld_ref if providing a matrix of LD values rather than a single row/column, reference variant identifier that is a column name in the LD matrix. This is the LD reference that you will want to plot the values from
#' @return \code{ld.df} 2 column dataframe that is input to the make_regional_plot function
#' @export
#' @importFrom data.table fread


load_ld <- function(file = NULL, df = NULL, ld_ref = NULL){
  # load the file
  if (!is.null(file)){
    ld.data <- fread(file, data.table = F, stringsAsFactors = F)
  } else if (!is.null(df)){
    ld.data <- df
  }

  # check the dimensions
  ld.dim <- dim(ld.data)

  # if a one row matrix
  if (ld.dim[1] == 1){
    # remove non-numeric values
    ld.matrix <- suppressWarnings(as.numeric(ld.data))

    # make the data frame
    ld.df <- data.frame(
      MarkerName = names(ld.data)[!is.na(ld.matrix)],
      ld = ld.matrix[!is.na(ld.matrix)],
      stringsAsFactors = F
    )
  } else if (ld.dim[1] != 1 & ld.dim[2] == 2){
    # only need to fix the names
    names(ld.df) <- c("MarkerName", "ld")
  } else if (ld.dim[1] > 1 & ld.dim[2] > 1){
    # fix the name
    names(ld.data)[1] <- "MarkerName"

    # make sure that marker names are in the right format
    ld.data[,"MarkerName"] <- gsub("[[:punct:]]", "_", ld.data[,"MarkerName"])
    names(ld.data)[2:ncol(ld.data)] <- gsub("[[:punct:]]", "_", names(ld.data)[2:ncol(ld.data)])
    ld_ref <- gsub("[[:punct:]]", "_", ld_ref)

    # fix chr encoding
    if (startsWith(ld_ref, "chr") | startsWith(ld.data[1,"MarkerName"], "chr")){
      ld_ref <- fix_chr(ld.data[1,"MarkerName"], ld_ref)
    }

    # make sure ld ref is in the matrix
    if (!(ld_ref %in% names(ld.data))){
      warning("LD reference variant not found in LD matrix, using first column.")
      ld_ref <- names(ld.data)[2]
    }

    # take only the ld ref column
    ld.df <- ld.data[, c("MarkerName", ld_ref)]
  }

  return(ld.df)
}





