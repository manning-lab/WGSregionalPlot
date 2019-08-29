#' Make a track for bed files
#' @param bed_data List of dataframes in bed format: chr start end label color
#' @param title title for sidebar
#' @param genome_build hg19 or hg38
#' @param background_color List of background colors for bed plot title bars
#' @param background_frame Should frames be added around each bed plot?
#' @return \code{bed.track} gviz track

#' @importFrom Gviz AnnotationTrack


make_bed_track <- function(bed_data, chr, start, end,
                         title = " ",
                         genome_build = "hg19",
                         background_color = "#56B4E9",
                         background_frame = T){
  # Assume that colors are stored in column 5
  # Need to make a better way of adding color to these annotation tracks

  # Make everything the same color if given a 4 column bed file
  if (ncol(bed_data) == 4){
    bed_data[,5] <- "#878787"
  }

  # fix chromosome encoding
  bed_data[,1] <- fix_chr(chr, bed_data[,1])

  # subset to only the range that we want
  bed_data <- bed_data[(bed_data[,1] == chr) & (bed_data[,2] <= end) & (bed_data[,3] >= start), ]

  bed.track <- AnnotationTrack(
    name = title,
    start = bed_data[,2],
    end = bed_data[,3],
    id = bed_data[,4],
    feature = bed_data[,5],
    chromosome = chr,
    genome = genome_build,
    stacking = "dense",
    col = "transparent",
    background.title = background_color,
    col.frame = background_color,
    frame = background_frame,
    stackHeight = 1,
    groupAnnotation = "feature",
    "255,255,255" = "#ffffff",
    "255,195,77" = "#ffc34d",
    "255,0,0" = "#ff0000",
    "255,69,0" = "#ff4500",
    "0,100,0" = "#006400",
    "255,255,0" = "#ffff00",
    "192,192,192" = "#c0c0c0",
    "205,92,92" = "#cd5c5c",
    "0,128,0" = "#008000",
    "128,128,128" = "#808080",
    "255,195,77" = "#ffc34d",
    "255,69,0" = "#ff4500",
    "194,225,5" = "#c2e105",
    "#ffffff" = "#ffffff",
    "#ffc34d" = "#ffc34d",
    "#ff0000" = "#ff0000",
    "#ff4500" = "#ff4500",
    "#006400" = "#006400",
    "#ffff00" = "#ffff00",
    "#c0c0c0" = "#c0c0c0",
    "#cd5c5c" = "#cd5c5c",
    "#008000" = "#008000",
    "#808080" = "#808080",
    "#ffc34d" = "#ffc34d",
    "#ff4500" = "#ff4500",
    "#c2e105" = "#c2e105"
  )
  return(bed.track)
}
