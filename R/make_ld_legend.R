#' Generate a regional plot from summary statistics for a desired interval
#' @param box.x x position of legend
#' @param box.ymin minimum y positions of legend
#' @param box.ymax maximum y positions of legend
#' @param box.width width of legend
#' @param box.colors vector of colors, one for each box in the legend
#' @param tic.width tic width
#' @param tic.text.offset.x x distance to offset tic label text
#' @param tic.text.offset.y y distance to offset tic label text
#' @param tic.text text to display next to each tic

#' @importFrom draw drawBox drawLine drawText

make_ld_legend <- function(box.x = 4,
                           box.ymin = 3.6,
                           box.ymax = 4.6,
                           box.width = NULL,
                           box.colors = c("#B8B8B8", "#357EBD", "#46B8DA", "#5CB85C", "#EEA236", "#D43F3A"),
                           tic.width = 0.1,
                           tic.text.offset.x = 0.05,
                           tic.text.offset.y = 0,
                           tic.text = c("NA", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0")){

  # calculate box height from y range
  box.height <- (box.ymax - box.ymin)/min(length(box.colors), (length(tic.text) - 1))

  # assign width if not provided
  if (is.null(box.width)) box.width <- box.height

  # generate box y positions
  box.ys <- seq(box.ymin, box.ymax, box.height)

  # loop through to make the legend
  for (bi in 1:length(box.colors)){
    drawBox(
      x = box.x,
      y = box.ys[bi],
      width = box.width,
      height = box.height,
      fillColor = box.colors[bi]
    )

    drawLine(
      x = c(box.x + tic.width, box.x + 2*tic.width),
      y = c(box.ys[bi] - box.height/2, box.ys[bi] - box.height/2)
    )

    drawText(
      x = box.x + 2*tic.width + tic.text.offset.x,
      y = box.ys[bi] - box.height/2 + tic.text.offset.y,
      text = tic.text[bi],
      just = "left"
    )
  }

  drawLine(
    x = c(box.x + tic.width, box.x + 2*tic.width),
    y = c(box.ys[bi] + box.height/2, box.ys[bi] + box.height/2)
  )

  drawText(
    x = box.x + 2*tic.width + tic.text.offset.x,
    y = box.ys[bi] + box.height/2 + tic.text.offset.y,
    text = tic.text[bi + 1],
    just = "left"
  )
}



