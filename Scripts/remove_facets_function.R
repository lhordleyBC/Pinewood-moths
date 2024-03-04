remove_facets <- function(plot, layout) {
  layout <- strsplit(layout, split = '\n')[[1]]
  layout <- lapply(layout, trimws)
  layout <- matrix(unlist(sapply(layout, strsplit, "")),
                   nrow = length(layout), byrow = T)
  layout <- which(layout == "#", arr.ind = TRUE)
  prm <- apply(layout,1,\(x) {
    c(glue::glue("panel-{x[1]}-{x[2]}"),
      glue::glue("strip-t-{x[2]}-{x[1]}"))
  })
  # https://stackoverflow.com/a/30372692/1296582
  g <- ggplot2::ggplotGrob(plot)
  rm_grobs <- g$layout$name %in% prm
  g$grobs[rm_grobs] <- NULL
  g$layout <- g$layout[!rm_grobs, ]
  ggpubr::as_ggplot(g)
}
