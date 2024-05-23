require(tidyverse)
#fix stat_ellipse
StatClipEllipse <- ggproto("StatClipEllipse", Stat,
                           required_aes = c("x", "y"),
                           compute_group = function(data, scales, type = "t", level = 0.95,
                                                    segments = 51, na.rm = FALSE) {
                             xx <- ggplot2:::calculate_ellipse(data = data, vars = c("x", "y"), type = type,
                                                               level = level, segments = segments)
                             xx %>% mutate(y=pmin(y, 1))
                           }
)
stat_clip_ellipse <- function(mapping = NULL, data = NULL,
                              geom = "path", position = "identity",
                              ...,
                              type = "t",
                              level = 0.95,
                              segments = 51,
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatClipEllipse,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      type = type,
      level = level,
      segments = segments,
      na.rm = na.rm,
      ...
    )
  )
}