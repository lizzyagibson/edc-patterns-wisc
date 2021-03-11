base_size = 12
base_family = 'serif'

col = as.character(color_scheme_get("mix-green-brightblue", i = c(3,4)))
fill = as.character(color_scheme_get("mix-green-brightblue", i=c(3,4)))
color_scheme_view("mix-green-brightblue")

options(
  ggplot2.discrete.color = col,
  ggplot2.discrete.fill = fill)

theme_set(theme_bw(
  base_family = base_family,
  base_size = base_size
) +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.4),
    axis.ticks = element_line(size = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(0.9)),
    strip.placement = "outside",
    panel.spacing = unit(1.5, "lines"),
    legend.position = "right",
    legend.background = element_blank(),
    legend.text = element_text(size = 13),
    legend.text.align = 0,
    legend.key = element_blank()
  ))
