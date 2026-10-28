source("scripts/import.R")

r12_families <- tf_table |> 
                dplyr::filter(Genome == "FW510-R12") |>
                dplyr::pull(TF_Name) |>
                unique()
                
b01_families <- tf_table |> 
  dplyr::filter(Genome == "FW104-10B01") |>
  dplyr::pull(TF_Name) |>
  unique()

r12_exclusive <- setdiff(r12_families, b01_families)
b01_exclusive <- setdiff(b01_families, r12_families)
intersection <- intersect(r12_families, b01_families)

# plotting
vennInfo <- data.frame(x = c(0.65, 1.5, 2.35),
                       y = rep(0, 3),
                       label = c(stringr::str_flatten(b01_exclusive, collapse = "\n"),
                                 stringr::str_flatten(intersection, collapse = "\n"),
                                 stringr::str_flatten(r12_exclusive, collapse = "\n")),
                       title = c("FW104-10B01", NA, "FW510-R12"),
                       x_title = c(0.4, NA, 2.5))

circles <- data.frame(x0 = c(1,2), y0 = c(0,0), r = c(1,1), groups = c("FW104-10B01", "FW510-R12"))

vennDiag <- ggplot(vennInfo) +
  ggforce::geom_circle(data = circles,
                       aes(x0 = x0, y0 = y0, r = r, fill = groups),
                       radius = 3, color = NA, alpha = 0.6, show.legend = FALSE) +
  geom_text(data = vennInfo,
            aes(x = x, y = y, label = label), size = 5) +
  geom_text(data = vennInfo,
            aes(x = x_title, y = y + 1.1, label = title), size = 8) +
  scale_fill_manual(values = strain_colors) +
  coord_fixed() +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(vennDiag, filename = "output/Figure 3_venn_diagram.svg", device = svg, width = 6, height = 6, units = "in", dpi = 300)

