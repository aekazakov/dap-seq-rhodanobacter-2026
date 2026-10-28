library("ggplot2")

growthData <- readr::read_csv("data/growth_data/FW104_B4.csv") |> 
  tidyr::pivot_longer(cols = dplyr::contains(c("10B01", "clp-")), names_to = "strain", values_to = "OD", names_prefix = "R. denitrificans ") |> 
  tidyr::separate(strain, sep = "\\.\\.\\.", into = c("strain", "replicate_code")) |> 
  dplyr::rename(time = `Time (h)`)

growthSummary <- growthData |> 
                  dplyr::group_by(strain, time) |> 
                  dplyr::summarise(OD_mean = mean(OD), OD_sd = sd(OD)) |> 
                  dplyr::mutate(
                    strain = dplyr::if_else(strain == "clp-", "tn::<i>clp</i>", "WT")
                  )

growth <- ggplot(growthSummary, aes(x = time, y = OD_mean, ymax = OD_mean+OD_sd, ymin = OD_mean-OD_sd, color = strain)) +
  geom_errorbar(linewidth = 1) +
  geom_point(size = 6) +
  labs(x = "Time (h)", y = "OD<sub>600</sub>") +
  theme_light(base_size = 20) +
  scale_color_manual(name = "Strain", values = c("gray20", "mediumpurple4")) +
  scale_x_continuous(breaks = c(12, 24, 36, 48, 60), expand = expansion(add = c(0.5, 1))) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4), limits = c(0, 0.4), expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  theme(
    
    axis.title.y = ggtext::element_markdown(),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 1),
    
    legend.title = element_text(size = 16),
    legend.text = ggtext::element_markdown(size = 16),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.85),
    legend.background = element_rect(color = "black", linewidth = 0.2),
    
    panel.border =  element_blank()
  )

ggsave(plot = growth, filename = "output/Figure S1_growth_curve.svg", device = "svg", width = 12, height = 7, units = "in")
