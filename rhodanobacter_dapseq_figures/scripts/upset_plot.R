source("scripts/import.R")

# organizing for the upset plot
r12_upset <- tf_table |> 
  dplyr::filter(Genome == "FW510-R12") |> 
  dplyr::select(TF_Name, `All genes of target operon`) |>
  tidyr::separate_longer_delim(delim = ",", cols = `All genes of target operon`) |> 
  table() |> 
  t() |> 
  as.data.frame.matrix()

b01_upset <- tf_table |> 
  dplyr::filter(Genome == "FW104-10B01") |> 
  dplyr::select(TF_Name, `All genes of target operon`) |>
  tidyr::separate_longer_delim(delim = ",", cols = `All genes of target operon`) |>  
  table() |> 
  t() |> 
  as.data.frame.matrix()

# plots
r12_plot <- UpSetR::upset(r12_upset, nsets = ncol(r12_upset),
              mainbar.y.label = "Number of targets in FW510-R12",
              sets.x.label = "Intersection size",
              mb.ratio = c(0.55, 0.45),
              main.bar.color = "gray10",
              sets.bar.color = "gray10",
              matrix.color = "gray10",
              text.scale = c(1.8, 1.8, 1.6, 1.8, 1.8, 1.8))

svg(filename="output/Figure 3_r12_upset.svg", width = 8.5, height = 7.5, pointsize = 12)
r12_plot
dev.off()


b01_plot <- UpSetR::upset(b01_upset, nsets = ncol(b01_upset),
              mainbar.y.label = "Number of targets in FW104-10B01",
              sets.x.label = "Intersection size",
              mb.ratio = c(0.55, 0.45),
              main.bar.color = "gray10",
              sets.bar.color = "gray10",
              matrix.color = "gray10",
              text.scale = c(1.8, 1.8, 1.6, 1.8, 1.8, 1.8))


svg(filename="output/Figure 3_b01_upset.svg", width = 8.5, height = 7.5, pointsize = 12)
b01_plot
dev.off()
