source("scripts/import.R")

# Making the circos plot for FW104-10B01 ----
b01_genomeInfo <- tibble::tibble(chr = "FW104-10B01", Start = 0, End = b01_size)
                           
b01_data <- makeGenomicDataset(b01_formatted, chrName = "FW104-10B01")
  
# preparing the B01 dataset for the genomic track with the symbols
b01_data_tfs_subset <- b01_data |> 
                       dplyr::filter(Gene %in% tf_lookup_table$TF) |> 
                       dplyr::distinct() |> 
                       dplyr::left_join(tf_pars)

# Preparing the B01 dataset for the links
 b01_interactions <- getTfTargetLinks(formattedDataset = b01_data, genomeSource = "FW104-10B01")
b01_links_targets <- purrr::map(b01_interactions, 1) |> dplyr::bind_rows()
    b01_links_tfs <- purrr::map(b01_interactions, 2) |>
                      dplyr::bind_rows() |> 
                      dplyr::left_join(tf_pars)

drawCircos(genomeInform = b01_genomeInfo,
           tfParameters = b01_data_tfs_subset,
           linksOrigins = b01_links_tfs,
           linksTargets = b01_links_targets,
           filename = "output/Figure 3_b01_circos.svg")


# Making the circos plot for FW510-R12 ----
r12_genomeInfo <- tibble::tibble(chr = "FW510-R12", Start = 0, End = r12_size)

r12_data <- makeGenomicDataset(r12_formatted, chrName = "FW510-R12")

# preparing the B01 dataset for the genomic track with the symbols
r12_data_tfs_subset <- r12_data |> 
  dplyr::filter(Gene %in% tf_lookup_table$TF) |> 
  dplyr::distinct() |> 
  dplyr::left_join(tf_pars)

# Preparing the B01 dataset for the links
r12_interactions <- getTfTargetLinks(formattedDataset = r12_data, genomeSource = "FW510-R12")
r12_links_targets <- purrr::map(r12_interactions, 1) |> dplyr::bind_rows()
r12_links_tfs <- purrr::map(r12_interactions, 2) |>
  dplyr::bind_rows() |> 
  dplyr::left_join(tf_pars)

drawCircos(genomeInform = r12_genomeInfo,
           tfParameters = r12_data_tfs_subset,
           linksOrigins = r12_links_tfs,
           linksTargets = r12_links_targets,
           filename = "output/Figure 3_r12_circos.svg")
