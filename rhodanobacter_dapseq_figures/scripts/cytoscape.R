source("scripts/import.R")

# make sure cytoscape is open before running this script...

# setting up the network for 10B01 ----
b01_cogs <- parseCOGforCytoscape(eggNOG = "data/eggNOG/FW104-10B01/out.emapper.annotations.xlsx") |> 
            tidyr::separate(col = query, into = c("query", "genome"), sep = "\\|") |> 
            dplyr::select(-genome)

b01_edges <- tf_table |> 
  dplyr::filter(Genome == "FW104-10B01") |> 
  dplyr::select(TF, `All genes of target operon`) |> 
  tidyr::separate_longer_delim(delim = ",", cols = `All genes of target operon`) |> 
  dplyr::rename(source = TF, target = `All genes of target operon`) |> 
  dplyr::mutate(interaction = "regulates") |> 
  dplyr::left_join(b01_cogs, by = c(target = "query"))

b01_nodes <- b01_edges |>
         dplyr::select(-interaction) |>
          tidyr::pivot_longer(cols = source:target, values_to = "id",names_to = "type") |> 
          dplyr::left_join(tf_lookup_table, by = c(id = "TF")) |> 
          dplyr::rename(label = TF_Name) |>
          dplyr::mutate(isTF = dplyr::if_else(!is.na(label), "yes", "no")) |> 
          dplyr::relocate(id) |> 
          dplyr::distinct()

plotCytoscape(nodes.df = b01_nodes,edges.df = b01_edges,  network.title = "10B01")


# setting up the network for R12 ----
r12_cogs <- parseCOGforCytoscape(eggNOG = "data/eggNOG/FW510-R12/out.emapper.annotations.xlsx") |> 
  tidyr::separate(col = query, into = c("query", "genome"), sep = "\\|") |> 
  dplyr::select(-genome)

r12_edges <- tf_table |> 
  dplyr::filter(Genome == "FW510-R12") |> 
  dplyr::select(TF, `All genes of target operon`) |> 
  tidyr::separate_longer_delim(delim = ",", cols = `All genes of target operon`) |> 
  dplyr::rename(source = TF, target = `All genes of target operon`) |> 
  dplyr::mutate(interaction = "regulates") |> 
  dplyr::left_join(r12_cogs, by = c(target = "query"))

r12_nodes <- r12_edges |>
  dplyr::select(-interaction) |>
  tidyr::pivot_longer(cols = source:target, values_to = "id",names_to = "type") |> 
  dplyr::left_join(tf_lookup_table, by = c(id = "TF")) |> 
  dplyr::rename(label = TF_Name) |>
  dplyr::mutate(isTF = dplyr::if_else(!is.na(label), "yes", "no")) |> 
  dplyr::relocate(id) |> 
  dplyr::distinct()

plotCytoscape(edges.df = r12_edges, nodes.df = r12_nodes, network.title = "R12")
