# "LRK53" --> R12 genes
# "LRK54" --> 10B01 genes

library(RCy3)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(UpSetR)
source("scripts/helpers.R")

# main tf/target data table ----
tf_table <- readxl::read_xlsx(path = "data/Supplementary_file_1_Supplementary_Tables.xlsx",
                                  sheet = 6, skip = 2) |> 
                dplyr::rename(TF = `TF gene`,
                              TF_Name = Regulator,
                              Target = `Gene associated with DAP-seq peak`)

tf_lookup_table <- tf_table |>
                   dplyr::select(TF_Name, TF) |>
                   dplyr::distinct()


# Preparing the data for the genome visualization ----
r12_gbff <- readLines("data/FW510-R12.gbff")
b01_gbff <- readLines("data/FW104-10B01.gbff")

# getting the proper lengths for the genome plots (from the gbff files)
r12_size <- r12_gbff[1] |> 
  stringr::str_extract(pattern = "\\d+(?=\\sbp)") |> 
  as.numeric()

b01_size <- b01_gbff[1] |> 
  stringr::str_extract(pattern = "\\d+(?=\\sbp)") |> 
  as.numeric()

# getting the genomic coordinates for the TFs
all_tfs <- unique(tf_table$TF)
all_targets <- unique(tf_table$Target)

b01_info <- NULL

for (tf in c(all_tfs, all_targets)) {
  matches <- stringr::str_which(b01_gbff, pattern = tf)
  
  # the first match is useful enough already
  out1 <- tibble::tibble(Gene = tf, coordinates = b01_gbff[matches[1]-1])
  out2 <- tibble::tibble(Gene = tf, coordinates = b01_gbff[matches[1]-2]) # in some cases there's an addditional line
  
  b01_info <- dplyr::bind_rows(b01_info, out1, out2)
}

r12_info <- NULL

for (tf in c(all_tfs, all_targets)) {
  matches <- stringr::str_which(r12_gbff, pattern = tf)
  
  # the first match is useful enough already
  out1 <- tibble::tibble(Gene = tf, coordinates = r12_gbff[matches[1]-1])
  out2 <- tibble::tibble(Gene = tf, coordinates = r12_gbff[matches[1]-2]) # in some cases there's an addditional line
  
  r12_info <- dplyr::bind_rows(r12_info, out1, out2)
}

# organizing the output

b01_formatted <- b01_info |>
  dplyr::filter(!is.na(coordinates),
                stringr::str_detect(coordinates, "\\sgene")) |> 
  dplyr::mutate(coordinates = stringr::str_extract(coordinates, "\\d+\\.\\.\\d+"),
                Start = stringr::str_extract(coordinates, "\\d+(?=\\.\\.)"),
                End = stringr::str_extract(coordinates, "(?<=\\.\\.)\\d+")) |> 
  dplyr::select(-coordinates) |> 
  dplyr::left_join(tf_lookup_table, by = c(Gene = "TF"))

r12_formatted <- r12_info |>
  dplyr::filter(!is.na(coordinates),
                stringr::str_detect(coordinates, "\\sgene")) |> 
  dplyr::mutate(coordinates = stringr::str_extract(coordinates, "\\d+\\.\\.\\d+"),
                Start = stringr::str_extract(coordinates, "\\d+(?=\\.\\.)"),
                End = stringr::str_extract(coordinates, "(?<=\\.\\.)\\d+")) |> 
  dplyr::select(-coordinates) |> 
  dplyr::left_join(tf_lookup_table, by = c(Gene = "TF"))

# organizing the COG key sheet from NCBI ----
cog_table <- readr::read_lines("data/cogs_keysheet.txt") |>
  stringr::str_subset(pattern = ".+", negate = FALSE) |> 
  stringr::str_squish() |>  purrr::map( ~{
    if (stringr::str_starts(.x, '[A-Z]')) {
      cogclass <- .x
      tibble::tibble(class = cogclass, name = NA)
    } else {
      name <- .x
      tibble::tibble(class = NA, name = .x)
    }
  }) |>
  dplyr::bind_rows() |>
  tidyr::fill(class, .direction = "down") |> 
  dplyr::filter(!is.na(name)) |>
  tidyr::separate_wider_regex(col = name, patterns = c(code = "\\[.\\]",
                                                       description = "\\s[A-Z].+$")) |> 
  dplyr::mutate(code = stringr::str_extract(code, pattern = "[A-Z]"),
                description = stringr::str_squish(description))



# Loading color vectors and everything else ----
# assigning colors for each strain for the Venn Diagrams (figure 2a) ----
strain_colors <- c(`FW104-10B01` = "#F6BB09", `FW510-R12` = "steelblue")

# Setting color parameters for the TFs plotted in the Circos tracks (figure 2b) ----
tf_target_count <- tf_table |> 
  dplyr::group_by(TF_Name) |> 
  dplyr::tally() |>  # we will prioritize each TF's # of targets for coloring
  dplyr::arrange(desc(n))

tf_names <- tf_target_count |>dplyr::pull(TF_Name)

custom_lwd <- scales::rescale(tf_target_count$n, to = c(0.1, 2))
custom_pch <- rep(c(23, 25, 21, 22, 24), length.out = length(tf_names))
custom_col <-rep(x = c("#C14040", "#DEAA10", "#125686", "#63B747", "#C180BF", "#FF6400"),
                 length.out = length(tf_names))

combinations <- dplyr::bind_cols(custom_col, custom_pch) |> dplyr::distinct()

tf_pars <- tibble::tibble(TF_Name = tf_names,
                          lwd = as.character(custom_lwd)) |> 
            dplyr::arrange(TF_Name) |>
            dplyr::mutate(col = as.character(combinations[[1]]),
                          pch = as.character(combinations[[2]]),
                          )


# Creating our own color assignments for each COG category (figure 3) ----
colors_cellProces <- c("#006CAD", "#02A82C", "#4E85BC", "#53B955", "#799FCB", "#7EC97B", "#A0B9DA", "#A5D9A0", "#C6D5E9", "#C9E8C5")
colors_informProc <- c("#63135F", "#B021A9", "#DE4FD6", "#EC9CE8", "#F4C3F1")
colors_metabolism <- c("#A80202", "#F0D116", "#BF472E", "#FAB460","#D37158", "#F09916", "#E49A85", "#F1C2B4")
colors_PoorlyChar <- c("#ADADAD", "#DBDBDB")

# Attributing each color to a cog category
all_colors <- c(colors_cellProces, colors_informProc, colors_metabolism, colors_PoorlyChar)
all_cogs <- cog_table |> 
  dplyr::arrange(class) |> 
  dplyr::bind_cols(color = all_colors) |> 
  dplyr::add_row(class = "POORLY CHARACTERIZED", code = "X", description = "Function unknown", color = "#DBDBDB")
