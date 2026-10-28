# generate the circos plot visualization for the genomes
drawCircos <- function(genomeInform = b01.genome, tfParameters = b01_data_tfs_subset, linksOrigins = link_tfs, linksTargets = link_targets, filename = "output/genome_test.svg") {
  
  svg(filename, width = 8, height = 6)
  
  # setting up graphical parameters of the output
  circle_size <- grid::unit(1, "snpc")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(x = 0,
                                    y = 0.5,
                                    width = circle_size,
                                    height = circle_size,
                                    just = c("left", "center")))
  suppressWarnings(
    par(omi = gridBase::gridOMI(), new = TRUE)
  )
  
  # initializing the {circlize} functions
  circos.clear()
  circos.par("start.degree" = 90, "track.height" = 0.05, "gap.after" = 0)
  circos.genomicInitialize(genomeInform, # uses the length information to draw the circle correctly
                           major.by = 0.5e6,
                           tickLabelsStartFromZero = TRUE,
                           labels = NULL,
                           axis.labels.cex = 0.8)
  
  
  # adding a track for the TFs using genomicPoints
  circos.genomicTrackPlotRegion(tfParameters, # dataframe with the genomic coordinates and other graphical parameters
                                stack = TRUE,
                                bg.border = NA,
                                panel.fun = function(region, value, ...) {
                                  
                                  circos.lines(get.cell.meta.data("xlim"), c(0, 0), col = "black", lwd = 1)
                                  circos.lines(get.cell.meta.data("xlim"), c(2,2), col = "black", lwd = 1)
                                  
                                  circos.genomicPoints(region, value, cex = 2, pch = as.numeric(value$pch),
                                                       bg = paste0(value$col, "AA"), # adds some transparency to the points
                                                       col = "black", ...)
                                  })
  
  # adding a track for the links between TFs and targets
  circos.genomicLink(              # the circos plot targets network expects 2 tables with the same nrow()
     linksOrigins,                 # 1 - all 'sources' (TFs)
     linksTargets,                 # 2 - all targets
     col = linksOrigins$col,    
     lwd = linksOrigins$lwd, 
     directional = 1,
     arr.type = "triangle",
     arr.length = 0.15,
     arr.lwd = 0.1,
     border = linksOrigins$col,
     h.ratio = 0.925)
  
  # Finally we have to 'manually' create a legend for the plot
  legend_data <- linksOrigins |> 
    dplyr::select(TF_Name, col, pch) |> 
    dplyr::distinct()
  
  # Creating a custom legend
  legend <-  ComplexHeatmap::Legend(at = legend_data$TF_Name,
                                    type = "points",
                                    pch = as.numeric(legend_data$pch),
                                    size = grid::unit(4.5, "mm"), 
                                    background = NA,
                                    row_gap = grid::unit(2, "mm"),
                                    legend_gp = grid::gpar(fill = legend_data$col),
                                    title_position = "topleft",
                                    title = "Transcription factors")
  
  # drawing the custom legend
  grid::upViewport()
  ComplexHeatmap::draw(legend, x = circle_size, y = circle_size /2, just = "left")
  
  # closing and saving the file
  dev.off()
}

# defines a function to generate the list of targets for each of the TFs to be used in the Circos plot
getTfTargetLinks <- function(fullDataset = tf_table, formattedDataset = b01_data, genomeSource = "FW104-10B01") {
  
  fullDataset |> 
    dplyr::filter(Genome == genomeSource) |>
    dplyr::group_split(TF_Name) |> 
    
    purrr::map( ~{
      targets <- dplyr::pull(.x, Target) |> unique()
      tfs   <- dplyr::pull(.x, TF) |> unique()
      
      targets_formatted <- dplyr::filter(formattedDataset,
                                         Gene %in% targets,
                                         !is.na(Start), # check out which is the gene causing trouble...
                                         !is.na(End))
      
      tf_formatted <- dplyr::filter(formattedDataset, Gene %in% tfs) |> 
                      dplyr::slice(rep(1, nrow(targets_formatted)))
      
      list(targets = targets_formatted, TF = tf_formatted)
    })
}

# defines a function to generate a simple table to initialize the Circos genome track
makeGenomicDataset <- function(dataset = b01_formatted, chrName = "FW104-10B01") {
  dataset |> 
    dplyr::mutate(chr = chrName, Start = as.numeric(Start), End = as.numeric(End)) |> 
    dplyr::relocate(chr, Start, End)
}

# defines a function for using edges and nodes to plot a cytoscape network
plotCytoscape <- function(edges.df = edges, nodes.df = nodes, network.title = "10B01", cogs.df = all_cogs) {
  
  # this resets the current cytoscape session
  commandsRun('session new destroyCurrentSession=true')
  
  # Creating the Cytoscape network live
  createNetworkFromDataFrames(nodes = nodes.df, edges = edges.df, title = network.title, collection = "basecollection") # removed the node argument and it worked --> investigate
  
  # preparing our custom style we will apply to the network
  style.name <- "go_terms"
  defaults <- list(EDGE_PAINT="black", NODE_LABEL_POSITION="C,C,C,0,0", # anchor points. Look the UI for the opts
                   EDGE_TARGET_ARROW_SHAPE = "Arrow", NODE_LABEL_FONT_SIZE = 24)
  
  label <- mapVisualProperty('node label', 'label', 'p')
  fill <- mapVisualProperty('node fill color', 'isTF','d',c("yes", "no"), c("#FFFFFF","#000000"))
  shape <- mapVisualProperty('node shape','isTF','d',c("yes", "no"), c("Rectangle","Circle"))
  size <- mapVisualProperty('node size','isTF','d',c("yes", "no"), c("75","40"))
  
  createVisualStyle(style.name, defaults, list(label, shape, size, fill))
  
  # Setting visual properties for nodes with pie charts
  piechart_columns <- nodes.df |>
    dplyr::select(dplyr::where(is.numeric)) |>
    colnames()
  
  cog_colors  <- cogs.df |> dplyr::pull(color)
  
  tfs <- nodes.df |>
    dplyr::filter(!is.na(label)) |>
    dplyr::pull(id) |>
    unique()
  
  not_assig <- nodes.df |>
    dplyr::filter(dplyr::if_any(dplyr::where(is.numeric), is.na),
                  !id %in% tfs) |>
    dplyr::pull(id) |>
    unique()
  
  targets <- setdiff(getAllNodes(), c(tfs,not_assig))
  
  # applying the style and overlaying the pie charts
  setVisualStyle(style.name)
  
  # removes piecharts for the TFs
  for (tf in tfs) {
    setNodePropertyBypass(node.names = tf, visual.property = "NODE_CUSTOMGRAPHICS_1", new.values = "null")
    setNodePropertyBypass(node.names = tf, visual.property = "NODE_BORDER_WIDTH", new.values = 5)
  }
  
  # formats targets that did not have any assignments
  for (node in not_assig) {
    setNodePropertyBypass(node.names = node, visual.property = "NODE_BORDER_STROKE", new.values = "EQUAL_DASH")
    setNodePropertyBypass(node.names = node, visual.property = "NODE_FILL_COLOR", new.values = "#FFFFFF")
    setNodePropertyBypass(node.names = node, visual.property ="NODE_BORDER_WIDTH", new.values = 2)
  }
  
  # formats the nodes underneath piecharts
  for (target in targets) {
    setNodePropertyBypass(node.names = target, visual.property ="NODE_BORDER_WIDTH", new.values = 0)
  }
  
  setNodeCustomPieChart(piechart_columns,
                        slot = 1,
                        colors = cog_colors,
                        startAngle = 45,
                        style.name = style.name)
  
  makeCOGlegend(nodes.df = nodes.df, title = network.title, cogs.list = cogs.df)
  
  message("DONE!")
  message("In Cytoscape, adjust the nodes as needed and export the network as an SVG file!")
}

# defines a function for parsing the the eggNog output excel files
parseCOGforCytoscape <- function(eggNOG = "data/eggNOG_10b01/out.emapper.annotations.xlsx") {
  
  # this function switches the "-" assignment to "X" because it seems like cytoscape doesn't like the dash
  # as a column/value name
  
  readxl::read_xlsx(eggNOG, skip = 2) |> 
    dplyr::select(query, COG_category) |>
    dplyr::filter(stringr::str_starts(query, pattern = "##", negate = TRUE)) |>
    dplyr::mutate(COG_category = stringr::str_replace_all(COG_category, "-", "X")) |> # updates the "-" assignment
    dplyr::cross_join(all_cogs) |> # helps identify multiple COG assignments per gene 
    dplyr::select(-class, -description, -color) |> 
    dplyr::mutate(match = dplyr::if_else(stringr::str_detect(COG_category, code), 1, NA)) |> 
    tidyr::pivot_wider(id_cols = query, names_from = code, values_from = match) |> 
    dplyr::group_split(query) |> 
    purrr::map(~ {
      
      transposed <- .x |> tidyr::pivot_longer(cols = dplyr::where(is.numeric), names_to = "code", values_to = "match")
      
      valid_cogs <- as.numeric(!is.na(transposed$match))
      prop_cogs <- valid_cogs/length(valid_cogs[valid_cogs > 0])
      
      transposed |>
        dplyr::mutate(match = prop_cogs) |> 
        tidyr::pivot_wider(names_from = code, values_from = match)
    }) |> 
    dplyr::bind_rows() 
}

# defines a function for creating a legend based on the current cytoscape network
makeCOGlegend <- function(nodes.df = nodes, out.dir = "output/cytoscape/", title = network.title, cogs.list = all_cogs) {
  cogs_in_plot <- nodes.df |>
    dplyr::select(dplyr::where(is.numeric)) |> 
    dplyr::filter(!dplyr::if_any(dplyr::where(is.numeric), is.na)) |>
    janitor::remove_constant() |> 
    colnames()
  
  cog_legend_table <- cogs.list |>
    dplyr::filter(code %in% cogs_in_plot) |> 
    dplyr::select(-code) |> 
    dplyr::distinct() |> 
    dplyr::add_row(class = "NONE", description = "No assignment", color = "#FFFFFF") |> 
    dplyr::mutate(
      lty = dplyr::if_else(class == "NONE", 2, 1),
      lwd = dplyr::if_else(class == "NONE", 0.8, 1),
    )
  
  graphical_pars <- cog_legend_table |> 
    dplyr::select(-class, -description) |>
    purrr::pmap(function(color, lty, lwd) {
      
      fun_call <- glue::glue("function(x, y, w, h, color, lty, lwd) {{
                                   grid::grid.circle(x, y, r = 0.02, gp = gpar(fill = '{color}',
                                                                               col = 'black',
                                                                               lty = {lty},
                                                                               lwd = {lwd}))}}")
      eval(parse(text = fun_call))
    })
  
  legend <- ComplexHeatmap::Legend(labels = cog_legend_table$description,
                                   legend_gp = gpar(fill = cog_legend_table$color),
                                   title = "COG Terms", background = NA,
                                   row_gap = grid::unit(2, "mm"),
                                   graphics = graphical_pars)
  
  dir.create(out.dir, showWarnings = FALSE)
  
  out_path <- glue::glue("./{out.dir}{title}_legend.svg")
  
  svg(out_path, width = 6, height = 8)
  ComplexHeatmap::draw(legend)
  dev.off()
  
  message(glue::glue("Legend created at {out.dir}"))
}
