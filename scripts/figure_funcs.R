rename_stage <- function(seurat_object, set_ident = TRUE) {
  levels(seurat_object@meta.data$stage) <- c(levels(seurat_object@meta.data$stage),
  	                                         "Cortico_medullary", "Ccl21a_high",
                                             "Early_Aire",
  	                                         "Aire_positive", "Late_Aire")

  seurat_object@meta.data$stage[seurat_object@meta.data$stage ==
                                "Immature"] <- "Ccl21a_high"
  seurat_object@meta.data$stage[seurat_object@meta.data$stage ==
                                "Intermediate"] <- "Early_Aire"
  seurat_object@meta.data$stage[seurat_object@meta.data$stage ==
                                "Mature"] <- "Aire_positive"
  seurat_object@meta.data$stage[seurat_object@meta.data$stage ==
                                "Late_mature"] <- "Late_Aire"

  seurat_object@meta.data$stage[seurat_object@meta.data$stage ==
                                "cTEC"] <- "Cortico_medullary"

  seurat_object@meta.data$stage <- factor(seurat_object@meta.data$stage,
  	                                      levels = c("Cortico_medullary",
                                                    "Ccl21a_high",
  	                                      	         "Early_Aire",
                                                     "Aire_positive",
  	                                      	         "Late_Aire",
                                                     "Tuft", 
                                                     "unknown"))

  if (set_ident){
    idents <- data.frame(seurat_object@ident)
    if (!identical(rownames(seurat_object@meta.data), rownames(idents))){
      seurat_object@meta.data <- seurat_object[match(rownames(idents),
                                                     rownames(mtec@meta.data))]
    }
    seurat_object <- Seurat::SetAllIdent(seurat_object, id = "stage")

    seurat_object@ident <- factor(seurat_object@ident,
  	                                        levels = c("Cortico_medullary",
                                                       "Ccl21a_high",
  	                                        	         "Early_Aire",
                                                       "Aire_positive",
  	                                        	         "Late_Aire",
                                                       "Tuft",
                                                       "unknown"))
  }
  return(seurat_object)
}

get_avg_exp <- function(mtec_obj, avg_expr_id = "stage") {
  idents <- data.frame(mtec_obj@ident)
  if (!identical(rownames(mtec_obj@meta.data), names(mtec_obj@ident))){
    mtec_obj@meta.data <- mtec_obj[match(rownames(idents),
                                                   rownames(mtec@meta.data))]
    }
  mtec_obj <- Seurat::SetAllIdent(mtec_obj, id = avg_expr_id)
  avg.expression <- log1p(Seurat::AverageExpression(mtec_obj))
  return(avg.expression)
}

plot_corr <- function(avg_expression_1, avg_expression_2, name_1,
                      name_2, color_df, density = FALSE) {
  stages_1 <- colnames(avg_expression_1)
  stages_2 <- colnames(avg_expression_2)
  stages <- intersect(stages_1, stages_2)
  
  cor_df <- NULL
  plots_list <- c()
  for (i in stages) {
    print(i)
    df_1 <- data.frame(row.names = rownames(avg_expression_1),
                       avg_exp = avg_expression_1[[i]])
    colnames(df_1)[1] <- name_1
    df_2 <- data.frame(row.names = rownames(avg_expression_2),
                       avg_exp = avg_expression_2[[i]])
    colnames(df_2)[1] <- name_2
    
    plot_df <- merge(df_1, df_2, by = "row.names")
    rownames(plot_df) <- plot_df$Row.names
    plot_df$Row.names <- NULL
    correlation <- cor(plot_df[[name_1]], plot_df[[name_2]])
    correlation_plot <- round(correlation, 2)
    correlation_plot <- paste0("r = ", correlation_plot)
    text_x <- max(plot_df[[name_1]]) - 1
    text_y <- max(plot_df[[name_2]]) / 2
    color_scale <- toString(color_df[[i]])
    
    if (density) {
      plot_df$density <- get_density(plot_df[[name_1]], plot_df[[name_2]])
      plot_1 <- ggplot2::ggplot(data = plot_df, ggplot2::aes_string(name_1,
                                                                    name_2)) +
        ggplot2::geom_point(ggplot2::aes(color = density)) +
        ggplot2::ggtitle(i) +
        ggplot2::scale_color_gradient(low = color_scale, high = "#A9A9A9") +
       # ggplot2::theme_classic() + 
        ggplot2::geom_text(x = text_x, y = text_y,
                                                      label = correlation_plot)
      #print(plot_1)
      
    } else {
      plot_1 <- ggplot2::ggplot(data = plot_df, ggplot2::aes_string(name_1,
                                                                    name_2)) +
        ggplot2::geom_point(color = color_scale) + ggplot2::ggtitle(i) +
        #ggplot2::theme_classic() + 
        ggplot2::geom_text(x = text_x, y = text_y,
                                                      label = correlation_plot)
      #print(plot_1)
    }
    if (i != "unknown") {
      plots_list[[i]] <- plot_1  
    }
    
    if (is.null(cor_df)) {
      cor_df <- data.frame(correlation)
      names(cor_df) <- i
    } else {
      cor_df[[i]] <- correlation
    }
  }
  nplots <- length(stages)
  if (nplots > 6) {
    print("Warning in plot_corr:")
    print(paste0("Works best if number of plots is less than 6. You have ",
                 nplots, " total plots"))
  }
  rows <- ceiling(nplots / 3)
  cols <- 3
  gridExtra::grid.arrange(grobs = plots_list, nrow = rows, ncol = cols)
  return(cor_df)
}

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

master_plot <- function(mtec_obj_1, name_1, mtec_obj_2,
                        name_2, stages_colors, density = FALSE) {
  avg_exp_1 <- get_avg_exp(mtec_obj_1)
  avg_exp_2 <- get_avg_exp(mtec_obj_2)
  cor_vals <- plot_corr(avg_exp_1, avg_exp_2, name_1, name_2, stages_colors)
  return(cor_vals)
}

populations_dfs <- function(seurat_object, sample_name, stage_df_all){
  stage_df <- data.frame(table(seurat_object@meta.data$stage))
  names(stage_df) <- c("stage", "count")
  stage_df$percent <- stage_df$count / sum(stage_df$count) * 100
  stage_df$sample <- sample_name
  if(is.null(stage_df_all)){
    stage_df_all <- stage_df
  } else {
    stage_df_all <- rbind(stage_df_all, stage_df)
  }
  return(stage_df_all)
}

populations_dfs_new <- function(seurat_object, sample_name, subsample = FALSE,
                                subsample_by = "exp", meta_data_col = "stage"){
  if (subsample) {
    cells_use <- rownames(seurat_object@meta.data)[
      seurat_object@meta.data[[subsample_by]] == sample_name]
    seurat_object <- Seurat::SubsetData(seurat_object, cells.use = cells_use)
  }
  stage_df <- data.frame(table(seurat_object@meta.data[[meta_data_col]]))
  names(stage_df) <- c("stage", "count")
  stage_df$percent <- stage_df$count / sum(stage_df$count) * 100
  stage_df$sample <- sample_name
  return(stage_df)
}

population_plots <- function(stage_df_all, color, save_plot = NULL){
  if (!(is.null(save_plot))){
    extension <- substr(save_plot, nchar(save_plot)-2, nchar(save_plot))
    if (extension == "pdf"){
      pdf(save_plot)
    } else if (extension == "png") {
      png(save_plot)
    } else {
      print("save plot must be .png or .pdf")
    }
  }
  plot_base <- ggplot2::ggplot(data = stage_df_all, ggplot2::aes_(x = ~sample,
                                                              y = ~percent,
                                                               fill = ~stage)) +
   # ggplot2::theme_classic() + 
    ggplot2::xlab("frequency")  +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = color, name = "stage")
  
  if (!(is.null(save_plot))){
    print(plot_base)
    dev.off()
  }
  return(plot_base)

}

gene_corr <- function(seurat_obj, gene_1, gene_2, stage_color) {
  seurat_df <- as.data.frame(as.matrix(seurat_obj@data))
  seurat_df <- as.data.frame(t(seurat_df[c(gene_1, gene_2), ]))
  genes <- sub("-", "_", c(gene_1, gene_2))
  gene_1 <- genes[1]
  gene_2 <- genes[2]
  seurat_obj@meta.data$stage <- factor (seurat_obj@meta.data$stage)
  names(seurat_df) <- sub("-", "_", names(seurat_df))
  if(!identical(rownames(seurat_df), rownames(seurat_obj@meta.data))){
    seurat_df <- seurat_df[match(rownames(seurat_obj@meta.data), 
                                 rownames(seurat_df)), ]
  }
  seurat_df$stage <- seurat_obj@meta.data$stage
  stages <- levels(seurat_df$stage)


  scatterPlot <- ggplot2::ggplot(data = seurat_df,
                                 ggplot2::aes_string(gene_1, gene_2)) +
    ggplot2::geom_point(ggplot2::aes(colour = stage)) +
    ggplot2::scale_color_manual(values = stage_color) + 
    #ggplot2::theme_classic() + 
    ggplot2::theme(legend.position= "none")
  
  xdensity <- ggplot2::ggplot(data = seurat_df, 
                              ggplot2::aes_string(gene_1)) +
    ggplot2::geom_density(ggplot2::aes(colour = stage)) +
    ggplot2::scale_color_manual(values = stage_color) + 
    #ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = "none")
  
  ydensity <- ggplot2::ggplot(data = seurat_df,
                              ggplot2::aes_string(gene_2)) +
    ggplot2::geom_density(ggplot2::aes(color = stage)) +
    ggplot2::scale_color_manual(values = stage_color) + 
    #ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = "none")
  
  ydensity <- ydensity + ggplot2::coord_flip()
  
  
  blankPlot <- ggplot2::ggplot() +
    ggplot2::geom_blank(ggplot2::aes(1,1))+
    ggplot2::theme(plot.background  = ggplot2::element_blank(), 
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), 
                   panel.border     = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.title.x     = ggplot2::element_blank(),
                   axis.title.y     = ggplot2::element_blank(),
                   axis.text.x      = ggplot2::element_blank(), 
                   axis.text.y      = ggplot2::element_blank(),
                   axis.ticks       = ggplot2::element_blank()
    )
  
  gridExtra::grid.arrange(xdensity, blankPlot, scatterPlot, ydensity,
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 5))
  
    
  #plot_1 <- ggplot2::ggplot(data = seurat_df, ggplot2::aes_string(gene_1,
   #                                                               gene_2)) +
  #  ggplot2::geom_point(ggplot2::aes(color = stage)) +
   # ggplot2::scale_color_manual(values = stage_color) + 
  #  ggplot2::theme_classic()
  #print(plot_1)
}

gene_exp_df <- function(combSubset, experiment, df_contents = "ranked") {

  # Find average expression of all genes at each stage
  avg.expression <- log1p(Seurat::AverageExpression(combSubset))
  
  if (df_contents == "ranked") {
  
    plot_df <- data.frame(genes = rownames(avg.expression))
  
    for (i in levels(combSubset@ident)) {
      print(i)
      # rank genes. Ties will all be given the lowest number. Ex.
      # c(1, 13, 5, 5, 7, 10)  wil be 1 6 2 2 4 5
      rank_expression <- rank(-avg.expression[[i]], ties.method = "min")
    
      plot_df[[i]] <- rank_expression
    }
    
    rownames(plot_df) <- plot_df$genes
    plot_df$genes <- NULL
    
  } else if ( df_contents == "expression") {
    plot_df <- avg.expression
  } else {
    print("df_contents must be either 'ranked' or 'expression'")
  }
  # Change the plot_df to be stages as rows and genes as columns
  plot_df <- as.data.frame(t(plot_df))
  plot_df$stage <- rownames(plot_df)
  plot_df$exp <- experiment
  return(plot_df)
}

plot_gene_exp <- function(plot_df_all, gene_name, low_lim = 0,
                          high_lim = 10000, col = NULL) {
  print(gene_name)
  stage <- "stage"
  experiment <- "exp"
  gene_plot <- ggplot2::ggplot(data = plot_df_all,
                               ggplot2::aes_string(experiment,
                                                   gene_name,
                                                   group = stage)) +
    ggplot2::geom_line(ggplot2::aes(color = stage)) +
    ggplot2::geom_point(ggplot2::aes(color = stage)) +
    ggplot2::ylim(low_lim, high_lim) 
    #ggplot2::theme_classic()

  if (is.null(col)) {
    gene_plot <- gene_plot + ggplot2::scale_color_brewer(palette = "Set1")
  } else {
    gene_plot <- gene_plot + ggplot2::scale_color_manual(values = col)
  }
    ggplot2::scale_color_brewer(palette = "Set1")
  gene_plot_2 <- ggplot2::ggplot(data = plot_df_all,
                                 ggplot2::aes_string(experiment,
                                                     gene_name,
                                                     group = stage)) +
    ggplot2::geom_line(ggplot2::aes(color = stage)) +
    ggplot2::geom_point(ggplot2::aes(color = stage)) 
    #ggplot2::theme_classic()
  if (is.null(col)) {
    gene_plot_2 <- gene_plot_2 + ggplot2::scale_color_brewer(palette = "Set1")
  } else {
    gene_plot_2 <- gene_plot_2 + ggplot2::scale_color_manual(values = col)
  }

  print(gene_plot)
  print(gene_plot_2)
}

plot_gene_set <- function(seurat_obj, gene_set, plot_name,
                          one_dataset = TRUE, data_set = NULL,
                          make_plot = TRUE, ...){
  print(head(gene_set))
  gene_set <- gene_set[gene_set %in%
                       rownames(seurat_obj@data)]
  mean_exp <- colMeans(as.matrix(seurat_obj@data[gene_set, ]), na.rm = TRUE)
  if (all(names(x = mean_exp) == rownames(x = seurat_obj@meta.data))) {
    print("Cell names order match in 'mean_exp' and 'object@meta.data': 
        adding gene set mean expression vaules in 'object@meta.data$gene.set.score'")
    seurat_obj@meta.data[[plot_name]] <- mean_exp
  }
  if (make_plot){
    if (one_dataset){
      print(tSNE_PCA(seurat_obj, plot_name, ...))
    } else {
      print(full_umap(seurat_obj, data_set, plot_name, ...))
    }
  }
  return(seurat_obj)
}

genes_per_group <- function(seurat_obj, gene_set, plot_name, group_by,
                                one_dataset = TRUE, data_set = NULL,
                                make_plot = FALSE, plot_group = NULL) {
  seurat_obj <- plot_gene_set(seurat_obj = seurat_obj, gene_set = gene_set,
                            plot_name = plot_name, one_dataset = one_dataset,
                            data_set = data_set, make_plot = make_plot)
  mean_all <- aggregate(seurat_obj@meta.data[[plot_name]],
                        list(seurat_obj@meta.data[[group_by]]),
                        mean)
  print(mean_all)
  names(mean_all) <- c("group", "average_expresion")
  if (!is.null(plot_group)) {
    mean_all <- mean_all[mean_all$group %in% plot_group, ]
  }
  print(mean_all)
  p <- ggplot2::ggplot(mean_all, ggplot2::aes(x = group, y = average_expresion,
                                              group = 1)) +
          ggplot2::geom_line() +
          #ggplot2::ylim(0, 1) +
          #ggplot2::theme_classic() +
          ggplot2::ggtitle(plot_name)
  print(p)
}

multiple_umap <- function(mtec, sample_list, col_by = "stage") {
  umap_list <- c()
  if (col_by == "stage") {
    new_umap <- lapply(sample_list, function(x) full_stage_umap(mtec, x))
  } else if (col_by %in% rownames(mtec@data) |
             col_by %in% colnames(mtec@meta.data)) {
    new_umap <- lapply(sample_list, function(x) full_gene_umap(mtec, x, col_by))
  }
  for (i in stage_list) {
    if (col_by == "stage") {
      new_umap <- full_stage_umap(mtec, i)
    } else if (col_by %in% rownames(mtec@data) |
               col_by %in% colnames(mtec@meta.data)) {
      new_umap <- full_gene_umap(mtec, i, col_by)
    }
    
    umap_list[[i]] <- new_umap
  }
  nplots <- length(sample_list)
  if (nplots > 6) {
    print("Warning in multiple_umap:")
    print(paste0("Works best if number of plots is less than 6. You have ",
          nplots, " total plots"))
  }
  rows <- ceiling(nplots / 2)
  cols <- 2
  gridExtra::grid.arrange(grobs = umap_list, nrow = rows, ncol = cols)
}

full_umap <- function(mtec, data_set, col_by, plot_type = "umap",
                      dims_use = NULL, meta_data_col = "exp", ...) {
  # Determine where in Seurat object to find variable to color by
  print(col_by)
  if (col_by %in% rownames(mtec@data)){
    col_by_data <- as.data.frame(mtec@data[col_by, ])
  }else if (col_by %in% colnames(mtec@meta.data)){
    col_by_data <- as.data.frame(mtec@meta.data[, col_by, drop = FALSE])
  }else if (col_by == "cluster" | col_by == "Cluster"){
    col_by_data <- as.data.frame(mtec@ident)
  }else {
    stop("col_by must be a gene, metric from meta data or 'cluster'")
  }

  # Make the name in the data frame the same regardless of what it was originally
  names(col_by_data) <- "colour_metric"
  
  col_by_data$all <- col_by_data$colour_metric
  if (is.null(dims_use)){
    dims_use <- c(1,2)
  }
  if (!identical(rownames(mtec@meta.data), rownames(col_by_data))) {
    print("must reorder cells")
    col_by_data <- col_by_data[match(rownames(mtec@meta.data),
                                     rownames(col_by_data)), , drop = FALSE]
  }
  col_by_data[[meta_data_col]] <- mtec@meta.data[[meta_data_col]]
  if (is.factor(col_by_data$all)){
    col_by_data$all <- factor(col_by_data$all,
      levels = c("all_samples", levels(col_by_data$all)))
  }
  col_by_data$all[!(col_by_data[[meta_data_col]] %in% data_set)] <- "all_samples"
  if (plot_type %in% names(mtec@dr)){
    plot_coord <- mtec@dr[[plot_type]]@cell.embeddings
    plot_names <- colnames(plot_coord)
    ndims <- length(plot_names)
    plot_cols <- lapply(dims_use, function(x){
      if (x > ndims) {
        stop("dims_use must be equal to or less than number of dimensions")
      } else {
        plot_col <- plot_names[x]
        return(plot_col)
      }
    })
    plot_cols <- unlist(plot_cols)
    plot_coord <- plot_coord[colnames(plot_coord) %in% plot_cols, ]
    axis_names <- colnames(plot_coord)
    colnames(plot_coord) <- c("dim1", "dim2")
    plot_df <- merge(plot_coord, col_by_data, by = "row.names")

  } else {
    stop("plot type must be a dimensional reduction in dr slot")
  }

   # Plot as discrete
  if (!is.numeric(col_by_data$colour_metric)){
    return_plot <- full_discrete_plots(data_set, plot_df, axis_names = axis_names,
      col_by = col_by, ...)
  # Plot as continuous
  }else{
    return_plot <- full_continuous_plots(data_set, plot_df, col_by = col_by, ...)
  }
  return(return_plot)

}

full_discrete_plots <- function(data_set, plot_df, col_by, axis_names = c("dim1", "dim2"),
                                color = NULL, save_plot = NULL, show_legend = TRUE) {
  if (!(is.null(save_plot))){
    extension <- substr(save_plot, nchar(save_plot)-2, nchar(save_plot))
    if (extension == "pdf"){
      pdf(save_plot)
    } else if (extension == "png") {
      png(save_plot)
    } else {
      print("save plot must be .png or .pdf")
    }
  }
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]
  
  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes_(~dim1,
                                                             ~dim2))
  
  base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                         ggplot2::aes_(~dim1, ~dim2), 
                                         color = "#DCDCDC",
                                         size = 1.5,
                                         show.legend = FALSE)
  base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                         ggplot2::aes_(~dim1, ~dim2,
                                                       color = ~all),
                                         size = 1.5,
                                         show.legend = show_legend)
  
  base_plot <- base_plot + #ggplot2::theme_classic() + 
    ggplot2::ggtitle(paste(data_set, collapse = "_")) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  if (is.null(color)) {
    nColors <- length(levels(factor(plot2$all)))
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
   } else {
    base_plot <- base_plot +
      ggplot2::scale_color_manual(values = color, name = col_by)
   }

  if (!(is.null(save_plot))){
    print(base_plot)
    dev.off()
  }
  return(base_plot)
}

full_continuous_plots <- function(data_set, plot_df, col_by, color = NULL,
                                  limits = NULL, axis_names = c("dim1", "dim2"),
                                  save_plot = NULL, show_legend = TRUE) {
  if (!(is.null(save_plot))){
    extension <- substr(save_plot, nchar(save_plot)-2, nchar(save_plot))
    if (extension == "pdf"){
      pdf(save_plot)
    } else if (extension == "png") {
      png(save_plot)
    } else {
      print("save plot must be .png or .pdf")
    }
  }
  plot_name_comb <- paste(data_set, collapse = "_")
  if (is.null(color)) {
    low <- "#00AFBB"
    high <- "#FC4E07"
  }
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]

  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes_(~dim1, ~dim2))
  
  base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                         ggplot2::aes_(~dim1, ~dim2), 
                                         color = "#DCDCDC",
                                         size = 1.5,
                                         show.legend = FALSE)
  base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                         ggplot2::aes_(~dim1, ~dim2,
                                                       color = ~colour_metric),
                                         size = 1.5,
                                         show.legend = show_legend)
  
  base_plot <- base_plot + #ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(plot_name_comb, " ", col_by)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if(is.null(limits)){
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low, high = high, 
                                                 name = col_by)
  } else {
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low, high = high, 
                                                 name = col_by, limits = limits)
  }

  if (!(is.null(save_plot))){
    print(base_plot)
    dev.off()
  }
  return(base_plot)
  
}

highlight_one_group <- function(seurat_object, meta_data_col, group, color_df = NULL,
                                ...){
  seurat_object@meta.data$highlight_group <- "other_cells"
  seurat_object@meta.data$highlight_group[
    seurat_object@meta.data[[meta_data_col]] == group] <- group
  if (!(is.null(color_df))){
    color_df <- color_df[group]
    color_df <- c(color_df, other_cells = "#DCDCDC")
  } else {
    color_df <- c(group = "#FF0000", other_cells = "#DCDCDC")
  }
  tSNE_PCA(seurat_object, "highlight_group", color = color_df, ...)
}

plot_avg_exp_genes <- function(seurat_object, gene_list, save_plot = NULL, ...){
  if (!(is.null(save_plot))){
    extension <- substr(save_plot, nchar(save_plot)-2, nchar(save_plot))
    if (extension == "pdf"){
      pdf(save_plot)
    } else if (extension == "png") {
      png(save_plot)
    } else {
      print("save plot must be .png or .pdf")
    }
  }
  avg_expression <- get_avg_exp(seurat_object, ...)
  avg_expression <- avg_expression[rownames(avg_expression) %in% average_gene_list, ]
  avg_expression$gene <- rownames(avg_expression)
  avg_expression_melt <- reshape2::melt(avg_expression)
  base_plot <- ggplot2::ggplot(avg_expression_melt, ggplot2::aes(x = variable, y = value,
                                                    group = gene)) +
    ggplot2::geom_line(ggplot2::aes(color = gene)) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::xlab("Experiment") +
    ggplot2::ylab("Average expression") 
    #ggplot2::theme_classic()

  if (!(is.null(save_plot))){
    print(base_plot)
    dev.off()
  }
}



percent_cycling_cells <- function(seurat_object, data_set, meta_data_col){
  cells_use <- rownames(seurat_object@meta.data)[
    seurat_object@meta.data[[meta_data_col]] == data_set]
  new_seurat <- Seurat::SubsetData(seurat_object, cells.use = cells_use)
  cycling_cells <- table(new_seurat@meta.data$cycle_phase)
  if (!("S" %in% cycling_cells)) {
    cycling_cells["S"] = 0 
  }
  cycling_percent <- (cycling_cells["G2M"] +
    cycling_cells["S"])/nrow(new_seurat@meta.data)
  names(cycling_percent) <- data_set
  return(cycling_percent)
}


# For this function, subset_by is best as batches
get_umi <- function(seurat_obj, subset_seurat = FALSE, subset_by = "exp",
  subset_val = "isoControlBeg"){
  if (subset_seurat){
    if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
      seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = subset_by)
      seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = subset_val,
        subset.raw = TRUE)
    } else {
      stop("ident and meta.data slots not in the same order")
    }
  }
  cell_matrix <- as.matrix(seurat_obj@raw.data)
  umi <- median(colSums(cell_matrix))
  return(umi)
}

percents_and_counts <- function(seurat_obj, gene_lists, downsample_UMI = FALSE,
  one_batch = FALSE, batch = "exp", batch_name = "all_cells",
  lowest_UMI = NULL, count = "genes"){
  # If not looking at all sample in a seurat object, than subset to the desired
  # batch
  if (one_batch){
    if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
      seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = batch)
      seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = batch_name,
        subset.raw = TRUE)
    } else {
      stop("ident and meta.data slots not in the same order")
    }
  }

  # Grab the raw.data slot from the seurat object
  cell_matrix <- as.matrix(seurat_obj@raw.data)

  # Downsample the UMI 
  if (downsample_UMI){
    # DropletUtils is required for this function
    if (!requireNamespace("DropletUtils", quietly = TRUE)){
      stop("Package \"DropletUtils\" needed for this function to work. Please install it.",
        call. = FALSE)
    }

    # Determine the median UMI for the dataset and determine factor based on the
    # given lowest UMI
    data_umi <- median(colSums(cell_matrix))
    if (is.null(lowest_UMI)) {
      stop("If downsampling, you must provide a value.
        This value will be the median number of UMI after downsampling")
    }
    factor <- lowest_UMI/data_umi
    
    # Use DropletUtils to downsample the raw matrix
    set.seed(0)
    cell_matrix <- DropletUtils::downsampleMatrix(cell_matrix, prop = factor)
  }
  
  return_list <- list()
  if ("genes" %in% count) {
    # Determine the number of genes in each gene set present in each
    # cell in the data set
    count_list <- lapply(names(gene_lists), function(x)
      gene_count_function(cell_matrix, gene_lists[[x]], x))
    count_df <- do.call(cbind, count_list)
    count_df$exp <- batch_name
    return_list$counts <- count_df
  }
  if ("umi" %in% count) {
    # Determine the number of UMIs in each gene set present in each
    # cell in the data set
    umi_list <- lapply(names(gene_lists), function(x)
      umi_count_function(cell_matrix, gene_lists[[x]], x))
    umi_df <- do.call(cbind, umi_list)
    umi_df$exp <- batch_name
    return_list$umi <- umi_df
  }
  if ("percent" %in% count) {
    # Determine the percent of genes in each set expressed in ANY cell in the data set
    # ie percent of genes seen in at least one cell.
    gene_percent_list <- sapply(names(gene_lists), function(x)
      percent_list(cell_matrix, gene_lists[[x]], x))
    return_list$percents <- gene_percent_list
  }

  # Return both and name based on the batch 
  return_list <- list(return_list)
  #names(return_list) <- batch_name
  return(return_list)
}

# percents_and_counts <- function(seurat_obj, downsample_UMI = FALSE,
#   one_batch = FALSE, batch = "exp", batch_name = "isoControlBeg"){
#   if (one_batch){
#     if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
#       seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = batch)
#       seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = batch_name,
#         subset.raw = TRUE)
#     } else {
#       stop("ident and meta.data slots not in the same order")
#     }
#   }
#   cell_matrix <- as.matrix(seurat_obj@raw.data)
#   if (downsample_UMI){
#     data_umi <- median(colSums(cell_matrix))

#     factor <- lowest_UMI/data_umi
    
#     set.seed(0)
#     cell_matrix <- DropletUtils::downsampleMatrix(cell_matrix, prop = factor)
#   }
#   # if (one_population) {
#   #   if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
#   #     seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = population)
#   #     seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = population_name)
#   #   } else {
#   #     stop("ident and meta.data slots not in the same order")
#   #   }
#   #   cell_matrix <- cell_matrix[ , colnames(cell_matrix) %in% colnames(seurat_obj@data)]
#   # }
#   gene_count_list <- lapply(names(gene_lists), function(x)
#     gene_count_function(cell_matrix, gene_lists[[x]], x))
#   gene_count_df <- do.call(cbind, gene_count_list)
#   gene_count_df$exp <- batch_name

#   gene_percent_list <- sapply(names(gene_lists), function(x)
#     percent_list(cell_matrix, gene_lists[[x]], x))
#   return_list <- list(list(counts = gene_count_df, percents = gene_percent_list))
#   names(return_list) <- batch_name
#   return(return_list)
# }

# # percents_and_counts <- function(seurat_obj, downsample_UMI = FALSE,
# #   one_batch = FALSE, batch = "exp", batch_name = "isoControlBeg",
# #   one_population = TRUE, population = "stage", population_name = "Aire_positive"){
# #   if (one_batch){
# #     if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
# #       seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = batch)
# #       seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = batch_name,
# #         subset.raw = TRUE)
# #     } else {
# #       stop("ident and meta.data slots not in the same order")
# #     }
# #   }
# #   cell_matrix <- as.matrix(seurat_obj@raw.data)
# #   if (downsample_UMI){
# #     data_umi <- median(colSums(cell_matrix))

# #     factor <- lowest_UMI/data_umi
    
# #     set.seed(0)
# #     cell_matrix <- DropletUtils::downsampleMatrix(cell_matrix, prop = factor)
# #   }
# #   if (one_population) {
# #     if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
# #       seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = population)
# #       seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = population_name)
# #     } else {
# #       stop("ident and meta.data slots not in the same order")
# #     }
# #     cell_matrix <- cell_matrix[ , colnames(cell_matrix) %in% colnames(seurat_obj@data)]
# #   }
# #   gene_count_list <- lapply(names(gene_lists), function(x)
# #     gene_count_function(cell_matrix, gene_lists[[x]], x))
# #   gene_count_df <- do.call(cbind, gene_count_list)
# #   gene_count_df$exp <- batch_name

# #   gene_percent_list <- sapply(names(gene_lists), function(x)
# #     percent_list(cell_matrix, gene_lists[[x]], x))
# #   return_list <- list(list(counts = gene_count_df, percents = gene_percent_list))
# #   names(return_list) <- batch_name
# #   return(return_list)
# # }

gene_count_function <- function(cell_matrix, gene_list, list_name){
  gene_matrix <- cell_matrix[rownames(cell_matrix) %in% gene_list, ]
  gene_count <- apply(gene_matrix, 2, function(x) sum(x > 0))
  gene_count <- data.frame(gene_count)
  names(gene_count) <- list_name
  return(gene_count)
}

umi_count_function <- function(cell_matrix, gene_list, list_name){
  gene_matrix <- cell_matrix[rownames(cell_matrix) %in% gene_list, ]
  umi_count <- apply(gene_matrix, 2, function(x) sum(x))
  umi_count <- data.frame(umi_count)
  names(umi_count) <- list_name
  return(umi_count)
}

percent_list <- function(cell_matrix, gene_list, gene_list_name){
  # Subset cell matrix to only be genes of interest
  gene_matrix <- cell_matrix[rownames(cell_matrix) %in% gene_list, ]
  # Counts number of cells expressing each gene.
  cell_count <- apply(gene_matrix, 1, function(x) sum(x > 0))
  # Counts number of genes expressed by at least 1 cell
  expr_genes <- cell_count[cell_count > 0]
  # Number of expressed genes
  n_expr_genes <- length(expr_genes)
  # percent of expressed genes
  percent <- n_expr_genes/length(gene_list)
  #names(percent) <- gene_list_name
  
  # Returns either the number or the updated list
  return(percent)
}


# This is a rather silly workaround. I'm sure there is a better way.
# get_perc_count <- function(percent_counts_list, list_slot, percent = FALSE,
#   count = FALSE){
#   percent_counts_one <- percent_counts_list[[list_slot]]
#   if (percent & count) {
#     stop("can only return percent or count")
#   } else if (!percent & !count){
#     stop("must bet either percent or count")
#   } else if (percent) {
#     return_val <- percent_counts_one$percents
#   } else if (count) {
#     return_val <- percent_counts_one$counts
#   }
#   return(return_val)
# }

get_perc_count <- function(percent_counts_list, list_slot, data_type = "counts"){
  percent_counts_one <- percent_counts_list[[list_slot]]
  if (data_type == "percent") {
    return_val <- percent_counts_one$percents
  } else if (data_type == "counts") {
    return_val <- percent_counts_one$counts
  } else if (data_type == "umi") {
    return_val <- percent_counts_one$umi
  } else {
    stop("must bet either percent, counts, or umi")
  }
  return(return_val)
}

##########################################################################
# New

plot_sling_pseudotime <- function(seurat_object, sling_object, y_val, col_by,
                                  pseudotime_curve, color = NULL,
                                  save_plot = NULL, range = NULL,
                                  plot_type = "dot_plot",
                                  height = 3, width = 7) {
  print(save_plot)
  pseudotime <- data.frame(slingshot::slingPseudotime(sling_object))
  pseudotime_df <- data.frame(pseudotime = pseudotime[[pseudotime_curve]],
    row.names = rownames(pseudotime))
  pseudotime_df <- pseudotime_df[!is.na(pseudotime_df$pseudotime), , drop = FALSE]
  plot_data <- make_plot_df(seurat_object = seurat_object, y_val = y_val,
                            x_val = col_by, col_by = col_by)
  plot_data <- merge(pseudotime_df, plot_data, by = "row.names", all = FALSE)
  if (is.null(color)){
    nColors <- length(levels(factor(plot_df$col_by)))
    color <- RColorBrewer::brewer.pal(nColors, "Set1")
  }
  base_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = pseudotime,
                                                       y = y_value))
  if (plot_type == "dot_plot"){
    base_plot <- base_plot + ggplot2::geom_point(ggplot2::aes(color = col_by)) +
      ggplot2::geom_smooth(se = FALSE, color = "black") +
      ggplot2::scale_color_manual(values = color, name = col_by)
  } else if (plot_type == "density"){
    base_plot <- base_plot + ggridges::geom_density_ridges(ggplot2::aes(fill = col_by)) +
      ggplot2::scale_fill_manual(values = color, name = col_by)
  } else {
    stop("plot_type must be dot_plot or density")
  }
  if (!(is.null(range))) {
    base_plot <- base_plot + ggplot2::xlim(range)
  }

  base_plot <- base_plot + ggplot2::ylab(y_val)
  
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot, height = height, width = width)
  }

  return(base_plot)
}

cumFreqFunc <- function(gene_matrix, gene_list, list_name){
  print(list_name)
  gene_matrix <- gene_matrix[rownames(gene_matrix) %in% gene_list, ]
  cellOrder <- order(colSums(gene_matrix > 0))
  cumFreqGenes <- sapply(seq_along(cellOrder), function(x){
    sum(!rowSums(gene_matrix[, cellOrder[seq_len(x)], drop = FALSE]) == 0)
    })
  cumFreqGenes <- data.frame(list_percent = cumFreqGenes/nrow(gene_matrix))
  cumFreqGenes$ID <- seq.int(nrow(cumFreqGenes))
  cumFreqGenes$gene_list <- list_name
  cumFreqGenes <- list(cumFreqGenes)
  return(cumFreqGenes)
}

get_cell_matrix <- function(seurat_obj, gene_list = NULL, subset = NULL,
                          run_dropout_percent = TRUE, downsample_matrix = FALSE,
                          lowest_UMI = NULL){
  gene_list <- unlist(gene_list)
  if (!(is.null(subset))){
    print(subset)
    cells_use <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$exp ==
                                                 subset]

    seurat_obj <- Seurat::SubsetData(seurat_obj, cells.use = cells_use,
                               subset.raw = TRUE)
  }
  raw_reads <- as.matrix(seurat_obj@raw.data)
  if (downsample_matrix){
    raw_reads <- get_downsampled_matrix(raw_reads, lowest_UMI)
  }
  if (run_dropout_percent){
    percents <- dropout_percent(raw_reads = raw_reads, gene_list = gene_list)
    return(percents)
  } else {
    return(raw_reads)
  }
}

get_downsampled_matrix <- function(raw_matrix, lowest_UMI) {
  data_umi <- median(colSums(as.matrix(raw_matrix)))
  
  factor <- lowest_UMI/data_umi

  set.seed(0)
  cell_matrix <- DropletUtils::downsampleMatrix(raw_matrix, prop = factor)
  cell_matrix <- as.matrix(cell_matrix)
  new_umi <- median(colSums(raw_matrix))
  return(cell_matrix)
}

dropout_percent <- function(raw_reads, gene_list = NULL){
  percents <- rowSums(raw_reads==0)/ncol(raw_reads)*100
  if (is.null(gene_list)){
    return(percents)
  } else {
    percents <- percents[gene_list]
    return(percents)
  }
}
