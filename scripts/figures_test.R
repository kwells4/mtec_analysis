library(mTEC.10x.pipeline)
library(dplyr)

###############
# This contains supplemental functions

###############
# Move all these functions to the package

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
      print(plot_1)
      
    } else {
      plot_1 <- ggplot2::ggplot(data = plot_df, ggplot2::aes_string(name_1,
                                                                    name_2)) +
        ggplot2::geom_point(color = color_scale) + ggplot2::ggtitle(i) +
        #ggplot2::theme_classic() + 
        ggplot2::geom_text(x = text_x, y = text_y,
                                                      label = correlation_plot)
      print(plot_1)
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
    ggplot2::geom_line(ggplot2::aes(linetype = gene))  
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

percents_and_counts <- function(seurat_obj, downsample_UMI = FALSE,
  one_batch = FALSE, batch = "exp", batch_name = "isoControlBeg"){
  if (one_batch){
    if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
      seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = batch)
      seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = batch_name,
        subset.raw = TRUE)
    } else {
      stop("ident and meta.data slots not in the same order")
    }
  }
  cell_matrix <- as.matrix(seurat_obj@raw.data)
  if (downsample_UMI){
    data_umi <- median(colSums(cell_matrix))

    factor <- lowest_UMI/data_umi
    
    set.seed(0)
    cell_matrix <- DropletUtils::downsampleMatrix(cell_matrix, prop = factor)
  }
  # if (one_population) {
  #   if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
  #     seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = population)
  #     seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = population_name)
  #   } else {
  #     stop("ident and meta.data slots not in the same order")
  #   }
  #   cell_matrix <- cell_matrix[ , colnames(cell_matrix) %in% colnames(seurat_obj@data)]
  # }
  gene_count_list <- lapply(names(gene_lists), function(x)
    gene_count_function(cell_matrix, gene_lists[[x]], x))
  gene_count_df <- do.call(cbind, gene_count_list)
  gene_count_df$exp <- batch_name

  gene_percent_list <- sapply(names(gene_lists), function(x)
    percent_list(cell_matrix, gene_lists[[x]], x))
  return_list <- list(list(counts = gene_count_df, percents = gene_percent_list))
  names(return_list) <- batch_name
  return(return_list)
}

# percents_and_counts <- function(seurat_obj, downsample_UMI = FALSE,
#   one_batch = FALSE, batch = "exp", batch_name = "isoControlBeg",
#   one_population = TRUE, population = "stage", population_name = "Aire_positive"){
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
#   if (one_population) {
#     if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
#       seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = population)
#       seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = population_name)
#     } else {
#       stop("ident and meta.data slots not in the same order")
#     }
#     cell_matrix <- cell_matrix[ , colnames(cell_matrix) %in% colnames(seurat_obj@data)]
#   }
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
get_perc_count <- function(percent_counts_list, list_slot, percent = FALSE,
  count = FALSE){
  percent_counts_one <- percent_counts_list[[list_slot]]
  if (percent & count) {
    stop("can only return percent or count")
  } else if (!percent & !count){
    stop("must bet either percent or count")
  } else if (percent) {
    return_val <- percent_counts_one$percents
  } else if (count) {
    return_val <- percent_counts_one$counts
  }
  return(return_val)
}

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

##################################################################################

# Get files from Snakemake
aireTrace <- snakemake@input[[1]]
controls <- snakemake@input[[2]]
allSamples <- snakemake@input[[3]]
controls_slingshot <- snakemake@input[[4]]
allSamples_slingshot <- snakemake@input[[5]]
early_aire_mtec <- snakemake@input[[6]]
sav_dir <- snakemake@output[[1]]
data_directory <- snakemake@params[[1]]
bootstrap_script <- snakemake@params[[2]]

# Colors for plotting
stage_color_df <- data.frame("Cortico_medullary" = "#CC6600", "Ccl21a_high" = "#009933",
                            "Early_Aire" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]

# Load in data
mtec <- get(load(aireTrace))

mtec_wt <- get(load(controls))

mtecCombined <- get(load(allSamples))

wt_slingshot <- get(load(controls_slingshot))

all_slingshot <- get(load(allSamples_slingshot))

progenitor_mtec <- get(load(early_aire_mtec))

TFs_all <- c("H2afz", "Top2a", "Hmgb1", "Hmgn1", "H2afx", as.character(TFs))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

############
# Figure 1 #
############
print("Figure 1")

# Make plots from aire trace

# Figure 1b
tSNE_PCA(mtec, "stage", color = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "figure_1b.pdf"))

# Figure 1c
# Violin plots of marker genes, pick a few more here
trio_plots(mtec, geneset = c("Ackr4", "Ccl21a", "Aire"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_1cI.pdf"))
trio_plots(mtec, geneset = c("Krt10", "Trpm5", "GFP"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_1cII.pdf"))

# Figure 1d
# Pseudotime of Aire and GFP


############
# Figure 2 #
############

print("Figure 2")
load(paste0(data_directory, "gene_lists.rda")


# Figure 2a
# Dot plot of all markers Change this to be most interesting markers
markers_to_plot_full <- c("Krt5", "Ccl21a", "Ascl1", "Hes1", "Hmgb2", "Hmgn2",
  "Hmgb1", "H2afx", "Stmn1", "Tubb5", "Mki67", "Ptma", "Aire", "Utf1", "Fezf2", 
  "Krt10", "Nupr1", "Cebpb", "Trpm5", "Pou2f3", "Dclk1")

pdf(paste0(save_dir, "figure_2aI.pdf"))
dot_plot <- Seurat::DotPlot(mtec, genes.plot = rev(markers_to_plot_full),
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = F, dot.scale = 8, do.return = T)


dev.off()

pdf(paste0(save_dir, "figure_2aII.pdf"))
dot_plot <- Seurat::DotPlot(mtec, genes.plot = rev(markers_to_plot_full),
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = T, dot.scale = 8, do.return = T)
dev.off()

# Figure 2b
# Jitter plots of chromatin modifiers overlayed with cell cycle state
trio_plots(mtec, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_2bI.pdf"), group_color = FALSE)
trio_plots(mtec, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_2bII.pdf"), group_color = FALSE)

# Figure 2c
# Flow of Mki76/Aire

# Figure 2d
# SC velocity

# Figure 2e
# UMAP of cell stage with slingshot overlay
##########
# This is my crummy workaround until I get my new package running. If plotDimRed
# is returning a ggplot object, I'm golden

umap_coords <- data.frame(mtec_wt@dr$umap@cell.embeddings)
umap_coords$stage <- mtec_wt@meta.data$stage
base_plot <- ggplot2::ggplot(umap_coords,
  ggplot2::aes(x = UMAP1, y = UMAP2, color = stage)) +
  ggplot2::geom_point() + 
  ggplot2::scale_color_manual(values = stage_color) +
  ggplot2::theme(legend.position = "none")

##############
# Keep this
c <- slingCurves(wt_slingshot)[[1]]
curve1_coord <- data.frame(c$s[c$ord, dims])
curve1_coord$stage <- "line"

# This line cuts off the long tail... Probably a better option is
# to just remove unknown before running slingshot.
curve1_coord <- curve1_coord[curve1_coord$UMAP2 > -5, ]
base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
  ggplot2::aes(UMAP1, UMAP2), color = "black", size = 1)

ggplot2::ggsave(paste0(save_dir, "figure_2e.pdf"), plot = base_plot)


gene_names <- c("Aire", "Fezf2")
plot_sets <- c("tra_fantom", "aire_genes", "fezf2_genes")

# Figure 2f
# Pseudotime of genes
plot_names <- list(protein_coding = paste0(save_dir, "figure_2fI.pdf"),
  tra_fantom = paste0(save_dir, "figure_2fII.pdf"),
  aire_genes = paste0(save_dir, "figure_2fIV.pdf"),
  Aire = paste0(save_dir, "figure_2fIII.pdf"),
  fezf2_genes = paste0(save_dir, "figure_2fVI.pdf"),
  Fezf2 = paste0(save_dir, "figure_2fV.pdf"))


for (gene_set in names(gene_lists)) {
  mtec_wt <- plot_gene_set(mtec_wt,
    gene_lists[[gene_set]], gene_set, make_plot = FALSE)
}


# Remove unknown cells
mtec_wt_plot <- mtec_wt
mtec_wt_plot@assay$DE <- NULL
mtec_wt_plot <- Seurat::SubsetData(mtec_wt_plot, ident.remove = "unknown")

# Plot each of the genes and gene sets in pseudotime, end at 16 because the
# is where the "unknown" cells are
plot_list <- lapply(names(plot_names), function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_plot, sling_object = wt_slingshot, y_val = x,
  col_by = "stage", pseudotime_curve = "curve1", color = stage_color,
  range = c(0, 16), save_plot = plot_names[[x]]))


# ############
# # Figure 3 #
# ############

print("Figure 3")

new_exp_names <- c(aireTrace = "Trace exp",
                  isoControlBeg = "Ctl wk 2",
                  isoControlEnd = "Ctl wk 10",
                  timepoint1 = "wk 2",
                  timepoint2 = "wk 4",
                  timepoint3 = "wk 6",
                  timepoint5 = "wk 10")

mtecCombined@meta.data$pub_exp <- new_exp_names[mtecCombined@meta.data$exp]
mtecCombined@meta.data$pub_exp <- factor(mtecCombined@meta.data$pub_exp,
                                         levels = unname(new_exp_names))

timecourse_color <- RColorBrewer::brewer.pal(8, "Set1")
timecourse_color <- c(timecourse_color[2:5], timecourse_color[7:8])

data_sets <- unique(mtecCombined@meta.data$pub_exp)
data_sets <- data_sets[data_sets != new_exp_names['aireTrace']]

# stage_color_df_2 <- data.frame("Cortico_medullary" = "#CC6600",
#                                "Ccl21a_high" = "#009933",
#                                "Early_Aire" = "#0066CC",
#                                "Aire_positive" = "#660099",
#                                "Late_Aire" = "#FF0000",
#                                "Tuft" = "#990000",
#                                "unknown" = "#D3D3D3")

# stage_color2 <- t(stage_color_df_2)[ , 1]

stage_color_df_3 <- data.frame("CorticoMedullary" = "#CC6600",
                               "Ccl21aHigh" = "#009933",
                               "EarlyAire" = "#0066CC",
                               "AirePositive" = "#660099",
                               "LateAire" = "#FF0000",
                               "Tuft" = "#990000",
                               "unknown" = "#D3D3D3")

stage_color3 <- t(stage_color_df_3)[ , 1]

# Figure 2B
# UMAP of all cells Either put the key on both or remove the key from both
tSNE_PCA(mtecCombined, "stage", save_plot = paste0(save_dir, "figure_3bI.pdf"),
  color = stage_color, show_legend = FALSE)

# UMAP highlighting aire_trace cells colored by aire_trace labels
full_umap(mtecCombined, "aireTrace", col_by = "at_stage", color = stage_color,
  save_plot = paste0(save_dir, "figure_3bII.pdf"), show_legend = FALSE)

# Figure 3C
# Barplots of recovery
mtecCombSub <- mtecCombined
mtecCombSub@assay$ablation_DE <- NULL

stage_list_all <- lapply(data_sets, function(x) populations_dfs_new(mtecCombSub,
                         x, subsample = TRUE, subsample_by = "pub_exp"))
stage_df_all <- do.call("rbind", stage_list_all)

stage_df_all$sample <- factor(stage_df_all$sample, levels = unname(new_exp_names))

population_plots(stage_df_all, color = stage_color,
  save_plot = paste0(save_dir, "figure_3c.pdf"))

# Figure 3D
# Umap with slingshot overlay
##########
# This is my crummy workaround until I get my new package running. If plotDimRed
# is returning a ggplot object, I'm golden

umap_coords <- data.frame(mtecCombined@dr$umap@cell.embeddings)
umap_coords$stage <- mtecCombined@meta.data$stage
base_plot <- ggplot2::ggplot(umap_coords,
  ggplot2::aes(x = UMAP1, y = UMAP2, color = stage)) +
  ggplot2::geom_point() + 
  ggplot2::scale_color_manual(values = stage_color_h) +
  ggplot2::theme(legend.position = "none")

##############
# Keep this
c <- slingCurves(allSamples_slingshot)[[3]]
curve1_coord <- data.frame(c$s[c$ord, dims])
curve1_coord$stage <- "line"

# This line cuts off the long tail... Probably a better option is
# to just remove unknown before running slingshot.
curve1_coord <- curve1_coord[curve1_coord$UMAP2 > -5, ]
base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
  ggplot2::aes(UMAP1, UMAP2), color = "black", size = 1)

ggplot2::ggsave(paste0(save_dir, "figure_3d.pdf"), plot = base_plot)


# Figure 3E
# Timecourse in pseudotime
mtec_no_at <- mtecCombSub
mtec_no_at <- Seurat::SetAllIdent(mtec_no_at, id = "exp")
mtec_no_at <- Seurat::SubsetData(mtec_no_at, ident.remove = "aireTrace")
plot_sling_pseudotime(
  seurat_object = mtec_no_at, sling_object = allSamples_slingshot, y_val = "pub_exp",
  col_by = "pub_exp", pseudotime_curve = "curve3", color = timecourse_color,
  plot_type = "density", width = 7, height = 7,
  save_plot = paste0(save_dir, "figure_3e.pdf"))


# ###############################################################################

# ############
# # Figure 4 #
# ############
print("Figure 4")


lowest_UMI_exp <- "timepoint2"

downsample_UMI <- TRUE

if (downsample_UMI) {
  
  lowest_UMI <- get_umi(mtecCombSub, subset_seurat = TRUE, subset_by = "stage_exp",
    subset_val = paste0("AirePositive_", lowest_UMI_exp))
}

# Figure 4a
# Start with full tSNE only coloring the Aire positive cluster
highlight_one_group(mtecCombined, meta_data_col = "stage", group = "Aire_positive",
  color_df = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "figure_4a.pdf"))

# Figure 4b
average_gene_list <- c("Aire", "Fezf2", "Gapdh", "Emc7")
mtec_aire_positive <- Seurat::SubsetData(mtecCombSub, ident.use = "AirePositive")
cells_use <- rownames(mtec_aire_positive@meta.data)[mtec_aire_positive@meta.data$exp !=
                                                 "aireTrace"]

no_at_mtec_aire <- Seurat::SubsetData(mtec_aire_positive, cells.use = cells_use)

# Make the lines thicker here
plot_avg_exp_genes(no_at_mtec_aire, average_gene_list,
                   save_plot = paste0(save_dir, "figure_4b.pdf"),
                   avg_expr_id = "pub_exp")



mtecCombined_all <- no_at_mtec_aire
# This is an okay place for a for loop (recursion) 
# http://adv-r.had.co.nz/Functionals.html
for (gene_set in names(gene_lists)) {
  mtecCombined_all <- plot_gene_set(mtecCombined_all,
    gene_lists[[gene_set]], gene_set, make_plot = FALSE)
}


# Figure 4c
gene_sets <- c("all_other_genes", "tra_fantom", "aire_genes")
trio_plots(mtecCombined_all, geneset = gene_sets,
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "figure_4c.pdf"))

# This could probably be made into it's own thing
# percents_counts_all <- sapply(data_sets, function(x) percents_and_counts(mtecCombSub,
#   batch_name = x, downsample_UMI = TRUE, one_batch = TRUE, batch = "pub_exp", 
#   one_population = TRUE, population = "stage", "population_name" = "Aire_positive"))

old_data_sets <- unique(mtecCombined@meta.data$exp)
old_data_sets <- old_data_sets[old_data_sets != "aireTrace"]

# I may have fixed it, run tomorrow
percents_counts_all <- sapply(old_data_sets, function(x) percents_and_counts(mtecCombSub,
  batch_name = paste0("AirePositive_", x), downsample_UMI = TRUE, one_batch = TRUE,
  batch = "stage_exp"))

# percents_counts_all <- sapply(old_data_sets, function(x) percents_and_counts(mtecCombSub,
#   batch_name = x, downsample_UMI = TRUE, one_batch = TRUE,
#   batch = "exp", one_population = TRUE, population_name = "Aire_positive",
#   population = "stage"))

percents <- sapply(names(percents_counts_all), function(x) 
  get_perc_count(percents_counts_all, x, percent = TRUE), USE.NAMES = TRUE)

counts <- lapply(names(percents_counts_all), function(x)
  get_perc_count(percents_counts_all, x, count = TRUE))

counts_df <- do.call(rbind, counts)
counts_df$exp <- sub("AirePositive_", "", counts_df$exp)
counts_df$pub_exp <- new_exp_names[counts_df$exp]
counts_df$pub_exp <- factor(counts_df$pub_exp,
                            levels = unname(new_exp_names))
counts_df_m <- reshape2::melt(counts_df, variable.name = "gene_list",
  value.name = "gene_count")

to_plot <- c("tra_fantom", "all_other_genes", "aire_genes", "fezf2_genes")
short_list <- c("tra_fantom", "aire_genes", "fezf2_genes")

counts_df_plot <- counts_df_m[counts_df_m$gene_list %in% to_plot, ]
counts_df_short <- counts_df_m[counts_df_m$gene_list %in% short_list, ]

percents_m <- reshape2::melt(percents)
names(percents_m) <- c("gene_list", "exp", "percent_of_genes")
percents_m$exp <- sub("AirePositive_", "", percents_m$exp)
percents_m$pub_exp <- new_exp_names[percents_m$exp]
percents_m$pub_exp <- factor(percents_m$pub_exp,
                             levels = unname(new_exp_names))
percents_plot <- percents_m[percents_m$gene_list %in% to_plot, ]

# Figure 4d
# Percent of gene lists
pdf(paste0(save_dir, "figure_4d.pdf"))
ggplot2::ggplot(percents_plot, ggplot2::aes(x = pub_exp, y = percent_of_genes,
                                         group = gene_list, color = gene_list)) +
  ggplot2::geom_line() +
  ggplot2::ylim(0,1)  
  #ggplot2::theme_classic()

dev.off()

# Figure 4e
# Number of genes per cell
full_plot <- ggplot2::ggplot(counts_df_plot, ggplot2::aes(x = gene_list,
                                                          y = gene_count,
                                                          fill = pub_exp)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = timecourse_color) +
  #ggplot2::theme_classic() +
  ggpubr::stat_compare_means(method = "anova", size = 2, label.y = 6150)

zoom_plot <- ggplot2::ggplot(counts_df_short, ggplot2::aes(x = gene_list,
                                                           y = gene_count,
                                                           fill = pub_exp)) +
  ggplot2::geom_boxplot(show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = timecourse_color) +
  #ggplot2::theme_classic() +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 panel.border = ggplot2::element_rect(color = "black",
                                                      fill = NA,
                                                      size = 1)) +
  ggpubr::stat_compare_means(method = "anova", size = 2, label.y = 300)

zoom_plot_g <- ggplot2::ggplotGrob(zoom_plot)

all_plots <- full_plot + ggplot2::annotation_custom(grob = zoom_plot_g,
                                                    xmin = 1.5,
                                                    xmax = Inf,
                                                    ymin = 1000,
                                                    ymax = Inf) +
  ggplot2::annotation_custom(grob = grid::rectGrob(gp = grid::gpar(fill = NA)),
                             xmin = 1.5,
                             xmax = Inf,
                             ymin = -Inf,
                             ymax = 500)
pdf(paste0(save_dir, "figure_4e.pdf"))
all_plots

dev.off()

# ############
# # Figure 5 #
# ############

print("Figure 5")
reanalysis_colors <- c("#603E95", "#009DA1", "#FAC22B", "#D7255D")

progenitor_mtec@meta.data$pub_exp <- new_exp_names[progenitor_mtec@meta.data$exp]
progenitor_mtec@meta.data$pub_exp <- factor(progenitor_mtec@meta.data$pub_exp,
                                         levels = unname(new_exp_names))
cells_use <- rownames(progenitor_mtec@meta.data)[progenitor_mtec@meta.data$exp !=
                                                 "aireTrace"]

no_at_mtec <- Seurat::SubsetData(progenitor_mtec, cells.use = cells_use)

# Figure 5a
# Highlight just the early aire cells 
highlight_one_group(mtecCombined, meta_data_col = "stage", group = "Early_Aire",
  color_df = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "figure_5a.pdf"))

# Figure 5b
# % of cells in G2
cell_cycle <- mtecCombined@meta.data[rownames(mtecCombined@meta.data) %in% 
  rownames(no_at_mtec@meta.data), ]
if (!identical(rownames(cell_cycle), rownames(no_at_mtec@meta.data))) {
  print("must reorder cells")
  cell_cycle <- cell_cycle[match(rownames(no_at_mtec@meta.data),
                                     rownames(cell_cycle)), , drop = FALSE]
}
no_at_mtec@meta.data$cycle_phase <- cell_cycle$cycle_phase
percent_cycling <- sapply(data_sets, USE.NAMES = TRUE,
  function(x) percent_cycling_cells(no_at_mtec,
  data_set = x, meta_data_col = "pub_exp"))
percent_cycling <- data.frame(percent_cycling)
percent_cycling$experiment <- rownames(percent_cycling)
percent_cycling$experiment <- factor(percent_cycling$experiment,
                                     levels = unname(new_exp_names))
pdf(paste0(save_dir, "figure_5b.pdf"))

ggplot2::ggplot(percent_cycling, ggplot2::aes(x = experiment,
                                              y = percent_cycling)) +
  ggplot2::geom_bar(ggplot2::aes(fill = experiment),
    stat = "identity") + 
  ggplot2::scale_fill_manual(values = timecourse_color)

dev.off()

# Figure 5c
trio_plots(no_at_mtec, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "figure_5cI.pdf"))

trio_plots(no_at_mtec, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "figure_5cII.pdf"))



# Figure 5d
# UMAP of reanalysis of early aire cells
tSNE_PCA(no_at_mtec, "cluster", color = reanalysis_colors,
  save_plot = paste0(save_dir, "figure_5d.pdf"))

# Figure 5e
trio_plots(no_at_mtec, geneset = c("Hmgb2", "Tubb5", "Stmn2"), cell_cycle = FALSE,
  plot_violin = TRUE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  sep_by = "pub_exp", color = timecourse_color,
  save_plot = paste0(save_dir, "figure_5eI.pdf"))

# Figure 5e
trio_plots(no_at_mtec, geneset = c("Aire", "Ccl21a", "Fezf2"), cell_cycle = FALSE,
  plot_violin = TRUE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  sep_by = "pub_exp", color = timecourse_color,
  save_plot = paste0(save_dir, "figure_5eII.pdf"))

# Figure 5f
clusters <- unique(no_at_mtec@meta.data$res.0.6)
percent_cycling <- sapply(clusters, USE.NAMES = TRUE,
  function(x) percent_cycling_cells(no_at_mtec,
  data_set = x, meta_data_col = "res.0.6"))
percent_cycling <- data.frame(percent_cycling)
percent_cycling$cluster <- sub("\\d\\.", "", rownames(percent_cycling))
pdf(paste0(save_dir, "figure_5f.pdf"))

ggplot2::ggplot(percent_cycling, ggplot2::aes(x = cluster,
                                              y = percent_cycling)) +
  ggplot2::geom_bar(ggplot2::aes(fill = cluster),
    stat = "identity") + 
  ggplot2::scale_fill_manual(values = reanalysis_colors)

dev.off()

# Figure 5g
stage_list_all <- lapply(data_sets, function(x) populations_dfs_new(no_at_mtec,
                         x, subsample = TRUE, subsample_by = "pub_exp",
                         meta_data_col = "res.0.6"))
stage_df_all <- do.call("rbind", stage_list_all)
stage_df_all$sample <- factor(stage_df_all$sample, levels = unname(new_exp_names))


population_plots(stage_df_all, color = reanalysis_colors,
  save_plot = paste0(save_dir, "figure_5g.pdf"))







############################################################################

########################
# Supplemental Figures #
########################

#########################
# Supplemental Figure 1 #
#########################

# S1a
# Heatmap of all TFs with interesting TFs highlighted
pdf(paste0(save_dir, "figure_s1a.pdf"))
plot_heatmap(mtec, subset_list = TFs_all,
  color_list = c("Cdx1", "Utf1", "Tcf7", "Spib", "Cdk4", "Ptma",
    "H2afx", "Hmgb1"),
  color_list2 = c("Aire", "Irf7", "Cited2", "Spib", "Hes1", "Pax1",
                "Relb", "Lmo4", "Pou2f3"),
  cell_color = stage_color)
def.off()

# S1b
# Marker genes on UMAP
tSNE_PCA(mtec, "Ackr4", save_plot = paste(save_dir, "figure_s2bI"))
tSNE_PCA(mtec, "Psmb11", save_plot = paste(save_dir, "figure_s2bII"))
tSNE_PCA(mtec, "Ccl21a", save_plot = paste(save_dir, "figure_s2bIII"))
tSNE_PCA(mtec, "Krt5", save_plot = paste(save_dir, "figure_s2bIV"))
tSNE_PCA(mtec, "Krt8", save_plot = paste(save_dir, "figure_s2bV"))
tSNE_PCA(mtec, "Ascl1", save_plot = paste(save_dir, "figure_s2bVI"))
tSNE_PCA(mtec, "Fezf2", save_plot = paste(save_dir, "figure_s2bVII"))
tSNE_PCA(mtec, "Aire", save_plot = paste(save_dir, "figure_s2bVIII"))
tSNE_PCA(mtec, "Tnfrsf11a", save_plot = paste(save_dir, "figure_s2bIX"))
tSNE_PCA(mtec, "Krt10", save_plot = paste(save_dir, "figure_s2bX"))
tSNE_PCA(mtec, "Trpm5", save_plot = paste(save_dir, "figure_s2bXI"))
tSNE_PCA(mtec, "Dclk1", save_plot = paste(save_dir, "figure_s2bXII"))
tSNE_PCA(mtec, "Pou2f3", save_plot = paste(save_dir, "figure_s2bXII"))
tSNE_PCA(mtec, "GFP", save_plot = paste(save_dir, "figure_s2bXIV"))

######################### 
# Supplemental Figure 2 #
#########################

# S2a
# Violin plots of marker genes for WT
trio_plots(mtec_wt, geneset = c("Ackr4", "Ccl21a", "Aire"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_s2aI.pdf"))
trio_plots(mtec_wt, geneset = c("Krt10", "Trpm5", "GFP"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_s2aII.pdf"))

# S2b
# Dotplot of WT
markers_to_plot_full <- c("Krt5", "Ccl21a", "Ascl1", "Hes1", "Hmgb2", "Hmgn2",
  "Hmgb1", "H2afx", "Stmn1", "Tubb5", "Mki67", "Ptma", "Aire", "Utf1", "Fezf2", 
  "Krt10", "Nupr1", "Cebpb", "Trpm5", "Pou2f3", "Dclk1")

pdf(paste0(save_dir, "figure_s2bI.pdf"))
dot_plot <- Seurat::DotPlot(mtec_wt, genes.plot = rev(markers_to_plot_full),
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = F, dot.scale = 8, do.return = T)


dev.off()

pdf(paste0(save_dir, "figure_s2bII.pdf"))
dot_plot <- Seurat::DotPlot(mtec_wt, genes.plot = rev(markers_to_plot_full),
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = T, dot.scale = 8, do.return = T)
dev.off()

# S2c
# Jitter plots of cycling with markers from AT
trio_plots(mtec_wt, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_s2cI.pdf"), group_color = FALSE)
trio_plots(mtec_wt, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_s2cII.pdf"), group_color = FALSE)

# S2d
# Correlation of WT with Aire Trace
pdf(paste0(save_dir, "figure_s2d"))
master_plot(mtec, "aire_trace", mtec_wt, "wt", stage_color_df)
dev.off()

#########################
# Supplemental Figure 3 #
#########################

# S3a
# Umap of stage recovery


# S3b
# Violin markers for all as in S2
trio_plots(mtecCombined, geneset = c("Ackr4", "Ccl21a", "Aire"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_s3bI.pdf"))
trio_plots(mtecCombined, geneset = c("Krt10", "Trpm5", "GFP"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "figure_s3bII.pdf"))

# S3c Correlation of TP5 with WT
mtecCombExp <- Seurat::SetAllIdent(mtecombSub, id = "exp")
mtec_tp5 <- Seurat::SubsetData(mtecCombExp, ident.use = "timepoint5")
mtec_end <- Seurat::SubsetData(mtecCombexp, ident.use = "isoControlEnd")
pdf(paste0(save_dir, "figure_s3c"))
master_plot(mtec_tp5, "wk 10", mtec_end, "ctl wk 10", stage_color_df)
dev.off()


#########################
# Supplemental Figure 4 #
#########################

# S4a
# Marker genes on UMAP
tSNE_PCA(mtecCombined, "Ackr4", save_plot = paste(save_dir, "figure_s4aI"))
tSNE_PCA(mtecCombined, "Psmb11", save_plot = paste(save_dir, "figure_s4aII"))
tSNE_PCA(mtecCombined, "Ccl21a", save_plot = paste(save_dir, "figure_s4aIII"))
tSNE_PCA(mtecCombined, "Krt5", save_plot = paste(save_dir, "figure_s4aIV"))
tSNE_PCA(mtecCombined, "Krt8", save_plot = paste(save_dir, "figure_s4aV"))
tSNE_PCA(mtecCombined, "Ascl1", save_plot = paste(save_dir, "figure_s4aVI"))
tSNE_PCA(mtecCombined, "Fezf2", save_plot = paste(save_dir, "figure_s4aVII"))
tSNE_PCA(mtecCombined, "Aire", save_plot = paste(save_dir, "figure_s4aVIII"))
tSNE_PCA(mtecCombined, "Tnfrsf11a", save_plot = paste(save_dir, "figure_s4aIX"))
tSNE_PCA(mtecCombined, "Krt10", save_plot = paste(save_dir, "figure_s4aX"))
tSNE_PCA(mtecCombined, "Trpm5", save_plot = paste(save_dir, "figure_s4aXI"))
tSNE_PCA(mtecCombined, "Dclk1", save_plot = paste(save_dir, "figure_s4aXII"))
tSNE_PCA(mtecCombined, "Pou2f3", save_plot = paste(save_dir, "figure_s4aXII"))
tSNE_PCA(mtecCombined, "GFP", save_plot = paste(save_dir, "figure_s4aXIV"))

#########################
# Supplemental Figure 5 #
#########################

# S5a 
# TRA recovery UMAPs
# TRAs in recovery
limit_list <- list(tra_fantom = c(0, 0.100),
                   aire_genes = c(0, 0.100),
                   fezf2_genes = c(0, 0.300))

plot_names_fig5 <- list(isoControlBeg = paste0(save_dir, "figure_s5aI.pdf"),
  isoControlEnd = paste0(save_dir, "figure_s5aII.pdf"),
  timepoint1 = paste0(save_dir, "figure_s5aIII.pdf"),
  timepoint2 = paste0(save_dir, "figure_s5aIV.pdf"),
  timepoint3 = paste0(save_dir, "figure_s5aV.pdf"),
  timepoint5 = paste0(save_dir, "figure_s5aVI.pdf"))

names(plot_names_fig5) <- new_exp_names[names(plot_names_fig5)]

lapply(names(plot_names_fig4), function(x) plot_gene_set(mtecCombined,
                                            gene_set = gene_lists[["tra_fantom"]],
                                            plot_name = "tra_fantom",
                                            one_dataset = FALSE,
                                            data_set = x,
                                            meta_data_col = "pub_exp",
                                            limits = limit_list[["tra_fantom"]],
                                            save_plot = plot_names_fig5[[x]]))

# S5b
# Bootstrap downsample plots
if (bootstrap) {
  source(bootstrap_script)

}

# S5c
# Number of protein coding genes seen in WT not downsampled

#########################
# Supplemental Figure 6 #
#########################
# S6a
# Ngene and nUMI before correction


# S6b
# nGene and nUMI after correction

# S6c
# Dropouts of different house keeping genes before and after





##################################################################################
# stage_to_plot <- c("CorticoMedullary", "Ccl21aHigh", "EarlyAire", "AirePositive",
#                    "LateAire", "Tuft")

# markers_to_plot <- c("Krt5", "Ccl21a", "Ascl1", "Hes1", "Hmgb2", "Hmgn2",
#                     "Hmgb1", "H2afx", "Ptma", "Aire", "Utf1",
#                     "Fezf2", "Krt10", "Trpm5", "Pou2f3", "Dclk1")

# sdp <- Seurat::SplitDotPlotGG(mtecCombined, genes.plot = rev(markers_to_plot),
#                               cols.use = RColorBrewer::brewer.pal(7, "Set1"),
#                               x.lab.rot = T, plot.legend = T, dot.scale = 8,
#                               do.return = T, grouping.var = "exp")

# for (i in stage_to_plot) {
#   mtecSubset <- Seurat::SubsetData(mtecCombined, ident.use = i)
#   sdp <- Seurat::SplitDotPlotGG(mtecSubset, genes.plot = rev(markers_to_plot),
#                                 cols.use = RColorBrewer::brewer.pal(7, "Set1"),
#                                 x.lab.rot = T, plot.legend = T, dot.scale = 8,
#                                 do.return = T, grouping.var = "exp")
# }

# # TRA in recovery


# my_pal <- RColorBrewer::brewer.pal(8, "Set1")
# my_pal <- c(my_pal[2:5], my_pal[7:8])
# avg_plots <- c("ranked", "expression")
# ranked_df_all <- NULL
# exp_df_all <- NULL
# mtecCombined@meta.data$exp <- factor(mtecCombined@meta.data$exp)

# # Some genes in recovery
# for (i in levels(mtecCombined@meta.data$exp)) {
#   print(i)
#   if (i != "aireTrace") {
#     cells_use <- rownames(mtecCombined@meta.data)[mtecCombined@meta.data$exp == i]
#     # Make a seurat object with only the data set cells. This object
#     # will maintain the tSNE locations from the combined analysis.
#     combSubset <- Seurat::SubsetData(mtecCombined, cells.use = cells_use)
#     if ("ranked" %in% avg_plots) {
#       ranked_df_one <- gene_exp_df(combSubset, i, df_contents = "ranked")
#       ranked_df_all <- rbind(ranked_df_all, ranked_df_one)
#     }
  
#     if ("expression" %in% avg_plots) {
#       exp_df_one <- gene_exp_df(combSubset, i, df_contents = "expression")
#       exp_df_all <- rbind(exp_df_all, exp_df_one)
#     }
#     new_umap <- full_stage_umap(mtecCombined, i)
#     print(new_umap)
#   }
# }


# gene_list <- c("Calcb", "Csn2", "Ubd", "H2afx", "Hes1", "Pou2f3", "Aire",
#                "Utf1", "Srgn", "Fezf2", "Anxa2", "Cyba", "Cxcl10",
#                "Lgals1", "Fgf21", "Plb1", "Skint10", "Rgs5", "Hagh", "Tmsb10",
#                "Fabp5", "Spib", "Hmgb2", "Hmgb1", "Hmgn2", "Ascl1", "Cldn7",
#                "Grap", "Cystm1", "Nos2", "Nupr1", "Cebpb", "Trpm5", "Dclk1",
#                "Cdkn1a", "Selm", "Cyp2a5", "Rptn", "Ahcyl2")

# for (gene in gene_list) {
#   if ("ranked" %in% avg_plots) {
#     ranked_df_all$stage <- factor(ranked_df_all$stage,
#                                   levels = stage_to_plot)
#     plot_gene_exp(ranked_df_all, gene, col = stage_color3)
#   }
#   if ("expression" %in% avg_plots) {
#     exp_df_all$stage <- factor(exp_df_all$stage,
#                                levels = stage_to_plot)
#     plot_gene_exp(exp_df_all, gene, low_lim = 0, high_lim = 4, col = stage_color3)
#   }
# }

# # Maybe recovery algorithm here...

# ############
# # Figure 6 #
# ############

# print("Figure 6")


# # Density plot of Ccl21a here
# meta_data <- mtecCombined@meta.data
# mtec_data <- t(as.matrix(mtecCombined@data))
# ccl21a_df <- as.data.frame(mtec_data[ , "Ccl21a"])
# aire_df <- as.data.frame(mtec_data[ , "Aire"])
# hmgb2_df <- as.data.frame(mtec_data[ , "Hmgb2"])
# names(ccl21a_df) <- "Ccl21a"
# names(aire_df) <- "Aire"
# names(hmgb2_df) <- "Hmgb2"
# meta_data <- merge(meta_data, ccl21a_df, by = "row.names")
# rownames(meta_data) <- meta_data$Row.names 
# meta_data$Row.names <- NULL

# meta_data <- merge(meta_data, aire_df, by = "row.names")
# rownames(meta_data) <- meta_data$Row.names 
# meta_data$Row.names <- NULL

# meta_data <- merge(meta_data, hmgb2_df, by = "row.names")
# rownames(meta_data) <- meta_data$Row.names 
# meta_data$Row.names <- NULL

# meta_data <- meta_data[meta_data$exp != "aireTrace", ]
# meta_data$exp <- factor(meta_data$exp)

# # Finish making the bar plot!
# meta_data$Ccl21a_high_expr <- meta_data$Ccl21a > 3 & meta_data$stage == "Early_Aire"
# meta_data$True_aire <- meta_data$stage == "Aire_positive"
# meta_data$True_progenitor <- meta_data$stage == "Early_Aire"

# percents_df <- data.frame(Timepoint = character(),
#                           Ccl21a_whole = double(),
#                           Progenitor_whole = double(),
#                           Aire_whole = double(),
#                           post_aire_whole = double())

# for (i in levels(meta_data$exp)) {
#   new_meta_data <- meta_data[meta_data$exp == i, ]
#   df_to_add <- data.frame(Timepoint = i,
#                           Ccl21a_whole = sum(new_meta_data$Ccl21a_high_expr) /
#                                          nrow(new_meta_data) * 100,
#                           Progenitor_whole = sum(new_meta_data$True_progenitor) /
#                                              nrow(new_meta_data) * 100,
#                           Aire_whole = sum(new_meta_data$True_aire) /
#                                            nrow(new_meta_data) * 100,
#                           Post_aire_whole = sum(new_meta_data$stage ==
#                                                 "Late_Aire") /
#                                             nrow(new_meta_data) * 100)
#   percents_df <- rbind(percents_df, df_to_add)
# }

# meta_data <- meta_data[meta_data$stage == "Early_Aire", ]



# p_density <- ggplot2::ggplot(meta_data, ggplot2::aes(Ccl21a, 
#       color = exp)) + ggplot2::geom_density() + 
#       ggplot2::scale_color_manual(values = my_pal, name = "timepoint") +
#       ggplot2::theme_classic() + ggplot2::xlab("Ccl21a")

# print(p_density)

# p_barplot_ccl21a <- ggplot2::ggplot(percents_df,
#                                     ggplot2::aes(Timepoint, Ccl21a_whole,
#                                                  fill = Timepoint)) +
#               ggplot2::geom_bar(stat = "identity") + 
#               ggplot2::scale_fill_manual(values = my_pal, name = "timepoint") +
#               ggplot2::theme_classic() + ggplot2::ylab("Ccl21a percent of whole")

# print(p_barplot_ccl21a)

# p_barplot_progenitor <- ggplot2::ggplot(percents_df,
#                                         ggplot2::aes(Timepoint, Progenitor_whole,
#                                                     fill = Timepoint)) +
#               ggplot2::geom_bar(stat = "identity") + 
#               ggplot2::scale_fill_manual(values = my_pal, name = "timepoint") +
#               ggplot2::theme_classic() + ggplot2::ylab("Early Aire percent of whole")

# print(p_barplot_progenitor)

# p_barplot_aire <- ggplot2::ggplot(percents_df,
#                                         ggplot2::aes(Timepoint, Aire_whole,
#                                                     fill = Timepoint)) +
#               ggplot2::geom_bar(stat = "identity") + 
#               ggplot2::scale_fill_manual(values = my_pal, name = "timepoint") +
#               ggplot2::theme_classic() + ggplot2::ylab("Aire Positive percent of whole")

# print(p_barplot_aire)

# p_barplot_post_aire <- ggplot2::ggplot(percents_df,
#                                         ggplot2::aes(Timepoint, Post_aire_whole,
#                                                     fill = Timepoint)) +
#               ggplot2::geom_bar(stat = "identity") + 
#               ggplot2::scale_fill_manual(values = my_pal, name = "timepoint") +
#               ggplot2::theme_classic() + ggplot2::ylab("Late Aire percent of whole")

# print(p_barplot_post_aire)

# # Focus on progenitor population

# ###############################################################################
# # Supplement

# #########################
# # Supplemental Figure 1 #
# #########################
# # UMAP of Aire trace
# tSNE_PCA(mtec, "stage", color = stage_color,
#   save_plot = paste0(save_dir, "figure_1b.pdf"))
# # Marker genes mapped on UMAP
# tSNE_PCA(mtec, "Ccl21a")
# tSNE_PCA(mtec, "Krt5")
# tSNE_PCA(mtec, "Trpm5")
# tSNE_PCA(mtec, "Krt10")
# tSNE_PCA(mtec, "Aire")
# tSNE_PCA(mtec, "GFP")

# # Put plots of UMI and gene counts here

# # Jitter and violin plots of genes that should be consistent
# trio_plots(mtec_wt, geneset = c("Actb", "Gapdh", "B2m"), cell_cycle = FALSE,
#          jitter_and_violin = TRUE, plot_jitter = FALSE, color = stage_color,
#          sep_by = "cluster")

# #########################
# # Supplemental Figure 2 #
# #########################

# # Heatmap of DE genes labeled by genes already known in Aire trace cells
# plot_heatmap(mtec, subset_list = TFs_all,
#         color_list = c("Cdx1", "Utf1", "Tcf7", "Spib", "Cdk4", "Ptma",
#                 "H2afx", "Hmgb1"),
#         color_list2 = c("Aire", "Irf7", "Cited2", "Spib", "Hes1", "Pax1",
#                 "Relb", "Lmo4", "Pou2f3"),
#         cell_color = stage_color)

# #########################
# # Supplemental Figure 3 #
# #########################


# # Violin plots of wt markers
# trio_plots(mtec_wt, geneset = c("Ackr4", "Ccl21a", "Aire"), cell_cycle = FALSE,
#           plot_violin = TRUE, jitter_and_violin = FALSE, plot_jitter = FALSE,
#           color = stage_color, sep_by = "cluster")
# trio_plots(mtec_wt, geneset = c("Krt10", "Trpm5", "Mki67"), cell_cycle = FALSE,
#            plot_violin = TRUE, jitter_and_violin = FALSE, plot_jitter = FALSE,
#            color = stage_color, sep_by = "cluster")

# trio_plots(mtec_wt, geneset = c("Hmgb2", "H2afx", "Hmgb1"), cell_cycle = TRUE,
#          jitter_and_violin = FALSE, plot_jitter = TRUE, plot_violin = FALSE,
#          color = stage_color, sep_by = "cluster")


# # Correlation plot
# master_plot(mtec, "aire_trace", mtec_wt, "wt", stage_color_df)

# # Dot plot of all markers
# markers_to_plot_full <- c("Krt5", "Ccl21a", "Ascl1", "Hes1", "Hmgb2", "Hmgn2",
#                           "Hmgb1", "H2afx", "Ptma", "Aire", "Utf1",
#                           "Fezf2", "Krt10", "Trpm5", "Pou2f3", "Dclk1")

# dot_plot <- Seurat::DotPlot(mtec_wt, genes.plot = rev(markers_to_plot_full),
#                             cols.use = c("blue", "red"), x.lab.rot = T,
#                             plot.legend = T, dot.scale = 8, do.return = T)

# # Dot plot of gene in WT population
# dot_plot <- Seurat::DotPlot(mtec_wt, genes.plot = rev(markers_to_plot),
#                               cols.use = c("blue", "red"), x.lab.rot = T,
#                               plot.legend = T, dot.scale = 8, do.return = T)

# # Heatmap of DE genes labeled by genes already known in WT cells
# plot_heatmap(mtec_wt, subset_list = TFs_all,
#         color_list = c("Cdx1", "Utf1", "Tcf7", "Spib", "Cdk4", "Ptma",
#                 "H2afx", "Hmgb1"),
#         color_list2 = c("Aire", "Irf7", "Cited2", "Spib", "Hes1", "Pax1",
#                 "Relb", "Lmo4", "Pou2f3"),
#         cell_color = stage_color)

# # Marker genes mapped on UMAP for wt population
# tSNE_PCA(mtec_wt, "Ccl21a")
# tSNE_PCA(mtec_wt, "Krt5")
# tSNE_PCA(mtec_wt, "Trpm5")
# tSNE_PCA(mtec_wt, "Krt10")
# tSNE_PCA(mtec_wt, "Aire")

# ########################
# # Supplemntal Figure 4 #
# ########################

# # UMAPs of gene expression at each stage, not including Aire trace

# dev.off()
