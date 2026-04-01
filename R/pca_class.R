#' @title S4 Class for PCA Data
#' @description An S4 class to store and process PLINK PCA results.
#' @slot eigenvec Data.frame containing principal components.
#' @slot eigenval Numeric vector of eigenvalues.
#' @slot variance_exp Numeric vector of percentage variance explained by each PC.
#' @slot group_info Data.frame containing Accession and Group assignments.
#' @export
setClass("PCAData",
         slots = c(
           eigenvec = "data.frame",
           eigenval = "numeric",
           variance_exp = "numeric",
           group_info = "data.frame"
         )
)

#' @title Read PLINK PCA Files
#' @description Reads PLINK .eigenvec and .eigenval files, calculates variance explained.
#' @param eigenvec_file Character. Path to the PLINK .eigenvec file.
#' @param eigenval_file Character. Path to the PLINK .eigenval file.
#' @param group_info Data.frame. Optional. A dataframe containing 'Accession' and 'Group' columns.
#' @return A PCAData S4 object.
#' @importFrom dplyr left_join
#' @export
read_plink_pca <- function(eigenvec_file, eigenval_file, group_info = NULL) {
  
  # Read eigenvalues and calculate percentage of variance explained
  evals <- read.table(eigenval_file, header = FALSE)$V1
  var_exp <- (evals / sum(evals)) * 100
  
  # Read eigenvectors (ignoring the first FID column)
  evecs <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)
  df_pca <- evecs[, -1] 
  colnames(df_pca) <- c("Accession", paste0("PC", 1:(ncol(df_pca)-1)))
  
  # Merge with grouping information if provided
  if (!is.null(group_info)) {
    if (!all(c("Accession", "Group") %in% colnames(group_info))) stop("group_info must contain 'Accession' and 'Group'.")
    df_pca <- dplyr::left_join(df_pca, group_info, by = "Accession")
    df_pca$Group[is.na(df_pca$Group)] <- "Unknown"
    # Ensure Group is a factor to maintain logical ordering
    df_pca$Group <- as.factor(df_pca$Group)
  } else {
    df_pca$Group <- as.factor("All")
    group_info <- data.frame()
  }
  
  new("PCAData", eigenvec = df_pca, eigenval = evals, variance_exp = var_exp, group_info = group_info)
}

#' @title Plot PCA Scree Plot (Elbow Method)
#' @param object A PCAData object.
#' @param n_pc Integer. Number of top principal components to display.
#' @return A ggplot object.
#' @import ggplot2
#' @export
setGeneric("plot_pca_scree", function(object, n_pc = 10) standardGeneric("plot_pca_scree"))

#' @rdname plot_pca_scree
#' @export
setMethod("plot_pca_scree", "PCAData", function(object, n_pc = 10) {
  n_pc <- min(n_pc, length(object@variance_exp))
  df_scree <- data.frame(PC = factor(paste0("PC", 1:n_pc), levels = paste0("PC", 1:n_pc)), Variance = object@variance_exp[1:n_pc])
  
  p <- ggplot2::ggplot(df_scree, ggplot2::aes(x = PC, y = Variance, group = 1)) +
    ggplot2::geom_line(color = "#2C3E50", linewidth = 1) +
    ggplot2::geom_point(color = "#E64B35", size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = "Principal Component", y = "Variance Explained (%)", title = "Scree Plot for PCA")
  return(p)
})

#' @title Plot 2D PCA Scatter Plot
#' @param object A PCAData object.
#' @param pc_x Integer. PC on X-axis.
#' @param pc_y Integer. PC on Y-axis.
#' @param colors Character. Palette name OR matched custom hex vector.
#' @param ellipse Logical. Draw confidence ellipse.
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom rlang .data
#' @export
setGeneric("plot_pca_2d", function(object, pc_x = 1, pc_y = 2, colors = "nature", ellipse = TRUE) standardGeneric("plot_pca_2d"))

#' @rdname plot_pca_2d
#' @export
setMethod("plot_pca_2d", "PCAData", function(object, pc_x = 1, pc_y = 2, colors = "nature", ellipse = TRUE) {
  df_plot <- object@eigenvec
  x_col <- paste0("PC", pc_x)
  y_col <- paste0("PC", pc_y)
  
  if (!all(c(x_col, y_col) %in% colnames(df_plot))) stop("Specified PCs not found.")
  
  x_label <- sprintf("%s (%.2f%%)", x_col, object@variance_exp[pc_x])
  y_label <- sprintf("%s (%.2f%%)", y_col, object@variance_exp[pc_y])
  
  # Use tidy evaluation (.data[[]]) to avoid aes_string deprecation warnings
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = Group, fill = Group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 14), panel.grid = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::labs(x = x_label, y = y_label)
  
  num_groups <- length(levels(df_plot$Group))
  
  # Add confidence ellipses if requested and applicable
  if (ellipse && num_groups > 1 && !("Unknown" %in% df_plot$Group)) {
    p <- p + ggplot2::stat_ellipse(geom = "polygon", alpha = 0.1, type = "norm", level = 0.95)
  }
  
  # Robust coloring logic
  if (exists("mytoolkits_palettes") && length(colors) == 1 && colors %in% names(mytoolkits_palettes)) {
    # Apply built-in package palettes
    p <- p + scale_color_mytoolkits(palette = colors) + scale_fill_mytoolkits(palette = colors)
  } else if (!is.null(names(colors))) {
    # If a named vector is provided, strictly map by name
    p <- p + ggplot2::scale_color_manual(values = colors) + ggplot2::scale_fill_manual(values = colors)
  } else {
    # Fallback to simple positional mapping
    custom_colors <- setNames(rep(colors, length.out = num_groups), levels(df_plot$Group))
    p <- p + ggplot2::scale_color_manual(values = custom_colors) + ggplot2::scale_fill_manual(values = custom_colors)
  }
  
  # Hide legend if no specific groups were defined
  if (num_groups == 1 && levels(df_plot$Group)[1] == "All") p <- p + ggplot2::theme(legend.position = "none")
  
  return(p)
})

#' @title Plot Lower Triangular PCA Pairs Plot
#' @param object A PCAData object.
#' @param n_pc Integer. Number of top PCs.
#' @param colors Character. Palette name OR matched custom hex vector.
#' @return A patchwork object.
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom rlang .data
#' @export
setGeneric("plot_pca_pairs", function(object, n_pc = 5, colors = "nature") standardGeneric("plot_pca_pairs"))

#' @rdname plot_pca_pairs
#' @export
setMethod("plot_pca_pairs", "PCAData", function(object, n_pc = 5, colors = "nature") {
  df_plot <- object@eigenvec
  max_pc <- min(n_pc, length(object@variance_exp))
  num_groups <- length(levels(df_plot$Group))
  p_list <- list()
  
  # Resolve colors before entering the loop
  if (exists("mytoolkits_palettes") && length(colors) == 1 && colors %in% names(mytoolkits_palettes)) {
    custom_colors <- setNames(rep(mytoolkits_palettes[[colors]], length.out = num_groups), levels(df_plot$Group))
  } else if (!is.null(names(colors))) {
    custom_colors <- colors
  } else {
    custom_colors <- setNames(rep(colors, length.out = num_groups), levels(df_plot$Group))
  }
  
  # Generate matrix of plots
  for (row in 1:max_pc) {
    for (col in 1:max_pc) {
      x_var <- paste0("PC", col)
      y_var <- paste0("PC", row)
      
      if (row == col) {
        # Diagonal: Density plot using tidy evaluation
        p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[[x_var]], fill = Group)) +
          ggplot2::geom_density(alpha = 0.5, color = NA) + ggplot2::scale_fill_manual(values = custom_colors) +
          ggplot2::theme_void() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5), plot.margin = ggplot2::margin(2, 2, 2, 2), legend.position = "none")
        p <- p + ggplot2::annotate("text", x = mean(range(df_plot[[x_var]])), y = 0, label = x_var, size = 5, face = "bold", vjust = -1)
        p_list[[length(p_list) + 1]] <- p
        
      } else if (col < row) {
        # Lower Triangle: Scatter plot using tidy evaluation
        p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], color = Group)) +
          ggplot2::geom_point(size = 0.8, alpha = 0.7) + ggplot2::scale_color_manual(values = custom_colors) +
          ggplot2::theme_bw() + ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 6), plot.margin = ggplot2::margin(2, 2, 2, 2), legend.position = "none")
        p_list[[length(p_list) + 1]] <- p
        
      } else {
        # Upper Triangle: Empty spacer
        p_list[[length(p_list) + 1]] <- patchwork::plot_spacer()
      }
    }
  }
  
  # Combine plots and collect legends
  final_plot <- patchwork::wrap_plots(p_list, ncol = max_pc, guides = "collect")
  if (!(num_groups == 1 && levels(df_plot$Group)[1] == "All")) final_plot <- final_plot & ggplot2::theme(legend.position = "right")
  
  return(final_plot)
})

#' @title Plot Interactive 3D PCA
#' @param object A PCAData object.
#' @param pc_x Integer. PC for X-axis.
#' @param pc_y Integer. PC for Y-axis.
#' @param pc_z Integer. PC for Z-axis.
#' @param colors Character. Palette name OR matched custom hex vector.
#' @return A plotly object.
#' @importFrom plotly plot_ly layout
#' @export
setGeneric("plot_pca_3d", function(object, pc_x = 1, pc_y = 2, pc_z = 3, colors = "nature") standardGeneric("plot_pca_3d"))

#' @rdname plot_pca_3d
#' @export
setMethod("plot_pca_3d", "PCAData", function(object, pc_x = 1, pc_y = 2, pc_z = 3, colors = "nature") {
  df_plot <- object@eigenvec
  x_col <- paste0("PC", pc_x); y_col <- paste0("PC", pc_y); z_col <- paste0("PC", pc_z)
  
  if (!all(c(x_col, y_col, z_col) %in% colnames(df_plot))) stop("Specified PCs not found.")
  
  num_groups <- length(levels(df_plot$Group))
  
  # Resolve colors
  if (exists("mytoolkits_palettes") && length(colors) == 1 && colors %in% names(mytoolkits_palettes)) {
    custom_colors <- setNames(rep(mytoolkits_palettes[[colors]], length.out = num_groups), levels(df_plot$Group))
  } else if (!is.null(names(colors))) {
    custom_colors <- colors
  } else {
    custom_colors <- setNames(rep(colors, length.out = num_groups), levels(df_plot$Group))
  }
  
  # Prepare hover text
  hover_text <- paste("<b>Accession:</b>", df_plot$Accession, "<br><b>Group:</b>", df_plot$Group)
  
  # Plotly uses as.formula, so it is unaffected by ggplot2 deprecation warnings
  fig <- plotly::plot_ly(data = df_plot, x = as.formula(paste0("~", x_col)), y = as.formula(paste0("~", y_col)), z = as.formula(paste0("~", z_col)),
                         color = ~Group, colors = custom_colors, type = "scatter3d", mode = "markers",
                         marker = list(size = 5, opacity = 0.8, line = list(color = 'white', width = 0.5)), text = hover_text, hoverinfo = "text")
  
  x_label <- sprintf("%s (%.2f%%)", x_col, object@variance_exp[pc_x])
  y_label <- sprintf("%s (%.2f%%)", y_col, object@variance_exp[pc_y])
  z_label <- sprintf("%s (%.2f%%)", z_col, object@variance_exp[pc_z])
  
  fig <- plotly::layout(fig, title = list(text = "Interactive 3D PCA Structure", font = list(size = 18)),
                        scene = list(xaxis = list(title = x_label), yaxis = list(title = y_label), zaxis = list(title = z_label)),
                        showlegend = ifelse(num_groups > 1, TRUE, FALSE))
  return(fig)
})