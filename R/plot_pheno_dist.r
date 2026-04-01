#' Phenotype Distribution Visualization (Histogram & Raincloud Plot)
#'
#' This function generates a combined plot for phenotype data visualization.
#' It handles standard GWAS phenotype tables. You can plot a single column (e.g., "Mean")
#' or multiple columns (e.g., replicates like "R1", "R2", "R3") to compare their distributions.
#'
#' @param data A dataframe containing the phenotype data.
#' @param val_cols Character or Vector. The column name(s) of the phenotype values.
#'        If a single string (e.g., "Mean"), plots a single distribution.
#'        If a vector (e.g., c("R1", "R2", "R3")), automatically compares them.
#' @param bins Integer. Number of bins for the histogram. Default is 30.
#' @param point_size Numeric. Size of the points in the plot. Default is 1.2.
#' @param fill_color Character or Vector. Custom color(s).
#'
#' @return A patchwork object combining two ggplot2 visualizations.
#' @export
#'
#' @importFrom rlang .data sym !!
#' @importFrom dplyr select all_of filter mutate
#' @importFrom tidyr pivot_longer everything
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_boxplot geom_jitter labs theme_classic theme element_blank element_text scale_fill_manual scale_color_manual after_stat position_nudge coord_flip element_line guides
#' @importFrom gghalves geom_half_violin
#' @importFrom patchwork plot_layout
#'
#' @examples
#' \dontrun{
#' # 1. Single distribution (Mean)
#' # p1 <- plot_pheno_dist(hemicellulose_trait_for_gwas, val_cols = "Mean")
#' # print(p1)
#'
#' # 2. Compare 3 replicates
#' # p2 <- plot_pheno_dist(
#' #   hemicellulose_trait_for_gwas,
#' #   val_cols = c("Hemicellulose-R1", "Hemicellulose-R2", "Hemicellulose-R3")
#' # )
#' # print(p2)
#' }
plot_pheno_dist <- function(data,
                            val_cols,
                            bins = 30,
                            point_size = 1.2,
                            fill_color = NULL) {

  # --- 1. Data Prep ---
  if (!all(val_cols %in% colnames(data))) {
    missing <- val_cols[!val_cols %in% colnames(data)]
    stop(paste("Column(s) not found in data:", paste(missing, collapse = ", ")))
  }

  if (length(val_cols) == 1) {
    plot_data <- data[!is.na(data[[val_cols]]), ]
    plot_data$Group <- "Phenotype"
    val_name <- val_cols
    group_col <- "Group"
    show_legend <- FALSE
    if (is.null(fill_color)) fill_color <- "#377EB8"

  } else {
    plot_data <- data %>%
      dplyr::select(dplyr::all_of(val_cols)) %>%
      tidyr::pivot_longer(cols = tidyr::everything(),
                          names_to = "Group",
                          values_to = "Value") %>%
      dplyr::filter(!is.na(Value))

    plot_data$Group <- factor(plot_data$Group, levels = val_cols)
    val_name <- "Value"
    group_col <- "Group"
    show_legend <- TRUE

    if (is.null(fill_color)) {
      fill_color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")[1:length(val_cols)]
    }
  }

  # --- 2. Plot A: Histogram + Density ---
  p_hist <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!rlang::sym(val_name))) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density), fill = !!rlang::sym(group_col)),
      bins = bins, color = "black", alpha = 0.5, position = "identity"
    ) +
    ggplot2::geom_density(
      ggplot2::aes(color = !!rlang::sym(group_col)),
      linewidth = 1, alpha = 0.2
    ) +
    ggplot2::scale_fill_manual(values = fill_color) +
    ggplot2::scale_color_manual(values = fill_color) +
    ggplot2::labs(x = NULL, y = "Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  # --- 3. Plot B: Raincloud Plot (Robust Version) ---
  p_rain <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!rlang::sym(group_col), y = !!rlang::sym(val_name), fill = !!rlang::sym(group_col))) +

    # 1. Left side (Rain): Pure ggplot2 Jitter (Bypasses gghalves bug)
    ggplot2::geom_jitter(
      width = 0.08, height = 0,
      alpha = 0.5, size = point_size, color = "grey30"
    ) +

    # 2. Middle: Narrow Boxplot (Shifted to x = 0.1)
    ggplot2::geom_boxplot(
      width = 0.08, outlier.shape = NA, alpha = 0.9, color = "black", fill = "white",
      position = ggplot2::position_nudge(x = 0.1, y = 0)
    ) +

    # 3. Right side (Cloud): Half Violin (Shifted to x = 0.15)
    gghalves::geom_half_violin(
      side = "r",
      position = ggplot2::position_nudge(x = 0.15, y = 0),
      alpha = 0.7, color = "black", trim = FALSE
    ) +

    ggplot2::scale_fill_manual(values = fill_color) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = if(length(val_cols) == 1) val_cols else "Phenotype Value") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = if(show_legend) "bottom" else "none",
      legend.title = ggplot2::element_blank(),
      axis.text.y = if(length(val_cols) == 1) ggplot2::element_blank() else ggplot2::element_text(color = "black", size = 10),
      axis.ticks.y = if(length(val_cols) == 1) ggplot2::element_blank() else ggplot2::element_line()
    )

  # --- 4. Combine Plots ---
  combined_plot <- p_hist / p_rain +
    patchwork::plot_layout(heights = c(1, 1.5))

  return(combined_plot)
}
