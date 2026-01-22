#' Plot Multi-Trait or Single-Trait Manhattan Plot from GWAS Results
#'
#' This function generates a Manhattan plot using ggplot2. It supports merging multiple
#' GWAS results into a single plot. It is designed to be flexible with input formats,
#' provided the correct column names are used.
#'
#' @param gwas_input A named list of dataframes (for multi-trait) or a single dataframe (for single-trait).
#' Each dataframe MUST contain at least three columns named exactly:
#' "CHR" (Chromosome), "BP" (Position), and "P" (P-value).
#' @param trait_name Character or Vector. Optional. The name(s) of the trait(s) to appear in the legend.
#' \itemize{
#'   \item If input is a single dataframe: Provide a single string (e.g., "Yield").
#'   \item If input is a list: Provide a vector of strings matching the list length. If NULL, list names are used.
#' }
#' @param threshold Numeric. The -log10(P) threshold line. SNPs above this line are considered significant. Default is 4.
#' @param sample_rate Numeric (0-1). The proportion of insignificant SNPs to retain for plotting to reduce file size and rendering time. Default is 0.5 (50 percent).
#' @param point_size_min Numeric. Point size for insignificant SNPs. Default is 0.3.
#' @param point_size_max Numeric. Point size for significant SNPs. Default is 1.0.
#' @param trait_colors Vector. Custom colors for significant SNPs (one color per trait). If NULL, uses a default high-contrast palette.
#' @param chr_colors Vector of length 2. Colors for alternating chromosomes (background points). Default is c("grey80", "lightblue").
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom purrr imap_dfr
#' @importFrom stringr str_sort
#' @importFrom stats runif setNames
#' @importFrom dplyr filter select mutate group_by summarise arrange left_join case_when n rename
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_size_manual
#'             scale_x_continuous scale_y_continuous geom_hline labs theme_classic
#'             theme element_blank element_text element_line
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' # ------------------------------------------------------------------
#' # Scenario 1: Using the package's helper function for EMMAX files
#' # ------------------------------------------------------------------
#' # If you have the raw .ps output from EMMAX, use read_emmax_result()
#' # to automatically handle the CHR:BP splitting.
#'
#' # Step 1: Read and clean data
#' df_catechin <- read_emmax_result("Catechin.ps")
#'
#' # Step 2: Plot
#' plot_emmax_manhattan(df_catechin, trait_name = "Catechin", threshold = 6)
#'
#' # ------------------------------------------------------------------
#' # Scenario 2: Plotting generic GWAS results (e.g., GAPIT, TASSEL)
#' # ------------------------------------------------------------------
#' # For other tools, you just need to ensure columns are named CHR, BP, P.
#'
#' df_other <- read.csv("other_gwas_result.csv", header = TRUE) %>%
#'   dplyr::rename(CHR = Chromosome, BP = Pos, P = P.value)
#'
#' plot_emmax_manhattan(df_other, trait_name = "Other_Trait", threshold = 6)
#'
#' # ------------------------------------------------------------------
#' # Scenario 3: Merging Multiple Traits
#' # ------------------------------------------------------------------
#'
#' my_list <- list("Catechin" = df_catechin, "Other" = df_other)
#'
#' plot_emmax_manhattan(my_list,
#'                      trait_name = c("Petal-Catechin", "Leaf-Trait"),
#'                      threshold = 5.13,
#'                      trait_colors = c("firebrick", "dodgerblue"))
#' }
plot_emmax_manhattan <- function(gwas_input,
                                 trait_name = NULL,
                                 threshold = 4,
                                 sample_rate = 0.5,
                                 point_size_min = 0.3,
                                 point_size_max = 1.0,
                                 trait_colors = NULL,
                                 chr_colors = c("grey80", "lightblue")) {

  # 1. Handle Input (Single Dataframe vs List) and Naming
  if (is.data.frame(gwas_input)) {
    # If single dataframe, set name from parameter or default to "Trait"
    final_name <- if (!is.null(trait_name)) trait_name[1] else "Trait"
    gwas_list <- stats::setNames(list(gwas_input), final_name)

  } else if (is.list(gwas_input)) {
    # If list, check if user wants to override names
    if (!is.null(trait_name)) {
      if (length(trait_name) != length(gwas_input)) {
        warning("Length of 'trait_name' does not match number of dataframes. Using original list names.")
        gwas_list <- gwas_input
      } else {
        gwas_list <- stats::setNames(gwas_input, trait_name)
      }
    } else {
      # No trait_name provided, check if list has names
      if (is.null(names(gwas_input))) {
        # List has no names and no trait_name provided
        warning("Input list has no names. Assigning Trait_1, Trait_2...")
        gwas_list <- stats::setNames(gwas_input, paste0("Trait_", seq_along(gwas_input)))
      } else {
        gwas_list <- gwas_input
      }
    }
  } else {
    stop("Input 'gwas_input' must be a dataframe or a named list of dataframes.")
  }

  n_traits <- length(gwas_list)

  # 2. Color Logic
  if (is.null(trait_colors)) {
    trait_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                      "#FFFF33", "#A65628", "#F781BF", "#999999", "#000000")
  }

  if (length(trait_colors) < n_traits) {
    warning("Provided colors are fewer than traits. Colors will be recycled.")
    trait_colors <- rep(trait_colors, length.out = n_traits)
  }

  my_trait_colors <- trait_colors[1:n_traits]
  names(my_trait_colors) <- names(gwas_list)

  message(">>> Processing ", n_traits, " trait(s): ", paste(names(gwas_list), collapse=", "))

  # 3. Data Processing & Merging
  # purrr::imap_dfr handles the loop and binding automatically
  combined_df <- imap_dfr(gwas_list, function(df, trait_name) {
    colnames(df) <- toupper(colnames(df))

    if (!all(c("CHR", "BP", "P") %in% colnames(df))) {
      stop(paste("Dataframe for", trait_name, "missing required columns: CHR, BP, P"))
    }

    df %>%
      select(CHR, BP, P) %>%
      mutate(
        trait = trait_name,
        logP = -log10(P),
        CHR = as.character(CHR) # Force character to avoid type mismatch
      )
  })

  # 4. Genomic Position Calculation
  message(">>> Calculating genomic positions...")

  # Natural sort for chromosomes (e.g., 1, 2, ... 10)
  chr_levels <- str_sort(unique(combined_df$CHR), numeric = TRUE)
  combined_df$CHR <- factor(combined_df$CHR, levels = chr_levels)

  chr_len_df <- combined_df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP), .groups = 'drop') %>%
    arrange(CHR) %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len)

  plot_data <- combined_df %>%
    left_join(chr_len_df, by = "CHR") %>%
    mutate(BPcum = BP + tot)

  # 5. Downsampling & Grouping
  message(">>> Downsampling insignificant SNPs to ", sample_rate * 100, "%...")

  plot_data_filtered <- plot_data %>%
    filter(logP >= threshold | runif(n()) <= sample_rate) %>%
    mutate(
      color_group = case_when(
        logP >= threshold ~ trait,
        as.numeric(CHR) %% 2 == 1 ~ "Odd_Chr",
        TRUE ~ "Even_Chr"
      ),
      size_group = ifelse(logP >= threshold, "Sig", "NonSig")
    )

  # Ensure factor levels for correct legend/color mapping
  color_levels <- c(names(gwas_list), "Odd_Chr", "Even_Chr")
  plot_data_filtered$color_group <- factor(plot_data_filtered$color_group, levels = color_levels)

  # 6. Plotting Parameters
  axisdf <- plot_data_filtered %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = 'drop')

  ymax <- ceiling(max(plot_data_filtered$logP, na.rm = TRUE))

  # Define final color mapping
  final_colors <- c(my_trait_colors, "Odd_Chr" = chr_colors[1], "Even_Chr" = chr_colors[2])

  # 7. Generate Plot
  message(">>> Rendering plot...")

  p <- ggplot(plot_data_filtered, aes(x = BPcum, y = logP)) +
    geom_point(aes(color = color_group, size = size_group), alpha = 0.8) +

    scale_color_manual(values = final_colors, breaks = names(gwas_list)) +
    scale_size_manual(values = c("NonSig" = point_size_min, "Sig" = point_size_max), guide = "none") +

    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax + 1)) +

    geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", linewidth = 0.5) +

    labs(x = NULL, y = expression(-log[10](italic(P)))) +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text.x = element_text(size = 10, angle = 0, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line = element_line(linewidth = 0.8),
      panel.grid = element_blank()
    )

  return(p)
}
