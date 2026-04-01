#' @title S4 Class for Admixture Data
#' @description An S4 class to store and process Admixture structure results.
#' @slot samples Character vector of sample IDs.
#' @slot q_list List of data.frames containing Q matrices for each K.
#' @slot k_range Numeric vector of K values.
#' @slot long_data Data.frame in long format for ggplot2.
#' @slot accid_order Character vector of ordered sample IDs for plotting.
#' @slot boundaries Numeric vector of X-axis positions for vertical group separator lines.
#' @slot group_info Data.frame containing the Accession and its assigned Group based on the reference K.
#' @export
setClass("AdmixData",
         slots = c(
           samples = "character",
           q_list = "list",
           k_range = "numeric",
           long_data = "data.frame",
           accid_order = "character",
           boundaries = "numeric",
           group_info = "data.frame"  
         )
)

#' @title Read Admixture Output Files
#' @description Reads a .fam file and corresponding .Q files to create an AdmixData object.
#' @param fam_file Character. Path to the PLINK .fam file.
#' @param q_prefix Character. Prefix of the .Q files (e.g., "independent_snp").
#' @param k_range Numeric vector. Range of K values to read (e.g., 2:4).
#' @return An AdmixData S4 object.
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select starts_with
#' @export
read_admixture <- function(fam_file, q_prefix, k_range) {
  # 1. Read sample IDs from the .fam file
  fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  samples <- fam$V1
  
  # 2. Iteratively read Q matrices into a list
  q_list <- list()
  for (k in k_range) {
    file_name <- paste0(q_prefix, ".", k, ".Q")
    if (!file.exists(file_name)) stop(paste("File not found:", file_name))
    q_df <- read.table(file_name, header = FALSE)
    colnames(q_df) <- paste0("K", k, "_", 1:k)
    q_list[[as.character(k)]] <- q_df
  }
  
  # 3. Combine Q matrices into a single wide dataframe
  df_combined <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  for (k in k_range) {
    df_combined <- cbind(df_combined, q_list[[as.character(k)]])
  }
  
  # 4. Convert wide dataframe to long format for ggplot2
  long_df <- df_combined %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Component", values_to = "Proportion") %>%
    dplyr::mutate(
      K_level = sub("_.*", "", Component),
      K_level = factor(K_level, levels = paste0("K", k_range), labels = paste0("K=", k_range))
    )
  
  # 5. Initialize the S4 object
  new("AdmixData",
      samples = samples,
      q_list = q_list,
      k_range = k_range,
      long_data = as.data.frame(long_df),
      accid_order = samples,
      boundaries = numeric(0),
      group_info = data.frame()  # <--- INITIALIZE EMPTY DATAFRAME
  )
}

#' @title Order Samples Based on a Specific K
#' @description Sorts the samples by grouping them based on their dominant ancestry component at a specific K. Generates boundary lines and assigns samples to sub-populations.
#' @param object An AdmixData object.
#' @param ordered_by_k Numeric. The K value used as the reference for sorting and sub-population assignment.
#' @return An updated AdmixData object.
#' @export
setGeneric("order_admix_accid", function(object, ordered_by_k) standardGeneric("order_admix_accid"))

#' @rdname order_admix_accid
#' @export
setMethod("order_admix_accid", "AdmixData", function(object, ordered_by_k) {
  # Validate if the specified K exists in the object
  if (!as.character(ordered_by_k) %in% names(object@q_list)) {
    stop("The specified ordered_by_k is not in the loaded K range.")
  }
  
  # Extract the reference Q matrix
  q_ref <- object@q_list[[as.character(ordered_by_k)]]
  
  # Find the dominant ancestral component index for each sample
  dom_comp <- apply(q_ref, 1, which.max)
  
  # Order groups based on their population size (smallest to largest)
  comp_counts <- table(dom_comp)
  group_order <- as.integer(names(sort(comp_counts)))
  
  # Construct a temporary dataframe for sorting
  df_sort <- data.frame(
    Sample = object@samples,
    Dom_Group = dom_comp,
    Dom_Prop = apply(q_ref, 1, max),
    stringsAsFactors = FALSE
  )
  
  # Map physical group indices to a ranked factor based on size
  df_sort$Group_Rank <- factor(df_sort$Dom_Group, levels = group_order)
  last_group <- tail(group_order, n = 1)
  
  # Sort logic: descending order for all groups EXCEPT the last/largest one (ascending)
  df_sort <- df_sort %>%
    dplyr::mutate(
      Sort_Score = ifelse(Dom_Group == last_group, Dom_Prop, -Dom_Prop)
    ) %>%
    dplyr::arrange(Group_Rank, Sort_Score)
  
  # Calculate X-axis intercept coordinates for vertical separator lines
  sorted_counts <- table(df_sort$Group_Rank)
  bounds <- cumsum(sorted_counts)[-length(sorted_counts)] + 0.5
  
  # ========================================================
  # NEW LOGIC: Construct the Group Information Dataframe
  # Map the sorted factor levels to sequential names: Group1, Group2...
  # ========================================================
  df_group <- data.frame(
    Accession = df_sort$Sample,
    Group = paste0("Group", as.numeric(df_sort$Group_Rank)),
    stringsAsFactors = FALSE
  )
  
  # Update all slots in the S4 object
  object@accid_order <- df_sort$Sample
  object@boundaries <- as.numeric(bounds)
  object@group_info <- df_group # Store the dataframe
  
  return(object)
})

#' @title Plot Admixture Structure
#' @description Generates a stacked barplot using ggplot2 with auto lineage-tracking across K values.
#' @param object An AdmixData object. Must run order_admix_accid() first.
#' @param colors Character. Either a MyToolkits palette name (e.g., "nature") OR a named character vector of custom hex colors.
#' @return A ggplot object.
#' @import ggplot2
#' @export
setGeneric("plot_admixture", function(object, colors = "nature") standardGeneric("plot_admixture"))

#' @rdname plot_admixture
#' @export
setMethod("plot_admixture", "AdmixData", function(object, colors = "nature") {
  
  df_plot <- object@long_data
  
  # Safety check for sample ordering
  if (length(object@accid_order) == 0) {
    stop("Sample order not found! Please run order_admix_accid() before plotting.")
  }
  
  # Lock the sample factor levels strictly to the calculated order
  df_plot$Sample <- factor(df_plot$Sample, levels = object@accid_order)
  
  # Initialize the base ggplot object
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Sample, y = Proportion, fill = Component)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = NA) +
    ggplot2::facet_grid(K_level ~ .) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.2, "lines"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 12, face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::labs(x = "Samples (Sorted)", y = "Ancestry Proportion")
  
  # Add vertical separator lines if boundaries exist
  if (length(object@boundaries) > 0) {
    p <- p + ggplot2::geom_vline(xintercept = object@boundaries, color = "black", linewidth = 0.8)
  }
  
  # ========================================================
  # Core Color Logic: Lineage Tracking & Auto-Alignment
  # ========================================================
  if (exists("mytoolkits_palettes") && length(colors) == 1 && colors %in% names(mytoolkits_palettes)) {
    
    base_palette <- mytoolkits_palettes[[colors]]
    k_vals <- sort(object@k_range)
    q_list <- object@q_list
    
    # 1. Assign base colors to the components of the smallest K
    k_min <- k_vals[1]
    col_map <- setNames(base_palette[1:k_min], colnames(q_list[[as.character(k_min)]]))
    next_col_idx <- k_min + 1
    
    # 2. Track lineage layer by layer to inherit colors downwards
    if (length(k_vals) > 1) {
      for (i in 2:length(k_vals)) {
        q_prev <- as.matrix(q_list[[as.character(k_vals[i-1])]])
        q_curr <- as.matrix(q_list[[as.character(k_vals[i])]])
        
        # Calculate similarity matrix between two K layers (Matrix Cross-product)
        sim_mat <- t(q_prev) %*% q_curr
        assigned_prev <- c()
        assigned_curr <- c()
        
        # Greedy matching: find highest overlapping populations to inherit colors
        for (step in 1:ncol(q_prev)) {
          valid_sim <- sim_mat
          if (length(assigned_prev) > 0) valid_sim[assigned_prev, ] <- -1
          if (length(assigned_curr) > 0) valid_sim[, assigned_curr] <- -1
          
          max_idx <- arrayInd(which.max(valid_sim), dim(valid_sim))
          p_idx <- max_idx[1, 1]
          c_idx <- max_idx[1, 2]
          
          assigned_prev <- c(assigned_prev, p_idx)
          assigned_curr <- c(assigned_curr, c_idx)
          
          # Inherit color from the matching parent component
          col_map[colnames(q_curr)[c_idx]] <- col_map[colnames(q_prev)[p_idx]]
        }
        
        # For newly differentiated populations, assign next available palette color
        unassigned <- setdiff(1:ncol(q_curr), assigned_curr)
        for (uc in unassigned) {
          # Prevent out-of-bounds error if the palette length is exceeded
          safe_idx <- ifelse(next_col_idx > length(base_palette), length(base_palette), next_col_idx)
          col_map[colnames(q_curr)[uc]] <- base_palette[safe_idx]
          next_col_idx <- next_col_idx + 1
        }
      }
    }
    
    # Apply the dynamically mapped lineage colors
    p <- p + ggplot2::scale_fill_manual(values = col_map)
    
  } else {
    # ========================================================
    # Custom Color Logic: Support user-provided named/unnamed vectors
    # ========================================================
    if (!is.null(names(colors))) {
      # Strictly apply named vector (e.g., c("K2_1" = "red", ...))
      p <- p + ggplot2::scale_fill_manual(values = colors)
    } else {
      # Fallback: Unnamed color vector mapped iteratively
      unique_comps <- unique(df_plot$Component)
      custom_colors <- setNames(rep(colors, length.out = length(unique_comps)), unique_comps)
      p <- p + ggplot2::scale_fill_manual(values = custom_colors)
    }
  }
  
  return(p)
})