#' @title S4 Class for Phylogenetic Tree Data
#' @description An S4 class to store and process phylogenetic tree results, integrating population group information.
#' @slot phylo A phylo object representing the raw or ladderized tree.
#' @slot group_info Data.frame containing Accession and Group assignments.
#' @slot grouped_tree A treedata object (from ggtree::groupOTU) with group mappings for ancestral branches.
#' @export
setClass("PhyloTreeData",
         slots = c(
           phylo = "ANY",          # Accepts ape::phylo
           group_info = "data.frame",
           grouped_tree = "ANY"    # Accepts ggtree::treedata
         )
)

#' @title Read and Process Phylogenetic Tree
#' @description Reads a .nwk file, automatically ladderizes it for better visual representation, and injects population group information into the tree nodes.
#' @param nwk_file Character. Path to the .nwk tree file.
#' @param group_info Data.frame. Optional. A dataframe containing 'Accession' and 'Group' columns (e.g., exported from AdmixData@group_info).
#' @param ladderize Logical. Whether to ladderize the tree to smooth out branch lengths (default: TRUE).
#' @return A PhyloTreeData S4 object.
#' @importFrom ape read.tree ladderize
#' @importFrom ggtree groupOTU
#' @export
read_phylo_tree <- function(nwk_file, group_info = NULL, ladderize = TRUE) {
  
  # 1. Read raw tree
  if (!file.exists(nwk_file)) stop(paste("File not found:", nwk_file))
  tree_raw <- ape::read.tree(nwk_file)
  
  # 2. Ladderize (sort nodes for smoother visualization)
  if (ladderize) {
    tree_processed <- ape::ladderize(tree_raw, right = FALSE)
  } else {
    tree_processed <- tree_raw
  }
  
  # 3. Inject group info (groupOTU)
  if (!is.null(group_info)) {
    if (!all(c("Accession", "Group") %in% colnames(group_info))) {
      stop("group_info must contain 'Accession' and 'Group' columns.")
    }
    # Ensure Group is treated as a factor
    group_info$Group <- as.factor(group_info$Group)
    
    # Split accessions into a list by group for ggtree parsing
    tip_groups <- split(group_info$Accession, group_info$Group)
    
    # Generate grouped treedata object
    tree_grouped <- ggtree::groupOTU(tree_processed, tip_groups)
  } else {
    tree_grouped <- tree_processed
    group_info <- data.frame()
  }
  
  # 4. Initialize S4 object
  new("PhyloTreeData",
      phylo = tree_processed,
      group_info = group_info,
      grouped_tree = tree_grouped
  )
}

#' @title Plot Phylogenetic Tree
#' @description Generates a highly customizable phylogenetic tree using ggtree, with automatic lineage color matching.
#' @param object A PhyloTreeData object.
#' @param layout Character. The layout of the tree (e.g., "equal_angle", "circular", "rectangular", "daylight"). Default: "equal_angle" (best for unrooted population trees).
#' @param colors Character. Either a MyToolkits palette name or a named character vector of custom hex colors.
#' @param linewidth Numeric. The thickness of the tree branches (default: 0.5).
#' @return A ggplot/ggtree object.
#' @import ggplot2
#' @import ggtree
#' @export
setGeneric("plot_phylo_tree", function(object, layout = "equal_angle", colors = "nature", linewidth = 0.5) standardGeneric("plot_phylo_tree"))

#' @rdname plot_phylo_tree
#' @export
setMethod("plot_phylo_tree", "PhyloTreeData", function(object, layout = "equal_angle", colors = "nature", linewidth = 0.5) {
  
  t_data <- object@grouped_tree
  has_groups <- nrow(object@group_info) > 0
  
  if (has_groups) {
    # Generate base plot mapping color to the injected 'group' metadata
    p <- ggtree::ggtree(t_data, layout = layout, linewidth = linewidth, ggplot2::aes(color = group))
    
    # Fetch exact group names
    group_levels <- levels(object@group_info$Group)
    if (is.null(group_levels)) group_levels <- unique(object@group_info$Group)
    num_groups <- length(group_levels)
    
    # Color Logic Resolution
    if (exists("mytoolkits_palettes") && length(colors) == 1 && colors %in% names(mytoolkits_palettes)) {
      pal <- mytoolkits_palettes[[colors]]
      custom_colors <- setNames(rep(pal, length.out = num_groups), group_levels)
    } else if (!is.null(names(colors))) {
      # Apply user-defined named vector (e.g., match_colors)
      custom_colors <- colors
    } else {
      # Fallback to simple mapping
      custom_colors <- setNames(rep(colors, length.out = num_groups), group_levels)
    }
    
    # Append "0" -> grey for unknown ancestral branches (lineage breaks)
    final_colors <- c(custom_colors, "0" = "grey70")
    
    # Apply scales
    p <- p + ggplot2::scale_color_manual(
      values = final_colors,
      breaks = names(custom_colors) # Strictly hide "0" from the legend
    ) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::labs(color = "Group")
    
  } else {
    # Fallback if no group_info is provided (simple black tree)
    p <- ggtree::ggtree(t_data, layout = layout, linewidth = linewidth, color = "black")
  }
  
  # Automatically add a scale bar for unrooted layouts to indicate genetic distance
  if (layout %in% c("equal_angle", "daylight", "unrooted")) {
    p <- p + ggtree::geom_treescale(x = 0, y = 0, color = "black", fontsize = 4, linesize = 0.5)
  }
  
  return(p)
})