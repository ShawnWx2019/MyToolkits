#' @title MyToolkits Color Palettes
#' @description A collection of high-quality, colorblind-friendly color palettes for academic publishing.
#' @export
mytoolkits_palettes <- list(
  nature = c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"),
  npg_classic = c("#3C5488", "#4DBBD5", "#E64B35", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148"),
  earthy = c("#3D5A80", "#74955B", "#A03530", "#E09F3E", "#E2D1F9", "#317773", "#8AAAE5", "#293241"),
  mp_cold = c("#1F618D", "#76D7C4", "#922B21", "#F1C40F", "#AF7AC5", "#5D6D7E", "#28B463", "#17202A"),
  cyber = c("#2C3E50", "#16A085", "#8E44AD", "#E67E22", "#F39C12", "#D35400", "#C0392B", "#7F8C8D")
)

#' @title Extract MyToolkits Palette Colors
#' @description Returns a color interpolator function for the chosen MyToolkits palette.
#' @param palette Character name of palette. Choices are: "nature", "npg_classic", "earthy", "mp_cold", "cyber".
#' @param reverse Boolean indicating whether the palette should be reversed.
#' @param ... Additional arguments passed to colorRampPalette.
#' @return A function that takes an integer argument and returns a character vector of colors.
#' @importFrom grDevices colorRampPalette
#' @export
pal_mytoolkits <- function(palette = "nature", reverse = FALSE, ...) {
  if (!palette %in% names(mytoolkits_palettes)) {
    stop(paste("Palette", palette, "not found! Available palettes are:", 
               paste(names(mytoolkits_palettes), collapse = ", ")))
  }
  
  pal <- mytoolkits_palettes[[palette]]
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}

#' @title MyToolkits Color Scale for ggplot2
#' @description Apply MyToolkits color palettes to ggplot2 color aesthetics.
#' @param palette Character name of palette.
#' @param discrete Boolean indicating whether color aesthetic is discrete or not.
#' @param reverse Boolean indicating whether the palette should be reversed.
#' @param ... Additional arguments passed to ggplot2::discrete_scale or ggplot2::scale_color_gradientn.
#' @importFrom ggplot2 discrete_scale scale_color_gradientn
#' @export
scale_color_mytoolkits <- function(palette = "nature", discrete = TRUE, reverse = FALSE, ...) {
  pal <- pal_mytoolkits(palette = palette, reverse = reverse)
  
  if (discrete) {
    ggplot2::discrete_scale("colour", paste0("mytoolkits_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_color_gradientn(colours = pal(256), ...)
  }
}

#' @title MyToolkits Fill Scale for ggplot2
#' @description Apply MyToolkits color palettes to ggplot2 fill aesthetics.
#' @param palette Character name of palette.
#' @param discrete Boolean indicating whether fill aesthetic is discrete or not.
#' @param reverse Boolean indicating whether the palette should be reversed.
#' @param ... Additional arguments passed to ggplot2::discrete_scale or ggplot2::scale_fill_gradientn.
#' @importFrom ggplot2 discrete_scale scale_fill_gradientn
#' @export
scale_fill_mytoolkits <- function(palette = "nature", discrete = TRUE, reverse = FALSE, ...) {
  pal <- pal_mytoolkits(palette = palette, reverse = reverse)
  
  if (discrete) {
    ggplot2::discrete_scale("fill", paste0("mytoolkits_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_fill_gradientn(colours = pal(256), ...)
  }
}