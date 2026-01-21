#' Extract Complete Data Table from BreedingValue Object
#'
#' @description
#' Converts the BreedingValue object into a wide-format data frame,
#' combining calculated values (BLUP/BLUE) with raw phenotypic data
#' from all environments.
#'
#' @param object An object of class \code{BreedingValue}.
#'
#' @return A data frame with Line ID, calculated trait value, and raw values for each environment.
#' @export
#'
#' @importFrom tidyr pivot_wider unite
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom methods is
#'
get_complete_table <- function(object) {
  if (!is(object, "BreedingValue")) {
    stop("Input must be a 'BreedingValue' object.")
  }

  # 1. Prepare Result (BLUP/BLUE)
  res_df <- object@result %>%
    rownames_to_column(var = "Line")

  # 2. Prepare Raw Data (Pivot to Wide)
  # Combine Year and Location into a single column header
  raw_wide <- object@raw_data %>%
    unite("Env", year, location, sep = "_", remove = TRUE) %>%
    pivot_wider(
      names_from = "Env",
      values_from = "value",
      names_prefix = "Raw_"
    )

  # 3. Merge
  final_df <- left_join(res_df, raw_wide, by = c("Line" = "line"))

  return(final_df)
}
