#' Phenotype Data Transformation
#'
#' A tidyverse-compatible vectorized function to transform numeric vectors.
#' Perfectly suited for use inside \code{dplyr::mutate()} to prepare phenotype
#' data for GWAS or other parametric statistical analyses.
#'
#' @param x A numeric vector.
#' @param method Character. The transformation method. Options are:
#'        "log2", "log10", "zscore", "arcsine_sqrt", "int" (Inverse Normal Transformation), or "logit".
#' @param offset Numeric. A constant to add to 'x' before logarithmic transformations
#'        to avoid log(0) or negative values. Default is 0.
#'
#' @return A numeric vector of the transformed values, of the same length as the input.
#' @export
#'
#' @importFrom stats sd qnorm
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' # Example dataset
#' df <- data.frame(
#'   Line = c("A1", "A2", "A3", "A4", "A5"),
#'   Yield = c(10, 20, 15, 0, NA),
#'   Disease_Rate = c(0.1, 0.5, 0.9, 0.2, 0.4)
#' )
#'
#' # Standard tidyverse pipeline usage
#' df_transformed <- df %>%
#'   mutate(
#'     # Add offset=1 to avoid log2(0)
#'     Yield_log2 = transform_pheno(Yield, "log2", offset = 1),
#'
#'     # Standard GWAS Inverse Normal Transformation
#'     Yield_INT = transform_pheno(Yield, "int"),
#'
#'     # Arcsine square root for proportions/percentages
#'     Disease_Arc = transform_pheno(Disease_Rate, "arcsine_sqrt")
#'   )
#'
#' # Batch transform multiple columns (e.g., Rep1, Rep2, Rep3)
#' # df %>% mutate(across(starts_with("Rep"), ~ transform_pheno(.x, "int")))
#' }
transform_pheno <- function(x,
                            method = c("log2", "log10", "zscore", "arcsine_sqrt", "int", "logit"),
                            offset = 0) {

  # Ensure the user picks exactly one valid method
  method <- match.arg(method)

  if (!is.numeric(x)) {
    stop("Input 'x' must be a numeric vector.")
  }

  res <- switch(method,

                "log2" = {
                  if (any(x + offset <= 0, na.rm = TRUE)) {
                    warning("log2: 'x + offset' contains values <= 0, which will result in NaNs/Inf. Consider adjusting the 'offset' parameter.")
                  }
                  log2(x + offset)
                },

                "log10" = {
                  if (any(x + offset <= 0, na.rm = TRUE)) {
                    warning("log10: 'x + offset' contains values <= 0, which will result in NaNs/Inf. Consider adjusting the 'offset' parameter.")
                  }
                  log10(x + offset)
                },

                "zscore" = {
                  # Standardize to Mean = 0, SD = 1
                  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
                },

                "arcsine_sqrt" = {
                  # Typically used for percentages/proportions (0 to 1)
                  if (any(x < 0 | x > 1, na.rm = TRUE)) {
                    warning("arcsine_sqrt: input should be proportions strictly between 0 and 1. Out of bounds values will return NaN.")
                  }
                  asin(sqrt(x))
                },

                "int" = {
                  # Inverse Normal Transformation (Rank-based)
                  # Using Blom's formula: (r - 3/8) / (n + 1/4)
                  # This is the gold standard for GWAS phenotype normalization
                  r <- rank(x, na.last = "keep", ties.method = "average")
                  n <- sum(!is.na(x))
                  stats::qnorm((r - 0.375) / (n + 0.25))
                },

                "logit" = {
                  # Log-odds typically used for probability values (0 to 1)
                  if (any(x <= 0 | x >= 1, na.rm = TRUE)) {
                    warning("logit: input should be strictly > 0 and < 1. Values at 0 or 1 will return Inf/-Inf.")
                  }
                  log(x / (1 - x))
                }
  )

  return(res)
}
