#' Show method for DescriptiveStat object
#'
#' @param object A DescriptiveStat object.
#' @export
#' @importFrom methods show
setMethod("show", "DescriptiveStat", function(object) {
  cat("\n=== Descriptive Statistics (Scientific) ===\n")
  cat("Variable : ", object@variable_name, "\n")
  cat("N        : ", object@statistics$num, "\n")

  # 格式化输出核心指标
  cat("\nCore Metrics:\n")
  core_stats <- object@statistics[, c("mean", "se", "sd", "cv", "min", "max")]
  # 统一保留 4 位有效数字打印
  print(format(core_stats, digits = 4), row.names = FALSE)

  cat("\nDistribution Shape (Unbiased):\n")
  cat("Skewness : ", round(object@statistics$skew, 3), "\n")
  cat("Kurtosis : ", round(object@statistics$kurt, 3), "\n")

  cat("\nTip: Use 'get_stat_table(object)' to get the dataframe.\n")
})


#' Extract Descriptive Statistics Table
#'
#' @description
#' Extracts the summary statistics data frame.
#'
#' @param object An object of class \code{DescriptiveStat}.
#'
#' @return A data frame containing all calculated statistics.
#' @export
#' @importFrom methods is
get_stat_table <- function(object) {
  if (!is(object, "DescriptiveStat")) {
    stop("Input must be a 'DescriptiveStat' object.")
  }
  return(object@statistics)
}
