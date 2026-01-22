#' Read EMMAX Output File (.ps) and Format for Manhattan Plot
#'
#' This function reads the standard output file (.ps) from EMMAX.
#' It assumes the first column is "CHR:BP" (e.g., "A01:2518"), the second is Beta, and the third is P-value.
#' It automatically splits the first column into CHR and BP.
#'
#' @param file Path to the EMMAX output file.
#' @param sep Character. The separator used in the SNP ID column. Default is ":".
#'
#' @return A dataframe with columns: SNP, CHR, BP, Beta, P.
#' @export
#' @importFrom dplyr rename mutate select
#' @importFrom tidyr separate
#' @importFrom utils read.table
#'
#' @examples
#' \dontrun{
#' # Read the raw file
#' df <- read_emmax_result("Catechin.ps")
#'
#' # Then plot it
#' plot_emmax_manhattan(df)
#' }
read_emmax_result <- function(file, sep = ":") {

  # 1. 读取数据 (EMMAX .ps 文件通常没有表头)
  raw_data <- utils::read.table(file, header = FALSE, stringsAsFactors = FALSE)

  # 2. 检查列数
  if (ncol(raw_data) < 3) {
    stop("EMMAX result file should have at least 3 columns (SNP, Beta, P-value).")
  }

  # 3. 处理数据
  # 假设格式: V1="A01:2518", V2=Beta, V3=P-value
  clean_data <- raw_data %>%
    # 重命名列
    dplyr::rename(SNP = V1, Beta = V2, P = V3) %>%
    # 拆分 SNP 列为 CHR 和 BP
    tidyr::separate(col = SNP, into = c("CHR", "BP"), sep = sep, remove = FALSE) %>%
    # 确保 BP 是数字类型 (这一点非常重要，否则画图排序会乱)
    dplyr::mutate(BP = as.numeric(BP)) %>%
    # 选择需要的列
    dplyr::select(SNP, CHR, BP, Beta, P)

  return(clean_data)
}
