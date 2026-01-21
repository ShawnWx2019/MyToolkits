#' Calculate Descriptive Statistics (Scientific Version)
#'
#' @description
#' Performs a comprehensive descriptive statistical analysis.
#' Includes unbiased estimates for Skewness and Kurtosis (SAS/SPSS standard),
#' and Standard Error (SE).
#'
#' @param lab_name A character string specifying the label/name of the variable.
#' @param x A numeric vector of data values.
#' @param na.omit Logical. If \code{TRUE}, NA values are removed before calculation.
#'
#' @return An object of class \code{\linkS4class{DescriptiveStat}}.
#' @export
#'
#' @importFrom stats sd var quantile
#' @importFrom methods new
#'
stat_descriptive <- function(lab_name, x, na.omit = FALSE) {

  # --- 1. 数据清洗与验证 ---
  if (!is.numeric(x)) {
    stop("Error: Input 'x' must be a numeric vector.")
  }

  if (na.omit) {
    x_clean <- x[!is.na(x)]
  } else {
    x_clean <- x
  }

  # 如果不移除 NA 但数据里有 NA，计算将返回 NA，这是正确的 R 行为
  # 但为了后续公式不报错，我们做个临时变量用于计算时刻
  x_calc <- x_clean[!is.na(x_clean)]
  n <- length(x_calc)

  # --- 2. 初始化结果 ---
  m <- NA_real_; s <- NA_real_; var_val <- NA_real_; cv_val <- NA_real_
  se <- NA_real_; skew <- NA_real_; kurt <- NA_real_
  quant <- rep(NA_real_, 5)

  # --- 3. 核心计算 ---
  if (n > 0) {
    # 基础指标
    m <- mean(x_calc)
    quant <- quantile(x_calc, probs = c(0, 0.25, 0.5, 0.75, 1))

    if (n > 1) {
      s <- sd(x_calc)
      var_val <- var(x_calc)
      se <- s / sqrt(n)  # 标准误 SE

      # 变异系数 CV (防除零)
      if (abs(m) > 1e-6) {
        cv_val <- (s / abs(m)) * 100
      } else {
        cv_val <- NA # 均值太接近0，CV无意义
      }

      # --- 更加科学的无偏估计 (Unbiased Skewness & Kurtosis) ---
      # 只有当 n >= 3 时才能计算无偏偏度
      if (n >= 3 && s > 0) {
        # 居中标准化
        z <- (x_calc - m)

        # 偏度公式 (Type 2, used by SPSS/Excel/SAS)
        # G1 = [n / ((n-1)(n-2))] * sum((xi-x_bar)^3) / s^3
        skew <- (n / ((n - 1) * (n - 2))) * sum(z^3) / (s^3)
      }

      # 只有当 n >= 4 时才能计算无偏峰度
      if (n >= 4 && s > 0) {
        # 峰度公式 (Type 2, Excess Kurtosis)
        # G2 = ... (复杂公式校正样本量)
        term1 <- (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))
        term2 <- sum(z^4) / (s^4)
        term3 <- (3 * (n - 1)^2) / ((n - 2) * (n - 3))
        kurt <- term1 * term2 - term3
      }
    }
  }

  # --- 4. 构建输出表格 ---
  tbl <- data.frame(
    name = as.character(lab_name),
    num = n,
    mean = m,
    se = se,          # 新增 SE
    sd = s,
    cv = cv_val,
    var = var_val,
    min = quant[1],
    `Quan_0.25` = quant[2],
    `Quan_0.50` = quant[3],
    `Quan_0.75` = quant[4],
    max = quant[5],
    skew = skew,
    kurt = kurt,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # --- 5. 返回 S4 对象 ---
  new("DescriptiveStat",
      statistics = tbl,
      raw_data = if(n > 0) x_calc else numeric(0),
      variable_name = as.character(lab_name)
  )
}
