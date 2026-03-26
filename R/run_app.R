#' Run the MyToolkits Shiny Dashboard
#'
#' @description
#' Launches the MyToolkits interactive dashboard.
#'
#' @param max_upload_size Numeric. Maximum upload file size in Megabytes (MB). Default is 200 MB.
#'
#' @return No return value, called for side effects.
#' @export
#' @importFrom shiny runApp
#'
run_mytoolkits_shiny <- function(max_upload_size = 200) {
  # 1. 寻找安装包内的 app 目录
  app_dir <- system.file("shiny", "dashboard", package = "MyToolkits")

  if (app_dir == "") {
    stop("Could not find the application directory. Try re-installing `MyToolkits`.", call. = FALSE)
  }

  if (!requireNamespace("bs4Dash", quietly = TRUE)) {
    stop("The 'bs4Dash' package is required. Please install it.", call. = FALSE)
  }

  # 2. 设置上传限制 (MB 转 Bytes)
  # 使用 on.exit 在 App 关闭后恢复默认设置，避免影响用户环境
  old_opt <- getOption("shiny.maxRequestSize")
  options(shiny.maxRequestSize = max_upload_size * 1024^2)
  on.exit(options(shiny.maxRequestSize = old_opt), add = TRUE)

  message(paste(">>> Maximum upload size set to:", max_upload_size, "MB"))

  # 3. 启动 App
  shiny::runApp(app_dir, display.mode = "normal")
}
