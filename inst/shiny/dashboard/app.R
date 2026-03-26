# =======================================================
# MyToolkits Dashboard
# Location: inst/shiny/dashboard/app.R
# =======================================================

library(shiny)
library(bs4Dash)
library(DT)
library(ggplot2)
library(dplyr)
library(readr)
library(shinycssloaders)
library(shinyalert)
library(MyToolkits)

# =======================================================
# UI Definition
# =======================================================
ui <- bs4DashPage(
  title = "MyToolkits Dashboard",
  # dark = NULL, # 保持默认，避免版本兼容问题

  # --- Header ---
  header = bs4DashNavbar(
    title = bs4DashBrand(
      title = "MyToolkits",
      color = "primary",
      href = "#"
    ),
    skin = "light",
    status = "white"
  ),

  # --- Sidebar ---
  sidebar = bs4DashSidebar(
    skin = "light",
    status = "primary",
    elevation = 3,

    bs4SidebarMenu(
      bs4SidebarMenuItem("Dashboard Home", tabName = "home", icon = icon("home")),

      bs4SidebarHeader("Genetic Analysis"),

      bs4SidebarMenuItem(
        "Population Genetics",
        icon = icon("dna"),
        startExpanded = TRUE,
        bs4SidebarMenuSubItem("BLUP/BLUE Calc", tabName = "blup_module", icon = icon("calculator")),
        bs4SidebarMenuSubItem("GWAS Visualizer", tabName = "gwas_module", icon = icon("chart-area"))
      ),

      bs4SidebarHeader("Omics Analysis"),
      bs4SidebarMenuItem("Transcriptomics", tabName = "rnaseq", icon = icon("chart-bar"), badgeLabel = "Dev", badgeColor = "secondary")
    )
  ),

  # --- Body ---
  body = bs4DashBody(
    includeCSS("www/custom.css"),
    useShinyalert(), # 初始化弹窗插件

    bs4TabItems(

      # ---------------- Home Tab ----------------
      bs4TabItem(
        tabName = "home",
        jumbotron(
          title = "Welcome to MyToolkits",
          lead = "A comprehensive platform for Plant Breeding & Multi-Omics Analysis.",
          status = "primary",
          btnName = "Start Analysis",
          href = "#"
        )
      ),

      # ---------------- BLUP Module ----------------
      bs4TabItem(
        tabName = "blup_module",
        fluidRow(
          # Left Panel
          column(
            width = 4,
            bs4Card(
              title = "Data Input & Parameters",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fileInput("blup_file", "Upload Phenotype Data (CSV)", accept = ".csv"),
              uiOutput("blup_col_selectors"),
              selectInput("analysis_method", "Method", choices = c("Random (BLUP)" = "blup", "Fixed (BLUE)" = "blue")),
              actionButton("run_blup", "Calculate", status = "success", icon = icon("play"), width = "100%")
            )
          ),
          # Right Panel
          column(
            width = 8,
            bs4Card(
              title = "Model Summary",
              width = 12,
              collapsible = TRUE,
              verbatimTextOutput("blup_summary_text")
            ),
            bs4Card(
              title = "Calculated Values",
              width = 12,
              downloadButton("download_blup", "Download Result CSV", class = "btn-sm btn-light"),
              hr(),
              withSpinner(DTOutput("blup_table"), type = 4, color = "#0073b7")
            )
          )
        )
      ),

      # ---------------- GWAS Module ----------------
      bs4TabItem(
        tabName = "gwas_module",
        fluidRow(
          # Left Panel: Settings
          column(
            width = 3,
            bs4Card(
              title = "Plot Settings",
              status = "danger",
              solidHeader = TRUE,
              width = 12,

              # 1. File Upload (允许多选)
              fileInput("gwas_file", "Upload GWAS Result(s)",
                        multiple = TRUE,
                        accept = c(".txt", ".csv", ".assoc", ".ps"),
                        placeholder = "Select one or more files"),

              helpText("Supported: Generic CSV (with Header) or EMMAX (No Header)."),

              # 2. Format Toggle
              checkboxInput("is_emmax", "Is EMMAX format?", value = TRUE),
              helpText("Checked: Use read_emmax_result(). Unchecked: Read header."),

              hr(),

              # 3. Parameters
              textInput("gwas_trait_name", "Legend Name (Single file only)", value = "Yield"),
              helpText("For multiple files, filenames are used automatically."),

              sliderInput("gwas_threshold", "-log10(P) Threshold", min = 0, max = 15, value = 4, step = 0.5),

              sliderInput("gwas_sample_rate", "Downsample Non-Sig SNPs", min = 0.01, max = 1, value = 0.5, step = 0.05),

              hr(),
              actionButton("plot_gwas", "Generate Plot", status = "danger", icon = icon("paint-brush"), width = "100%")
            )
          ),

          # Right Panel: Plot
          column(
            width = 9,
            bs4Card(
              title = "Manhattan Plot",
              width = 12,
              height = "800px",
              maximizable = TRUE,
              withSpinner(
                plotOutput("manhattan_plot", height = "700px"),
                type = 4,
                color = "#dc3545",
                size = 1
              )
            )
          )
        )
      )
    )
  )
)

# =======================================================
# Server Logic
# =======================================================
server <- function(input, output, session) {

  # -----------------------------------------------------
  # 1. BLUP / BLUE Logic
  # -----------------------------------------------------
  blup_data <- reactive({
    req(input$blup_file)
    tryCatch({
      read.csv(input$blup_file$datapath)
    }, error = function(e) {
      shinyalert::shinyalert("Load Error", paste("Failed to read CSV:", e$message), type = "error")
      return(NULL)
    })
  })

  output$blup_col_selectors <- renderUI({
    req(blup_data())
    cols <- colnames(blup_data())
    tagList(
      selectInput("col_trait", "Select Trait Column", choices = cols, selected = cols[grep("trait|pheno", tolower(cols))[1]]),
      selectInput("col_line", "Select Line/Genotype Column", choices = cols, selected = cols[grep("line|geno|taxa", tolower(cols))[1]]),
      selectInput("col_year", "Select Year Column", choices = cols, selected = cols[grep("year|env", tolower(cols))[1]]),
      selectInput("col_loc", "Select Location Column", choices = cols, selected = cols[grep("loc", tolower(cols))[1]]),
      selectInput("col_val", "Select Value Column", choices = cols, selected = cols[grep("val|y", tolower(cols))[1]]),
      textInput("target_trait_level", "Target Trait Level (e.g. Yield)", value = "Yield")
    )
  })

  blup_result_obj <- eventReactive(input$run_blup, {
    req(blup_data(), input$col_trait)

    # 简单的进度提示
    withProgress(message = 'Calculating...', value = 0.5, {
      tryCatch({
        MyToolkits::get_unbiase_by_year_loc(
          x = blup_data(),
          Tag = input$target_trait_level,
          method = input$analysis_method,
          year.name = input$col_year,
          loc.name = input$col_loc,
          trait.name = input$col_trait,
          line.name = input$col_line,
          value.name = input$col_val
        )
      }, error = function(e) {
        shinyalert::shinyalert("Model Failed", e$message, type = "error")
        return(NULL)
      })
    })
  })

  output$blup_summary_text <- renderPrint({
    req(blup_result_obj())
    obj <- blup_result_obj()
    cat("Method:", toupper(obj@method), "\nTrait:", obj@trait, "\nHeritability:", round(obj@heritability, 4), "\nInfo:", obj@model_info$message)
  })

  final_blup_df <- reactive({
    req(blup_result_obj())
    MyToolkits::get_complete_table(blup_result_obj())
  })

  output$blup_table <- renderDT({
    req(final_blup_df())
    datatable(final_blup_df(), options = list(pageLength = 10, scrollX = TRUE, dom = 'tip'))
  })

  output$download_blup <- downloadHandler(
    filename = function() { paste0("Result_", Sys.Date(), ".csv") },
    content = function(file) { write.csv(final_blup_df(), file, row.names = FALSE) }
  )


  # -----------------------------------------------------
  # 2. GWAS Visualization Logic (Multi-Trait Support)
  # -----------------------------------------------------

  gwas_data_list <- reactive({
    req(input$gwas_file)
    files <- input$gwas_file # 这里是 data.frame: name, size, datapath

    # 大文件提示
    if(sum(files$size) > 500 * 1024^2) {
      shinyalert::shinyalert("Large Data", "Total file size > 500MB. Plotting might take a while.", type = "info")
    }

    result_list <- list()

    tryCatch({

      # 循环读取所有上传的文件
      for(i in 1:nrow(files)) {

        f_path <- files$datapath[i]
        f_name <- files$name[i]
        # 使用不带后缀的文件名作为 Trait Name (e.g., "Yield_2023.ps" -> "Yield_2023")
        clean_name <- tools::file_path_sans_ext(f_name)

        df <- NULL

        # --- 读取逻辑分支 ---
        if (input$is_emmax) {
          # EMMAX 模式：使用 read_emmax_result
          df <- MyToolkits::read_emmax_result(f_path)

          # 防御性转换：确保列名为大写
          names(df) <- toupper(names(df))

        } else {
          # 通用模式：带表头读取
          df <- read.table(f_path, header = TRUE, stringsAsFactors = FALSE)
          names(df) <- toupper(names(df))

          # 智能匹配 CHR/BP/P
          if (!all(c("CHR", "BP", "P") %in% names(df))) {
            # 尝试寻找常见变体
            chr_idx <- grep("CHR|CHROM", names(df))[1]
            bp_idx  <- grep("BP|POS", names(df))[1]
            p_idx   <- grep("PVAL|P.VAL|P_VAL|P", names(df))[1]

            if (!is.na(chr_idx) && !is.na(bp_idx) && !is.na(p_idx)) {
              names(df)[c(chr_idx, bp_idx, p_idx)] <- c("CHR", "BP", "P")
            } else {
              stop(paste0("File '", f_name, "' is missing CHR, BP, or P columns."))
            }
          }
        }

        # 简单校验
        if(!is.numeric(df$P)) stop(paste0("P-value column in '", f_name, "' must be numeric."))

        # 存入 List
        result_list[[clean_name]] <- df
      }

      return(result_list)

    }, error = function(e) {
      shinyalert::shinyalert("Data Read Error", e$message, type = "error")
      return(NULL)
    })
  })

  output$manhattan_plot <- renderPlot({

    # 获取数据 List
    data_list <- gwas_data_list()
    req(data_list)

    # 获取参数 (隔离，防止滑动滑块直接重绘)
    thresh <- isolate(input$gwas_threshold)
    s_rate <- isolate(input$gwas_sample_rate)
    user_name <- isolate(input$gwas_trait_name)

    message(">>> Starting GWAS Plot rendering...")

    # 准备绘图输入
    final_input <- NULL
    final_trait_name <- NULL

    # 如果只上传了1个文件，且用户在输入框里填了名字，则使用用户填的名字
    if (length(data_list) == 1 && user_name != "") {
      final_input <- data_list[[1]] # 取出单个 Dataframe
      final_trait_name <- user_name
    } else {
      # 多个文件，或用户没填名字 -> 直接传 List，函数会自动用 List 的 names 做图例
      final_input <- data_list
      final_trait_name <- NULL
    }

    tryCatch({
      MyToolkits::plot_emmax_manhattan(
        gwas_input = final_input,
        trait_name = final_trait_name,
        threshold = thresh,
        sample_rate = s_rate
      )
    }, error = function(e) {
      shinyalert::shinyalert("Plotting Error", e$message, type = "error")
      return(NULL)
    })

  }) %>%
    bindEvent(input$plot_gwas) # 必须点击按钮才触发

}

shinyApp(ui, server)
