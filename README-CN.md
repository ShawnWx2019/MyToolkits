[![Read in English](https://img.shields.io/badge/README-English-blue.svg)](./README.md)
# MyToolkits <img src="man/figures/logo.png" align="right" height="139" alt="" />

MyToolkits 是一个专门为植物育种、数量遗传学以及多组学研究数据分析设计的 R 包。基于稳健的 S4 面向对象类构建，它提供了从多环境试验（BLUP/BLUE）、科学描述性统计、GWAS 可视化、表型数据转换到群体结构分析的完整工作流程。

## 功能特点

- **多环境试验分析**：使用 `lme4` 计算 BLUP（最佳线性无偏预测）和 BLUE（最佳线性无偏估计）值。自动估计广义遗传力（$H^2$）。
- **科学描述性统计**：生成可直接用于学术发表的统计摘要，包括标准误（SE）以及采用 SAS/SPSS 标准的无偏偏度和峰度估计。
- **表型数据处理**：提供兼容 Tidyverse 的函数来对性状进行标准化处理（例如逆正态变换、Z-score），并生成云雨分布图以在下游分析前评估数据质量。
- **灵活的 GWAS 可视化**：生成高度兼容 EMMAX、GAPIT 和 TASSEL 输出的曼哈顿图。支持多性状叠加对比和大规模数据的自动降采样。
- **群体结构可视化**：基于综合的 S4 类工具，使用内置的学术级且对色盲友好的调色板，对 PLINK PCA、ADMIXTURE Q 矩阵（支持自动谱系颜色追踪）以及系统发育树进行高质量可视化。

## 安装

您可以从 GitHub 安装开发版的 MyToolkits：
*我们强烈建议在安装时设置 `build_vignettes = TRUE`，以便获取并编译包内置的详细离线教程。*

```r
# 安装 devtools 包（如未安装）
# install.packages("devtools")
devtools::install_github("ShawnWx2019/MyToolkits", build_vignettes = TRUE)
```

## 使用示例

### 1. 育种值分析（BLUP & BLUE）

计算基因型在多年份和多地点的育种值。该函数可自动处理基因型与环境的互作效应（GxE）。

```r
library(MyToolkits)

# --- 1. 计算 BLUPs 并提取遗传力 ---
blup_res <- get_unbiase_by_year_loc(
  x = my_pheno_data,
  Tag = "Yield",           # 要分析的具体性状水平
  method = "blup",         # "blup" 或 "blue"
  trait.name = "TraitID",  
  year.name = "Year",      
  loc.name = "Loc",        
  line.name = "Genotype",  
  value.name = "Value"     
)

print(blup_res) # 打印查看 S4 对象摘要及模型状态

# --- 2. 导出完整数据表 ---
# 将计算的 BLUPs 与原始数据合并
final_table <- get_complete_table(blup_res)
```

### 2. 描述性统计

计算详细的统计指标，使用科学的无偏估计（特别适合小样本量）。

```r
phenotypes <- c(10.5, 12.1, 9.8, 11.0, 10.2, 13.5, NA, 10.8)

stats_obj <- stat_descriptive(lab_name = "Plant_Height", x = phenotypes, na.omit = TRUE)

# 访问并提取结果表格
get_stat_table(stats_obj)
```

### 3. 表型数据转换与分布

在进行 GWAS 或下游统计建模之前，清洗并可视化您的表型数据。

📖 **教程文档 (Vignettes)**：有关为何以及如何标准化表型数据的全面指南，请通过 `vignette("description_pheno_data_CN", package = "MyToolkits")` 在 R 中阅读内置教程，或直接在 GitHub 上查看源码：
- [表型数据标准化说明 (中文版)](./vignettes/description_pheno_data_CN.qmd)
- [Phenotypic Data Description (English)](./vignettes/description_pheno_data.qmd)

```r
library(dplyr)

# 3.1 数据转换（兼容 Tidyverse）
# 使用 GWAS 标准的逆正态变换 (INT) 处理异常值
norm_data <- my_pheno_data %>%
  mutate(Yield_INT = transform_pheno(Yield, method = "int"))

# 3.2 可视化（直方图与云雨图联合绘制）
# 单个性状分布
plot_pheno_dist(norm_data, val_cols = "Yield_INT", fill_color = "#377EB8")

# 比较多个生物学重复
plot_pheno_dist(norm_data, val_cols = c("Rep1", "Rep2", "Rep3"))
```

### 4. GWAS 可视化（曼哈顿图）

高度灵活的曼哈顿图绘制工具，支持多性状叠加以及提升渲染速度的降采样功能。

```r
# 场景 A：单性状分析（使用针对 EMMAX 格式的辅助读取函数）
df_gwas <- read_emmax_result("emmax_result.ps")
plot_emmax_manhattan(df_gwas, trait_name = "Yield", threshold = 5)

# 场景 B：多性状叠加比较
gwas_list <- list("Trait1" = df_trait1, "Trait2" = df_trait2)
plot_emmax_manhattan(
  gwas_input = gwas_list,
  threshold = 6, 
  sample_rate = 0.5, # 对非显著 SNP 降采样 50%
  trait_colors = c("firebrick", "dodgerblue")
)
```

### 5. 群体结构可视化

用于可视化群体遗传结构的完整工具集，原生支持 PLINK 和 ADMIXTURE 的输出格式。
*注：MyToolkits 内置了多款学术级配色方案，如 "nature", "npg_classic", "earthy", "mp_cold" 及 "cyber" 等。*

#### 5.1 ADMIXTURE 谱系图

解析多个 `.Q` 文件，并自动在不同的 K 值层级间追踪和继承亚群颜色。

```r
# 读取 K=2 到 K=5 的 Admixture 结果文件
admix_data <- read_admixture(fam_file = "pop.fam", q_prefix = "pop_structure", k_range = 2:5)

# 基于 K=3 的核心亚群对样本进行排序分组
admix_ordered <- order_admix_accid(admix_data, ordered_by_k = 3)

# 绘制堆叠柱状图
plot_admixture(admix_ordered, colors = "nature")
```

#### 5.2 主成分分析 (PCA)

支持从 PLINK 输出结果生成碎石图、2D 散点图、多维矩阵图以及 3D 交互图。

```r
# 读取 PLINK PCA 结果并与分组元数据合并
pca_data <- read_plink_pca(eigenvec_file = "plink.eigenvec", eigenval_file = "plink.eigenval", group_info = sample_info)

# 绘制带有置信椭圆的 2D 散点图，程序会自动修正内部因子的颜色匹配
plot_pca_2d(pca_data, pc_x = 1, pc_y = 2, colors = "npg_classic", ellipse = TRUE)

# 其他可用的 PCA 绘图功能：
# plot_pca_scree(pca_data, n_pc = 10)
# plot_pca_pairs(pca_data, n_pc = 5)
# plot_pca_3d(pca_data, pc_x = 1, pc_y = 2, pc_z = 3)
```

#### 5.3 系统发育树

基于 `ggtree` 构建，支持导入元数据并实现高度定制化的系统树绘制（自动进行节点排序平滑处理）。

```r
# 读取 Newick 格式树文件并注入群体分组信息
tree_data <- read_phylo_tree(nwk_file = "population.nwk", group_info = sample_info, ladderize = TRUE)

# 绘制带有分支颜色的环形树
plot_phylo_tree(tree_data, layout = "circular", colors = "earthy", linewidth = 0.5)
```

## 依赖包

- **统计与模型**：`lme4`, `emmeans`, `stats`
- **数据处理**：`dplyr`, `tidyr`, `tibble`, `magrittr`, `purrr`, `stringr`
- **可视化**：`ggplot2`, `ggtree`, `gghalves`, `patchwork`, `plotly`
- **系统发育**：`ape`

## 许可证

MIT 许可证