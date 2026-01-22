[![Read in English](https://img.shields.io/badge/README-English-blue.svg)](./README.md)
# MyToolkits

MyToolkits 是一个专门为植物育种和数量遗传学数据分析设计的 R 包。它提供了从多环境试验计算育种值、执行全面的描述性统计，以及高度灵活地可视化 GWAS 结果的完整工作流程。

## 功能特点

- **多环境试验分析**：使用 lme4 计算 BLUP（最佳线性无偏预测）和 BLUE（最佳线性无偏估计）值。自动估计广义遗传力（$H^2$）。

- **科学描述性统计**：生成可直接用于发表的统计摘要，包括标准误（SE）以及采用 SAS/SPSS 标准的无偏偏度和峰度估计。

- **灵活的 GWAS 可视化**：生成兼容 EMMAX、GAPIT 和 TASSEL 输出的曼哈顿图。支持多性状可视化和大规模数据的自动降采样。

- **S4 面向对象设计**：使用稳健的 S4 类（BreedingValue, DescriptiveStat）确保数据完整性。

## 安装

您可以从 GitHub 安装开发版 MyToolkits：

```{r, eval=FALSE}
# 安装 devtools 包（如未安装）
# install.packages("devtools")
devtools::install_github("ShawnWx2019/MyToolkits")
```

## 使用示例

### 1. 育种值分析（BLUP & BLUE）

计算基因型在多年份和多地点的育种值。该函数自动处理基因型与环境互作效应。

```{r, eval=FALSE}
library(MyToolkits)

# --- 1. 加载多环境试验数据 ---
# 数据必须包含：Line, Year, Location, Trait, Value
data("my_pheno_data")

# --- 2. 计算 BLUPs（随机效应模型）---
# 拟合模型：Value ~ (1|Line) + (1|Year) + (1|Loc) + (1|Line:Year) + ...
blup_res <- get_unbiase_by_year_loc(
  x = my_pheno_data,
  Tag = "Yield",           # 要分析的具体性状水平
  method = "blup",         # "blup" 或 "blue"
  trait.name = "TraitID",  # 性状ID的列名
  year.name = "Year",      # 年份列名
  loc.name = "Loc",        # 地点列名
  line.name = "Genotype",  # 基因型/品系列名
  value.name = "Value"     # 表型值列名
)

# 查看 S4 对象摘要
print(blup_res)
# 输出包含：模型收敛状态和遗传力（H2）

# --- 3. 导出完整数据表 ---
# 将计算的 BLUPs 与原始数据合并（按环境横向展开）
final_table <- get_complete_table(blup_res)
head(final_table)
```

### 2. 描述性统计

计算详细的统计指标，使用无偏估计（特别适合小样本量）。

```{r, eval=FALSE}
# 创建数值向量
phenotypes <- c(10.5, 12.1, 9.8, 11.0, 10.2, 13.5, NA, 10.8)

# 计算统计量
stats_obj <- stat_descriptive(
  lab_name = "Plant_Height",
  x = phenotypes,
  na.omit = TRUE
)

# 访问结果表格
print(stats_obj@statistics)
# 返回：样本量、均值、标准误、标准差、变异系数、最小值/最大值、分位数、偏度、峰度
```

### 3. GWAS 可视化（曼哈顿图）

plot_emmax_manhattan 函数具有高度兼容性，可接受 EMMAX（无表头）或 GAPIT（有表头）等工具的原始输出。

#### 场景 A：单性状分析（EMMAX 格式）

```{r, eval=FALSE}
# EMMAX 输出通常没有表头
# 必须设置列名为："CHR", "BP", "P"
# 对于Emmax软件输出的.ps文件，我们提供了read_emmax_result函数直接进行读取和数据转换
df_gwas_res <- read_emmax_result("emmax_result.ps")
# 对于其他软件导入的结果：
df_gwas_res <- read.table("gwas_result.txt", header = FALSE)
colnames(df_gwas_res) <- c("CHR", "BP", "P")
# 使用默认阈值（-log10(P) = 4）绘图
plot_emmax_manhattan(df_emmax, trait_name = "Fiber_Length", threshold = 5)
```

#### 场景 B：多性状比较

可通过命名列表传递多个数据框，在同一图表中比较多个性状（通过颜色区分）。

```{r, eval=FALSE}
# 准备数据框列表
gwas_list <- list(
  "Petal_Color" = df_petal,  # 数据框必须包含 CHR, BP, P 列
  "Leaf_Shape"  = df_leaf
)

plot_emmax_manhattan(
  gwas_input = gwas_list,
  threshold = 6,
  sample_rate = 0.5,       # 对非显著SNP降采样50%以提高绘图速度
  trait_colors = c("firebrick", "dodgerblue")
)
```

## 依赖包

- lme4：用于混合线性模型
- emmeans：用于最小二乘均值（BLUE）
- ggplot2：用于可视化
- dplyr, tidyr, tibble, magrittr：用于数据整理

## 许可证

MIT 许可证
