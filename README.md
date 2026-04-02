[![Read in Chinese](https://img.shields.io/badge/README-%E4%B8%AD%E6%96%87-red.svg)](./README-CN.md)
# MyToolkits <img src="man/figures/logo.png" align="right" height="139" alt="" />

MyToolkits is an R package designed to streamline data analysis in plant breeding, quantitative genetics, and multi-omics studies. Built on robust S4 object-oriented classes, it provides a comprehensive workflow for multi-environment trials (BLUP/BLUE), scientific descriptive statistics, GWAS visualization, phenotypic data transformation, and population structure analysis.

## Features

- **Multi-Environment Trial Analysis**: Calculate BLUP and BLUE values using `lme4`. Automatically estimates Broad-sense Heritability ($H^2$).
- **Scientific Descriptive Statistics**: Generates publication-ready summary statistics including Standard Error (SE) and unbiased estimates for Skewness and Kurtosis (SAS/SPSS standard).
- **Phenotypic Data Processing**: Tidyverse-compatible functions to normalize traits (e.g., Inverse Normal Transformation, z-score) and generate raincloud distribution plots.
- **Flexible GWAS Visualization**: Generate Manhattan plots compatible with EMMAX, GAPIT, and TASSEL outputs, supporting multi-trait overlay and downsampling.
- **Population Structure Visualization**: Comprehensive tools built on S4 classes to visualize PLINK PCA, ADMIXTURE Q-matrices (with automatic lineage-tracking), and Phylogenetic trees using custom, colorblind-friendly palettes.

## Installation

You can install the development version of MyToolkits from GitHub. 
*We highly recommend setting `build_vignettes = TRUE` to access the detailed offline tutorials included in the package.*

```r
# install.packages("devtools")
devtools::install_github("ShawnWx2019/MyToolkits", build_vignettes = TRUE)
```

## Usage Examples

### 1. Breeding Value Analysis (BLUP & BLUE)

Calculate breeding values across multiple years and locations, handling genotype-by-environment (GxE) interactions.

```r
library(MyToolkits)

# Calculate BLUPs and extract Heritability
blup_res <- get_unbiase_by_year_loc(
  x = my_pheno_data,
  Tag = "Yield",           # The specific trait to analyze
  method = "blup",         # "blup" or "blue"
  trait.name = "TraitID",
  year.name = "Year",
  loc.name = "Loc",
  line.name = "Genotype",
  value.name = "Value"
)

print(blup_res) # Shows S4 object summary

# Export Complete Data Table
final_table <- get_complete_table(blup_res)
```

### 2. Descriptive Statistics

Calculate detailed statistics with unbiased estimators (ideal for small sample sizes).

```r
phenotypes <- c(10.5, 12.1, 9.8, 11.0, 10.2, 13.5, NA, 10.8)

stats_obj <- stat_descriptive(lab_name = "Plant_Height", x = phenotypes, na.omit = TRUE)

# Access the results table
get_stat_table(stats_obj)
```

### 3. Phenotypic Data Transformation & Distribution

Clean and visualize your phenotypic data before performing GWAS or downstream modeling. 

📖 **Tutorial Vignettes:** For a comprehensive guide on why and how to standardize phenotypic data, please read our vignette inside R via `vignette("description_pheno_data", package = "MyToolkits")`, or view the source files directly on GitHub:
- [Phenotypic Data Description (English)](./vignettes/description_pheno_data.qmd)
- [表型数据标准化说明 (中文版)](./vignettes/description_pheno_data_CN.qmd)

```r
library(dplyr)

# 3.1 Data Transformation (Tidyverse-compatible)
# Standard GWAS Inverse Normal Transformation (INT)
norm_data <- my_pheno_data %>%
  mutate(Yield_INT = transform_pheno(Yield, method = "int"))

# 3.2 Visualization (Histogram & Raincloud Plot)
# Single Trait
plot_pheno_dist(norm_data, val_cols = "Yield_INT", fill_color = "#377EB8")

# Compare Multiple Replicates
plot_pheno_dist(norm_data, val_cols = c("Rep1", "Rep2", "Rep3"))
```

### 4. GWAS Visualization (Manhattan Plots)

Highly flexible Manhattan plotting tool supporting multiple traits and downsampling for speed.

```r
# Scenario A: Single Trait (EMMAX format)
df_gwas <- read_emmax_result("emmax_result.ps")
plot_emmax_manhattan(df_gwas, trait_name = "Yield", threshold = 5)

# Scenario B: Multi-Trait Overlay
gwas_list <- list("Trait1" = df_trait1, "Trait2" = df_trait2)
plot_emmax_manhattan(
  gwas_input = gwas_list,
  threshold = 6, 
  sample_rate = 0.5, # Downsample non-significant SNPs by 50%
  trait_colors = c("firebrick", "dodgerblue")
)
```

### 5. Population Structure Visualization

Tools to visualize genetic structure, natively supporting PLINK and ADMIXTURE formats. *Note: MyToolkits comes with built-in academic palettes like "nature", "npg_classic", "earthy", "mp_cold", and "cyber".*

#### 5.1 ADMIXTURE Lineage

Parses multiple `.Q` files and automatically tracks color lineages across different K values.

```r
# Read Admixture results across K=2 to K=5
admix_data <- read_admixture(fam_file = "pop.fam", q_prefix = "pop_structure", k_range = 2:5)

# Order samples based on sub-populations at K=3
admix_ordered <- order_admix_accid(admix_data, ordered_by_k = 3)

# Plot stacked bar chart
plot_admixture(admix_ordered, colors = "nature")
```

#### 5.2 Principal Component Analysis (PCA)

Generate Scree, 2D, Pairs, and 3D PCA plots from PLINK outputs.

```r
# Read PLINK PCA results and merge with grouping metadata
pca_data <- read_plink_pca(eigenvec_file = "plink.eigenvec", eigenval_file = "plink.eigenval", group_info = sample_info)

# 2D Scatter Plot with confidence ellipses
plot_pca_2d(pca_data, pc_x = 1, pc_y = 2, colors = "npg_classic", ellipse = TRUE)

# Other available PCA plots:
# plot_pca_scree(pca_data, n_pc = 10)
# plot_pca_pairs(pca_data, n_pc = 5)
# plot_pca_3d(pca_data, pc_x = 1, pc_y = 2, pc_z = 3)
```

#### 5.3 Phylogenetic Tree

Customizable tree plotting using `ggtree` with auto-ladderization.

```r
# Read Newick tree and inject group information
tree_data <- read_phylo_tree(nwk_file = "population.nwk", group_info = sample_info, ladderize = TRUE)

# Plot circular tree
plot_phylo_tree(tree_data, layout = "circular", colors = "earthy", linewidth = 0.5)
```

## Dependencies

- **Statistics & Models:** `lme4`, `emmeans`, `stats`
- **Data Manipulation:** `dplyr`, `tidyr`, `tibble`, `magrittr`, `purrr`, `stringr`
- **Visualization:** `ggplot2`, `ggtree`, `gghalves`, `patchwork`, `plotly`
- **Phylogenetics:** `ape`

## License

MIT