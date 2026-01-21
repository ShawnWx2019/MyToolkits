
# MyToolkits

MyToolkits is an R package designed to streamline data analysis in plant breeding and quantitative genetics. It provides a robust workflow for calculating breeding values (BLUP/BLUE) from multi-environment trials (MET), performing comprehensive descriptive statistics, and visualizing GWAS results with high-flexibility Manhattan plots.

## Features

- **Multi-Environment Trial Analysis**: Calculate BLUP (Best Linear Unbiased Prediction) and BLUE (Best Linear Unbiased Estimation) values using lme4. Automatically estimates Broad-sense Heritability ($H^2$).

- **Scientific Descriptive Statistics**: Generates publication-ready summary statistics including Standard Error (SE) and unbiased estimates for Skewness and Kurtosis (SAS/SPSS standard).

- **Flexible GWAS Visualization**: Generate Manhattan plots compatible with EMMAX, GAPIT, and TASSEL outputs. Supports multi-trait visualization and automatic downsampling for large datasets.

- **S4 Object Oriented**: Uses robust S4 classes (BreedingValue, DescriptiveStat) to ensure data integrity.

## Installation

You can install the development version of MyToolkits from GitHub:

```{r, eval=FALSE}
# install.packages("devtools")
devtools::install_github("ShawnWx2019/MyToolkits")
```

## Usage Examples

### 1. Breeding Value Analysis (BLUP & BLUE)

Calculate breeding values for genotypes across multiple years and locations. The function handles genotype-by-environment (GxE) interactions automatically.

```{r, eval=FALSE}
library(MyToolkits)

# --- 1. Load your Multi-Environment Trial data ---
# Data must contain: Line, Year, Location, Trait, Value
data("my_pheno_data") 

# --- 2. Calculate BLUPs (Random Effects Model) ---
# This fits: Value ~ (1|Line) + (1|Year) + (1|Loc) + (1|Line:Year) + ...
blup_res <- get_unbiase_by_year_loc(
  x = my_pheno_data,
  Tag = "Yield",           # The specific trait level to analyze
  method = "blup",         # "blup" or "blue"
  trait.name = "TraitID",  # Column name for trait IDs
  year.name = "Year",      # Column name for Year
  loc.name = "Loc",        # Column name for Location
  line.name = "Genotype",  # Column name for Line/Genotype
  value.name = "Value"     # Column name for phenotype values
)

# View the S4 Object Summary
print(blup_res)
# Output includes: Model convergence status and Heritability (H2)

# --- 3. Export Complete Data Table ---
# Merges calculated BLUPs with raw data (pivoted wide by environment)
final_table <- get_complete_table(blup_res)
head(final_table)
```

### 2. Descriptive Statistics

Calculate detailed statistics with unbiased estimators (ideal for small sample sizes).

```{r, eval=FALSE}
# Create a numeric vector
phenotypes <- c(10.5, 12.1, 9.8, 11.0, 10.2, 13.5, NA, 10.8)

# Calculate statistics
stats_obj <- stat_descriptive(
  lab_name = "Plant_Height", 
  x = phenotypes, 
  na.omit = TRUE
)

# Access the results table
print(stats_obj@statistics)
# Returns: N, Mean, SE, SD, CV, Min/Max, Quantiles, Skewness, Kurtosis
```

### 3. GWAS Visualization (Manhattan Plots)

The plot_emmax_manhattan function is designed to be highly compatible. It accepts raw outputs from EMMAX (no headers) or standard tools like GAPIT (with headers).

#### Scenario A: Single Trait (EMMAX style)

```{r, eval=FALSE}
# EMMAX usually has no header. 
# You MUST set column names to: "CHR", "BP", "P"
df_emmax <- read.table("emmax_output.ps", header = FALSE)
colnames(df_emmax) <- c("CHR", "BP", "P")

# Plot with default threshold (-log10(P) = 4)
plot_emmax_manhattan(df_emmax, trait_name = "Fiber_Length", threshold = 5)
```

#### Scenario B: Multi-Trait Comparison

You can pass a named list of dataframes to plot multiple traits on the same track (distinguished by color).

```{r, eval=FALSE}
# Prepare list of dataframes
gwas_list <- list(
  "Petal_Color" = df_petal,  # df must have CHR, BP, P
  "Leaf_Shape"  = df_leaf
)

plot_emmax_manhattan(
  gwas_input = gwas_list,
  threshold = 6, 
  sample_rate = 0.5,       # Downsample non-significant SNPs to 50% for speed
  trait_colors = c("firebrick", "dodgerblue")
)
```

## Dependencies

- lme4: For mixed linear models.
- emmeans: For Least-Squares Means (BLUE).
- ggplot2: For visualization.
- dplyr, tidyr, tibble, magrittr: For data manipulation.

## License

MIT
