# MyToolkits 0.1.0

## 新增功能 (New Features)
* **表型数据处理**：新增 `transform_pheno()` 函数，支持逆正态变换等标准化方法；新增 `plot_pheno_dist()` 用于绘制云雨分布图。
* **群体结构可视化**：引入了全新的 `AdmixData` S4 类，新增 `read_admixture()` 和 `plot_admixture()` 函数，支持 K 值谱系的自动颜色追踪。
* **主成分分析 (PCA)**：新增对 PLINK 结果的读取和多种可视化支持，包括 `plot_pca_2d()`, `plot_pca_scree()`, `plot_pca_pairs()` 和 `plot_pca_3d()`。
* **系统发育树**：新增 `plot_phylo_tree()` 函数，支持结合元数据绘制高定制化的进化树。

## 改进 (Improvements)
* 优化了 `plot_emmax_manhattan()` 的运行逻辑。