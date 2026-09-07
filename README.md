# scexpr

`scexpr` is an R package for exploring and analyzing single-cell RNA-sequencing data stored in Seurat objects. It collects practical helpers for visualization, quality control, clustering, marker analysis, gene-set enrichment, label transfer, and related workflows.

The plotting functions are designed for exploratory analysis: they accept familiar Seurat inputs, provide useful defaults, and return standard `ggplot2` or `patchwork` objects that can be customized or saved with the usual R plotting tools.

## Installation

Install the development version from GitHub:

```r
install.packages("pak")
pak::pak("Close-your-eyes/scexpr")
```

To install a local checkout instead, run this from the directory containing the package:

```r
pak::pak("./scexpr")
```

Some functions use optional packages. `scexpr` checks for the packages needed by a function when it runs and attempts to install any that are missing. The first call may therefore take longer and may require internet access.

## Quick start

The package includes a small PBMC Seurat object that can be used to try the examples:

```r
library(scexpr)

pbmc <- readRDS(system.file(
  "extdata",
  "SO_5k_pbmc_v3_RNA_none_1_800_12_small.rds",
  package = "scexpr"
))
```

## Plot features on an embedding with `feature_plot2()`

`feature_plot2()` displays gene expression or Seurat metadata on a two-dimensional reduction such as UMAP, t-SNE, or PCA.

```r
feature_plot2(
  SO = pbmc,
  features = c("MS4A1", "CD3E", "CD68"),
  reduction = "umap",
  assay = "RNA",
  ncol_combine = 3
)
```

For several features, the plots are combined into a `patchwork` layout by default. Use `combine = FALSE` to receive a list of individual `ggplot2` objects.

```r
plots <- feature_plot2(
  SO = pbmc,
  features = c("MS4A1", "CD3E"),
  combine = FALSE
)

plots[[1]] + ggplot2::labs(subtitle = "B-cell marker")
```

### Why use it?

- **One interface for genes and metadata.** Plot expression values, cluster assignments, sample labels, or QC columns without switching functions.
- **Useful multi-panel behavior.** Plot several features at once, control the layout, and optionally collect guides.
- **Clear handling of sparse expression.** Non-expressing cells remain visible, expressing cells can be emphasized, and `freq_plot` can annotate the fraction of expressing cells.
- **Flexible color scales.** Use continuous, stepped, quantile-based, or binary expression scales. For example, `col_steps = "quartiles0"` derives quartiles from values above zero.
- **Built-in biological context.** Split panels by metadata, label clusters, add density contours, highlight selected cells, or filter cells by another feature.
- **Works across objects.** A named list of Seurat objects can be plotted through the same call, which is useful for comparing samples or processing variants.

Here is a more annotated example:

```r
feature_plot2(
  SO = pbmc,
  features = "MS4A1",
  reduction = "umap",
  label_feature = "pca12_rna800_snn_res_0.1",
  contour_feature = "pca12_rna800_snn_res_0.1",
  col_steps = "quartiles0",
  pt_size = 0.4,
  pt_size_fct = 1.5,
  title = "B-cell marker expression"
)
```

Common arguments include:

| Argument | Purpose |
|---|---|
| `features` | Gene names or metadata columns to plot |
| `reduction`, `dims` | Reduction and dimensions to display |
| `cells` | Cells to emphasize; other cells remain as background |
| `split_feature` | Metadata column used to create facets |
| `label_feature` | Metadata column used for labels |
| `contour_feature` | Metadata column used for density contours |
| `col_binary` | Show detected versus non-detected expression |
| `col_steps` | Control continuous or binned color scales |
| `combine`, `ncol_combine`, `nrow_combine` | Control multi-feature output |

## Compare groups with `feature_plot_stat()`

An embedding shows *where* expression occurs, but it can be difficult to judge group-level differences from color alone. `feature_plot_stat()` complements `feature_plot2()` by plotting expression distributions across a metadata grouping variable.

```r
feature_plot_stat(
  SO = pbmc,
  features = c("MS4A1", "CD3E"),
  meta_col = "pca12_rna800_snn_res_0.1",
  geom1 = "jitter",
  geom2 = "violin",
  expr_freq_pct = TRUE
)
```

### Why use it?

- **Distribution and detection in one view.** Combine cell-level points with violin or box plots and annotate the percentage of cells with expression above zero.
- **Faceted comparison of several features.** Each feature gets its own panel, with free y-scales by default.
- **Control over zeros.** Keep non-expressing cells to show dropout and detection frequency, or set `plot_non_expr = FALSE` to focus on expressing cells.
- **Flexible grouping and coloring.** Group cells with `meta_col` and, when useful, color them with a different metadata column via `color_by`.
- **Optional statistical annotations.** Set `add_pwc = TRUE` to add pairwise comparisons through `ggpubr`; customize them with `pwc_args`.

For a compact box-plot summary:

```r
feature_plot_stat(
  SO = pbmc,
  features = c("MS4A1", "CD3E", "CD68"),
  meta_col = "pca12_rna800_snn_res_0.1",
  geom1 = "none",
  geom2 = "boxplot",
  expr_freq_pct = TRUE,
  pt_size = 0.3
)
```

Both functions also accept a list of Seurat objects. Naming the list makes the object-of-origin legend easier to interpret:

```r
pbmc_8k <- readRDS(system.file(
  "extdata",
  "SO_pbmc8k_RNA_none_1_900_14_small.rds",
  package = "scexpr"
))

feature_plot2(
  SO = list(pbmc_5k = pbmc, pbmc_8k = pbmc_8k),
  features = "MS4A1",
  reduction = "umap"
)
```

## Getting help

Use the function documentation for the full set of options:

```r
?feature_plot2
?feature_plot_stat
vignette("feature_plot2", package = "scexpr")
```

## License

`scexpr` is licensed under GPL (>= 3).
