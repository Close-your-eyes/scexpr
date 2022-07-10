.onLoad <- function(libname = find.package("scexpr"), pkgname = "scexpr") {

  # https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when

  # CRAN Note avoidance
  if(getRversion() >= "2.15.1")
    utils::globalVariables(

      #.calc_vd
      c("an",
        #expr_jitter
        "Feature", "expr", "group", "max.feat.expr", "pct.expr", "log2.fc",
        #.plot_vp
        "xpos", "ypos", "n.genes",
        #cluster_correlation_matrix
        "Var1", "Var2", "value",
        #convert_gene_identifier
        "SYMBOL", "ENTREZID",
        #feature_plot
        "feature", "xmin", "xmax", "ymin", "ymax",
        #
        "avgExpr", "pct_in", "padj", "group", "cluster", "norm_avgexpr", "y", "GENENAME",
        #hla_typing
        "unique_explained_reads", "double_explained_reads", "total_explained_reads", "total_explained_reads_rank", "unique_explained_reads_rank",
        "allele.1.2", "combined.rank", "explained_reads_allele.1", "explained_reads_allele.2", "p_group.1", "p_group.2",
        "allele.1", "allele.2", "allele_group.1", "allele_group.2", "p_group.1.2", "row.number.rank",
        "allele_group", "plot.color", "g_group.1", "g_group.2",
        #qc_params_meta_cols
        "level", "value", "qc_param", "meta.col",
        #misc
        "ID", "abs_diff", "count", "data", "dbl_score_log", "end_angle", "frac", "meta_UMAP_1",
        "meta_UMAP_2", "mid_angle", "nCount_RNA", "nCount_RNA_log", "nFeature_RNA",
        "nFeature_RNA_log", "n_expr_corrected", "n_expr_uncorrected", "pct_mt",
        "pct_mt_log", "r", "ref_feature", "residuals", "start_angle", "stat", "x",
        "Alt.symbol", "Alt.symbol.HS", "Gene.HS", "Gene.MM", "Gene.stable.ID.HS",
        "Gene.stable.ID.MM",
        #qc_diagnostics
        "logFC", "statistic", "pval", "pct_out")


    )
}
