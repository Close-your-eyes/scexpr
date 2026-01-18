#df <- SeuratData::AvailableData()
so2pt <- subset(so2, celltype == "PT")

Seurat::Idents(so2pt) <- so2pt$orig.ident2
# Subset to cells of interest
cells_A <- Seurat::WhichCells(so2pt, idents = "DTX3")
cells_B <- Seurat::WhichCells(so2pt, idents = "Muto")
cells_use <- c(cells_A, cells_B)
so2pt <- subset(so2pt, cells = cells_use)
so2pt@meta.data$idents <- as.factor(Seurat::Idents(so2pt))


# create single cell exp
sca <- MAST::FromMatrix(
  exprsArray = scexpr::get_layer(so2pt, as = "dense"),
  cData = so2pt@meta.data,
  fData = data.frame(primerid = rownames(scexpr::get_layer(so2pt)))
)
options("mc.cores")
# fit hurdle model
options("mc.cores" = 8)
zlm_fit <- MAST::zlm(
  ~ idents + nFeature_RNA,
  sca
)
# Likelihood Ratio Test
lrt <- MAST::lrTest(zlm_fit, "condition")
# extract results
lrt_table <- lrt[, , "Pr(>Chisq)"]
pvals <- lrt_table[, "hurdle"]


coef_table <- summary(zlm_fit)$coef
logFC <- coef_table["conditionB", , "logFC"]

results <- data.frame(
  gene = rownames(expr),
  avg_logFC = logFC,
  p_val = pvals,
  p_val_adj = p.adjust(pvals, method = "BH")
)

results <- results[order(results$p_val), ]

