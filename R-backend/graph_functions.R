

#=========================================================
# FUNCTION: plot_deg_heatmap()
#=========================================================
plot_deg_pheatmap <- function(expr_mat,
                              group_vector,
                              deg_results,
                              up_n = 25,
                              down_n = 25,
                              scale_rows = TRUE,
                              main_title = "DEGs Heatmap") {
  
  # -----------------------------
  # 1. Identify top up- and down-regulated genes
  # -----------------------------
  deg_results$Gene <- rownames(deg_results)
  
  up_degs <- deg_results[deg_results$logFC > 0, ]
  down_degs <- deg_results[deg_results$logFC < 0, ]
  
  top_up_genes <- head(up_degs[order(up_degs$adj.P.Val), ], up_n)$Gene
  top_down_genes <- head(down_degs[order(down_degs$adj.P.Val), ], down_n)$Gene
  
  top_genes <- c(top_up_genes, top_down_genes)
  
  # -----------------------------
  # 2. Subset expression matrix
  # -----------------------------
  expr_subset <- expr_mat[top_genes, intersect(colnames(expr_mat), names(group_vector)), drop = FALSE]
  
  # Set matrix name to control the colorbar legend
  attr(expr_subset, "name") <- "Expression"
  
  # -----------------------------
  # 3. Scale rows (optional)
  # -----------------------------
  if(scale_rows) {
    expr_subset <- t(scale(t(expr_subset)))
  }
  
  # -----------------------------
  # 4. Column annotation
  # -----------------------------
  annotation_col <- data.frame(Group = factor(group_vector[colnames(expr_subset)]))
  rownames(annotation_col) <- colnames(expr_subset)
  
  annotation_colors <- list(
    Group = c(NC = "#60C9CC", RA = "#FA8787")
  )
  
  # -----------------------------
  # 5. Color palette
  # -----------------------------
  heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
  
  # -----------------------------
  # 6. Plot pheatmap
  # -----------------------------
  pheatmap(
    expr_subset,
    color = heatmap_colors,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    fontsize_row = 5,
    scale = "none",
    border_color = "white",
    main = main_title
  )
}

#==============================
# FUNCTION: plot_gsva_heatmap()
#==============================
plot_gsva_heatmap <- function(gsva_scores,
                              group_vector,
                              top_n = 20,
                              scale_rows = TRUE,
                              main_title = "GSVA Pathway Heatmap") {
  
  # -----------------------------
  # 1️⃣ Differential analysis using limma
  # -----------------------------
  group_vector <- factor(group_vector)
  design <- model.matrix(~0 + group_vector)
  colnames(design) <- levels(group_vector)
  contrast <- makeContrasts(RA - NC, levels = design)
  
  fit <- lmFit(gsva_scores, design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Select top pathways by absolute t-stat
  top_pathways <- topTable(fit2, number = nrow(gsva_scores), sort.by = "t")
  top_pathways <- top_pathways[order(-abs(top_pathways$t)), ]
  top_pathway_names <- rownames(top_pathways)[1:top_n]
  
  # -----------------------------
  # 2️⃣ Subset GSVA scores
  # -----------------------------
  expr_subset <- gsva_scores[top_pathway_names, names(group_vector), drop = FALSE]
  
  # -----------------------------
  # 3️⃣ Scale rows if needed
  # -----------------------------
  if (scale_rows) {
    expr_subset <- t(scale(t(expr_subset)))
  }
  
  # -----------------------------
  # 4️⃣ Column annotation
  # -----------------------------
  annotation_col <- data.frame(Group = factor(group_vector[colnames(expr_subset)]))
  rownames(annotation_col) <- colnames(expr_subset)
  
  annotation_colors <- list(
    Group = c(NC = "#006400", RA = "#66C266")  # NC dark green, RA light green
  )
  
  # -----------------------------
  # 5️⃣ Color palette: blue-white-red
  # -----------------------------
  heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # -----------------------------
  # 6️⃣ Plot heatmap
  # -----------------------------
  pheatmap(
    expr_subset,
    color = heatmap_colors,
    cluster_rows = TRUE,       # dendrogram on rows only
    cluster_cols = FALSE,      # no clustering on columns
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    fontsize_row = 6,
    scale = "none",
    border_color = "white",
    main = main_title
  )
}




#=========================================================
# FUNCTION: plot_volcano()
#=========================================================
plot_volcano <- function(deg_results,
                         top_n = 10,
                         logFC_cutoff = 1,
                         pval_cutoff = 0.05,
                         main_title = "Volcano Plot") {
  
  # Ensure gene names column exists
  deg_results$Gene <- rownames(deg_results)
  
  # Categorize genes
  deg_results$Category <- "Unchanged"
  deg_results$Category[deg_results$logFC > logFC_cutoff & deg_results$adj.P.Val < pval_cutoff] <- "Up"
  deg_results$Category[deg_results$logFC < -logFC_cutoff & deg_results$adj.P.Val < pval_cutoff] <- "Down"
  
  # Identify top N up and down genes
  up_genes <- head(deg_results[deg_results$Category == "Up", ][order(deg_results$adj.P.Val[deg_results$Category=="Up"]), ], top_n)$Gene
  down_genes <- head(deg_results[deg_results$Category == "Down", ][order(deg_results$adj.P.Val[deg_results$Category=="Down"]), ], top_n)$Gene
  
  deg_results$Label <- ""  # no label for others
  
  # Subset only top genes for labeling
  label_df <- deg_results[deg_results$Gene %in% c(up_genes, down_genes), ]
  label_df$Label <- label_df$Gene  # labels only for top genes
  
  # Plot
  ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Category)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Up" = "blue", "Down" = "orange", "Unchanged" = "grey70")) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black") +
    geom_label_repel(
      data = label_df,   # only top genes
      aes(x = logFC, y = -log10(adj.P.Val), label = Label),
      box.padding = 0.3,
      segment.color = "grey50",
      size = 3,
      fill = "white",
      max.overlaps = 20,
      show.legend = FALSE   # ✅ THIS LINE
    ) +
    theme_minimal(base_size = 14) +
    labs(title = main_title, x = "logFC", y = "-log10(P-value)") +
    theme(legend.title = element_blank())
}


#=========================================================
# FUNCTION: plot_soft_threshold()
#=========================================================
plot_soft_threshold <- function(sft_data, r2_cut = 0.8) {
  
  a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.05) +
    geom_hline(yintercept = r2_cut, color = "red") +
    labs(
      x = "",
      y = "Scale Free Topology Model Fit (R²)"
    )
  
  a2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.05) +
    labs(
      x = "Soft Threshold (power)",
      y = "Mean connectivity"
    ) +
    theme_classic()
  
  grid.arrange(a1, a2, nrow = 2)
}

#=========================================================
# FUNCTION: plot_gene_dendrogram()
#=========================================================
plot_gene_dendrogram <- function(datExprMAD, softPower, dissTOM) {

  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  plot(
    geneTree,
    labels = FALSE,
    hang = 0.03,
    main = "Gene clustering of top 25% MAD genes"
  )
  
  invisible(geneTree)
}

#=========================================================
# FUNCTION: plot_module_dendrogram()
#=========================================================
plot_module_dendrogram <- function(
    geneTree,
    dynamicColors,
    mergedColors
) {
  plotDendroAndColors(
    geneTree,
    cbind(dynamicColors, mergedColors),
    c("Dynamic Tree Cut", "Merged modules"),
    dendroLabels = FALSE,
    hang = 0.03,
    main = "Gene dendogram and module colors",
    addGuide = TRUE,
    guideHang = 0.05
  )
}

#=========================================================
# FUNCTION: plot_module_trait_heatmap()
#=========================================================
plot_module_trait_heatmap <- function(
    moduleTraitCor,
    moduleTraitPvalue,
    mergedMEs,
    trait_name = colnames(moduleTraitCor),
    width = 6,
    height = 8,
    sigThreshold = 0.05
) {
  
  corRounded <- formatC(moduleTraitCor, format = "f", digits = 2)
  
  textMatrix <- matrix(nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
  for (i in 1:nrow(moduleTraitCor)) {
    for (j in 1:ncol(moduleTraitCor)) {
      if (moduleTraitPvalue[i, j] < sigThreshold) {
        textMatrix[i, j] <- paste0(corRounded[i, j], "*")
      } else {
        textMatrix[i, j] <- corRounded[i, j]
      }
    }
  }
  
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = trait_name,
    yLabels = names(mergedMEs),
    ySymbols = names(mergedMEs),
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    main = "Module–trait relationships",
    cex.text = 0.7,
    zlim = c(-1, 1)
  )
}

#=========================================================
# FUNCTION: plot_GS_MM()
#=========================================================
plot_GS_MM <- function(
    geneInfo_mod,
    sigModule,
    sigColor
) {
  plot(
    geneInfo_mod$MM,
    geneInfo_mod$GS,
    xlab = paste("Module Membership in", sigModule, "module"),
    ylab = "Gene Significance for proliferating",
    main = paste("Module Membership vs Gene Significance"),
    col = sigColor,
    pch = 16
  )
  
  abline(
    lm(GS ~ MM, data = geneInfo_mod),
    col = "black",
    lwd = 2
  )
  
  cor_val <- cor(geneInfo_mod$MM, geneInfo_mod$GS, use = "p")
  
  legend(
    "topleft",
    legend = paste("Pearson r =", round(cor_val, 2)),
    bty = "n"
  )
}

#=========================================================
# FUNCTION: plot_deg_wgcna_venn()
#=========================================================
plot_deg_wgcna_venn <- function(
    deg_genes,
    wgcna_hub_genes,
    labels = c("DEGs", "WGCNA"),
    fill_colors = c("red", "blue"),
    alpha = 0.5,
    cex = 1.5,
    cat_cex = 2.2
) {
  
  venn.plot <- venn.diagram(
    x = list(
      DEGs = deg_genes,
      WGCNA = wgcna_hub_genes
    ),
    category.names = labels,
    filename = NULL,
    fill = fill_colors,
    alpha = alpha,
    cex = cex,
    cat.cex = cat_cex
  )
  
  grid.newpage()
  grid.draw(venn.plot)
  
  invisible(venn.plot)
}

#=========================================================
# FUNCTION: plot_ml_venn()
#=========================================================
plot_ml_venn <- function(lasso_genes, svm_genes, rf_genes) {
  
  venn.plot <- venn.diagram(
    x = list(
      LASSO = lasso_genes,
      SVM_RFE = svm_genes,
      RF = rf_genes
    ),
    filename = NULL,
    fill = c("red", "blue", "green"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.2
  )
  
  grid.newpage()
  grid.draw(venn.plot)
}
