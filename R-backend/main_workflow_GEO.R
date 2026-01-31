
#=========================================================
# Install all required packages
#=========================================================
#install.packages("R.utils")
#install.packages("gridExtra")
#install.packages("VennDiagram")
#install.packages("glmnet")
#install.packages("caret")
#install.packages("randomForest")
#install.packages("rms")
#install.packages("rmda")
#install.packages("Matrix")
#install.pckages("cowplot")
#install.packages(c("ggplot2", "pheatmap", "openxlsx"))
#install.packages("hdf5r")
#install.packages("WGCNA")
#install.packages("flashClust")
#install.packages("tcltk")
#install.packages("BiocManager")

#BiocManager::install(c(
#  "GEOquery", "Biobase", "affy", "affyio", "AnnotationDbi",
#  "sva", "limma", "annotate", "hgu133a.db", "hgu133plus2.db",
#  "arrayQualityMetrics", "rhdf5", "GO.db", "impute",
#  "oligo", "pd.hugene.1.0.st.v1", "hugene10sttranscriptcluster.db",
#  "clusterProfiler","enrichplot","DOSE","GSVA"
#), ask = FALSE, update = FALSE)

# Load necessary libraries
library(openxlsx)
library(R.utils)
library(tcltk)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(sva)
library(dplyr)
library(oligo)
library(GEOquery)
library(affy)
library(affyio)
library(hgu133a.db)
library(hgu133plus2.db)
library(AnnotationDbi)
library(grid)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(DOSE)
library(WGCNA)
library(flashClust)
library(VennDiagram)
library(glmnet)
library(caret)
library(e1071)
library(randomForest)
library(rms)
library(pROC)
library(rmda)
library(cowplot)
library(patchwork)


#import custom-made script with functions
source("D:/Master Thesis Bioinformatics/Code/GEO_lib_functions.R")
source("D:/Master Thesis Bioinformatics/Code/graph_functions.R")

#------------------1.DOWNLOAD DATA AND METADATA----------------------#

#download .CEL and series_matrix files
gse_ids <- c("GSE77298", "GSE55235","GSE12021")
downloads <- download_GEO(gse_ids, type = c("series_matrix", "CEL"))


#Get the metadata of the files and find groups of disease-non disease
group_vector <- get_geo_groups()

#Validation GSE datasets
validation_sets <- c("GSE55457")

#------------------2.DATA PREPROCESSING----------------------#

#create xlsx files from .CEL and RMA normalize them
cel_matrices <- gene_mat_preprocess(input_type = "CEL")
#create xlsx files from series_matrix.txt
#series_matrices <- gene_mat_preprocess(input_type = "series_matrix")
combat_final_rma_mat <- merge_and_combat(group_vector)

rowMeans(combat_final_rma_mat)  # are there very low-expressed genes?
apply(combat_final_rma_mat, 1, sd)  # check variance

# Keep only genes with mean >= 5
gene_means <- rowMeans(combat_final_rma_mat)
combat_final_rma_mat <- combat_final_rma_mat[gene_means >= 5, ]

cat("Number of genes after filtering:", nrow(combat_final_rma_mat), "\n")

#call the DEG analysis function
limma_out <- run_limma_deg(
  expr_mat = combat_final_rma_mat,
  group_vector = group_vector,
  group_case = "RA",
  group_control = "NC",
  logfc_cutoff = 1,
  p_cutoff = 0.013,
  adjust_method = "BH"
)


deg_results <- limma_out$deg_results
degs <- limma_out$degs


plot_deg_pheatmap(
  expr_mat = combat_final_rma_mat,
  group_vector = group_vector,
  deg_results = limma_out$deg_results,
  up_n = 25,
  down_n = 25,
  scale_rows = TRUE,
  main_title = "Top DEGs Heatmap"
)

plot_volcano(
  deg_results = limma_out$deg_results,
  top_n = 10,
  logFC_cutoff = 1,
  pval_cutoff = 0.0017,
  main_title = "Top DEGs Volcano Plot"
)


#=========================================================
# GO enrichment
#=========================================================
# Assume `deg_genes` is a vector of significant gene symbols
deg_genes <- rownames(limma_out$degs)

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#Run GO enrichment (ALL)
ego <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",       # BP + CC + MF
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Convert enrichGO object to data.frame
ego_df <- as.data.frame(ego)

# Keep top 5 terms PER ontology (as in paper)
ego_df <- ego_df %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup()

# Convert GeneRatio from "x/y" → numeric
ego_df$GeneRatio <- sapply(ego_df$GeneRatio, function(x) {
  eval(parse(text = x))
})

# Order GO terms within each ontology
ego_df$Description <- factor(
  ego_df$Description,
  levels = rev(unique(ego_df$Description))
)

#-----------------------------
# 5. Bubble plot (Figure 3A)
#-----------------------------
ggplot(
  ego_df,
  aes(
    x     = Count,
    y     = Description,
    size  = Count,
    color = -log10(p.adjust)
  )
) +
  geom_point(alpha = 0.9) +
  facet_grid(
    ONTOLOGY ~ .,
    scales = "free_y",
    space  = "free_y"
  ) +
  scale_color_gradient(low = "pink", high = "red") +
  theme_bw(base_size = 14) +
  labs(
    title = "GO Enrichment",
    x     = "Gene Number",
    y     = "Gene Ratio",
    color = "-log10(adj.P)",
    size  = "Count"
  )




#=========================================================
# KEGG enrichment
#=========================================================
ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = 'hsa',     # human
  pvalueCutoff = 0.05
)

# Convert Entrez IDs back to gene symbols
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

#convert to dataframe
ekegg_df <- as.data.frame(ekegg)

# Keep top 10 by adjusted p-value
ekegg_df <- ekegg_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)

# Compute enrichment factor (GeneRatio as numeric)
#ekegg_df$RichFactor <- sapply(ekegg_df$GeneRatio, function(x) eval(parse(text = x)))

# Order pathways for plotting (largest on top)
ekegg_df$Description <- factor(ekegg_df$Description, levels = rev(ekegg_df$Description))

ggplot(
  ekegg_df,
  aes(
    x     = RichFactor,        # enrichment factor
    y     = Description,
    size  = Count,             # number of DEGs in pathway
    color = -log10(p.adjust)  # significance
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "pink", high = "red") +
  theme_bw(base_size = 14) +
  labs(
    title = "KEGG Enrichment",
    x     = "Enrichment Factor",
    y     = "KEGG term",
    color = "-log10(adj.P)",
    size  = "Count"
  )


#=========================================================
# KEGG GSEA – multi-pathway composite enrichment plot
#=========================================================

#---------------------------------------------------------
# 1️⃣ Prepare ranked gene list (ALL genes)
#---------------------------------------------------------
deg_results <- limma_out$deg_results
deg_results$SYMBOL <- rownames(deg_results)

# Map gene symbols to Entrez IDs
gene_map <- bitr(deg_results$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#remove duplicates found in db  
gene_map <- gene_map[!duplicated(gene_map$SYMBOL), ]

deg_results_mapped <- merge(deg_results, gene_map, by = "SYMBOL")

gene_list <- deg_results_mapped$logFC
names(gene_list) <- as.character(deg_results_mapped$ENTREZID)
gene_list <- sort(gene_list, decreasing = TRUE)

#---------------------------------------------------------
# 2️⃣ KEGG gene sets (MSigDB C2)
#---------------------------------------------------------
kegg_sets <- msigdbr(
  species  = "Homo sapiens",
  category = "C2"
) %>%
  filter(grepl("^KEGG_", gs_name))

term2gene <- kegg_sets %>%
  select(gs_name, entrez_gene) %>%
  distinct() %>%
  rename(geneSet = gs_name, gene = entrez_gene)

#---------------------------------------------------------
# 3️⃣ Run GSEA
#---------------------------------------------------------
gsea_res <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

#---------------------------------------------------------
# 4️⃣ Select top 5 KEGG pathways
#---------------------------------------------------------
gsea_tbl <- as.data.frame(gsea_res) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)

gsea_tbl$Description <- gsea_tbl$Description %>%
  gsub("^KEGG_", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

top_ids    <- gsea_tbl$ID
top_labels <- gsea_tbl$Description
n_sets     <- length(top_ids)

#---------------------------------------------------------
# 5️⃣ Running enrichment scores
#---------------------------------------------------------
get_running_es <- function(gsea_obj, geneSetID) {
  geneset <- gsea_obj@geneSets[[geneSetID]]
  metric  <- gsea_obj@geneList
  ranks   <- names(metric)
  
  hits <- ranks %in% geneset
  Nh   <- sum(hits)
  Nm   <- length(metric) - Nh
  
  runningES <- cumsum(
    ifelse(
      hits,
      abs(metric) / sum(abs(metric[hits])),
      -1 / Nm
    )
  )
  
  data.frame(
    Rank    = seq_along(runningES),
    ES      = runningES,
    Pathway = geneSetID
  )
}

es_df <- bind_rows(
  lapply(top_ids, function(id) get_running_es(gsea_res, id))
)

es_df$Pathway <- factor(
  es_df$Pathway,
  levels = top_ids,
  labels = top_labels
)

#---------------------------------------------------------
# 6️⃣ Correct stacked-color barcode data
#---------------------------------------------------------
barcode_df <- bind_rows(lapply(seq_along(top_ids), function(i) {
  geneset <- gsea_res@geneSets[[top_ids[i]]]
  hits <- which(names(gene_list) %in% geneset)
  
  data.frame(
    Rank = hits,
    ymin = (i - 1) / n_sets,
    ymax = i / n_sets,
    Pathway = top_labels[i]
  )
}))

barcode_df$Pathway <- factor(barcode_df$Pathway, levels = top_labels)

#---------------------------------------------------------
# 7️⃣ Plot A: Running Enrichment Scores
#---------------------------------------------------------
p_es <- ggplot(es_df, aes(Rank, ES, color = Pathway)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Running Enrichment Score",
    color = NULL
  ) +
  theme(
    legend.position = "right",
    axis.title.y = element_text(face = "bold")
  )

#---------------------------------------------------------
# 8️⃣ Plot B: Stacked-color gene hit barcode (FIXED)
#---------------------------------------------------------
p_hits <- ggplot(barcode_df, aes(x = Rank, color = Pathway)) +
  geom_segment(
    aes(xend = Rank, y = ymin, yend = ymax),
    linewidth = 0.9
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "none"
  )

#---------------------------------------------------------
# 9️⃣ Plot C: Ranked List Metric
#---------------------------------------------------------
metric_df <- data.frame(
  Rank   = seq_along(gene_list),
  Metric = gene_list
)

p_metric <- ggplot(metric_df, aes(x = Rank)) +
  geom_ribbon(
    aes(
      ymin = pmin(Metric, 0),
      ymax = pmax(Metric, 0)
    ),
    fill = "grey80"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Rank in Ordered Dataset",
    y = "Ranked List Metric"
  ) +
  theme(
    axis.title = element_text(face = "bold")
  )

#---------------------------------------------------------
# 🔟 Stack ALL THREE panels (FINAL FIGURE)
#---------------------------------------------------------
final_plot <- p_es / p_hits / p_metric +
  plot_layout(heights = c(3, 0.7, 1.2))

final_plot



#=========================================================
# GSVA (ssGSEA) for KEGG pathways
#=========================================================

# 1️⃣ Prepare expression matrix
expr_mat <- as.matrix(combat_final_rma_mat)
mode(expr_mat) <- "numeric"

# 2️⃣ Get KEGG gene sets from MSigDB
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
kegg_sets <- subset(c2_sets, grepl("KEGG", gs_name))
kegg_list <- split(kegg_sets$gene_symbol, kegg_sets$gs_name)

# 3️⃣ Create ssGSEA parameter object
ssgsea_params <- ssgseaParam(
  exprData = expr_mat,
  geneSets = kegg_list,
  alpha = 0.25,      # weighting exponent
  normalize = TRUE
)

# 4️⃣ Run GSVA
gsva_scores <- gsva(ssgsea_params, verbose = TRUE)

# gsva_scores: rows = pathways, cols = samples

# 5️⃣ Differential pathway analysis using limma
group_vector <- factor(group_vector)  # ensure factor
design <- model.matrix(~ 0 + group_vector)
colnames(design) <- levels(group_vector)

contrast <- makeContrasts(RA - NC, levels = design)

fit <- lmFit(gsva_scores, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Top 20 enriched pathways (by p-value)
top_pathways <- topTable(fit2, number = 20, sort.by = "P")
top_pathways$Pathway <- factor(rownames(top_pathways), levels = rev(rownames(top_pathways)))


plot_gsva_heatmap(
  gsva_scores = gsva_scores,
  group_vector = group_vector,
  top_n = 20,            # show top 20 pathways
  scale_rows = TRUE,
  main_title = "GSVA Analysis"
)


#========================WGCNA ANALYIS=======================#

# --------------Prepare expression data-----------------------
datExpr <- t(combat_final_rma_mat)
datExpr <- as.data.frame(datExpr)
datExpr[] <- lapply(datExpr, as.numeric)


#choose scale-free topology index power beta coeff
powers <- 1:20
sft <- pickSoftThreshold(datExpr,powerVector = powers,verbose = 5)
sft_data <- sft$fitIndices
softPower <- 8

#Run WGCNA
wgcna_res <- run_wgcna(datExpr,group_vector)

dissTOM<-wgcna_res$dissTOM
wgcna_res$sigModule
hubGenes<-(wgcna_res$hubGenes)

#---------------------INTERESECT DEGs and Hub Genes---------------------#
# DEG gene list
deg_genes <- rownames(degs)

# WGCNA hub gene list
wgcna_hub_genes <- hubGenes

# Intersection
candidate_genes <- intersect(deg_genes, wgcna_hub_genes)

plot_deg_wgcna_venn(deg_genes = deg_genes,wgcna_hub_genes = wgcna_hub_genes)


#===============ML Preps=====================#
expr_ml <- t(combat_final_rma_mat[candidate_genes, ])
expr_ml <- as.data.frame(expr_ml)
expr_ml[] <- lapply(expr_ml, as.numeric)
# Binary outcome (0/1)
y <- as.numeric(group_vector == "RA")


#========================LASSO 10-FOLD CV=========================#
lasso_out <- run_lasso_cv(expr_ml, group_vector)
lasso_genes <- lasso_out$selected_genes
length(lasso_genes)

#========================SVM-RFE 5-FOLD CV=========================#
svm_out <- run_svm_rfe(expr_ml, group_vector)
svm_genes <- svm_out$selected_genes
length(svm_genes)

#=========================Random-Forest============================#
rf_out <- run_random_forest(expr_ml, group_vector)
rf_genes <- rf_out$selected_genes
length(rf_genes)

#=============Intersection of 3 ML methods==================#
final_hub_genes <- Reduce(
  intersect,
  list(lasso_genes, svm_genes, rf_genes)
)

final_hub_genes
length(final_hub_genes)

plot_ml_venn(lasso_genes, svm_genes, rf_genes)


#panel B style A
cv <- lasso_out$cv_model
plot(cv)
abline(v = log(cv$lambda.min), col = "red", lty = 2)
abline(v = log(cv$lambda.1se), col = "blue", lty = 2)
legend("topright",
       legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"),
       lty = 2)


#panel B style B
cv <- lasso_out$cv_model
log_lambda <- log(cv$lambda)
cvm   <- cv$cvm    # mean cross-validated error
cvsd  <- cv$cvsd   # standard deviation
nzero <- cv$nzero  # number of non-zero coefficients
plot(log_lambda, cvm,
     pch = 16,
     col="red",
     ylim = range(cvm + cvsd, cvm - cvsd),
     xlab = "Log(λ)",
     ylab = "Misclassification Error",
     main = "LASSO CV Error")

# error bars
arrows(log_lambda,
       cvm - cvsd,
       log_lambda,
       cvm + cvsd,
       angle = 90,
       code = 3,
       length = 0.05,
       col = "grey")

# lambda.min and lambda.1se
abline(v = log(cv$lambda.min), lty = 2)
abline(v = log(cv$lambda.1se), lty = 2)

# numbers on top (non-zero features)
text(log_lambda,
     max(cvm + cvsd) * 1.05,
     labels = nzero,
     cex = 0.2)


#Panel D
svm_res <- svm_out$rfe_model
df <- svm_res$results

ggplot(df, aes(x = Variables, y = Accuracy)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 2, shape = 1) +  # open circles
  theme_bw() +
  labs(
    x = "Variables",
    y = "Accuracy (Cross-Validation)"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

# Panel E
plot(rf_out$rf_model,
     main = "ERROR & TREES")
legend("topright",
       legend = colnames(rf_out$rf_model$err.rate),
       col = 1:3, lty = 1)

#Panel F
# Sort all genes by importance
imp <- rf_out$importance
imp <- imp[order(imp[, "MeanDecreaseGini"], decreasing = TRUE), ]
all_genes <- names(imp)  # 172 genes

cv_error <- numeric(length(all_genes))

# Compute OOB error for top 1,2,...,172 genes
for(i in seq_along(all_genes)){
  genes_subset <- all_genes[1:i]
  
  rf_tmp <- randomForest(
    x = expr_ml[, genes_subset, drop = FALSE],
    y = factor(group_vector),
    ntree = 500
  )
  
  cv_error[i] <- rf_tmp$err.rate[nrow(rf_tmp$err.rate), "OOB"]
}

# Number of genes RF actually selected
opt_n <- length(rf_out$selected_genes)  # 43 in your case

# Plot
plot(1:length(all_genes), cv_error,
     type = "l",           # blue line
     lwd = 2,
     col = "blue",
     xlab = "Number of OTUs",
     ylab = "Cross-Validation Error",
     main = "Random Forest Feature Selection")
abline(v = opt_n, col = "black", lty = 1, lwd = 1)  # vertical line at 43 genes
grid()





#===================nomogram on the final hub genes================#

# Expression matrix for hub genes
expr_hub <- t(combat_final_rma_mat[final_hub_genes, ])

expr_hub <- as.data.frame(expr_hub)
expr_hub$RA <- ifelse(group_vector == "RA", 1, 0)

dd <- datadist(expr_hub)
options(datadist = "dd")

nomogram_model <- lrm(
  RA ~ CRYBG1 + `IGKV1OR1-1` + SKAP2 + BTN2A2 + CXCL13 + QPCT,
  data = expr_hub,
  x = TRUE,
  y = TRUE
)

nom <- nomogram(
  nomogram_model,
  fun = plogis,
  funlabel = "Risk of RA",
  lp = FALSE
)

plot(nom, xfrac = 0.4)


cal <- calibrate(
  nomogram_model,
  method = "boot",
  B = 1000
)

plot(
  cal,
  xlab = "Predicted probability",
  ylab = "Observed probability",
  main = "Calibration curve"
)

dca_model <- decision_curve(
  RA ~ CRYBG1 + `IGKV1OR1-1` + SKAP2 + BTN2A2 + CXCL13 + QPCT,
  data = expr_hub,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = FALSE
)

plot_decision_curve(
  dca_model,
  curve.names = "Nomogram",
  xlab = "Threshold probability",
  ylab = "Net benefit",
  legend.position = "topright"
)

nom_score <- predict(
  nomogram_model,
  type = "fitted"
)

roc_nom <- roc(expr_hub$RA, nom_score)

plot(
  roc_nom,
  col = "red",
  lwd = 2,
  main = "ROC curve of nomogram model"
)

auc(roc_nom)

par(mfrow = c(2,3))

for (gene in final_hub_genes) {
  roc_gene <- roc(expr_hub$RA, expr_hub[[gene]])
  
  plot(
    roc_gene,
    main = paste("ROC:", gene),
    col = "blue",
    lwd = 2
  )
  
  text(
    0.6, 0.2,
    paste("AUC =", round(auc(roc_gene), 3))
  )
}
