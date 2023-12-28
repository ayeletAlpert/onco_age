# load packages -----------------------------------------------------------
rm(list = ls())
memory.limit(size = 100000)
# Load the readr library
library(readr)
library(Matrix)
library(Seurat)
require(Biobase)
library(AnnotationDbi)
library(org.Hs.eg.db)
require(stringr)
require(ggplot2)
library(BayesPrism)
library(devtools)
library(matrixStats)
library(SummarizedExperiment)
library(InstaPrism)
library(survival)
library(ComplexHeatmap)
library(clusterProfiler)
theme_set(theme_bw())

# Read data structures -------------------------------------------------------------
bp.res <- readRDS("~/MDPhD/oncology/BayesPrismLUAD/Data/bp_res_high_res_cells.rds")
LUAD_RnaseqSE <- readRDS(file = "C:\\Users\\ShenorLab\\Documents\\MDPhD\\oncology\\RNAseqTCGACOAD\\data\\LUAD-tcga_data.rds")
plot(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")), 
     bp.res@Post.ini.cs@theta["Epithelial cells",], xlab = "tumor purity by TCGA", 
     ylab = "fraction of tumor cells by BayesPrism", main = "Correlation between tumor purity and estimated tumor fraction")
summary(lm(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")) ~
             bp.res@Post.ini.cs@theta["Epithelial cells",]))
cor(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")),
    bp.res@Post.ini.cs@theta["Epithelial cells",], use = "complete.obs", method = "spearman")

#process data purity from TCGA from the manuscript: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4671203/#S1
TCGA_purity <- read.csv(file = "D:/Ayelet/bulk_TCGA/estimate_purity_TCGA_samples.csv", row.names = 1)
cor(TCGA_purity[gsub("^((?:[^-]+-){3}[^-]+).*", "\\1", colnames(LUAD_RnaseqSE), perl = TRUE),"ESTIMATE"], 
    bp.res@Post.ini.cs@theta["Epithelial cells",], use = "complete.obs", method = "spearman")

plot(TCGA_purity[gsub("^((?:[^-]+-){3}[^-]+).*", "\\1", colnames(LUAD_RnaseqSE), perl = TRUE),"ESTIMATE"], 
     bp.res@Post.ini.cs@theta["Epithelial cells",], xlab = "tumor purity by TCGA", 
     ylab = "fraction of tumor cells by BayesPrism", main = "Correlation between tumor purity and estimated tumor fraction")


# correlation of cell type abundance with age ----------------------------------------------------------------
# cell_types_immune <- c("T lymphocytes", "Myeloid cells", "B lymphocytes", "NK cells")
Kim_Lee_cell_annot <- readRDS( file = "~/MDPhD/oncology/BayesPrismLUAD/Data/Kim_Lee_cell_annot.rds")

cell_types_immune <- data.frame(cell_type = c("Follicular B cells", "CD4+ Th", "CD8 low T", "Monocytes", "Exhausted CD8+ T", "Cytotoxic CD8+ T",
                                              "CD8+/CD4+ Mixed Th", "mo-Mac", "Treg","Naive CD4+ T","Naive CD8+ T", "MALT B cells", "NK","CD1c+ DCs",
                                              "Exhausted Tfh", "Plasma cells", "Activated DCs", "CD141+ DCs", "pDCs"))
cell_types_immune$parent_cell_type <- sapply(cell_types_immune$cell_type, function(cell_type){
  return(unique(Kim_Lee_cell_annot$Cell_type.refined[Kim_Lee_cell_annot$Cell_subtype == cell_type]))})

cell_types_immune$parent_cell_type[cell_types_immune$cell_type %in% c("CD8 low T", "Exhausted CD8+ T", "Cytotoxic CD8+ T", "Naive CD8+ T")] <- "CD8+T"
cell_types_immune$parent_cell_type[cell_types_immune$cell_type %in% c("CD4+ Th", "Exhausted Tfh", "Naive CD4+ T", "Treg", "CD8+/CD4+ Mixed Th")] <- "CD4+T"
cell_types_immune$parent_cell_type[cell_types_immune$cell_type %in% c("Monocytes", "mo-Mac")] <- "Monocytes"
cell_types_immune$parent_cell_type[cell_types_immune$cell_type %in% c("pDCs", "CD141+ DCs", "Activated DCs", "CD1c+ DCs")] <- "DC"
cell_types_immune$parent_cell_type[cell_types_immune$cell_type %in% c("NK")] <- "NK"

#calculate the proportion of each immune cell out of all:
immune_cell_types_norm <- t(apply(bp.res@Post.ini.cs@theta[cell_types_immune$cell_type,],2,function(x){return(x/sum(x))}))

#calculate the relative abundance of each immune cell out of the parent population:
immune_cell_types_norm_to_parent <- do.call('cbind', lapply(unique(cell_types_immune$parent_cell_type), function(parent_cell_type){
  print(parent_cell_type)
  child_cell_type <- cell_types_immune$cell_type[cell_types_immune$parent_cell_type == parent_cell_type]
  if(length(child_cell_type) == 1){return(NULL)}
  t(apply(bp.res@Post.ini.cs@theta[child_cell_type,],2,function(x){return(x/sum(x))}))
}))
  
#correlate the cell type abundance (relative to all cells) with age:
p_cell_type_abun_age_total <- do.call('rbind', lapply(colnames(immune_cell_types_norm), function(cell_type){
  print(cell_type)
  
  lmres_tumor <- summary(lm(immune_cell_types_norm[colnames(LUAD_RnaseqSE)[LUAD_RnaseqSE$definition == "Primary solid Tumor"],cell_type] ~
                              LUAD_RnaseqSE$age_at_diagnosis[LUAD_RnaseqSE$definition == "Primary solid Tumor"]))
  lmres_healthy <- summary(lm(immune_cell_types_norm[colnames(LUAD_RnaseqSE)[LUAD_RnaseqSE$definition == "Solid Tissue Normal"],cell_type] ~
                                LUAD_RnaseqSE$age_at_diagnosis[LUAD_RnaseqSE$definition == "Solid Tissue Normal"]))
  
  return(data.frame(cell_type = cell_type, 
                    p_age_tumor = coef(lmres_tumor)[2,4], dir_age_tumor = sign(coef(lmres_tumor)[2,1]), 
                    p_age_healthy = coef(lmres_healthy)[2,4], dir_age_healthy = sign(coef(lmres_healthy)[2,1])))
}))

#correlate the cell type abundance (relative to parent) with age:
p_cell_type_abun_age_parent <- do.call('rbind', lapply(colnames(immune_cell_types_norm_to_parent), function(cell_type){
  print(cell_type)
  
  lmres_tumor_parent <- summary(lm(immune_cell_types_norm_to_parent[colnames(LUAD_RnaseqSE)[LUAD_RnaseqSE$definition == "Primary solid Tumor"],cell_type] ~
                                     LUAD_RnaseqSE$age_at_diagnosis[LUAD_RnaseqSE$definition == "Primary solid Tumor"]))
  lmres_healthy_parent <- summary(lm(immune_cell_types_norm_to_parent[colnames(LUAD_RnaseqSE)[LUAD_RnaseqSE$definition == "Solid Tissue Normal"],cell_type] ~
                                       LUAD_RnaseqSE$age_at_diagnosis[LUAD_RnaseqSE$definition == "Solid Tissue Normal"]))
  
  return(data.frame(cell_type = cell_type, 
                    p_age_tumor_parent = coef(lmres_tumor_parent)[2,4], dir_age_tumor_parent = sign(coef(lmres_tumor_parent)[2,1]), 
                    p_age_healthy_parent = coef(lmres_healthy_parent)[2,4], dir_age_healthy_parent = sign(coef(lmres_healthy_parent)[2,1])))
}))

hist(p_cell_type_abun_age_parent$p_age_tumor_parent,10)
hist(p_cell_type_abun_age_total$p_age_tumor,10)

#plot as circular heatmap:
p_cell_type_abun_age_parent_4_hm <- data.frame(tumor_parent = -1*p_cell_type_abun_age_parent$dir_age_tumor_parent*log(p_cell_type_abun_age_parent$p_age_tumor_parent),
                                               healthy_parent = -1*p_cell_type_abun_age_parent$dir_age_healthy_parent*log(p_cell_type_abun_age_parent$p_age_healthy_parent), 
                                               row.names = p_cell_type_abun_age_parent$cell_type)
pheatmap(p_cell_type_abun_age_parent_4_hm)

#read mutation data:
mutation_data <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/mutation_data.rds")
samples_with_driver_mutation <- unique(mutation_data$Tumor_Sample_Barcode[mutation_data$Hugo_Symbol %in% c("EGFR", "ROS1", "ALK", "BRAF")])
samples_with_driver_mutation <- sub("(-[^-]+){3}$", "", samples_with_driver_mutation)
LUAD_RnaseqSE$driver_mut <- sapply(LUAD_RnaseqSE$sample, function(sample){sample %in% samples_with_driver_mutation})

#plot the relation between the cell type frequency and age
df <- data.frame(age = LUAD_RnaseqSE$age_at_diagnosis[LUAD_RnaseqSE$definition == "Primary solid Tumor"],
                 immune_cell_types_norm_to_parent[colnames(LUAD_RnaseqSE)[LUAD_RnaseqSE$definition == "Primary solid Tumor"],],
                 mut = LUAD_RnaseqSE$driver_mut[LUAD_RnaseqSE$definition == "Primary solid Tumor"], check.names = T)
ggplot(df, aes(x = age, y = Cytotoxic.CD8..T, color = mut)) + geom_point() + geom_smooth(method = "lm")

# cell-specific expression ------------------------------------------------
cell_types <- rownames(bp.res@Post.ini.cs@theta)
cell_specific_exp_dir <- "~/MDPhD/oncology/BayesPrismLUAD/Data/cell_specific_expression"

# cell_type_specific_expression_age_p <- do.call('rbind', lapply(cell_types, function(cell_type){
#   print(cell_type)
#   cell_type_specific_expression <- readRDS(file = file.path(cell_specific_exp_dir, paste0(str_replace(cell_type, "/",""), "_specific_expression_high_res_cells.rds")))
#   
#   #filter-out those low-abundant genes with expression > 0 in < 10% of cells
#   cell_type_specific_expression <- cell_type_specific_expression[!rowSums(cell_type_specific_expression > 0) < 0.1*ncol(cell_type_specific_expression),]
#   
#   #filter-out the healthy samples:
#   cell_type_specific_expression <- cell_type_specific_expression[,colnames(LUAD_RnaseqSE)[colData(LUAD_RnaseqSE)[,"definition"] == "Primary solid Tumor"]]
#   
#   p_age_gene <- do.call('rbind', lapply(rownames(cell_type_specific_expression), function(gene_name){
#     print(gene_name)
#     lmres_gene_age <- summary(lm(log(cell_type_specific_expression[gene_name,] + 1) ~ 
#                                    colData(LUAD_RnaseqSE)[colnames(cell_type_specific_expression),"age_at_diagnosis"]))
#     # plot(log(cell_type_specific_expression[gene_name,] + 1),
#     #      colData(LUAD_RnaseqSE)[colnames(cell_type_specific_expression),"age_at_diagnosis"])
#     
#     return(data.frame(gene = gene_name, p_age = coef(lmres_gene_age)[2,4]))
#   }))
#   p_age_gene$p_adj_BH <- p.adjust(p_age_gene$p_age, method = "BH")
#   return(data.frame(cell_type = cell_type, p_age_gene))
# }))
# saveRDS(cell_type_specific_expression_age_p, file = file.path(cell_specific_exp_dir, 'cell_type_specific_expression_high_res_cells_age_p.rds'))
cell_type_specific_expression_age_p <- readRDS(file = file.path(cell_specific_exp_dir, 'cell_type_specific_expression_high_res_cells_age_p.rds'))

ggplot(cell_type_specific_expression_age_p, aes(x = p_age)) + geom_histogram() + facet_wrap(~cell_type, scales = "free") +
  ggtitle("p value of linear regression - expression VS age - per cell type")
ggplot(cell_type_specific_expression_age_p, aes(x = p_adj_BH)) + geom_histogram() + facet_wrap(~cell_type, scales = "free") +
  ggtitle("adjusted p value of linear regression - expression VS age - per cell type")

# correlation of cell type abundance with survival ------------------------
surv_data <- cbind(colData(LUAD_RnaseqSE)[,c("days_to_last_follow_up", "vital_status", "days_to_death", "ajcc_pathologic_stage", "definition", "age_at_diagnosis")],
                   immune_cell_types_norm_to_parent[colnames(LUAD_RnaseqSE),])
surv_data <- surv_data[surv_data$definition == "Primary solid Tumor",]
surv_data$stage_short <- str_extract(surv_data$ajcc_pathologic_stage, "[I]+V*")
surv_data$surv_state <- surv_data$vital_status == "Dead"
surv_data$TTE <- apply(surv_data,1,function(x){return(max(as.numeric(x['days_to_death']), as.numeric(x['days_to_last_follow_up']), na.rm = T))})

sum(is.infinite(surv_data[,c("surv_state")]))
surv_data <- surv_data[!is.na(surv_data$surv_state) & !is.na(surv_data$TTE) & !is.na(surv_data$stage_short) & !is.infinite(surv_data$TTE),]

p_abundance_cell_types_survival <- do.call('rbind', lapply(colnames(immune_cell_types_norm_to_parent), function(cell_type){
  cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$stage_short + surv_data$age_at_diagnosis + surv_data[,cell_type])
  cox_res <- summary(cox)
  return(data.frame(cell_type = cell_type, p = coef(cox_res)[nrow(coef(cox_res)),5], dir = sign(coef(cox_res)[nrow(coef(cox_res)),1])))
}))

# cell-specific expression ------------------------------------------------
cell_types <- rownames(bp.res@Post.ini.cs@theta)
cell_specific_exp_dir <- "~/MDPhD/oncology/BayesPrismLUAD/Data/cell_specific_expression"

cell_type_specific_expression_age_p <- do.call('rbind', lapply(cell_types, function(cell_type){
  print(cell_type)
  cell_type_specific_expression <- readRDS(file = file.path(cell_specific_exp_dir, paste0(cell_type, "_specific_expression.rds")))
  
  #filter-out those low-abundant genes with expression > 0 in < 10% of cells
  cell_type_specific_expression <- cell_type_specific_expression[!rowSums(cell_type_specific_expression > 0) < 0.1*ncol(cell_type_specific_expression),]
  
  #filter-out the healthy samples:
  cell_type_specific_expression <- cell_type_specific_expression[,colnames(LUAD_RnaseqSE)[colData(LUAD_RnaseqSE)[,"definition"] == "Primary solid Tumor"]]
  
  p_age_gene <- do.call('rbind', lapply(rownames(cell_type_specific_expression), function(gene_name){
    print(gene_name)
    lmres_gene_age <- summary(lm(log(cell_type_specific_expression[gene_name,] + 1) ~ 
                                   colData(LUAD_RnaseqSE)[colnames(cell_type_specific_expression),"age_at_diagnosis"]))
    # plot(log(cell_type_specific_expression[gene_name,] + 1),
    #      colData(LUAD_RnaseqSE)[colnames(cell_type_specific_expression),"age_at_diagnosis"])
    
    return(data.frame(gene = gene_name, p_age = coef(lmres_gene_age)[2,4]))
  }))
  p_age_gene$p_adj_BH <- p.adjust(p_age_gene$p_age, method = "BH")
  return(data.frame(cell_type = cell_type, p_age_gene))
}))
ggplot(cell_type_specific_expression_age_p, aes(x = p_age)) + geom_histogram() + facet_wrap(~cell_type, scales = "free") +
  ggtitle("p value of linear regression - expression VS age - per cell type")
ggplot(cell_type_specific_expression_age_p, aes(x = p_adj_BH)) + geom_histogram() + facet_wrap(~cell_type, scales = "free") +
  ggtitle("adjusted p value of linear regression - expression VS age - per cell type")


# Association of genes with tumor progression -----------------------------
LUAD_RnaseqSE <- readRDS(file = "C:\\Users\\ShenorLab\\Documents\\MDPhD\\oncology\\RNAseqTCGACOAD\\data\\LUAD-tcga_data.rds")
relevant_stages <- c("Stage IA", "Stage IB", "Stage IIA", "Stage IIB", "Stage IIIA", "Stage IIIB", "Stage IV")
rel_samples_tumor <- colnames(LUAD_RnaseqSE)[(LUAD_RnaseqSE$ajcc_pathologic_stage %in% relevant_stages) &
                                               (LUAD_RnaseqSE$definition == "Primary solid Tumor")]
rel_samples_healthy <- colnames(LUAD_RnaseqSE)[(LUAD_RnaseqSE$ajcc_pathologic_stage %in% relevant_stages) &
                                                 (LUAD_RnaseqSE$definition == "Solid Tissue Normal")]

LUAD_RnaseqSE_stage <- sapply(rel_samples_tumor, function(sample){return(which(relevant_stages == colData(LUAD_RnaseqSE)[sample,"ajcc_pathologic_stage"]))})
LUAD_RnaseqSE_meta_stage <- sapply(rel_samples_tumor, function(sample){return(str_extract(colData(LUAD_RnaseqSE)[sample,"ajcc_pathologic_stage"], "[I]+V*"))})
  
#convert rownames (emsemble id) to gene symbols:
rownames(LUAD_RnaseqSE) <- str_replace(rownames(LUAD_RnaseqSE), "\\.[0-9]+","")
rownames(count_data) <- ensembl2sym(rownames(count_data))

#  FUNCTIONS FOR MAPPING GENE SYMBOLS TO EMSEBLE IDS ------------------------
res <- mapIds(org.Hs.eg.db, keys <- rownames(rowData(LUAD_RnaseqSE)), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
ensembl2sym <- function(ensembl){return(res[ensembl])}
res_op_sym <- names(res)
names(res_op_sym) <- res
sym2ensembl <- function(sym){return(res_op_sym[sym])}

#  FUNCTIONS FOR MAPPING ENTREZ_ID TO EMSEBLE IDS ------------------------
res <- mapIds(org.Hs.eg.db, keys <- rownames(rowData(LUAD_RnaseqSE)), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
ensembl2entrez <- function(ensembl){return(res[ensembl])}
res_op_entrez <- names(res)
names(res_op_entrez) <- res
entrez2ensembl <- function(sym){return(res_op_entrez[sym])}

p_vals_genes_stage_tumor = do.call('rbind', lapply(rownames(LUAD_RnaseqSE)[!is.na(rowData(LUAD_RnaseqSE)[,"hgnc_id"])], function(gene_ensembl){
  print(gene_ensembl)
  lm_res <- summary(lm(assays(LUAD_RnaseqSE)$tpm_unstrand[gene_ensembl,rel_samples_tumor] ~ LUAD_RnaseqSE_stage))
  return(data.frame(gene = gene_ensembl, gene_symbol = rowData(LUAD_RnaseqSE)[str_replace(gene_ensembl, "\\.[0-9]+",""),"gene_name"], 
                    p_res = coefficients(lm_res)[2,4], coef = coefficients(lm_res)[2,1]))
}))

gene_ensembl <- "ENSG00000004478"
df_gene <- data.frame(exp = assays(LUAD_RnaseqSE)$tpm_unstrand[gene_ensembl,rel_samples_tumor], stage = LUAD_RnaseqSE_stage)
ggplot(df_gene, aes(x = factor(stage), y = exp)) + geom_boxplot()
t.test(df_gene$exp[df_gene$stage %in% c(1,2)], df_gene$exp[df_gene$stage %in% c(5,6,7)])

p_vals_genes_stage_tumor$ENTREZ <- ensembl2entrez(p_vals_genes_stage_tumor$gene)
hist(p_vals_genes_stage_tumor$p_res, main = "p values histogram of DE genes", xlab = "p value")
p_vals_genes_stage_tumor$bh_p <- p.adjust(p_vals_genes_stage_tumor$p_res)
hist(p_vals_genes_stage_tumor$bh_p)
p_vals_genes_stage_tumor[!is.na(p_vals_genes_stage_tumor$bh_p) & p_vals_genes_stage_tumor$bh_p < 0.1,]
sum(p_vals_genes_stage_tumor$bh_p < 0.1, na.rm = T)

#functional enrichment analysis:
P_THRESH <- 0.05
sig_genes <- p_vals_genes_stage_tumor$ENTREZ[!is.na(p_vals_genes_stage_tumor$bh_p) & p_vals_genes_stage_tumor$bh_p < P_THRESH]
sig_genes_ensemble <- p_vals_genes_stage_tumor$gene[!is.na(p_vals_genes_stage_tumor$bh_p) & p_vals_genes_stage_tumor$bh_p < P_THRESH]
length(sig_genes)

#over representation analysis:
enrich_res_stage <- enrichKEGG(sig_genes, organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                               universe = p_vals_genes_stage_tumor$ENTREZ[!is.na(p_vals_genes_stage_tumor$bh_p)], 
                               minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
barplot(enrich_res_stage, showCategory=40) 
enrich_res_stage@result

pca_res <- data.frame(prcomp(t(assays(LUAD_RnaseqSE)$tpm_unstrand[sig_genes_ensemble,rel_samples_tumor]), scale. = T)$x[,c("PC1","PC2")])
pca_res$stage <- LUAD_RnaseqSE_stage
pca_res$metastage <-LUAD_RnaseqSE_meta_stage
ggplot(pca_res, aes(x = PC1, y = PC2, color = stage)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = metastage)) + geom_point()

ggplot(pca_res, aes(x = factor(stage), y = PC1)) + geom_boxplot() + ggtitle("PC1 distribution of advancing stages")
ggplot(pca_res, aes(x = factor(metastage), y = PC1)) + geom_boxplot()

ggplot(pca_res, aes(x = PC1, y = ..density.., fill = factor(stage))) + geom_density(alpha = 0.2)
ggplot(pca_res, aes(x = PC1, y = ..density.., fill = factor(metastage))) + geom_density(alpha = 0.2)

#check the difference in outcome in class 1 and 2
rel_stage = c(1:4)
QUANT = 0.25
low_PC1_stage <- rownames(pca_res)[pca_res$stage %in% rel_stage & pca_res$PC1 < quantile(pca_res$PC1[pca_res$stage %in% rel_stage], QUANT)]
high_PC1_stage <- rownames(pca_res)[pca_res$stage %in% rel_stage & pca_res$PC1 > quantile(pca_res$PC1[pca_res$stage %in% rel_stage], (1-QUANT))]

surv_data <- cbind(colData(LUAD_RnaseqSE)[c(low_PC1_stage, high_PC1_stage),
                                          c("days_to_last_follow_up", "vital_status", "days_to_death", "ajcc_pathologic_stage", "definition", "age_at_diagnosis")],
                   PC1_class = c(low_PC1_stage, high_PC1_stage) %in% low_PC1_stage, 
                   PC1 = pca_res[c(low_PC1_stage, high_PC1_stage),"PC1"])
surv_data$surv_state <- surv_data$vital_status == "Dead"
surv_data$TTE <- apply(surv_data,1,function(x){return(max(as.numeric(x['days_to_death']), as.numeric(x['days_to_last_follow_up']), na.rm = T))})

sum(is.infinite(surv_data[,c("surv_state")]))
surv_data <- surv_data[!is.na(surv_data$surv_state) & !is.na(surv_data$TTE) & !is.infinite(surv_data$TTE),]

cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$PC1_class + surv_data$age_at_diagnosis + surv_data$ajcc_pathologic_stage)
cox_res <- summary(cox)
cox_res

cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$PC1 + surv_data$age_at_diagnosis + surv_data$ajcc_pathologic_stage)
cox_res <- summary(cox)
cox_res


time <- surv_data$TTE
event <- surv_data$surv_state
fit <- survfit(Surv(TTE, surv_state) ~ as.factor(PC1_class) + age_at_diagnosis + ajcc_pathologic_stage, data = surv_data)

# Plot Kaplan-Meier curve
ggsurvplot(fit, data = surv_data, 
           risk.table = TRUE, 
           pval = TRUE,
           surv.median.line = "hv",
           palette = "jco")

#illustrate the dynamics of different gene modules along PC1:
#hallmark geneset:
gene_sets <- read.gmt(gmtfile = "~/MDPhD/oncology/BayesPrismLUAD/Data/hallmark_genesets_gene_symbols.gmt")

scaled_exp_matrix <- do.call('rbind', lapply(as.character(unique(gene_sets$term)), function(rel_term){
  print(rel_term)
  term_gene_symbols <- gene_sets$gene[gene_sets$term == rel_term]
  term_gene_ensmbl <- sym2ensembl(term_gene_symbols)
  exp_matrix <- assays(LUAD_RnaseqSE)$tpm_unstrand[term_gene_ensmbl[term_gene_ensmbl %in% rownames(LUAD_RnaseqSE)], 
                                                   rownames(pca_res)[order(pca_res$PC1)]]
  exp_matrix <- exp_matrix[apply(exp_matrix,1,var) > 0,]
  scaled_exp_matrix <- t(apply(exp_matrix,1,scale))
  colnames(scaled_exp_matrix) <- colnames(exp_matrix)
  
  return(apply(scaled_exp_matrix,2,median))
}))
rownames(scaled_exp_matrix) <- as.character(unique(gene_sets$term))

#calculate p value of linear regression against the trajectory
p_vals_lm_PC1 <- apply(scaled_exp_matrix,1,function(x){
  lm_res_term <- summary(lm(x ~ pca_res$PC1[order(pca_res$PC1)]))
  return(coefficients(lm_res_term)[2,4])
})

ordered_scaled_exp_matrix <- scaled_exp_matrix[order(p_vals_lm_PC1),]
ordered_scaled_exp_matrix[ordered_scaled_exp_matrix > 2] <- 2
ordered_scaled_exp_matrix[ordered_scaled_exp_matrix < -2] <- -2
pheatmap(ordered_scaled_exp_matrix, cluster_rows = F, cluster_cols = F, show_colnames = F)
  
rel_terms <- c("HALLMARK_DNA_REPAIR", "HALLMARK_P53_PATHWAY")
df_plot_functions_trajectory <- data.frame(PC1 = pca_res$PC1[order(pca_res$PC1)], sample = colnames(ordered_scaled_exp_matrix),
                                           term = rep(rel_terms, each = ncol(ordered_scaled_exp_matrix)),
                                           median_exp = c(t(ordered_scaled_exp_matrix[rel_terms,])))
ggplot(df_plot_functions_trajectory, aes(x = PC1, y = median_exp, color = term)) + geom_point(alpha = 0.2) + geom_smooth()

#correlate with mutation data:
mutation_data <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/mutation_data.rds")
TP53_mutation <- mutation_data$Tumor_Sample_Barcode[mutation_data$Hugo_Symbol == "TP53"]
TP53_mutation <- sub("(-[^-]+){3}$", "", TP53_mutation)
p53_data <- df_plot_functions_trajectory[df_plot_functions_trajectory$term == "HALLMARK_P53_PATHWAY",]
p53_data$sample <- sub("(-[^-]+){3}$", "", p53_data$sample)
p53_data$mut_p53 <- sapply(p53_data$sample, function(x){return(x %in% TP53_mutation)})
ggplot(p53_data, aes(x = PC1, y = median_exp, color = mut_p53)) + geom_point() + geom_smooth()
ggplot(p53_data, aes(x = mut_p53, y = PC1)) + geom_boxplot()
t.test(p53_data$PC1[p53_data$mut_p53], p53_data$PC1[!p53_data$mut_p53])

#apply BayesPrism on the gene expression to reflect changes in cell abundances along the trajectory
bp.res <- readRDS("~/MDPhD/oncology/BayesPrismLUAD/Data/bp_res_high_res_cells.rds")

scaled_cellular_abundance_matrix <- t(apply(bp.res@Post.ini.cs@theta[,rel_samples_tumor],1,scale))
colnames(scaled_cellular_abundance_matrix) <- rel_samples_tumor

p_vals_lm_PC1_cells <- apply(scaled_cellular_abundance_matrix,1,function(x){
  lm_res_cell <- summary(lm(x ~ pca_res[colnames(scaled_cellular_abundance_matrix),"PC1"]))
  return(coefficients(lm_res_cell)[2,4])
})

#normalize the scaled value of outliers from both sides (extremely high or extremely low)
NORM = T
if(NORM){
  scaled_cellular_abundance_matrix[scaled_cellular_abundance_matrix > 2] <- 2
  scaled_cellular_abundance_matrix[scaled_cellular_abundance_matrix < -2] <- -2
}

rel_cells <- c("Cytotoxic CD8+ T", "Tumor ECs")
df_plot_cells_trajectory <- data.frame(PC1 = pca_res[colnames(scaled_cellular_abundance_matrix),"PC1"], 
                                       sample = colnames(scaled_cellular_abundance_matrix),
                                       cell = rep(rel_cells, each = ncol(scaled_cellular_abundance_matrix)),
                                       scaled_abundance = c(t(scaled_cellular_abundance_matrix[rel_cells,])))
ggplot(df_plot_cells_trajectory, aes(x = PC1, y = scaled_abundance, color = cell)) + geom_point(alpha = 0.2) + geom_smooth()

#====================================================================







# Convert the covariate matrix to a data frame
cov_matrix <- data.frame(cbind(surv_data$stage_short,surv_data$age_at_diagnosis,surv_data[,cell_types_immune]))
cov_matrix_df <- as.data.frame(cov_matrix)

# Combine the survival data and the covariate data frame
coxph_data <- data.frame(
  Time = surv_data$TTE,
  Status = surv_data$surv_state,
  cov_matrix_df
)

# Run the Cox proportional hazards model
cox_model <- coxph(Surv(Time, Status) ~ ., data = coxph_data)
summary(cox_model)

cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$stage_short + surv_data$age_at_diagnosis + surv_data[,cell_types_immune])
cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ cov_matrix)
cox <- coxph(Surv(TTE, surv_state) ~ cov_matrix, data = data.frame(TTE = surv_data$TTE, surv_state = surv_data$surv_state, cov_matrix))

cox_model <- coxph(Surv(Time, Status) ~ covariate_matrix, data = data.frame(Time = surv_obj$time, Status = surv_obj$status, covariate_matrix))


 


cox <- coxph(Surv(TTE, surv_state) ~ stage_short, data = surv_data)
cox <- coxph(formula = Surv(TTE, surv_state) ~ stage_short, data = surv_data)
summary(cox)

summary(cox)

autoplot(km_fit)

summary(km_fit)


bp.res@Post.ini.ct@theta

#Association of cell type abundance with age:
pca_res = data.frame(prcomp(t(bp.res@Post.ini.cs@theta), scale. = T)$x[,c("PC1","PC2")], colData(LUAD_RnaseqSE),
                     t(bp.res@Post.ini.cs@theta))
ggplot(pca_res, aes(x = PC1, y = PC2, color = sample_type)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Epithelial.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = T.lymphocytes)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Fibroblasts)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Myeloid.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = B.lymphocytes)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = NK.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Endothelial.cells)) + geom_point()
cell_types <- c("Epithelial.cells", "T.lymphocytes", "Fibroblasts", "Myeloid.cells", "B.lymphocytes", "NK.cells", "Endothelial.cells")

annot_rows <- data.frame(pca_res[,c("definition", "paper_Age.at.diagnosis")])
annot_rows$paper_Age.at.diagnosis <- as.numeric(as.character(annot_rows$paper_Age.at.diagnosis))
pheatmap(pca_res[,cell_types], annotation_row = annot_rows, show_rownames = F)

cell_types_immune <- c("T.lymphocytes", "Myeloid.cells", "B.lymphocytes", "NK.cells")
normalized_cell_type_immune <- t(apply(pca_res[,cell_types_immune],1,function(x){return(x/sum(x))}))
colnames(normalized_cell_type_immune) <- paste0("norm_",cell_types_immune)
pca_res <- cbind(pca_res, normalized_cell_type_immune)

sapply(colnames(normalized_cell_type_immune), function(cell_type){
  lmres <- summary(lm(pca_res[pca_res$definition == "Primary solid Tumor",cell_type] ~ pca_res[pca_res$definition == "Primary solid Tumor",
                                                                                               "age_at_diagnosis"]))
  return(data.frame(cell_type = cell_type, p = coef(lmres)[2,4]))
})
abun_cell_type_age <- melt(pca_res[pca_res$definition == "Primary solid Tumor", c(colnames(normalized_cell_type_immune), 'age_at_diagnosis')],
                           id.vars = 'age_at_diagnosis')

ggplot(abun_cell_type_age, aes(x = age_at_diagnosis, y = value)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~ variable, scales = "free")

sapply(colnames(normalized_cell_type_immune), function(cell_type){
  lmres <- summary(lm(pca_res[pca_res$definition == "Solid Tissue Normal",cell_type] ~ pca_res[pca_res$definition == "Solid Tissue Normal",
                                                                                               "age_at_diagnosis"]))
  return(data.frame(cell_type = cell_type, p = coef(lmres)[2,4]))
})
abun_cell_type_age_H <- melt(pca_res[pca_res$definition == "Solid Tissue Normal", c(colnames(normalized_cell_type_immune), 'age_at_diagnosis')],
                           id.vars = 'age_at_diagnosis')

p_values <- abun_cell_type_age_H %>%
  group_by(variable) %>%
  summarise(p_value = summary(lm(value ~ age_at_diagnosis))$coefficients[2, 4])

ggplot(abun_cell_type_age_H, aes(x = age_at_diagnosis, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ variable, scales = "free")+
  
  # Add p-values as text annotations to each facet
  geom_text(data = p_values, aes(label = sprintf("p = %.6f", p_value)),
            x = Inf, y = -Inf, hjust = 3.8, vjust = -25, size = 3)

# plot-david --------------------------------------------------------------

p_values <- abun_cell_type_age %>%
  group_by(variable) %>%
  summarise(p_value = summary(lm(value ~ age_at_diagnosis))$coefficients[2, 4])

# Create the ggplot with facets
ggplot(abun_cell_type_age, aes(x = age_at_diagnosis, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ variable, scales = "free") +
  
  # Add p-values as text annotations to each facet
  geom_text(data = p_values, aes(label = sprintf("p = %.6f", p_value)),
            x = Inf, y = -Inf, hjust = 3.8, vjust = -25, size = 3)


# end of plot david ---------------



plot(pca_res[pca_res$definition == "Primary solid Tumor","norm_NK.cells"] ~ pca_res[pca_res$definition == "Primary solid Tumor",
                                                                              "age_at_diagnosis"])



#survival analysis with T lymphocytes:
surv_data <- cbind(pca_res[,c("days_to_last_follow_up", "vital_status", "days_to_death", "ajcc_pathologic_stage", cell_types, "sample_type", "age_at_diagnosis")],
                   normalized_cell_type_immune)
surv_data <- surv_data[surv_data$sample_type == "Primary Tumor",]
surv_data$stage_short <- str_extract(surv_data$ajcc_pathologic_stage, "[I]+V*")
surv_data$surv_state <- surv_data$vital_status == "Dead"
surv_data$TTE <- apply(surv_data,1,function(x){return(max(as.numeric(x['days_to_death']), as.numeric(x['days_to_last_follow_up']), na.rm = T))})

sum(is.infinite(surv_data[,c("surv_state")]))
surv_data <- surv_data[!is.na(surv_data$surv_state) & !is.na(surv_data$TTE) & !is.na(surv_data$stage_short) & !is.infinite(surv_data$TTE),]

cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$stage_short + surv_data$age_at_diagnosis)
summary(cox)

cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$stage_short + surv_data$norm_Myeloid.cells)
summary(cox)

cox <- coxph(Surv(TTE, surv_state) ~ stage_short, data = surv_data)
cox <- coxph(formula = Surv(TTE, surv_state) ~ stage_short, data = surv_data)
summary(cox)

summary(cox)

autoplot(km_fit)

summary(km_fit)

apply(bp.res@Post.ini.cs@theta, 2, function(x){return(x/sum())})

pca_res_LUAD = data.frame(prcomp(t(subset(bp.res@Post.ini.cs@theta), LUAD_RnaseqSE$), scale. = T)$x[,c("PC1","PC2")], colData(LUAD_RnaseqSE),
                          t(bp.res@Post.ini.cs@theta))
ggplot(pca_res, aes(x = PC1, y = PC2, color = sample_type)) + geom_point()
ggplot(pca_res, aes(x = age_at_index, y = PC2, color = sample_type)) + geom_point()
cor(pca_res$PC2[pca_res$sample_type == "Primary Tumor"], 
    pca_res$age_at_diagnosis[pca_res$sample_type == "Primary Tumor"], use = "complete.obs")

ggplot(pca_res, aes(x = PC1, y = PC2, color = age_at_index)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = NK.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Epithelial.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = age_at_index)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = paper_expression_subtype)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = sample_type)) + geom_point()

ggplot(pca_res, aes(x = age_at_index, y = PC1)) + geom_point()
ggplot(pca_res, aes(x = ajcc_pathologic_stage, y = PC1)) + geom_boxplot()
ggplot(pca_res, aes(x = pack_years_smoked, y = PC1)) + geom_point()
ggplot(pca_res, aes(x = age_at_index, y = PC1)) + geom_point()

#extract only the samples for which we have > 100 cells:
big_samples = unique(scRNAseq_data_CD8$sample)[table(scRNAseq_data_CD8$sample) >= 100]
df_age = unique(data.frame(age = scRNAseq_data_CD8$age, sample = scRNAseq_data_CD8$sample))
hist(unique(data.frame(age = scRNAseq_data$age, sample = scRNAseq_data$sample))['age',])

exp_data_dense <- as.matrix(exp_data)
# S


# ANALYSIS HIGH RES DECONVOLUTION -----------------------------------------
bp.res_high_res <- read_rds(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/bp.res_high_res.rds")
dim(bp.res_high_res@Post.ini.cs@theta)
bp.res_high_res@Post.ini.cs@theta[1:3,1:3]
LUAD_RnaseqSE <- readRDS(file = "C:\\Users\\ShenorLab\\Documents\\MDPhD\\oncology\\RNAseqTCGACOAD\\data\\LUAD-tcga_data.rds")

bp.res@Post.ini.ct@theta
plot(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")), 
     bp.res_high_res@Post.ini.cs@theta["tumor",], xlab = "tumor purity by TCGA", 
     ylab = "fraction of tumor cells by BayesPrism", main = "Correlation between tumor purity and estimated tumor fraction")
plot(bp.res@Post.ini.cs@theta["Epithelial cells",], 
     bp.res_high_res@Post.ini.cs@theta["tumor",], xlab = "fraction of tumor cells low res", 
     ylab = "fraction of tumor cells high res", main = "Correlation between high and low res deconvolution")

summary(lm(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")) ~
             bp.res_high_res@Post.ini.cs@theta["tumor",]))
cor(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")),
    bp.res_high_res@Post.ini.cs@theta["tumor",], use="complete.obs", method = "spearman")

cell_types <- rownames(bp.res_high_res@Post.ini.cs@theta)

annot_rows <- data.frame(pca_res[,c("definition", "paper_Age.at.diagnosis")])
annot_rows$paper_Age.at.diagnosis <- as.numeric(as.character(annot_rows$paper_Age.at.diagnosis))
pheatmap(pca_res[,cell_types], annotation_row = annot_rows, show_rownames = F)

pca_res = data.frame(prcomp(t(bp.res_high_res@Post.ini.cs@theta), scale. = T)$x[,c("PC1","PC2")], colData(LUAD_RnaseqSE),
                     t(bp.res_high_res@Post.ini.cs@theta), check.names = F)
ggplot(pca_res, aes(x = PC1, y = PC2, color = sample_type)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = tumor)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = T lymphocytes)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Fibroblasts)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Myeloid.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = B.lymphocytes)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = NK.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Endothelial.cells)) + geom_point()

cell_types_immune <- c("Follicular B cells", "CD4+ Th", "CD8 low T", "Monocytes", "MAST", "Exhausted CD8+ T", "Cytotoxic CD8+ T",
                       "CD8+/CD4+ Mixed Th", "mo-Mac", "Treg","Naive CD4+ T","Naive CD8+ T", "MALT B cells", "NK","CD1c+ DCs",
                       "Exhausted Tfh", "Plasma cells", "Activated DCs", "CD141+ DCs", "pDCs")
CD8_T_cells <- c("Cytotoxic CD8+ T", "Naive CD8+ T", "Exhausted CD8+ T", "CD8 low T")
normalized_cell_type_immune <- t(apply(pca_res[,CD8_T_cells],1,function(x){return(x/sum(x))}))
colnames(normalized_cell_type_immune) <- paste0("norm_",CD8_T_cells)
pca_res <- cbind(pca_res, normalized_cell_type_immune)

sapply(colnames(normalized_cell_type_immune), function(cell_type){
  lmres <- summary(lm(pca_res[pca_res$definition == "Primary solid Tumor",cell_type] ~ pca_res[pca_res$definition == "Primary solid Tumor",
                                                                                               "age_at_diagnosis"]))
  return(data.frame(cell_type = cell_type, p = coef(lmres)[2,4]))
})

abun_cell_type_age <- melt(pca_res[pca_res$definition == "Primary solid Tumor", c(colnames(normalized_cell_type_immune), 'age_at_diagnosis')],
                           id.vars = 'age_at_diagnosis')

ggplot(abun_cell_type_age, aes(x = age_at_diagnosis, y = value)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~ variable, scales = "free")

