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
theme_set(theme_bw())

# Read data structures -------------------------------------------------------------
bp.res <- readRDS("~/MDPhD/oncology/BayesPrismLUAD/Data/bp_res.rds")
LUAD_RnaseqSE <- readRDS(file = "C:\\Users\\ShenorLab\\Documents\\MDPhD\\oncology\\RNAseqTCGACOAD\\data\\LUAD-tcga_data.rds")
plot(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")), 
     bp.res@Post.ini.cs@theta["Epithelial cells",], xlab = "tumor purity by TCGA", 
     ylab = "fraction of tumor cells by BayesPrism", main = "Correlation between tumor purity and estimated tumor fraction")
summary(lm(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")) ~
             bp.res@Post.ini.cs@theta["Epithelial cells",]))
cor(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")),
    bp.res@Post.ini.cs@theta["Epithelial cells",], use = "complete.obs", method = "spearman")

# correlation of cell type abundance with age ----------------------------------------------------------------
cell_types_immune <- c("T lymphocytes", "Myeloid cells", "B lymphocytes", "NK cells")
immune_cell_types_norm <- t(apply(bp.res@Post.ini.cs@theta[cell_types_immune,],2,function(x){return(x/sum(x))}))

p_cell_type_abun_age <- do.call('rbind', lapply(colnames(immune_cell_types_norm), function(cell_type){
  lmres <- summary(lm(immune_cell_types_norm[colnames(LUAD_RnaseqSE)[LUAD_RnaseqSE$definition == "Primary solid Tumor"],cell_type] ~
                        LUAD_RnaseqSE$age_at_diagnosis[LUAD_RnaseqSE$definition == "Primary solid Tumor"]))
  return(data.frame(cell_type = cell_type, p = coef(lmres)[2,4], dir = sign(coef(lmres)[2,1])))
}))


# correlation of cell type abundance with survival ------------------------
surv_data <- cbind(colData(LUAD_RnaseqSE)[,c("days_to_last_follow_up", "vital_status", "days_to_death", "ajcc_pathologic_stage", "sample_type", "age_at_diagnosis")],
                   immune_cell_types_norm[colnames(LUAD_RnaseqSE),])
surv_data <- surv_data[surv_data$sample_type == "Primary Tumor",]
surv_data$stage_short <- str_extract(surv_data$ajcc_pathologic_stage, "[I]+V*")
surv_data$surv_state <- surv_data$vital_status == "Dead"
surv_data$TTE <- apply(surv_data,1,function(x){return(max(as.numeric(x['days_to_death']), as.numeric(x['days_to_last_follow_up']), na.rm = T))})

sum(is.infinite(surv_data[,c("surv_state")]))
surv_data <- surv_data[!is.na(surv_data$surv_state) & !is.na(surv_data$TTE) & !is.na(surv_data$stage_short) & !is.infinite(surv_data$TTE),]

p_abundance_cell_types_survival <- do.call('rbind', lapply(cell_types_immune, function(cell_type){
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

