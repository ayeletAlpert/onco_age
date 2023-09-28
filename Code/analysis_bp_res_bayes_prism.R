
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
library("devtools")
library(matrixStats)
library(SummarizedExperiment)
library(InstaPrism)
library(survival)
theme_set(theme_bw())


# Read bp_res -------------------------------------------------------------
bp.res <- readRDS("~/MDPhD/oncology/BayesPrismLUAD/Data/bp_res.rds")

# analysis ----------------------------------------------------------------
hist(bp.res@Post.ini.cs@theta)
dim(bp.res@Post.ini.cs@theta)
bp.res@Post.ini.cs@theta[1:3,1:3]
plot(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")), 
     bp.res@Post.ini.cs@theta["Epithelial cells",], xlab = "tumor purity by TCGA", 
     ylab = "fraction of tumor cells by BayesPrism", main = "Correlation between tumor purity and estimated tumor fraction")
summary(lm(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")) ~
             bp.res@Post.ini.cs@theta["Epithelial cells",]))
cor(as.numeric(str_replace(as.character(LUAD_RnaseqSE$paper_Purity.ABSOLUTE.calls),"\\,","\\.")),
    bp.res@Post.ini.cs@theta["Epithelial cells",], use="complete.obs", method = "spearman")

#Association of cell type abundance with age:
pca_res = data.frame(prcomp(t(bp.res@Post.ini.cs@theta), scale. = T)$x[,c("PC1","PC2")], colData(LUAD_RnaseqSE),
                     t(bp.res@Post.ini.cs@theta))
ggplot(pca_res, aes(x = PC1, y = PC2, color = sample_type)) + geom_point()
pca_res+m <- melt(pca_res[c('PC1','PC2'),])
ggplot(pca_res, aes(x = PC1, y = PC2, color = Epithelial.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = T.lymphocytes)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Fibroblasts)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Myeloid.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = B.lymphocytes)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = NK.cells)) + geom_point()
ggplot(pca_res, aes(x = PC1, y = PC2, color = Endothelial.cells)) + geom_point()
cell_types <- c("Epithelial.cells", "T.lymphocytes", "Fibroblasts", "Myeloid.cells", "B.lymphocytes", "NK.cells", "Endothelial.cells")
#survival analysis with T lymphocytes:
surv_data <- pca_res[,c("days_to_last_follow_up", "vital_status", "days_to_death", "ajcc_pathologic_stage", cell_types, "sample_type")]
surv_data <- surv_data[surv_data$sample_type == "Primary Tumor",]
surv_data$stage_short <- str_extract(surv_data$ajcc_pathologic_stage, "[I]+V*")
surv_data$surv_state <- surv_data$vital_status == "Dead"
surv_data$TTE <- apply(surv_data,1,function(x){return(max(as.numeric(x['days_to_death']), as.numeric(x['days_to_last_follow_up']), na.rm = T))})

surv_data[,c("surv_state", "TTE",'stage_short')]
sum(is.infinite(surv_data[,c("surv_state")]))
surv_data <- surv_data[!is.na(surv_data$surv_state) & !is.na(surv_data$TTE) & !is.na(surv_data$stage_short) & !is.infinite(surv_data$TTE),]

# surv_data <- surv_data[!is.na(surv_data$surv_state) & !is.na(surv_data$TTE) & !is.na(surv_data$stage_short) & !is.na(surv_data$T.lymphocytes),]
surv_data$t_under_median <- as.character(surv_data$T.lymphocytes < median(surv_data$T.lymphocytes))
df <- surv_data[,c("TTE", "surv_state", "")]
cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$stage_short + surv_data$Endothelial.cells)
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
