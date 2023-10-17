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
data_dir <- "~/MDPhD/oncology/scRNAseqLUAD/Data"


# CELLULAR ANNOTATIONS ----------------------------------------------------
Kim_Lee_cell_annot <- read.table(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/GSE131907_Lung_Cancer_cell_annotation.txt", skip = 1, sep = "\t")
colnames(Kim_Lee_cell_annot) <- as.character(read.table(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/GSE131907_Lung_Cancer_cell_annotation.txt", nrows = 1))
Kim_Lee_cell_annot <- Kim_Lee_cell_annot[Kim_Lee_cell_annot$Index %in% colnames(scRNAseq_data_Kim_Lee),]
rownames(Kim_Lee_cell_annot) <- Kim_Lee_cell_annot$Index
Kim_Lee_cell_annot <- Kim_Lee_cell_annot[colnames(scRNAseq_data_Kim_Lee),]
Kim_Lee_cell_annot <- Kim_Lee_cell_annot[!is.na(Kim_Lee_cell_annot$Cell_subtype),]
saveRDS(Kim_Lee_cell_annot, file = "~/MDPhD/oncology/BayesPrismLUAD/Data/Kim_Lee_cell_annot.rds")

Kim_Lee_cell_annot= readRDS( file = "~/MDPhD/oncology/BayesPrismLUAD/Data/Kim_Lee_cell_annot.rds")



# scRNAseq data - LUAD - normalized data ----------------------------------
# scRNAseq_data <- readRDS(file = file.path(data_dir, "core_atlas.rds"))
# scRNAseq_data_Kim_Lee <- scRNAseq_data[,scRNAseq_data$dataset == "Kim_Lee_2020"]
# lung_samples <- paste0("Kim_Lee_2020_",c("LUNG_T06", "LUNG_T08", "LUNG_T09", "LUNG_T18", "LUNG_T19", "LUNG_T20", "LUNG_T25", "LUNG_T28", "LUNG_T30", "LUNG_T31", "LUNG_T34"))
# scRNAseq_data_Kim_Lee <- scRNAseq_data_Kim_Lee[,as.character(scRNAseq_data_Kim_Lee$sample) %in% lung_samples]
# scRNAseq_data_Kim_Lee <- FindVariableFeatures(scRNAseq_data_Kim_Lee, selection.method = "vst", nfeatures = 2000)
# colnames(scRNAseq_data_Kim_Lee) <- str_replace(colnames(scRNAseq_data_Kim_Lee), "\\-5", "")
# Kim_Lee_cell_annot <- readRDS( file = "~/MDPhD/oncology/BayesPrismLUAD/Data/Kim_Lee_cell_annot.rds")
# scRNAseq_data_Kim_Lee$cell_type <- Kim_Lee_cell_annot[,"Cell_type"]
# scRNAseq_data_Kim_Lee$cell_subtype <- Kim_Lee_cell_annot[,"Cell_subtype"]
# all.genes <- rownames(scRNAseq_data_Kim_Lee)
# scRNAseq_data_Kim_Lee <- ScaleData(scRNAseq_data_Kim_Lee, features = all.genes)
# scRNAseq_data_Kim_Lee <- RunPCA(scRNAseq_data_Kim_Lee, features = VariableFeatures(object = scRNAseq_data_Kim_Lee))
# scRNAseq_data_Kim_Lee <- FindNeighbors(scRNAseq_data_Kim_Lee, dims = 1:10)
# scRNAseq_data_Kim_Lee <- FindClusters(scRNAseq_data_Kim_Lee, resolution = 0.5)
# DimPlot(scRNAseq_data_Kim_Lee, reduction = "umap", label = T)

# saveRDS(scRNAseq_data_Kim_Lee, "~/MDPhD/oncology/scRNAseqLUAD/savedData/scRNAseq_data_Kim_Lee.rds")


#  FUNCTIONS FOR MAPPING GENE SYMBOLS TO EMSEBLE IDS ------------------------
scRNAseq_data_Kim_Lee <- readRDS(file = "~/MDPhD/oncology/scRNAseqLUAD/savedData/scRNAseq_data_Kim_Lee.rds")
res <- mapIds(org.Hs.eg.db, keys <- row.names(scRNAseq_data_Kim_Lee), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
ensembl2sym <- function(ensembl){return(res[ensembl])}
res_op <- names(res)
names(res_op) <- res
sym2ensembl <- function(sym){return(res_op[sym])}


# bulk LUAD data ----------------------------------------------------------
LUAD_RnaseqSE <- readRDS(file = "C:\\Users\\ShenorLab\\Documents\\MDPhD\\oncology\\RNAseqTCGACOAD\\data\\LUAD-tcga_data.rds")
count_data <- assays(LUAD_RnaseqSE)$unstranded
rownames(count_data) <- str_replace(rownames(count_data), "\\.[0-9]+","")
rownames(count_data) <- ensembl2sym(rownames(count_data))
saveRDS(count_data, file = "~/MDPhD/oncology/scRNAseqLUAD/savedData/bulk_count_data.rds")
count_data <- readRDS(file = "~/MDPhD/oncology/scRNAseqLUAD/savedData/bulk_count_data.rds")

# single-cell RNA seq LUAD data - counts ----------------------------------
count_filer <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
count_file_split <- count_filer[,str_detect(colnames(count_filer), "[A,C,G,T]+\\_T.")]
row_sums <- rowSums(count_file_split[, -1])
sc_dat <- count_file_split[row_sums  > 0, ]
sc_dat <- t(sc_dat)
saveRDS(sc_dat, file = "~/MDPhD/oncology/BayesPrismLUAD/Data/GSE131907_Lung_Cancer_raw_UMI_matrix_filter.rds")

sc_dat <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/GSE131907_Lung_Cancer_raw_UMI_matrix_filter.rds")

# BAYES PRISM -------------------------------------------------------------
sc_dat <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/GSE131907_Lung_Cancer_raw_UMI_matrix_filter.rds")
count_data <- readRDS(file = "~/MDPhD/oncology/scRNAseqLUAD/savedData/bulk_count_data.rds")
Kim_Lee_cell_annot <- readRDS( file = "~/MDPhD/oncology/BayesPrismLUAD/Data/Kim_Lee_cell_annot.rds")

shared_genes <- intersect(colnames(sc_dat), rownames(count_data))
count_data <- count_data[shared_genes,]
count_data <- t(count_data)
sc_dat <- sc_dat[,shared_genes]
sc_dat <- sc_dat[rownames(Kim_Lee_cell_annot),]

cell.type.labels <- Kim_Lee_cell_annot[rownames(sc_dat),'Cell_type']
cell.state.labels <- Kim_Lee_cell_annot[rownames(sc_dat),'Cell_type']

sc.dat.filtered <- cleanup.genes (input=sc_dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = count_data)

sc.dat.filtered.pc <- select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")

myPrism <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=count_data,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key="Epithelial cells",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# saveRDS(myPrism, file = "~/MDPhD/oncology/BayesPrismLUAD/Data/my_prism_luad.rds")

# BAYES PRISM - HIGHER RESOLUTION CELL TYPES -------------------------------------------------------------
sc_dat <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/GSE131907_Lung_Cancer_raw_UMI_matrix_filter.rds")
count_data <- readRDS(file = "~/MDPhD/oncology/scRNAseqLUAD/savedData/bulk_count_data.rds")
Kim_Lee_cell_annot <- readRDS( file = "~/MDPhD/oncology/BayesPrismLUAD/Data/Kim_Lee_cell_annot.rds")

shared_genes <- intersect(colnames(sc_dat), rownames(count_data))
count_data <- count_data[shared_genes,]
count_data <- t(count_data)
sc_dat <- sc_dat[,shared_genes]
sc_dat <- sc_dat[rownames(Kim_Lee_cell_annot),]

cell.type.labels <- Kim_Lee_cell_annot[rownames(sc_dat),'Cell_type']
cell.state.labels <- Kim_Lee_cell_annot[rownames(sc_dat),'Cell_type']

sc.dat.filtered <- cleanup.genes (input=sc_dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = count_data)

sc.dat.filtered.pc <- select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")

myPrism <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=count_data,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key="Epithelial cells",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# saveRDS(myPrism, file = "~/MDPhD/oncology/BayesPrismLUAD/Data/my_prism_luad.rds")

# RUN FAST BAYESPRISM -----------------------------------------------------
myPrism <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/my_prism_luad.rds")
bp.res <- InstaPrism(input_type = 'prism',prismObj = myPrism)
#bp.res <- run.prism(prism = myPrism, n.cores=8)
saveRDS(bp.res,file = "~/MDPhD/oncology/BayesPrismLUAD/Data/bp_res.rds")


# CELL SPECIFIC EXPRESSION ------------------------------------------------
bp.res <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/bp_res.rds")
cell_types <- rownames(bp.res@Post.ini.cs@theta)
for(cell_type in cell_types){
  saveRDS(reconstruct_Z_ct_initial(InstaPrism_obj = bp.res, cell.type.of.interest = cell_type), 
          file = paste0("~/MDPhD/oncology/BayesPrismLUAD/Data/cell_specific_expression/", cell_type, "_specific_expression.rds"))
  print(0)
}


# RUN BAYESPRISM WITH HIGHER RESOLUTION -----------------------------------
sc_dat <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/GSE131907_Lung_Cancer_raw_UMI_matrix_filter.rds")
count_data <- readRDS(file = "~/MDPhD/oncology/scRNAseqLUAD/savedData/bulk_count_data.rds")
Kim_Lee_cell_annot <- readRDS( file = "~/MDPhD/oncology/BayesPrismLUAD/Data/Kim_Lee_cell_annot.rds")

shared_genes <- intersect(colnames(sc_dat), rownames(count_data))
count_data <- count_data[shared_genes,]
count_data <- t(count_data)
sc_dat <- sc_dat[,shared_genes]
sc_dat <- sc_dat[rownames(Kim_Lee_cell_annot),]

cell.type.labels <- Kim_Lee_cell_annot[rownames(sc_dat),'Cell_subtype']
cell.state.labels <- Kim_Lee_cell_annot[rownames(sc_dat),'Cell_subtype']

tumor_cells <- c("tS1","tS2","tS3")
cell.type.labels[cell.type.labels %in% tumor_cells] <- "tumor"
cell.state.labels[cell.state.labels %in% tumor_cells] <- "tumor"

sc.dat.filtered <- cleanup.genes (input=sc_dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = count_data)

sc.dat.filtered.pc <- select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")

myPrism <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=count_data,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
saveRDS(myPrism, file = "~/MDPhD/oncology/BayesPrismLUAD/Data/my_prism_luad_high_res.rds")

# RUN FAST BAYESPRISM -----------------------------------------------------
myPrism <- readRDS(file = "~/MDPhD/oncology/BayesPrismLUAD/Data/my_prism_luad_high_res.rds")
bp.res_high_res <- InstaPrism(input_type = 'prism',prismObj = myPrism)
saveRDS(bp.res_high_res,file = "~/MDPhD/oncology/BayesPrismLUAD/Data/bp.res_high_res.rds")
