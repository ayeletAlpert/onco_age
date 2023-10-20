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
library(TCGAbiolinks)
theme_set(theme_bw())

CANCER_TYPE <- "COAD"
# process bulk TCGA data --------------------------------------------------
# Define your data download directory
data_dir_bulk <- "Z:/Ayelet/bulk_TCGA"  # Change this to your preferred directory
#data_dir_bulk <- file.path(data_dir_bulk,CANCER_TYPE)
# Specify the project and other query parameters

project <- paste0("TCGA-", CANCER_TYPE)

# Create a query specifying the RNA-Seq data with the "STAR - Counts" workflow
query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification")

# Download the data with the correct query
GDCdownload(query,
  method = "api",
  directory = data_dir_bulk,
  files.per.chunk = 10)

# Define the directory for the specific case (e.g., TCGA-COAD)
case_dir <- file.path(data_dir_bulk, "TCGA", CANCER_TYPE)

# Use GDCprepare to prepare the data for the specified case and data type
RnaseqSE <- GDCprepare(query = query, directory = case_dir)

saveRDS(RnaseqSE, file = file.path(data_dir_bulk, paste0(CANCER_TYPE, "_TCGA.rds")))

count_data <- assays(RnaseqSE)$unstranded
rownames(count_data) <- str_replace(rownames(count_data), "\\.[0-9]+","")
saveRDS(count_data, file = file.path(data_dir_bulk, paste0(CANCER_TYPE, "_count_data.rds")))

# read scRNAseq data ------------------------------------------------------
data_dir_sc_data <- file.path("Z:/Ayelet/scRNASeq_Oncology/scRNAseq", CANCER_TYPE, "unzipped_files")

# Check if the directory exists
if (!dir.exists(data_dir_sc_data)) {
  # If it doesn't exist, create the directory
  dir.create(data_dir_sc_data, recursive = TRUE)
  cat(paste("Directory", data_dir_sc_data, "created.\n"))
} else {
  cat(paste("Directory", data_dir_sc_data, "already exists.\n"))
}

#unzip scRNAseq files:
zipped_files <- list.files(file.path("Z:/Ayelet/scRNASeq_Oncology/scRNAseq", CANCER_TYPE), pattern = "\\.zip$", full.names = TRUE)

for(dataset in zipped_files){
  zipped_files_name <- str_extract(dataset, "[A-Za-z]+[0-9]+")
  dir.create(file.path(data_dir_sc_data, zipped_files_name), recursive = TRUE)
  
  unzip(dataset, exdir = file.path(data_dir_sc_data, zipped_files_name))
}

# get the relevant scRNAseq datasets --------------------------------------
datasets <- list.files(data_dir_sc_data, full.names = TRUE)
for(dataset in datasets){
  if(sum(str_detect(list.files(dataset), "counts")) == 0){
    cat("no count data")
  }else{
    expression_matrix <- readMM(file.path(dataset, "Exp_data_UMIcounts.mtx"))
    
    cellular_annotations <- read.csv(file = file.path(dataset, "Cells.csv"))
    cellular_annotations <- cellular_annotations[cellular_annotations$cell_type != "",]
    
    #assign cell names to columns:
    colnames(expression_matrix) <- cellular_annotations[,"cell_name"]
    expression_matrix <- expression_matrix[,cellular_annotations[,"cell_name"]]
    
    #assign gene names to rows:
    gene_names <- read.table(file = file.path(dataset, "Genes.txt"))
    rownames(expression_matrix) <- gene_names[,1]
    
    #filter the genes and the cells that are completely 0:
    expression_matrix <- expression_matrix[rowSums(expression_matrix) > 0,]
    expression_matrix <- expression_matrix[,colSums(expression_matrix) > 0]
    
    #functions for translation gene symbols:
    res <- mapIds(org.Hs.eg.db, keys <- row.names(expression_matrix), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    ensembl2sym <- function(ensembl){return(res[ensembl])}
    res_op <- names(res)
    names(res_op) <- res
    sym2ensembl <- function(sym){return(res_op[sym])}
    
    #translate the ensemble genes to gene symbols:
    count_data <- readRDS(file = file.path(data_dir_bulk, paste0(CANCER_TYPE, "_count_data.rds")))
    rownames(count_data) <- ensembl2sym(rownames(count_data))
    
    #get the shared genes in the bulk and scRNAseq:
    shared_genes <- intersect(colnames(expression_matrix), rownames(count_data))
    count_data <- count_data[shared_genes,]
    expression_matrix <- expression_matrix[,shared_genes]
    
    #get annotations:
    cell.type.labels <- cellular_annotations[,'cell_type']
    cell.state.labels <- cellular_annotations[,'cell_type']
    
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
      key="Malignant",
      outlier.cut=0.01,
      outlier.fraction=0.1,
    )
    
    bp_res <- InstaPrism(input_type = 'prism',prismObj = myPrism)

    saveRDS(bp.res,file = file.path(dataset, "bp_res.rds"))
  }
}
