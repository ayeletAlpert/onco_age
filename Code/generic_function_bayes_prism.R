rm(list = ls())
memory.limit(size = 100000)
# Load the required libraries
library(readr)
library(Matrix)
library(Seurat)
library(Biobase)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
library(BayesPrism)
library(devtools)
library(matrixStats)
library(SummarizedExperiment)
library(InstaPrism)
library(survival)
library(TCGAbiolinks)
library(estimate)
theme_set(theme_bw())

CANCER_TYPES <- c("BRCA","CHOL","COAD","HNSC" ,"KIRC" ,"OV" ,  "PAAD" ,"PRAD" ,"SARC")
# CANCER_TYPES <- c("PRAD" ,"SARC")

for (CANCER_TYPE in CANCER_TYPES) {
  # process bulk TCGA data --------------------------------------------------
  # Define your data download directory
 
  # data_dir_bulk <- file.path("C:/Users/ShenorLab/Documents/MDPhD/oncology/temp",CANCER_TYPE)
  data_dir_bulk <- file.path("D:/Ayelet/bulk_TCGA",CANCER_TYPE)
  
  # Check if the directory exists
  if (!dir.exists(data_dir_bulk)) {
    # If it doesn't exist, create the directory
    dir.create(data_dir_bulk, recursive = TRUE)
    cat(paste("Directory", data_dir_bulk, "created.\n"))
  } else {
    cat(paste("Directory", data_dir_bulk, "already exists.\n"))
  }
  files_in_dir <- list.files(data_dir_bulk, full.names = TRUE)
  
  if (length(files_in_dir)<3)
  {
  # Specify the project and other query parameters
  project <- paste0("TCGA-", CANCER_TYPE)
  
  # Create a query specifying the RNA-Seq data with the "STAR - Counts" workflow
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification")
  
  # Check if download_dir exists
  if (dir.exists(data_dir_bulk)) {
    # List files in the directory
    files_in_dir <- list.files(data_dir_bulk, full.names = TRUE)
    
    # Check if the directory contains at least one file
    if (length(files_in_dir) == 0) {
      # Run GDCdownload
      GDCdownload(query, method = "api", directory = data_dir_bulk, files.per.chunk = 50)
    } else {
      cat("The directory is already exist\n")
    }
  } else {
    cat("The directory does not exist.\n")
  }
  

    # Use GDCprepare to prepare the data for the specified case and data type
    RnaseqSE <- GDCprepare(query = query, directory = data_dir_bulk)
    
    saveRDS(RnaseqSE, file = file.path(data_dir_bulk, paste0(CANCER_TYPE, "_TCGA.rds")))
    
    count_data <- assays(RnaseqSE)$unstranded
    rownames(count_data) <- str_replace(rownames(count_data), "\\.[0-9]+","")
    saveRDS(count_data, file = file.path(data_dir_bulk, paste0(CANCER_TYPE, "_count_data.rds")))  
  }
  
  # read scRNAseq data ------------------------------------------------------
  data_dir_sc_data <- file.path("D:/Ayelet/scRNASeq_Oncology/scRNAseq", CANCER_TYPE, "unzipped_files")
  
  # Check if the directory exists
  if (!dir.exists(data_dir_sc_data)) {
    # If it doesn't exist, create the directory
    dir.create(data_dir_sc_data, recursive = TRUE)
    cat(paste("Directory", data_dir_sc_data, "created.\n"))
  } else {
    cat(paste("Directory", data_dir_sc_data, "already exists.\n"))
  }
  
  #unzip scRNAseq files:
  zipped_files <- list.files(file.path("D:/Ayelet/scRNASeq_Oncology/scRNAseq", CANCER_TYPE), pattern = "\\.zip$", full.names = TRUE)
  
  for(dataset in zipped_files){
    zipped_files_name <- str_extract(dataset, "[A-Za-z]+[0-9]+")
    dir.create(file.path(data_dir_sc_data, zipped_files_name), recursive = TRUE)
    
    unzip(dataset, exdir = file.path(data_dir_sc_data, zipped_files_name))
  }
  
  # get the relevant scRNAseq datasets --------------------------------------
  datasets <- list.files(data_dir_sc_data, full.names = TRUE)
  for(dataset in datasets){
    cellular_annotations <- read.csv(file = file.path(dataset, "Cells.csv"))
    
    if(sum(str_detect(list.files(dataset), "counts")) == 0 | sum(cellular_annotations$cell_type == "Malignant") ==  0){
      cat("no count data or malignant cell ID")
    }else{
      
      expression_matrix <- readMM(list.files(dataset, full.names = T)[str_detect(list.files(dataset), "counts")])
      # expression_matrix <- readMM(file.path(dataset, "Exp_data_UMIcounts.mtx"))
      
      #assign cell names to columns:
      colnames(expression_matrix) <- cellular_annotations[,"cell_name"]
      
      cellular_annotations <- cellular_annotations[cellular_annotations$cell_type != "",]
      expression_matrix <- expression_matrix[,cellular_annotations[,"cell_name"]]
      
      #assign gene names to rows:
      gene_names <- read.table(file = file.path(dataset, "Genes.txt"))
      rownames(expression_matrix) <- gene_names[,1]
      
      #filter the genes and the cells that are completely 0:
      expression_matrix <- expression_matrix[rowSums(expression_matrix) > 0,]
      expression_matrix <- expression_matrix[,colSums(expression_matrix) > 0]
      
      #functions for translation gene symbols:
      res <- mapIds(org.Hs.eg.db, keys <- row.names(expression_matrix), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
      sym2ensembl <- function(sym){return(res[sym])}
      res_op <- names(res)
      names(res_op) <- res
      ensembl2sym <- function(ensembl){return(res_op[ensembl])}
      
      #translate the ensemble genes to gene symbols:
      count_data <- readRDS(file = file.path(data_dir_bulk, paste0(CANCER_TYPE, "_count_data.rds")))
      rownames(count_data) <- ensembl2sym(rownames(count_data))
      count_data <- count_data[!is.na(rownames(count_data)),]
      
      #get the shared genes in the bulk and scRNAseq:
      shared_genes <- intersect(rownames(expression_matrix), rownames(count_data))
      count_data <- count_data[shared_genes,]
      expression_matrix <- as.matrix(expression_matrix[shared_genes,])
      
      #transpose expression data of both bulk and sc data:
      count_data <- t(count_data)
      expression_matrix <- t(expression_matrix)
      
      #get annotations:
      cell.type.labels <- cellular_annotations[,'cell_type']
      cell.state.labels <- cellular_annotations[,'cell_type']
      
      sc.dat.filtered <- cleanup.genes (input=expression_matrix,
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
      
      bp_res <- InstaPrism(input_type = 'prism', prismObj = myPrism, convergence.plot = TRUE, n.iter = 1000)
      
      # bp_res <- InstaPrism(input_type = 'prism',prismObj = myPrism)
      
      # Check if the directory exists
      if (!dir.exists(dataset)) {
        # If it doesn't exist, create the directory
        dir.create(dataset, recursive = TRUE)
        cat(paste("Directory", dataset, "created.\n"))
      } else {
        cat(paste("Directory", dataset, "already exists.\n"))
      }
      
      saveRDS(bp_res,file = file.path(dataset, "bp_res.rds"))
    }
  }
}

# -----------------------------------------------------
# CHECK THE CORREALTION BETWEEN THE PERCENTAGE OF MALIGNANT CELLS AND TUMOR PURITY AS DEFINED BY TCGA DATABASE
# -----------------------------------------------------
parent_dir_scRNAseq <- "D:/Ayelet/scRNASeq_Oncology/scRNAseq"
parent_dir_bulk <- "D:/Ayelet/bulk_TCGA"
#process data purity from TCGA from the manuscript: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4671203/#S1
TCGA_purity <- read.csv(file = "D:/Ayelet/bulk_TCGA/estimate_purity_TCGA_samples.csv")
rownames(TCGA_purity) <- TCGA_purity$Sample.ID

#correlate the tumor purity obtained by different methodologies with the one obtained by bayes prism:
bpres_files <- list.files(parent_dir_scRNAseq, pattern = "bp_res.rds", recursive = TRUE, full.names = TRUE)
cor_methodologies_malignant_cells <- do.call('rbind', lapply(bpres_files, function(bp_res_file){
  print(bp_res_file)
  
  #extract the relevant cancer type:
  path_elements <- unlist(str_split(str_replace(bp_res_file, parent_dir_scRNAseq, ""), "/"))
  cancer_type <- path_elements[which(path_elements != "")[1]]
  data_set <- path_elements[which(path_elements != "")[3]]
  
  #exclude sarcoma samples:
  if(cancer_type == "SARC"){return(NULL)}
  
  #get the bp_res data:
  bp_res <- readRDS(file = bp_res_file)
  percentage_malig_cells <- bp_res@Post.ini.cs@theta["Malignant",]
  
  if(sum(gsub("^((?:[^-]+-){3}[^-]+).*", "\\1", names(percentage_malig_cells), perl = TRUE) %in% TCGA_purity[,"Sample.ID"]) == 0){
    ESTIMATE_purity <- read.table(file = file.path(parent_dir_bulk, paste0("ESTIMATE_SCORE_", cancer_type, ".txt")), header = T, row.names = 1)
    rownames(ESTIMATE_purity) <- paste0(rownames(ESTIMATE_purity), "A")
    return(data.frame(cancer_type = cancer_type,
                      data_set = data_set,
                      num_cells = nrow(bp_res@Post.ini.cs@theta),
                      method = "ESTIMATE",
                      corr_prop = cor(ESTIMATE_purity[gsub("^((?:[^-]+-){3}[^-]+).*", "\\1", names(percentage_malig_cells), perl = TRUE), "ESTIMATE_score"], 
                                      bp_res@Post.ini.cs@theta["Malignant",], use = "complete.obs", method = "spearman")))
    
  }else{
    #get the tumor purity obtained by the different methodologies:
    methodologies <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC")
    cor_methodologies <- do.call('rbind', lapply(methodologies, function(method){
      print(method)
      if(sum(!is.na(TCGA_purity[gsub("^((?:[^-]+-){3}[^-]+).*", "\\1", names(percentage_malig_cells), perl = TRUE),method])) == 0){return(NULL)}
      return(data.frame(cancer_type = cancer_type,
                        data_set = data_set,
                        num_cells = nrow(bp_res@Post.ini.cs@theta),
                        method = method,
                        corr_prop = cor(TCGA_purity[gsub("^((?:[^-]+-){3}[^-]+).*", "\\1", names(percentage_malig_cells), perl = TRUE),method], 
                                        bp_res@Post.ini.cs@theta["Malignant",], use = "complete.obs", method = "spearman")))
    }))
  }
}))

#choose per cancer type the optimal deconvolution scRNAseq dataset:
cancer_types <- unique(cor_methodologies_malignant_cells$cancer_type)
optimal_datasets_per_cancer_type <- do.call('rbind', lapply(cancer_types, function(cancer){
  subset_cor_values <- subset(cor_methodologies_malignant_cells, method == "ESTIMATE" & cancer_type == cancer)
  return(subset_cor_values[which.max(abs(subset_cor_values$corr_prop)),])
}))

#get the cell types chosen per optimal dataset:
cell_types_per_cancer <- lapply(cancer_types, function(cancer){
  chosen_dataset <- optimal_datasets_per_cancer_type$data_set[optimal_datasets_per_cancer_type$cancer_type == cancer]
  
  #get the bp_res file:
  bp_res <- readRDS(file = file.path(parent_dir_scRNAseq, cancer, "unzipped_files", chosen_dataset, "bp_res.rds"))
  
  return(rownames(bp_res@Post.ini.cs@theta))
})
names(cell_types_per_cancer) <- cancer_types

#association with age & outcome:
cancer_types <- c("COAD","HNSC","KIRC","OV")

p_imm_cell_age_surv <- do.call('rbind', lapply(cancer_types, function(cancer){
  
  print(cancer)
  
  #read bulk data:
  bulk_data <- readRDS(RnaseqSE, file = file.path("D:/Ayelet/bulk_TCGA", cancer, paste0(cancer, "_TCGA.rds")))
  
  #extract only the relevant samples (non-healthy, MSS):
  if(sum(colnames(colData(bulk_data)) == "paper_MSI_status") > 0){
    rel_samples <- colnames(bulk_data)[bulk_data$definition != "Solid Tissue Normal" & !(is.na(bulk_data$paper_MSI_status) | bulk_data$paper_MSI_status == "MSI-H") ]
  }else{
    rel_samples <- colnames(bulk_data)[bulk_data$definition != "Solid Tissue Normal"]
  }
  
  #read bp_res:
  chosen_dataset <- optimal_datasets_per_cancer_type$data_set[optimal_datasets_per_cancer_type$cancer_type == cancer]
  bp_res <- readRDS(file = file.path(parent_dir_scRNAseq, cancer, "unzipped_files", chosen_dataset, "bp_res.rds"))
  
  #calculate the frequency of immune cells alone:
  IMM_CELLS <- c("T_cell", "Macrophage", "B_cell", "NK_cell", "Myeloid", "Monocyte", "Dendritic", "Lymphoid")
  freq_all_cells <- bp_res@Post.ini.cs@theta[,rel_samples]
  freq_imm_cells <- bp_res@Post.ini.cs@theta[rownames(bp_res@Post.ini.cs@theta) %in% IMM_CELLS, rel_samples]
  freq_imm_cells <- apply(freq_imm_cells,2,function(x){return(x/sum(x))})
  
  #ratios between cell types:
  imm_cells_comb <- t(combn(rownames(freq_imm_cells),2))
  freq_imm_cells_ratio <- matrix(NA, nrow(imm_cells_comb), ncol(freq_imm_cells))
  for(i in 1:nrow(imm_cells_comb)){
    freq_imm_cells_ratio[i,] <- freq_imm_cells[imm_cells_comb[i,1],]/freq_imm_cells[imm_cells_comb[i,2],]
  }
  rownames(freq_imm_cells_ratio) <- paste0(imm_cells_comb[,1], "_", imm_cells_comb[,2])
  
  #association of immune cell types with age:
  p_vals_age_imm <- apply(freq_imm_cells, 1, function(x){
    lm_res <- summary(lm(log(x) ~ colData(bulk_data)[rel_samples,'age_at_diagnosis']))
    return(coef(lm_res)[2,4])
  })
  
  p_vals_age_ratio <- apply(freq_imm_cells_ratio, 1, function(x){
    lm_res <- summary(lm(x ~ colData(bulk_data)[rel_samples,'age_at_diagnosis']))
    return(coef(lm_res)[2,4])
  })
  
  # correlation of cell type abundance with survival ------------------------
  # stage_col_name <- colnames(colData(bulk_data))[str_detect(colnames(colData(bulk_data)), "stage") & !(str_detect(colnames(colData(bulk_data)), "paper"))]
  stage_col_name <- colnames(colData(bulk_data))[str_detect(colnames(colData(bulk_data)), "ajcc_pathologic_stage|figo_stage|primary_gleason_grade")]
  surv_data <- cbind(colData(bulk_data)[rel_samples, c("days_to_last_follow_up", "vital_status", "days_to_death", stage_col_name, "age_at_diagnosis")],
                     t(freq_imm_cells), t(freq_imm_cells_ratio))
  surv_data$stage_short <- str_extract(surv_data[,stage_col_name], "[I]+V*")
  surv_data$surv_state <- surv_data$vital_status == "Dead"
  surv_data$TTE <- apply(surv_data,1,function(x){return(max(as.numeric(x['days_to_death']), as.numeric(x['days_to_last_follow_up']), na.rm = T))})
  
  sum(is.infinite(surv_data[,c("surv_state")]))
  surv_data <- surv_data[!is.na(surv_data$surv_state) & !is.na(surv_data$TTE) & !is.na(surv_data$stage_short) & !is.infinite(surv_data$TTE),]
  
  p_abundance_cell_types_survival <- do.call('rbind', lapply(names(p_vals_age_imm), function(cell_type){
    print(cell_type)
    cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$stage_short + surv_data$age_at_diagnosis + surv_data[,cell_type])
    cox_res <- summary(cox)
    return(data.frame(cell_type = cell_type, p_surv = coef(cox_res)[nrow(coef(cox_res)),5], dir = sign(coef(cox_res)[nrow(coef(cox_res)),1])))
  }))
  
  p_abundance_cell_types_survival$p_age <- p_vals_age_imm
  
  #ratios:
  p_abundance_cell_types_ratio_survival <- do.call('rbind', lapply(names(p_vals_age_ratio), function(cell_type_ratio){
    print(cell_type)
    cox <- coxph(Surv(surv_data$TTE, surv_data$surv_state) ~ surv_data$stage_short + surv_data$age_at_diagnosis + surv_data[,cell_type_ratio])
    cox_res <- summary(cox)
    return(data.frame(cell_type = cell_type_ratio, p_surv = coef(cox_res)[nrow(coef(cox_res)),5], dir = sign(coef(cox_res)[nrow(coef(cox_res)),1])))
  }))
  
  p_abundance_cell_types_ratio_survival$p_age <- p_vals_age_ratio
  
  return(data.frame(cancer_type = cancer, rbind(p_abundance_cell_types_survival, p_abundance_cell_types_ratio_survival)))
}))
