library("devtools")
remotes::install_github("Danko-Lab/BayesPrism/BayesPrism")
library(BayesPrism)
memory.limit(size = 1000000)

data_dir = "~/MDPhD/oncology/BayesPrismLUAD/Data"
load(file.path(data_dir, "tutorial.gbm.rdata"))
ls()
dim(bk.dat)
head(rownames(bk.dat))
head(colnames(bk.dat))

dim(sc.dat)
head(rownames(sc.dat))
head(colnames(sc.dat))

## TCGA: GENCODE v22

table(cell.type.labels)
table(cell.state.labels)

table(cbind.data.frame(cell.state.labels, cell.type.labels))

sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)
dim(sc.dat.filtered)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)

sc.dat.filtered.pc <- select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")
dim(sc.dat.filtered.pc)

myPrism <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=bk.dat[1:5,],
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores=50)
slotNames(bp.res)

# extract posterior mean of cell type fraction theta
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")

head(theta)
# extract coefficient of variation (CV) of cell type fraction
theta.cv <- bp.res@posterior.theta_f@theta.cv

head(theta.cv)

# extract posterior mean of cell type-specific gene expression count matrix Z  
Z.tumor <- get.exp (bp=bp.res,
                    state.or.type="type",
                    cell.name="tumor")

head(t(Z.tumor[1:5,]))

#Try the Lung carcinoma file:
#load the clstering file:
clustering_file <- read.table(file = "C:\\Users\\ShenorLab\\Downloads/E-MTAB-6653.clusters.tsv")
