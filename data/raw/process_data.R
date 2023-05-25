# Data reading, quality check and writing in R's convenient rda format
#  - otu_HMP2_16S.csv = abundances table of OTUs
#  - meta_HMP2.csv = samples meta info
#  - taxonomy_HMP2_16S.csv = taxonomy classification of OTUs

# Read tables
otu <- read.table("data/raw/otu_HMP2_16S.csv", header=TRUE, 
                  sep=",", row.names=1)
meta <- read.table("data/raw/meta_HMP2.csv", header=TRUE,
                   sep=",", row.names=1)
taxa <- read.table("data/raw/taxonomy_HMP2_16S.csv", header=TRUE,
                   sep=",", row.names=1)

# Checks 
stopifnot(all(rownames(otu)==rownames(meta)))
stopifnot(all(colnames(otu)==rownames(taxa)))
stopifnot(all(apply(otu,c(1,2),is.numeric)))
stopifnot(all(otu>=0))

# Write Data
saveRDS(otu, "data/otu_HMP2.rds")
saveRDS(meta, "data/meta_HMP2.rds")
saveRDS(taxa, "data/taxonomy.rds")