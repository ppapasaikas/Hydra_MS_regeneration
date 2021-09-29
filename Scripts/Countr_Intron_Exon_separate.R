library(eisaR)
library(GenomicFeatures)
library(QuasR)

# Preparing Annotation:
## Create TxDB for Hydra from gtf
txdb <- makeTxDbFromGFF("Hydra_vulgaris_r102/Hydra_vulgaris_ncbi_r102_processed.gtf", dataSource="ncbi", organism="Hydra vulgaris", taxonomyId=6087)
## Extract regions:
reg <- getRegionsFromTxDb(txdb = txdb, strandedData = FALSE)

## Define QuasR sample and gneome files: 
sampleFile <-  paste0(base_dir,"QuasR_SampleSheet.txt") 
genomeFile <- "/work/gbioinfo/DB/genomes/Hydra_vulgaris_r102/Hydra_vulgaris_ncbi_r102_processed.fa"
proj <- qAlign(sampleFile, genomeFile,paired='no')

## Count alignments in exons and gene bodies
clObj <- makeCluster(48)
cntEx <- qCount(proj, reg$exons, orientation = "any", clObj=clObj, collapseBySample=TRUE)
cntGb <- qCount(proj, reg$genebodies, orientation = "any", clObj=clObj, collapseBySample=TRUE)
cntIn <- cntGb - cntEx
head(cntEx)
head(cntIn)
stopCluster(clObj)

## Save Count Tables
saveRDS(cntEx,paste0(base_dir,"Analysis/RDS_files/rawcounts_exonic.rds") )
saveRDS(cntGb,paste0(base_dir,"Analysis/RDS_files/rawcounts_genebody.rds") )
saveRDS(cntIn,paste0(base_dir,"Analysis/RDS_files/rawcounts_intronic.rds") )


















