setwd('/s/project/sys_gen_students/2023_2024/outrider')
# setwd('C:/Users/tiant/Documents/4-ws2023/Computational modelling for system genetics/Project_outlier_calling')

print("Start outrider R script")
library(saseR)
library(OUTRIDER)
library(dplyr)
# library(SummarizedExperiment)

register(MulticoreParam(workers = 30L))
###################### set method, file name, and user id #######################
method <- 'saseR'

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No file name provided", call. = FALSE)
}
file_name <- args[1]
# file_name <- './counts/Geuvadis.tsv'

# Extract the file base name without the extension (dataset name)
file_base_name <- tools::file_path_sans_ext(basename(file_name))

# Use the user's system login name as a unique identifier or a timestamp
unique_id <- Sys.getenv("LOGNAME")

# Create output filenames that include the input base name and the unique identifier
aberrant_genes_filename <- paste0('./results/aberrant_genes_', 
                                  method, 
                                  '_', file_base_name, 
                                  '_', unique_id, 
                                  '.csv')

########################## get filtered count ##################################
# get start time
# moved here so time_taken is comparable to OutSingle python script
time_start <- Sys.time()

# get raw data
ctsTable <- read.table(file = file_name, check.names=FALSE, header=TRUE)

# make gene name column as rowname of count matrix
rownames(ctsTable) <- ctsTable[, 1]
ctsTable <- ctsTable[, -1]
col_nb <- dim(ctsTable)[2]

# transform dataset format to RangedSummarizedExperiment
ods <- OutriderDataSet(countData=ctsTable)
# below line is quicker? & no need to call OUTRIDER library, better choice if no need to filter dataset
# ods <- SummarizedExperiment(assays=list(counts=as.matrix(ctsTable)))

# dataset pre-filtered already, only purpose of this line is to make time & memory comparable
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)

########################### run analysis #######################################
# no known control variable to be added to model
metadata(ods)$design <- ~1

# calculate offsets, either method is ok
ods <- calculateOffsets(ods, method='geommean')  # geometric mean by DESeq2
# ods <- calculateOffsets(ods, method='TMM')  # trimmed mean of M-values by edgeR

# calculate optimal latent factor with SVD OHT method
ods <- saseRfindEncodingDim(ods, method = "GD", analysis = 'AE')  # aberrant expression
best_q <- metadata(ods)$optimalEncDim

# fit final model with fast parameter estimation
# uses a overdispersed quadratic mean-variance relationship
ods <- saseRfit(ods, analysis = "AE", padjust = "BY", fit = "fast")

############################### get results ####################################
# export number of aberrant genes per sample
count_padj <- function(padj){
  count_outliers <- sum(padj < 0.05)
  return(count_outliers)
}
aberrant_genes_per_sample <- apply(assays(ods)$pValueAdjust, 2, count_padj)
aberrant_genes_per_sample <- arrange(data.frame(aberrant_per_sample), 
                               aberrant_per_sample)  # sort by count
geneID <- rownames(aberrant_genes_per_sample)
aberrant_genes_per_sample <- cbind(geneID, aberrant_genes_per_sample)

write.csv(aberrant_genes_per_sample, 
          aberrant_genes_filename, 
          row.names=F)

# export ourlier gene results for every dataset
# merge P values & adjusted P values & sort by adjust P values
pvals <- assays(ods)$pValue
pvals_melted <- reshape2::melt(pvals)
colnames(pvals_melted) <- c('geneID', 'sampleID', 'pValue')

pvals_adj <- assays(ods)$pValueAdjust
pvals_adj_melted <- reshape2::melt(pvals_adj)
colnames(pvals_adj_melted) <- c('geneID', 'sampleID', 'padjust')

gene_results <- merge(pvals_melted, pvals_adj_melted, by=c('geneID', 'sampleID'))
gene_results <- arrange(gene_results, padjust)

# add aberrant column TRUE/FALSE
gene_results <- gene_results %>% 
  mutate(aberrant = ifelse(padjust<0.05, TRUE, FALSE))

write.csv(gene_results, 
          paste0('./results/gene_results_', method, 
                 '_', file_base_name, 
                 '_', unique_id, 
                 '.csv'))

# get execution time in minutes
# moved end time here so time_taken is comparable to OutSingle's python script
time_taken <- round(difftime(Sys.time(), time_start, units='mins'), 2)

# export summary_stats
summary <- c(file_base_name, method, time_taken, best_q, col_nb)
print(summary)
# summary_stats <- read.csv(file = 'results/summary_stats.csv', 
#                           check.names=FALSE, header=TRUE)
# summary_stats[nrow(summary_stats) + 1,] <- c(file_base_name, 
#                                              method, 
#                                              time_taken, 
#                                              best_q, 
#                                              col_nb)
# write.csv(summary_stats, 'results/summary_stats.csv', row.names=F)
