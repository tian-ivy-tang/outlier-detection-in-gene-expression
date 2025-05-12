#Test run OUTRIDER vignette
Sys.setenv(LANGUAGE="en")
setwd('C:/Users/tiant/Documents/4-ws2023/Computational modelling for system genetics')

#install Packages
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("org.Hs.eg.db")
#install.packages('RMariaDB')
#BiocManager::install("AnnotationDbi")

#load Packages
library(OUTRIDER)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(RMariaDB)
library(AnnotationDbi)

#get data
ctsFile <- system.file('extdata', 'KremerNBaderSmall.tsv',
                       package = 'OUTRIDER')
ctsTable <- read.table(ctsFile, check.names=FALSE)
ods <- OutriderDataSet(countData = ctsTable)

#filter out non expressed genes
ods <- filterExpression(ods, minCounts=TRUE, 
                        filterGenes=TRUE)

#run full OUTRIDER pipeline (control, fit model, calculate P-values)
ods <- OUTRIDER(ods)
#Time difference of 1.269528 mins, 15 Final nb-AE loss: 4.14779369576384"

#results (only significant)
res <- results(ods)
head(res)

#example of QQ plot for most significant outlier
plotQQ(ods, res[1, geneID])

################################################
#OUTRIDER analysis in detail
################################################
#4.1 OutriderDataSet
#small testing dataset
odsSmall <- makeExampleOutriderDataSet(dataset = 'Kremer')

#full dataset from Kremer et al.
baseURL <- paste0("https://static-content.springer.com/esm/",
                  "art%3A10.1038%2Fncomms15824/MediaObjects/")
count_URL <- paste0(baseURL, "41467_2017_BFncomms15824_MOESM390_ESM.txt")
anno_URL <- paste0(baseURL, "41467_2017_BFncomms15824_MOESM397_ESM.txt")
ctsTable <- read.table(count_URL, sep="\t")
annoTable <- read.table(anno_URL, sep="\t", header=TRUE)
annoTable$sampleID <- annoTable$RNA_ID

#create OutriderDataSet object
ods <- OutriderDataSet(countData = ctsTable, 
                       colData = annoTable)

#4.2 Preprocessing
# get annotation
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
map <- select(org.Hs.eg.db, keys = keys(txdb, keytype = 'GENEID'),
              keytype = 'ENTREZID', columns = c('SYMBOL'))

try({
  library(RMariaDB)
  library(AnnotationDbi)
  con <- dbConnect(MariaDB(), host='genome-mysql.cse.ucsc.edu',
                   dbname="hg19", user='genome')
  map <- dbGetQuery(con, 'select kgId AS TXNAME, geneSymbol from kgXref')
  txdbUrl <- paste0("https://cmm.in.tum.de/public/",
                    "paper/mitoMultiOmics/ucsc.knownGenes.db")
  download.file(txdbUrl, "ucsc.knownGenes.db")
  txdb <- loadDb("ucsc.knownGenes.db")
})

# calculate FPKM values and label not expressed genes
ods <- filterExpression(ods, txdb, mapping=map,
                        filterGenes=FALSE, savefpkm=TRUE)

# display the FPKM distribution of counts.
plotFPKM(ods)

# display gene filter summary statistics
plotExpressedGenes(ods)

# do the actual subsetting based on filtering labels
ods <- ods[mcols(ods)$passedFilter,]

#4.3 Control for confounders
# heatmap of the sample correlation
# it can also annotate the clusters resulting from dendrogram
ods <- plotCountCorHeatmap(ods, 
                           colGroup=c('SEX', 'RNA_HOX_GROUP'),
                           normalized=FALSE,
                           nRowCluster=4)

# heatmap of the gene/sample expression
ods <- plotCountGeneSampleHeatmap(ods,
                            colGroup=c('SEX', 'RNA_HOX_GROUP'),
                            normalized=FALSE,
                            nRowCluster=4)

# automatically control for confounders
# we use only 3 iterations to make the vignette faster. The default is 15.
ods <- estimateSizeFactors(ods)
ods <- controlForConfounders(ods, q=21, iterations=3)

# heatmap of the sample correlation after controlling
ods <- plotCountCorHeatmap(ods,
                           normalized=TRUE,
                           colGroups=c('SEX', 'RNA_HOX_GROUP'))

#4.4 Find the right encoding dimension q
# find the optimal encoding dimension q
ods <- findEncodingDim(ods)

# visualize the hyper parameter optimization
plotEncDimSearch(ods)

#4.4.1 Exclude samples from the autoencoder fit
# set exclusion mask
sampleExclusionMask(ods) <- FALSE
sampleExclusionMask(ods[, 'MUC1365']) <- TRUE

# check which samples are excluded from the autoencoder fit
sampleExclusionMask(ods)

#4.5 Fit the negative binomial (NB) model
# fit the model when alternative methods were used in the control step
ods <- fit(ods)
hist(theta(ods))

#4.6 P-value calculation
# compute P-values (nominal & adjusted)
ods <- computePvalues(ods, 
                      alternative='two.sided',
                      method='BY')

#4.7 Z-score calculation
# compute Z-scores
ods <- computeZscores(ods)

################################################
#Results
################################################
# get result (default only significant, padj < 0.05)
res <- results(ods)
head(res)

dim(res)

# set a different significance level and filter by Z-scores
res <- results(ods,
               padjCutoff = 0.1,
               zScoreCutoff=2)
head(res)

dim(res)

################################################
#5.2 Number of aberrant genes per sample
# number of aberrant genes per sample
tail(sort(aberrant(ods, by="sample")))

tail(sort(aberrant(ods, by="gene", zScoreCutoff=1)))

# plot the aberrant events per sample
plotAberrantPerSample(ods, padjCutoff=0.05)


################################################
#5.3 Volcano plots
# MUC1344 is a diagnosed sample from Kremer et al.
plotVolcano(ods, "MUC1344", basePlot=TRUE)

#5.4 Gene level plots
# expression rank of a gene with outlier events
plotExpressionRank(ods, 'TIMMDC1', basePlot = TRUE)

# QQ plot for a given gene
plotQQ(ods, 'TIMMDC1')

# observed vs. expected gene expression
plotExpectedVsObservedCounts(ods,
                             'TIMMDC1',
                             basePlot = TRUE)
