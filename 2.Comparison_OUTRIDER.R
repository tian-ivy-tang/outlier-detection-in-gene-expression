list.files('/s/project/sys_gen_students/2023_2024/outrider/counts')
#list.dirs()

setwd('/s/project/sys_gen_students/2023_2024/outrider')
library(OUTRIDER)

# define empty vector to record data, only run ONCE at beginning of a session
file_name_vector <- c()
time_vector <- c()
q_vector <- c()
col_nb_vector <- c()

##############get filtered count####################
# get raw data
file_name <- './counts/Geuvadis.tsv'
ctsTable <- read.table(file = file_name, check.names=FALSE, header=TRUE)

# make gene name column as rowname
rownames(ctsTable) <- ctsTable[, 1]
ctsTable <- ctsTable[, -1]
col_nb <- dim(ctsTable)[2]

# sanity check if there is no 0 count of gene -> corresponds to filterExpression()
sum(ctsTable=0)

# transform dataset format to RangedSummarizedExperiment
ods <- OutriderDataSet(countData=ctsTable)


ods <- makeExampleOutriderDataSet(dataset="Kremer")


# filter out non expressed genes
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)

#############run analysis###########################
# get start time
time_start <- Sys.time()

# find encoding dimension q first with autoencoder
implementation <- 'autoencoder'
ods <- findEncodingDim(ods, 
                       implementation = implementation)
getBestQ(ods)

# run OUTRIDER pipeline
ods <- OUTRIDER(ods, controlData=TRUE)

# get execution time in minutes
time_taken <- difftime(Sys.time(), time_start, units = 'mins')

##########get comparison stats#######################
# get best q
best_q <- getBestQ(ods)

# number of aberrant genes per sample
aberrant_genes_per_sample <- tail(sort(aberrant(ods, by="sample")))

# append file name, time, q, column number to their vectors
file_name_vector <- c(file_name_vector, substr(file_name,1,nchar(file_name)-4))
time_vector <- c(time_vector, time_taken)
q_vector <- c(q_vector, best_q)
col_nb_vector <- c(col_nb_vector, col_nb)

########export comparison stats######################
# combine vectors into a dataframe
summary_stats <- data.frame(File = file_name_vector,
                            Time = time_vector,
                            Q = q_vector,
                            Column_nb = col_nb_vector)
print(summary_stats)

# !!!add your name to csv file so you don't overwrite other people's results!!!
#write.csv(summary_stats, './results/summary_stats.csv')

# export aberrant genes per sample to csv file
write.csv(aberrant_genes_per_sample, 
          paste0(substr(file_name,1,nchar(file_name)-4),'.csv'))
