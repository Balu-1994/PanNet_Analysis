# The following code essentially follows the vignette:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Change the working directory to a local directory that contains data.
setwd("C:/Users/br23n469/Desktop/PanNet_RNASeq_Analysis")

# Import library DESeq2
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(gplots)

# Setting up color profiles from color brewer
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_1 <- "T.1"
sample_2 <- "T.2"
sample_3 <- "T.3"
sample_4 <- "T.4"


## Create Samples Data-frame
condition_column <- c(rep(sample_1, 8),
                      rep(sample_2, 8),
                      rep(sample_3, 6),
                      rep(sample_4, 3))

run_column <- c(paste(sample_1, "1", sep = "_"), paste(sample_1, "2", sep = "_"), 
                paste(sample_1, "3", sep = "_"), paste(sample_1, "4", sep = "_"), 
                paste(sample_1, "5", sep = "_"), paste(sample_1, "6", sep = "_"), 
                paste(sample_1, "7", sep = "_"), paste(sample_1, "8", sep = "_"), 
                paste(sample_2, "1", sep = "_"), paste(sample_2, "2", sep = "_"), 
                paste(sample_2, "3", sep = "_"), paste(sample_2, "4", sep = "_"), 
                paste(sample_2, "5", sep = "_"), paste(sample_2, "6", sep = "_"), 
                paste(sample_2, "7", sep = "_"), paste(sample_2, "8", sep = "_"),   
                paste(sample_3, "1", sep = "_"), paste(sample_3, "2", sep = "_"), 
                paste(sample_3, "3", sep = "_"), paste(sample_3, "4", sep = "_"), 
                paste(sample_3, "5", sep = "_"), paste(sample_3, "6", sep = "_"), 
                paste(sample_4, "1", sep = "_"), paste(sample_4, "2", sep = "_"), 
                paste(sample_4, "3", sep = "_"))

rep_column <- c("1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "2.1", 
                "2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8", "3.1", "3.2", 
                "3.3", "3.4", "3.5", "3.6", "4.1", "4.2", "4.3")


samples_df <- data.frame(condition_column,
                         run_column,
                         rep_column)

colnames(samples_df) <- c("condition", "run", "rep")

rownames(samples_df) <- samples_df$run

samples_df$condition <- factor(rep(c(rep(sample_1, 8), rep(sample_2, 8),
                                     rep(sample_3, 6), rep(sample_4, 3))))


# The table contains values for each CDS included in the sub_read_count_data  variable.
sub_read_count_data <- read.table("subreadcounts.txt", header = TRUE, row.names = 1)

# Reorder to make the order consistent with samples$run
sub_read_count_data <- sub_read_count_data[ , -c(3, 4, 1, 2,5)]

# Change colnames
colnames(sub_read_count_data) <- rownames(samples_df)

# Import as DESeq-DataSet (dds)
dds <- DESeqDataSetFromMatrix(countData = sub_read_count_data,
                              colData = samples_df,
                              design = ~ condition)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Differential expression analysis
dds <- DESeq(dds)

colData(dds) # to check whether names are correct

# Log transformation for data quality assessment
rld <- rlog(dds, blind = FALSE)


# PCA plot
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
ggplot(pcaData, aes(PC1, PC2, color = condition)) + geom_point()









