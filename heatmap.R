library(readr)
xp1 <- read_delim("expression/8_dmso_vs_8_rot50_12.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
xp2 <- read_delim("expression/8_dmso_vs_8_rot50_24.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
xp3 <- read_delim("expression/8_dmso_vs_8_rot100_24.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
xp4 <- read_delim("expression/15_dmso_vs_15_rot50_12.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
xp5<- read_delim("expression/15_dmso_vs_15_rot50_24.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
xp6 <- read_delim("expression/15_dmso_vs_15_rot100_24.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
  #Loading GEO2R data sets

install.packages("gdata")
library(gdata)
  #gdata has the cbindX function, allowing cbind with differing number of rows (rows are extended and filled with NA)
df <- cbindX( xp1[,"GENE_SYMBOL"],xp1[,"logFC"], xp2[,"logFC"], xp3[,"logFC"], xp4[,"logFC"], xp5[,"logFC"], xp6[,"logFC"])
colnames(df) <- c("Gene symbol", "xp1", "xp2", "xp3", "xp4", "xp5", "xp6")
  #Creating a df containing gene symbols and the logFC of all samples

rownames(df) <- make.names(df[,1], unique = TRUE)
xpmatrix <- as.matrix(df[2:7])
  #Heatmap functions require matrix input, the previous 2 lines convert the df to a matrix with the gene names as row names
heatmap(xpmatrix, Rowv=NA, Colv=NA)

