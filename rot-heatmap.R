library(readr)
setwd("~/GitHub/rotenone-geo2r-heatmap")
xp1 <- read_delim("expression/8_dmso_vs_8_rot50_12.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
xp2 <- read_delim("expression/8_dmso_vs_8_rot50_24.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
xp3 <- read_delim("expression/8_dmso_vs_8_rot100_24.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
  #Loading expression datasets 

install.packages("reshape")
install.packages("gdata")
library(reshape)
library(gdata)
#gdata has the cbindX function, allowing cbind with differing number of rows (rows are extended and filled with NA)
df <- cbindX( xp1[,"GENE_SYMBOL"],xp1[,"logFC"], xp2[,"logFC"], xp3[,"logFC"])
colnames(df) <- c("Genesymbol", "xp1", "xp2", "xp3")
#Creating a df containing gene symbols and the logFC of all samples
rownames(df) <- make.names(df[,1], unique = TRUE)
xpmatrix <- as.matrix(df[2:4])
df<-melt(xpmatrix)
colnames(df) <- c("Gene","dataset","logFC")
  #Making a dataframe that can be used with ggplot2
library(ggplot2)
rotheatmap <- ggplot(df, aes(x = dataset, y = Gene, fill = logFC)) +
  geom_tile()+
  scale_fill_gradient2(low = "#0000ff",
                       mid = "#ffffff",
                       high = "#ff0000",
                       limits=c(-0.26,0.26),
                       breaks=seq(-0.26,0.26,by=0.13))
rotheatmap=rotheatmap
install.packages("svglite")
library(svglite)
ggsave(file="rotheatmap.svg", plot=rotheatmap, width=10,height=10)

