library(ggplot2)
library(corrplot)
library(reshape2)


setwd("../result/Anti/eigenspecies/")
args <- commandArgs(trailingOnly = TRUE)

prefix <- args[1]
g1 <- strsplit(prefix, '\\.')[[1]][1]
g2 <- strsplit(prefix, '\\.')[[1]][2]


fill_matrix <- function(df, clu_order) {
  missing_names <- setdiff(clu_order, rownames(df))
  if (length(missing_names) > 0) {
    df <- rbind(df, matrix(0, nrow = length(missing_names), ncol = ncol(df), dimnames = list(missing_names, colnames(df))))
    df <- cbind(df, matrix(0, nrow = nrow(df), ncol = length(missing_names), dimnames = list(rownames(df), missing_names)))
  }
  return(df[clu_order, clu_order])
}

clu_order<-c('S1_C1','S1_C2','S1_C3','S1_C4','S1_C5','S1_C6','S1_C8','S1_C9','S1_C10','S1_C14','S1_C15','S1_C16','S1_C17','S1_C20','S1_C22','S1_C24','S2_C3','S2_C4','S2_C5','S2_C7','S3_C1','S3_C3','S5_C2','S6_C1','S6_C3')

eigen1<-read.table(paste(prefix, ".eigenspecies_cor.", g1,".tsv", sep = ""),header=T,row.names = 1,sep = "\t")
eigen2<-read.table(paste(prefix, ".eigenspecies_cor.", g2,".tsv", sep = ""),header=T,row.names = 1,sep = "\t")
eigen1<-as.matrix(eigen1)
eigen2<-as.matrix(eigen2)
eigen1_filled <- fill_matrix(eigen1, clu_order)
eigen2_filled <- fill_matrix(eigen2, clu_order)


pdf(paste(prefix, ".eigenspecies_cor.",g1,".pdf", sep = ""),width = 6,height = 5)
corrplot(eigen1_filled ,method = 'square', order = 'AOE',  tl.pos = "lt", tl.col = "black", tl.cex = 0.8,  number.cex=0.8, col = COL2("PRGn",10))
dev.off()

pdf(paste(prefix, ".eigenspecies_cor.",g2,".pdf", sep = ""),width = 6,height = 5)
corrplot(eigen2_filled ,method = 'square', order = 'AOE',  tl.pos = "lt", tl.col = "black", tl.cex = 0.8,  number.cex=0.8, col = COL2("PRGn",10))
dev.off()

mat<-read.table(paste(prefix, ".preserv_matrix.tsv", sep = ""),header=T,row.names = 1,sep = "\t")
mat<-as.matrix(mat)

corr_matrix_filled <- fill_matrix(mat, clu_order)
pdf(paste(prefix, ".preserv_matrix.pdf", sep = ""),width = 6,height = 5)
corrplot.mixed(corr_matrix_filled, upper = 'pie', lower = "ellipse", tl.pos = "lt", tl.col = "black", tl.cex = 0.8,  number.cex=0.8,upper.col = COL2("PRGn",10),lower.col = COL2("PRGn"))
dev.off()



mat_long<-melt(mat,varnames = c('Cluster','to_cluster'),value.name = 'Preservation')
mat_long <- mat_long[mat_long$Cluster != mat_long$to_cluster, ]
mat_long$Cluster<-factor(mat_long$Cluster,levels = clu_order)

colors_map<-colorRampPalette(c("#85D696",  "#FFC68A","#CC99DB"))(length(clu_order))


pdf(paste(prefix, ".preserv_density.pdf", sep = ""),width = 6,height = 5)

n <- nrow(mat)
diag(mat) <- 0
density <- sum(mat) / (n * (n - 1))

ggplot(mat_long, aes(x=Cluster, y=Preservation,fill=Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_minimal() +
  scale_fill_manual(values=colors_map) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(sprintf("Preservation Network Density: D = %.4f", density))

dev.off()
print(paste(prefix,density,sep = ":"))
