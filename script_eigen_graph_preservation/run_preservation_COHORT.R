library(ggplot2)
library(corrplot)
library(reshape2)
library(svglite)

args <- commandArgs(trailingOnly = TRUE)

#g1 <- 'EB_0'
#g2 <- 'EB_7'

g1 <- args[1]
g2 <- args[2]
new_path <- args[3]
cluster_order <- readLines("cluster.order")
setwd(new_path)
prefix<-paste(g1, g2, sep = ".")
# svg(paste(prefix, ".eigengene.svg", sep = ""),width = 24,height = 5)

#prefix <- args[1]
#g1 <- strsplit(prefix, '\\.')[[1]][1]
#g2 <- strsplit(prefix, '\\.')[[1]][2]
#paste(g1, g2, sep = ".")


fill_matrix <- function(df, clu_order) {
  # 补齐行列
  missing_names <- setdiff(clu_order, rownames(df))
  if (length(missing_names) > 0) {
    df <- rbind(df, matrix(0, nrow = length(missing_names), ncol = ncol(df), dimnames = list(missing_names, colnames(df))))
    df <- cbind(df, matrix(0, nrow = nrow(df), ncol = length(missing_names), dimnames = list(rownames(df), missing_names)))
  }
  return(df[clu_order, clu_order])
}

#clu_order<-c('S1_C1','S1_C2','S1_C3','S1_C4','S1_C5','S1_C8','S1_C9','S1_C10','S1_C14','S1_C15','S1_C16','S1_C17','S1_C20','S1_C24','S2_C4','S2_C5','S3_C1','S3_C3','S5_C2','S6_C1','S6_C3')




eigen1<-read.table(paste(g1, ".eigengene_cor.tsv", sep = ""),header=T,row.names = 1,sep = "\t")
eigen2<-read.table(paste(g2, ".eigengene_cor.tsv", sep = ""),header=T,row.names = 1,sep = "\t")
eigen1<-as.matrix(eigen1)
eigen2<-as.matrix(eigen2)



all_modules <- sort(unique(rownames(eigen1), rownames(eigen2)))
clu_order <- c(cluster_order[cluster_order %in% all_modules],
               all_modules[!all_modules %in% cluster_order])


eigen1_filled <- fill_matrix(eigen1, clu_order)
eigen2_filled <- fill_matrix(eigen2, clu_order)


svg(paste(g1, ".eigengene_cor.svg", sep = ""),width = 6,height = 5)
corrplot(eigen1_filled ,method = 'square', order = 'AOE',  tl.pos = "lt", tl.col = "black", tl.cex = 1,  number.cex=0.8, col = COL2("PRGn",10))
dev.off()


svg(paste(g2, ".eigengene_cor.svg", sep = ""),width = 6,height = 5)
corrplot(eigen2_filled ,method = 'square', order = 'AOE',  tl.pos = "lt", tl.col = "black", tl.cex = 1,  number.cex=0.8, col = COL2("PRGn",10))
dev.off()

mat<-read.table(paste(prefix, ".preserv_matrix.tsv", sep = ""),header=T,row.names = 1,sep = "\t")
mat<-as.matrix(mat)

corr_matrix_filled <- fill_matrix(mat, clu_order)
svg(paste(prefix, ".preserv_matrix.svg", sep = ""),width = 6,height = 5)
corrplot.mixed(corr_matrix_filled, upper = 'pie', lower = "ellipse", tl.pos = "lt", tl.col = "black", tl.cex = 1,  number.cex=0.8,upper.col = COL2("PRGn",10),lower.col = COL2("PRGn"))
dev.off()



mat_long<-melt(mat,varnames = c('Cluster','to_cluster'),value.name = 'Preservation')
mat_long <- mat_long[mat_long$Cluster != mat_long$to_cluster, ]
mat_long$Cluster<-factor(mat_long$Cluster,levels = clu_order)

colors_map<-colorRampPalette(c("#85D696",  "#FFC68A","#CC99DB"))(length(clu_order))

#ggplot(mat_long, aes(x= Preservation, y=Cluster,fill=Cluster)) + 
#  geom_density_ridges(alpha = 0.8,jittered_points = TRUE, point_shape = "|", point_size = 2,  position = position_points_jitter(height = 0), scale = 1.2)+ #theme_minimal() + scale_fill_manual(values=colors_map)+
#  theme(legend.position = "none")+  theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlim(-1,1)



svg(paste(prefix, ".preserv_density.svg", sep = ""),width = 6,height = 5)


n <- nrow(mat)
density <- sum(mat) / (n * (n - 1))

#print(density)
ggplot(mat_long, aes(x=Cluster, y=Preservation,fill=Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_minimal() +
  scale_fill_manual(values=colors_map) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) +
  ggtitle(sprintf("Preservation Network Density: D = %.4f", density))

dev.off()
print(paste(prefix,density,sep = ":"))