#setwd("D:/FR/result_20240417")
library(effsize)

args <- commandArgs(trailingOnly = TRUE)


genome_module_file <- args[1]
cluster_species_file <- args[2]
out_file <- args[3]
indir <- args[4]
setwd(indir)


genome_module <- read.table(genome_module_file, header=T, row.names = 1)
cluster_species <- read.table(cluster_species_file, header=F)  

cluster_species$cluster <- cluster_species[[1]]
cluster_species$species <- cluster_species[[2]]


# 结果表
results <- data.frame(
  module=character(), 
  cluster=character(),
  n_c=numeric(),
  com_c=numeric(),
  mean_c=numeric(),
  var_c=numeric(), 
  n_not_c=numeric(),
  com_not_c=numeric(),
  mean_not_c=numeric(),
  var_not_c=numeric(),
  rank_sum=numeric(),
  effect_size=numeric(),
  pvalue=numeric(),  
  padj=numeric(),
  stringsAsFactors=F
)



# 计算效应大小



# 循环每个模块
for(m in colnames(genome_module)){
  
  # 循环每个集群
  for(c in unique(cluster_species$cluster)){
    
    # 获取这一集群下的genome
    c_genomes <- cluster_species[cluster_species$cluster==c,]$species 
    
    # 获取其他 genome
    not_c_genomes <- setdiff(rownames(genome_module), c_genomes) 
    
    # 分组获取完整度
    c_completeness <- genome_module[rownames(genome_module) %in% c_genomes, m]
    not_c_completeness <- genome_module[rownames(genome_module) %in% not_c_genomes, m]
    
    #c_completeness <- c_completeness[!is.na(c_completeness)]
    #not_c_completeness <- not_c_completeness[!is.na(not_c_completeness)]
    
    #if(sum(is.na(c_completeness)) == length(c_completeness)){
    #  print(m)
    #  print(c)
    #  next
    #}
    
    n_c <- length(c_completeness)
    n_not_c <- length(not_c_completeness)
    
    com_c <- length(c_completeness[c_completeness!=0])/n_c
    com_not_c <- length(not_c_completeness[not_c_completeness!=0])/n_not_c
    
    # Wilcoxon检验
    wt <- wilcox.test(c_completeness, not_c_completeness)
    effect_size <- cliff.delta(c_completeness, not_c_completeness)$estimate
    
    # 记录结果
    
    results <- rbind(results, data.frame(
      module=m,
      cluster=c,
      n_c, 
      com_c,
      mean_c = mean(c_completeness),
      var_c = var(c_completeness),
      n_not_c,
      com_not_c,
      mean_not_c = mean(not_c_completeness), 
      var_not_c = var(not_c_completeness),
      rank_sum = wt$statistic,
      effect_size,
      pvalue=wt$p.value  
    ))
    
  }
  
}

# FDR校正
results$padj <- p.adjust(results$pvalue, method = "BH")

# 输出

write.table(results, out_file, quote=F, row.name=F,sep = "\t")


