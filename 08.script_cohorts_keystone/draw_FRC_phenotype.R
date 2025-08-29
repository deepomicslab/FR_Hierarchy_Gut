df <- read.table("age_FRC.corrlation.tsv",sep= "\t",header = T)
df$phenotype <- factor(df$phenotype)
df$FRC <- factor(df$FRC)
ggplot(df, aes(FRC, phenotype)) +
  geom_point(aes(size = -log(pvalue), color =  correlation)) + 
  geom_point(aes(size = -log(pvalue)),shape=21) +
  scale_size_continuous(range = c(2,8))+
  scale_colour_gradientn(
    colours = colorRampPalette(c("#1f77b4", "grey", "#ff8419"))(50)
  ) +  theme_bw() + theme(axis.text.x = element_text(angle = 45 , hjust = 1)) 
  


df <- read.table("bmi_FRC.corrlation.tsv",sep= "\t",header = T)
df$phenotype <- factor(df$phenotype)
df$FRC <- factor(df$FRC)
ggplot(df, aes(FRC, phenotype)) +
  geom_point(aes(size = -log(pvalue), color =  correlation)) + 
  geom_point(aes(size = -log(pvalue)),shape=21) +
  scale_size_continuous(range = c(2,8))+
  scale_colour_gradientn(
    colours = colorRampPalette(c("#1f77b4", "grey", "#ff8419"))(50)
  ) +  theme_bw() + theme(axis.text.x = element_text(angle = 45 , hjust = 1)) 


df <- read.table("other_FRC.corrlation.tsv",sep= "\t",header = T)
df$phenotype <- factor(df$phenotype, levels = sort(unique(df$phenotype)))
df$FRC <- factor(df$FRC)
ggplot(df, aes(FRC, phenotype)) +
  geom_point(aes(size = -log(pvalue), color =  correlation)) + 
  geom_point(aes(size = -log(pvalue)),shape=21) +
  scale_size_continuous(range = c(2,8))+
  scale_colour_gradientn(
    colours = colorRampPalette(c("#1f77b4", "grey", "#ff8419"))(50)
  ) +  theme_bw() + theme(axis.text.x = element_text(angle = 45 , hjust = 1)) 
