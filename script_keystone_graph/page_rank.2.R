
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

prefix <- args[1]
df <- read.table(paste(prefix, ".group_PR.tsv", sep = ""), sep = "\t", header = TRUE)
sig <- read.table(paste(prefix, ".diff.tsv", sep = ""),header = T,sep = "\t")



# set Taxa as factor，PR_case order descending
df$Taxa <- factor(df$Taxa, levels=df$Taxa[order(df$PR_case, decreasing = TRUE)])
sig$Taxa<-factor(sig$Taxa,levels = levels(df$Taxa))

# Find the range for y
y_min <- min(c(df$PR_control, df$PR_case))
y_max <- max(c(df$PR_control, df$PR_case))

df_long <- df %>%
  pivot_longer(cols = c(PR_control, PR_case), names_to = "Condition", values_to = "Value") %>%
  mutate(Taxa = as.factor(Taxa))

# add mark
df_long <- df_long %>%
  mutate(
    Mark = case_when(
      eigen_control == 1 & Condition == "PR_control" ~ "*",
      eigen_control == 2 & Condition == "PR_control" ~ "**",
      eigen_case == 1 & Condition == "PR_case" ~ "*",
      eigen_case == 2 & Condition == "PR_case" ~ "**",
      TRUE ~ ""
    ),
    MarkColor = case_when(
      Condition == "PR_control" ~ "#3784bb",
      Condition == "PR_case" ~ "#c84940",
      TRUE ~ NA_character_
    )
  )


p1<-ggplot(df_long, aes(x = Taxa, y = Value, fill = Condition, group = Condition)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
  coord_cartesian(ylim = c(y_min - (y_max - y_min) * 0.05, y_max*1.05)) +
  geom_text(aes(label = Mark, color = MarkColor), vjust = 0, alpha = 0.7) +
  scale_color_identity() +
  scale_fill_manual(values = c("PR_control" = "#3784bb", "PR_case" = "#c84940")) + 
  theme_bw() + labs(x = "Taxa", y = "Page Rank Score", fill = "Condition") +
  ggtitle(paste("Comparison of Keystone species and cluster in", prefix))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())  + theme(legend.position = "none")


keystone_tax_list <- readLines(paste(prefix, ".keystone.species.list", sep = ""))
df$IsKeystone <- ifelse(df$Taxa %in% keystone_tax_list, TRUE, FALSE)




res_counts <- sig %>% 
  filter(Feature == "Cluster") %>%
  count(Res) %>%
  filter(n > 5) %>%
  pull(Res)

color_palette <- colorRampPalette(c("#176B87","#86B6F6", "#B4D4FF"))(length(res_counts))

color_map <- c("control" = "#3784bb", "case" = "#c84940", "NS" = "grey", "Firmicutes"="#85D696","Proteobacteria"="#7EC9CE","Actinobacteria"="#FFC68A","Bacteroidetes"="#999FDB","Fusobacteria"="#C0E481","Euryarchaeota"="#CC99DB")
color_map[res_counts] <- color_palette

color_map['Other'] = 'grey'

sig$Feature<-factor(sig$Feature,levels = rev(c("Abundance","Phylum","Cluster")))


p2<-sig %>%
  mutate(Res_mapped = ifelse(Feature == "Cluster" & !Res %in% res_counts, "Other", Res)) %>% ggplot( aes(Taxa, Feature, colour = Res_mapped)) +
  geom_point(size = 3) +
  scale_color_manual(values = color_map,breaks = names(color_map), labels = names(color_map)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(colour = guide_legend(override.aes = list(size = 3), title = "Res")) + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8))                   #face = ifelse(sig$Taxa %in% keystone_tax_list, "bold", "plain"), color = ifelse(sig$Taxa %in% keystone_tax_list, "red", "black"))) 


combined_plot<-ggarrange(p1,p2, heights = c(1, 1),
          ncol = 1, nrow = 2, align = "v")
#combined_plot

ggsave(filename = paste(prefix,".PR.svg",sep=""),
       plot = combined_plot,
       width = 12,
       height = 4,
       units = "in", 
       dpi = 300)

