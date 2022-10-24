library(ggplot2)                    
library(cowplot)
library(tidyverse)  

### Sequencing depth
data <- read.csv("data_summary.csv")
pdf("sample_depth.pdf",width = 4.5, height = 3)
ggplot(data.depth.merge, aes(x=group, y=depth)) + 
  geom_violin(aes(color=group,fill=group)) +
  labs(x="", y="Sequecing depth") + #theme_cowplot(font_size = 18, font_family = "", line_size = 1) +
  scale_color_brewer(palette = "Reds", direction = -1) + scale_fill_brewer(palette = "Reds", direction = -1)+
  stat_summary(fun.data =  data_summary, colour = "black", size = 0.7) +
  scale_y_continuous(breaks=seq(0, 70000000, 10000000), labels =  c("0M",paste0(1:7,"0M")))
dev.off()

### useful
pdf("sample_useful.pdf",width = 4.5, height = 3)
ggplot(data.ratio, aes(x=group, y=useful)) + 
  geom_violin(aes(color=group,fill=group)) + 
  labs(x="", y="Useful reads percent") +
  scale_color_brewer(palette = "Reds", direction = -1) + scale_fill_brewer(palette = "Reds", direction = -1)+expand_limits(y=c(0,1))+
  stat_summary(fun.data =  data_summary, colour = "black", size = 0.7) +
  scale_y_continuous(labels = percent, breaks = seq(0,1,0.2))
dev.off()

mean_ci(data$useful)
# 0.8236428 0.807032 0.8402536

### mapped
pdf("sample_mapping.pdf",width = 4.5, height = 3)
ggplot(data, aes(x=group, y=map)) + 
  geom_violin(aes(color=group, fill=group)) +
  labs(x="", y="Mapped reads percent") + 
  scale_color_brewer(palette = "Reds", direction = -1) + scale_fill_brewer(palette = "Reds", direction = -1)+expand_limits(y=c(0,1))+
  stat_summary(fun.data =  data_summary, colour = "black", size = 0.7) +
  scale_y_continuous(labels = percent, breaks = seq(0,1,0.2))
dev.off()

mean_ci(data$map)
# 0.6524456 0.6315517 0.6733396

### miRNA kinds
miRNAkinds <- read.csv("miRNA_species_group.csv")
pdf("sample_miRNAspecies.pdf",width = 4.5, height = 3)
ggplot(miRNAkinds, aes(x=group, y=counts)) + 
  geom_violin(aes(color=group, fill=group)) +
  labs(x="", y="miRNA species") + 
  scale_color_brewer(palette = "Blues", direction = -1) + scale_fill_brewer(palette = "Blues", direction = -1)+
  stat_summary(fun.data =  data_summary, colour = "black", size = 0.7) +
  coord_cartesian(ylim = c(0, 600))
dev.off()


### Length distribution
length.dis <- read.csv("lenth.dis.csv")
length.dis <- reshape2::melt(length.dis, id.vars="X")
length.dis.sum <- sum(length.dis$value)
pdf("length_dis_percent.pdf", width = 5, height = 3)
ggplot(length.dis, aes(x=X, y=value/length.dis.sum, fill=variable, color=NULL)) + labs(x="", y="Percent")+
  geom_bar(position = position_stack(reverse = T), stat = "identity", width = 1)+
  scale_color_brewer(palette = "Set1") +scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent) + 
  coord_cartesian(x=c(17,50)) +
  cowplot::theme_cowplot()
dev.off()

### smallRNA ratio
smallRNAratio <- read.csv("smallRNA_ratio.csv")
smallRNAratio.mean <- apply(smallRNAratio[,2:9],2,mean)

df <- data.frame(value=as.numeric(smallRNAratio.mean), Group=names(smallRNAratio.mean)) %>%
  mutate(Group = factor(Group, levels = names(smallRNAratio.mean)),
         cumulative = cumsum(value),
         midpoint = cumulative - value / 2,
         label = paste0(Group, " ", round(value / sum(value) * 100, 1), "%"))
  
pdf("smallRNA_ratio_pie.pdf")
ggplot(df, aes(x = 1, weight = value, fill = Group)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") + theme_nothing() +
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "right")+
  geom_text(aes(x = 1.3, y = midpoint, label = label)) 
dev.off()


### Depth, miRNA kinds
data.depth.kinds <- read.csv("Depth_miRNA-kinds.csv")

pdf("depth_miRNAkinds.pdf", width = 5, height = 4)
ggplot(data.depth.kinds, aes(x=depth/1000000, y=kinds))+geom_point(fill ="#990033", color="#990033") + 
  geom_smooth(formula = y ~ log(x), se = T) + 
  coord_cartesian(ylim = c(0,500)) + theme_cowplot(font_size = 15, font_family = "", line_size = 1) +
  scale_x_continuous(breaks = seq(0,80,by=10)) + 
  labs(x="Sequencing depth(million reads)", y="Number of high condifence\ndeteced miRNA species") +
  geom_vline(aes(xintercept=10), linetype="dashed", size=1)
dev.off()
