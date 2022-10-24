
library(ggplot2)                    
library(cowplot)
library(pheatmap)

### Boxplot of ysRNA ratio
data <- <- read.csv("Cas9-sgRNA_treated_ratio.csv")
p1 <- ggplot(data, aes(x=group, y=yrna, fill=group, color=group)) + geom_boxplot()+
  coord_cartesian(ylim = c(0,40)) + labs(x="", y="ysRNA percent (%)") +
  scale_fill_brewer(palette =  "Set1", direction = 1)+ scale_color_brewer(palette =  "Set1")+
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  theme_cowplot() +theme(legend.position = "none")

# Average 27.97%, 11.86%
pdf("yrna_percent.pdf", width = 3, height = 4)
p1
dev.off()

### scatterplot plot
library(scales)
tmp <- miRNA[rowSums(miRNA[,1:2])>2,]
p2<- ggplot(data=tmp, aes(x=Ctrl1, y=Treat1)) +geom_point(color="#B2182B", size=3, alpha=0.8)+ #background_grid(minor='none')+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  labs(x="Ctrl (log10,RPM)", y="Treat (log10,RPM)")+  theme_half_open() + background_grid(size.major = 0.2)  +
  geom_smooth(method = MASS::rlm,se=F, color= "#B2182B")+expand_limits(x=c(0,100000),y=c(0,100000)) +
  annotate("text", x = 100, y = 100000, label =  "italic(R) ^ 2 == 0.9932", parse = TRUE, color="blue", size=6)

pdf("point_plot.pdf", width = 4.5, height = 4)
p2
dev.off()

# Adjusted R-squared:  0.9932
# linear regression
fit <- lm(Ctrl1~Treat1, miRNA)
summary(fit)


### heatmap correlation
library(pheatmap)
miRNA <- read.csv("Cas9-sgRNA_treated_miRNA_expression.txt", sep = "\t", row.names = 1)

miRNA.corr <- data.frame(cor(miRNA, method = "spearman"))
pheatmap(miRNA.corr, cluster_rows = F, cluster_cols = F, scale = "none", fontsize_row = 8,
         display_numbers=T, number_color = "black",  fontsize_number = 10,
         filename = "heatmap_cor_spearman.pdf", height = 4, width = 4.5, border_color = "white", show_rownames = T, 
         breaks=c(seq(0,1,length=101)),
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))

miRNA.corr2 <- data.frame(cor(miRNA, method = "pearson"))
pheatmap(miRNA.corr2, cluster_rows = F, cluster_cols = F, scale = "none", fontsize_row = 8,
         display_numbers=T, number_color = "black",  fontsize_number = 10,
         filename = "heatmap_cor_pearson.pdf", height = 4, width = 4.5, border_color = "white", show_rownames = T, 
         breaks=c(seq(0,1,length=101)),
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))

