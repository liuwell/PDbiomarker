library(ggplot2)
library(cowplot)
library(tidyverse)
library(stringr)
library(Mfuzz)

### Expression pattern
data <- read.csv("miRNA_nomalized_mean.txt", sep = "\t", header = T)
data <- data[data$kruskal.p<0.05, ]

df <- data
row.names(df)<- df[,1]
df2 <- df[, 2:4]
df3a<-as.matrix(df2)
df3Ex<- ExpressionSet(assayData = df3a)
df3F <- filter.NA(df3Ex,thres = 0.25)
df3F <- fill.NA(df3F,mode = 'mean')
df3F <- standardise(df3F)

set.seed(2022)
c <- 8
m <- mestimate(df3F) #
cl <- mfuzz(df3F, c = c, m = m) 

cl$size
cl$cluster[cl$cluster == 2]
cl$membership 

cl.thres <- acore(df3F, cl, min.acore=0) 
unlist(lapply(cl.thres, nrow))
lapply(cl.thres, head)

dir.create(path="mfuzz",recursive = TRUE)
for(i in 1:8){
  #potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl.thres[i],paste0("mfuzz","/mfuzz_",i,".csv"))
}

color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
pdf("mfuzz_200_c8.5.pdf", width = 10, height = 5)
mfuzz.plot2(df3F, cl=cl, mfrow=c(2,4),centre=TRUE, x11=F, centre.lwd=1, ylim.set = c(-1.5,1.5),
            time.labels = c("Healthy", "iRBD", "PD"), xlab = "Stage",
            colo = "fancy", min.mem = 0)
#X11(w=1.5,h=5);
#par(mar=c(1,1,1,5))
#mfuzzColorBar(col = "fancy", horizontal = F, main="Membership value")
dev.off()

### Healthy-iRBD-PD hierarchy
# increase
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

p.tab <- read.csv("miRNA_nomalized_mean.txt", sep = "\t", header = T)
p.tab.filter1 <- p.tab %>% as_tibble() %>% filter(kruskal.q<0.01) %>% filter(Healthy_mean<RBD_mean) %>%filter(RBD_mean<PD_mean)
miRNA.hierarchy.up <- c("let-7i-5p",  "miR-4433b-5p", "miR-766-3p", "miR-335-3p", "miR-130b-5p", "miR-744-5p",  "miR-1301-3p",
                        "miR-4433a-3p", "miR-889-3p", "miR-4433a-5p")
boxp <- limmaBatch[rownames(limmaBatch) %in% miRNA.hierarchy.up,1:169]
boxp$gene <- factor(rownames(boxp),levels = rownames(boxp))
boxp.melt <- melt(boxp, id.vars="gene")
boxp.melt$group <- factor(c(rep("Healthy", 10*60), rep("iRBD", 10*56), rep("PD",10*53)),levels = c("Healthy","iRBD", "PD"))
my_comparisons <- list( c("Healthy", "iRBD"), c("iRBD", "PD"), c("Healthy", "PD") )
pdf("boxplot_hierarchy_up2.pdf", width = 12,height = 5)
ggplot(boxp.melt, aes(x=group, y=value, color=group, fill=group)) + #geom_boxplot(outlier.shape = 18,outlier.size = 1)+ 
  geom_violin() +
  labs(x="", y="log2(normalized read count)") + #theme_cowplot(font_size = 18, font_family = "", line_size = 1) +
  scale_color_brewer(palette = "OrRd") + scale_fill_brewer(palette = "OrRd")+
  #geom_jitter(aes(color=group,fill=group), width = 0.25, shape=17)+
  facet_wrap(.~gene, scales = "free", ncol = 5) +theme(legend.position="none")+
  stat_compare_means(size=3) + stat_summary(fun.data =  data_summary, colour = "black", size = 0.7) 
  #stat_compare_means( comparisons = my_comparisons, size=3, method = "wilcox.test")
dev.off()

### decrease
p.tab.filter2 <- p.tab %>% as_tibble() %>% filter(kruskal.q<0.05) %>% filter(Healthy_mean>RBD_mean) %>%filter(RBD_mean>PD_mean)
miRNA.hierarchy.down <- c("miR-10b-5p",  "miR-150-5p",  "miR-342-3p",  "miR-192-5p",  "miR-361-3p",
                          "miR-155-5p",  "miR-3615",  "miR-150-3p",  "miR-500a-3p", "miR-186-5p")
boxp <- limmaBatch[rownames(limmaBatch) %in% miRNA.hierarchy.down,1:169]
boxp$gene <- factor(rownames(boxp),levels = rownames(boxp))
boxp.melt <- melt(boxp, id.vars="gene")
boxp.melt$group <- factor(c(rep("Healthy", 10*60), rep("iRBD", 10*56), rep("PD",10*53)),levels = c("Healthy","iRBD", "PD"))
my_comparisons <- list( c("Healthy", "iRBD"), c("iRBD", "PD"), c("Healthy", "PD") )

pdf("boxplot_hierarchy_down.pdf", width = 12,height = 5)
ggplot(boxp.melt, aes(x=group, y=value, color=group, fill=group)) + #geom_boxplot(outlier.shape = 18,outlier.size = 1)+ 
  geom_violin() +
  labs(x="", y="log2(normalized read count)") + #theme_cowplot(font_size = 18, font_family = "", line_size = 1) +
  scale_color_brewer(palette = "PuBu", direction = -1) + scale_fill_brewer(palette = "PuBu", direction = -1)+
  #geom_jitter(aes(color=group,fill=group), width = 0.25, shape=17)+
  facet_wrap(.~gene, scales = "free", ncol = 5) +theme(legend.position="none")+
  stat_compare_means( comparisons = my_comparisons, size=3, method = "wilcox.test")
  stat_compare_means(size=3) + stat_summary(fun.data =  data_summary, colour = "black", size = 0.7) 
  #stat_compare_means( comparisons = my_comparisons, size=3, method = "wilcox.test")
dev.off()

