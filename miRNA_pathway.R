library(ggplot)
library(tidyverse)

### miRPathDB
### the pathway of miRNA targets
path.miR.7 <- read.csv("miRPathDB_miR.7.5p.csv")
path.miR.7 <- path.miR.7[1:20, ]
path.miR.7$Pathway <- factor(path.miR.7$Pathway, levels = path.miR.7$Pathway)
pdf("Path_miR_7_5p.pdf", height = 9, width = 5)
ggplot(path.miR.7[1:20,], aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95)) 
dev.off()

path.miR.27b <- read.csv("miRPathDB_miR.27b.3p.csv")
path.miR.27b <- path.miR.27b[1:20, ]
path.miR.27b$Pathway <- factor(path.miR.27b$Pathway, levels = path.miR.27b$Pathway)
pdf("Path_miR_27b_3p.pdf", height = 9, width = 5)
ggplot(path.miR.27b, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95)) 
dev.off()

path.miR.182 <- read.csv("miRPathDB_miR.182.5p.csv")
path.miR.182 <- path.miR.182[-7, ]
path.miR.182 <- path.miR.182[1:20, ]
path.miR.182$Pathway <- factor(path.miR.182$Pathway, levels = path.miR.182$Pathway)
pdf("Path_miR_182_5p.pdf", height = 9, width = 5)
ggplot(path.miR.182, aes(x=Pathway, y=-log10(P.value))) + geom_col(color="#A50F15", fill="#A50F15", width=0.85) +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust = 0.95)) 
dev.off()


### miRNA expression in human tissue
### PD-specifc, iRBD-specific
biomarks3 <- c("miR-199a-5p", "miR-4433b-5p", "miR-27b-3p", "miR-151a-3p", "miR-584-5p",  "miR-130b-5p", "miR-197-3p", "miR-96-5p",
              "miR-150-3p", "miR-4433a-3p", "miR-619-5p", "miR-155-5p", "miR-150-5p", "miR-3615",  "miR-889-3p", "miR-182-5p", "miR-7-5p")

human.tissue.miRNA <- read.csv("miRNA_tissue_matrix_quantile.txt", sep = "\t")
rownames(human.tissue.miRNA) <- gsub("hsa-", "", rownames(human.tissue.miRNA))
human.tissue.miRNA.mark <- human.tissue.miRNA[rownames(human.tissue.miRNA)%in% biomarks3 , ] # lack miR-3615, miR-4433b-3p
human.tissue.miRNA.mark <- as.data.frame(t(human.tissue.miRNA.mark))
human.tissue.miRNA.mark$group <- rownames(human.tissue.miRNA.mark)

x <- strsplit(rownames(human.tissue.miRNA.mark), "[.]")
tissueName <- NULL
for (i in x) {
 tissueName <- c(tissueName, i[1])  
}

human.tissue.miRNA.mark$group <- tissueName
write.table(human.tissue.miRNA.mark, file = "human.tissue.miRNA.markers2.txt", sep="\t", quote = F)

#human.tissue.miRNA.marker <- read.csv("human.tissue.miRNA.markers.txt", sep = "\t", header = T)
human.tissue.miRNA.marker.melt <- reshape2::melt(human.tissue.miRNA.mark, ids.var="group")
human.tissue.miRNA.marker2 <-  human.tissue.miRNA.marker.melt %>% as_tibble()  %>% 
  group_by(group, variable) %>% dplyr::summarise(value=mean(value)) 

pdf("human.tissue.miRNA.biomarkers5.pdf", width = 12, height = 6)
ggplot(human.tissue.miRNA.marker2, aes(x=factor(group), y=variable, size=value)) +
  #geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", )
  geom_point(color="#A50F15") + #scale_color_gradient(low = "#FEE5D9", high = "#A50F15") + scale_fill_gradient(low = "#FEE5D9", high = "#A50F15") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95)) +
  labs(y=NULL, x="Human tissue")
dev.off()


