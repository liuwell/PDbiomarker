library(ggplot2)
library(cowplot)
library(tidyverse)
library(RColorBrewer)
library(limma)
library(Boruta)
library(e1071)
library(ggrepel)

### filter low frequency miRNAs
sample.freq <- function(data){ 
  a <- NULL
  for(i in 1:nrow(data)){
    x <- 0
    for(j in 1:ncol(data)){
      if(data[i,j]>0){ 
        x <- x+1
      }
    }
    a <- append(a,x)
  }
  return(a)
}

### normalzied by size factor
calc_sf <- function (expr_mat){
  geomeans <- exp(rowMeans(log(expr_mat)))
  SF <- function(cnts){
    median((cnts/geomeans)[(is.finite(geomeans) & geomeans >0)])
  }
  norm_factor <- apply(expr_mat,2,SF)
  return(t(t(expr_mat)/norm_factor))
}

### PCA
PCA.point <- function(data, g1,g2,g3,g4, angle=45, output){
  hd<-t(data)
  pca<-prcomp(hd, center = T, scale. =T)
  
  pca2 <- as.data.frame(pca$x)
  pca2$group <- factor(rep(c(g1, g2),c(g3, g4)), levels = c(g1, g2))
  
  labx <- paste("PC1", sprintf('(%0.1f%%)', 100 * pca$sdev[1]^2/sum(pca$sdev^2)))
  laby <- paste("PC2", sprintf('(%0.1f%%)', 100 * pca$sdev[2]^2/sum(pca$sdev^2)))
  
  #output.2d <- paste(substitute(data), "PCA.2d.pdf", sep = ".")
  output.2d <- paste(output, "PCA.2d.pdf", sep = ".")
  p <- ggplot(pca2, aes(x=PC1, y=PC2, color=group, shape=group)) +geom_point(size=3, alpha=0.8) +labs(x=labx, y=laby) + 
    scale_color_brewer(palette ="Set1") + theme_cowplot(font_size = 23, font_family = "", line_size = 1)
  
  pdf(output.2d, width = 6.5, height = 4)
  print(p)
  dev.off()
  
  # PCA scatterplot 3D
  labx <- paste("PC1", sprintf('(%0.1f%%)', 100 * pca$sdev[1]^2/sum(pca$sdev^2)))
  laby <- paste("PC2", sprintf('(%0.1f%%)', 100 * pca$sdev[2]^2/sum(pca$sdev^2)))
  labz <- paste("PC3", sprintf('(%0.1f%%)', 100 * pca$sdev[3]^2/sum(pca$sdev^2)))
  
  shapes = c(16, 17, 18) 
  shapes <- shapes[as.numeric(pca2$group)]
  
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A","#984EA3")
  colors <- colors[as.numeric(pca2$group)]
  output.3d <- paste(output, "PCA.3d.pdf", sep = ".")
  pdf(output.3d, width = 4, height = 4)
  scatterplot3d(pca2[,c(1,2,3)], pch = shapes, color = colors, angle = angle, xlab = labx, ylab = laby, zlab = labz)
  dev.off()
}

### feature selection by Boruta
feature.select <- function(data,a,b){
  data.t <- as.data.frame(t(data))
  data.t$class <- factor(c(rep(1, a), rep(0, b)), levels = c(1, 0))
  
  set.seed(2022)
  Boruta.ivf.extended <- Boruta(class~.,data=data.t, pValue=0.01, maxRuns=200)
  print(Boruta.ivf.extended)
  
  output <- paste(substitute(data), "boruta.pdf", sep = ".")
  pdf(output)
  plot(Boruta.ivf.extended)
  dev.off()
  
  Boruta.stats <- attStats(Boruta.ivf.extended)
  Boruta.confirmed <- Boruta.stats[which(Boruta.stats$decision=="Confirmed"),]
  Boruta.confirmed <- Boruta.confirmed[order(Boruta.confirmed$meanImp, decreasing = T),]
  Boruta.confirmed$gene <- gsub("`", "", rownames(Boruta.confirmed))
  return(Boruta.confirmed)
}

### svm model
acc.dis.singleModle <- function(data, g1, g2, seed=100, info=FALSE){ 
  
  tmp <- as.data.frame(t(data))
  tmp$class <- factor(c(rep("ctrl", g1), rep("treat", g2)), levels = c("ctrl", "treat"))
  
  train.ac <- NULL; train.sp <- NULL; train.se <- NULL
  validate.ac <- NULL; validate.sp <- NULL; validate.se <- NULL
  
  n1 <- nrow(tmp) # total samples
  n2 <- n1*0.6    # 60% samples
  wts <- 100 / table(tmp$class)
  
  set.seed(seed)
  train <- sample(n1, n2)
  
  #train_svm <- train ### preserve training dataset
  tmp.train <- tmp[train,]
  tmp.validate <- tmp[-train,]
  
  # set best parameters of gamma and cost
  set.seed(seed)
  svm.tuned <- tune.svm(class~., data = tmp.train, gamma = 10^(-3:3), cost = 10^(-3:3), class.weights = wts) 
  set.seed(seed)
  fit.svm <- svm(class~., data = tmp.train, gamma=svm.tuned$best.parameters[1], cost=svm.tuned$best.parameters[2], cross=5,
                 class.weights = wts)
  
  # for train
  pred_data <- tmp.train
  svm.pred <- predict(fit.svm, na.omit(pred_data))
  svm.perf <- table(na.omit(pred_data)$class, svm.pred, dnn=c("Actual", "Predict"))
  #print(svm.perf)
  
  tp = svm.perf[1,1];fn = svm.perf[1,2];fp = svm.perf[2,1];tn = svm.perf[2,2]
  se = round(tp/(tp+fn), 4)
  sp = round(tn/(tn+fp), 4)
  ac = round((tp+tn)/(tp+tn+fp+fn), 4)
  train_score <- c(se,sp,ac)
  # for validate 
  pred_data <- tmp.validate
  svm.pred2 <- predict(fit.svm, na.omit(pred_data))
  svm.perf2 <- table(na.omit(pred_data)$class, svm.pred2, dnn=c("Actual", "Predict"))
  #print(svm.perf2)
  
  tp = svm.perf2[1,1];fn = svm.perf2[1,2];fp = svm.perf2[2,1];tn = svm.perf2[2,2]
  se = round(tp/(tp+fn), 4)
  sp = round(tn/(tn+fp), 4)
  ac = round((tp+tn)/(tp+tn+fp+fn), 4)
  test_score <- c(se,sp,ac)
  
  result <- list(train_table=svm.perf, test_table=svm.perf2, train_score=train_score, 
                    test_score=test_score, train_pred=svm.pred, test_pred=svm.pred2
                 )
  #return(fit.svm)
  if(info){ 
  print(svm.perf)
  print(svm.perf2)
  print(result$train_score)
  print(result$test_score)
  print(paste0(seed, "-----------------------"))
  return(result)
  }else{
    return(fit.svm)}
}

acc.dis.singleModle.2 <- function(data, g1, g2, seed=100){ 
  
  tmp <- as.data.frame(t(data))
  tmp$class <- c(rep(1, g1), rep(0, g2))
  
  train.ac <- NULL; train.sp <- NULL; train.se <- NULL
  validate.ac <- NULL; validate.sp <- NULL; validate.se <- NULL
  
  n1 <- nrow(tmp) # total samples
  n2 <- n1*0.6    # 60% samples
  wts <- 100 / table(tmp$class)
  
  set.seed(seed)
  train <- sample(n1, n2)
  
  #train_svm <- train ### preserve training dataset
  tmp.train <- tmp[train,]
  tmp.validate <- tmp[-train,]
  
  # set best parameters of gamma and cost
  set.seed(seed)
  svm.tuned <- tune.svm(class~., data = tmp.train, gamma = 10^(-3:3), cost = 10^(-3:3))
  set.seed(seed)
  fit.svm <- svm(class~., data = tmp.train, gamma=svm.tuned$best.parameters[1], cost=svm.tuned$best.parameters[2], cross=5,
                 class.weights = wts)
  
  # for train
  pred_data <- tmp.train
  svm.pred <- predict(fit.svm, na.omit(pred_data))
  svm.perf <- table(na.omit(pred_data)$class, svm.pred, dnn=c("Actual", "Predict"))
  #print(svm.perf)
  
  tp = svm.perf[1,1];fn = svm.perf[1,2];fp = svm.perf[2,1];tn = svm.perf[2,2]
  se = round(tp/(tp+fn), 4)
  sp = round(tn/(tn+fp), 4)
  ac = round((tp+tn)/(tp+tn+fp+fn), 4)
  train_score <- c(se,sp,ac)
  # for validate 
  pred_data <- tmp.validate
  svm.pred <- predict(fit.svm, na.omit(pred_data))
  svm.perf2 <- table(na.omit(pred_data)$class, svm.pred, dnn=c("Actual", "Predict"))
  #print(svm.perf2)
  
  tp = svm.perf2[1,1];fn = svm.perf2[1,2];fp = svm.perf2[2,1];tn = svm.perf2[2,2]
  se = round(tp/(tp+fn), 4)
  sp = round(tn/(tn+fp), 4)
  ac = round((tp+tn)/(tp+tn+fp+fn), 4)
  test_score <- c(se,sp,ac)
  
  result <- list(train_table=svm.perf, test_table=svm.perf2, train_score=train_score, test_score=test_score)
  return(fit.svm)
}

### ROC
roc.plot3 <- function(data, g1, g2, output, r=100, svm.model){
  tmp <- as.data.frame(t(data))
  tmp$class <- c(rep(1, g1), rep(0, g2))
  
  auc.train.set <- NULL; auc.test.set <- NULL
  train.tpr.set <- NULL; train.fpr.set <- NULL;test.tpr.set <- NULL; test.fpr.set <- NULL
  
  n1 <- nrow(tmp) # total samples
  n2 <- n1*0.6    # 60% samples
  
  for(i in 1:r){
    train <- sample(n1, n2)
    
    tmp.train <- tmp[train,]
    tmp.validate <- tmp[-train,]
    
    # train dataset
    svm.prob <- predict(svm.model, tmp.train, type="response")
    svm.pred <- prediction(svm.prob, tmp.train$class)
    svm.perf <- performance(svm.pred, measure = "tpr", x.measure = "fpr")
    auc <- performance(svm.pred, measure = "auc")
    auc.train <- round(auc@y.values[[1]], 3)
    auc.train.set <- append(auc.train.set, auc.train)
    
    train.fpr=unlist(svm.perf@x.values)
    train.tpr=unlist(svm.perf@y.values)
    #roc.train <- data.frame(fpr=fpr, tpr=tpr, model="train")
    train.fpr.set <- append(train.fpr.set, train.fpr)
    train.tpr.set <- append(train.tpr.set, train.tpr)
    
    # test dataset
    svm.prob <- predict(svm.model, tmp.validate, type="response")
    svm.pred <- prediction(svm.prob, tmp.validate$class)
    svm.perf <- performance(svm.pred, measure = "tpr", x.measure = "fpr")
    auc <- performance(svm.pred, measure = "auc")
    auc.test <- round(auc@y.values[[1]], 3)
    auc.test.set <- append(auc.test.set, auc.test)
    
    test.fpr=unlist(svm.perf@x.values)
    test.tpr=unlist(svm.perf@y.values)
    #roc.test <- data.frame(fpr=fpr,tpr=tpr,model="test")
    test.fpr.set <- append(test.fpr.set, test.fpr)
    test.tpr.set <- append(test.tpr.set, test.tpr)
    
    #print(i)
  }
  
  ### plot train
  train.fpr.matrix <- matrix(train.fpr.set, ncol = r)
  train.tpr.matrix <- matrix(train.tpr.set, ncol = r)
  
  roc.train <- data.frame(fpr=train.fpr, tpr=train.tpr, model="train")
  m.ci <- round(mean_ci(auc.train.set),3)
  m.ci.label <- paste0("AUC=", m.ci[[1]], ",","CI=[", m.ci[[2]],",",m.ci[[3]],"]")
  
  ### plot test
  test.fpr.matrix <- matrix(test.fpr.set, ncol = r)
  test.tpr.matrix <- matrix(test.tpr.set, ncol = r)
  
  roc.test <- data.frame(fpr=test.fpr, tpr=test.tpr, model="validation")
  m.ci2 <- round(mean_ci(auc.test.set),3)
  m.ci2.label <- paste0("AUC=", m.ci2[[1]], ",","CI=[", m.ci2[[2]],",",m.ci2[[3]],"]")
  
  for(i in 1:length(auc.test.set)){
    x1 <- round(auc.test.set[i],3)
    x2 <- round(m.ci2[[1]],3)
    if(x1==x2){
      roc.train <- data.frame(fpr=train.fpr.matrix[,i], tpr=train.tpr.matrix[,i], model="train")
      roc.test <- data.frame(fpr=test.fpr.matrix[,i], tpr=test.tpr.matrix[,i], model="validation")
      roc.combine <- rbind(roc.train, roc.test)
      p1 <- ggplot(roc.combine, aes(x=fpr, ymin=0, ymax=tpr, group=model, color=model)) +
        labs(x="False positive rate",y="True posotive rate", title="")+
        theme_cowplot(font_size = 16, font_family = "", line_size = 1)+
        geom_abline(intercept = 0, slope = 1,size=0.5)+
        geom_line(aes(y=tpr), size=1) + scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
        geom_label(aes(x=0.7,y=0.35,label= m.ci.label),size=3, family="serif", fontface="plain", color="#E41A1C") + 
        geom_label(aes(x=0.7,y=0.25,label= m.ci2.label),size=3, family="serif", fontface="plain", color="#377EB8") + 
        coord_cartesian(xlim = c(0.03,1.0),ylim = c(0.03,1.0)) +
        theme(legend.position = "none")
      
      pdf(paste0(output, ".",i,".roc.pdf"), width = 4,height = 3.5)
      print(p1)
      dev.off()
    }
  }
}

### data input
data.PD <- read.csv("PD_miRNA_counts.csv", sep = "\t", header = T, row.names = 1)
data.PD <- data.PD[, sort(names(data.PD))]

a <- sample.freq(data.PD)
data.PD.freq <- data.PD[a>42,]
data.PD.sf <- as.data.frame(calc_sf(data.PD.freq))

### for batch effect and gender bias
PD.sex.age <- read.csv("Samples_Sex_Age.csv")

batch <- rep(c(1,2,1,2,1,2), c(33,27, 25,28, 32,24 ))
sex <- PD.sex.age$Sex

data.PD.sf <- data.PD.sf[order(rowSums(data.PD.sf), decreasing = T),]
data.PD.sf2 <- log2(as.matrix(data.PD.sf+1))
limmaBatch <- as.data.frame(removeBatchEffect(data.PD.sf2, batch, sex))

limmaBatch <- limmaBatch[, c(1:60, 114:169, 61:113)]
write.csv(limmaBatch, file = "PD_miRNA_normalzied.csv")

###################################################################################
# Healthy vs PD
condition <-factor(c(rep("C",60),rep("PD",53)))
design <- model.matrix(~condition)
fit <- lmFit(limmaBatch[,c(1:60, 117:169)], design)
fit <- eBayes(fit)
topTable(fit, coef = 2)

DE.miRNA.CP <- topTable(fit, coef = 2,sort.by = "P", number = nrow(limmaBatch))
DE.miRNA.CP$miRNA <- rownames(DE.miRNA.CP)
DE.miRNA.CP.sig <- DE.miRNA.CP　%>%　as_tibble() %>% filter(P.Value<0.05) %>% filter(abs(logFC)>0.5) ### nrow=127, up=89, down=38
DE.miRNA.CP.expr <-  limmaBatch[rownames(limmaBatch)%in% DE.miRNA.CP.sig$miRNA, c(1:60, 117:169)]

# heatmap
df <- data.frame(Group=rep(c("Healthy", "PD"),c(60, 53)), Gender=PD.sex.age[c(1:60, 117:169),2], Age=PD.sex.age[c(1:60, 117:169),3])
rownames(df) <- colnames(DE.miRNA.CP.expr)

nc<-brewer.pal(9,"Set1")
nc2 <-  brewer.pal(9,"Greys")
anno_color = list(Group = c("Healthy" = nc[1], "PD" = nc[2]), Gender=c("Female"=nc[3], "Male"=nc[4]),
                  Age = c("35-50"=nc2[3], "51-60"=nc2[5], "61-70"=nc2[7], "71-80"=nc2[9]))

pheatmap(DE.miRNA.CP.expr, cluster_rows = T, cluster_cols = T, scale = "row", fontsize_row = 5, 
         display_numbers=F, number_color = "black",  fontsize_col = 4, 
         annotation_col = df, annotation_colors = anno_color, #clustering_method = "mcquitty",
         filename = "DE.miRNA.CP.expr.pdf", height = 6, width = 5, show_rownames = F, show_colnames = F,border_color = NA,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(21))

# PCA
PCA.point(DE.miRNA.CP.expr, "Healthy", "PD", 60, 53, output = "DE.miRNA.CP", angle = 45)

# MA-plot
dMA <- DE.miRNA.CP
dMA$change <- as.factor(ifelse(dMA$P.Value < 0.05 & abs(dMA$logFC) > 0.5,ifelse(dMA$logFC > 0.5,'Up','Down'),'Stable'))
table(dMA$change) # down 39, up 87,
dMA$gene <- rownames(dMA)

pdf("DE.miRNA.CP.MA.pdf", width = 6,height = 4)
ggplot(dMA, aes(x=AveExpr, y=logFC, color=change)) + geom_point(alpha=0.8, size = 3)+
  #scale_color_manual(values =c("#E41A1C","#377EB8", "#984EA3"))+
  scale_color_manual(values =c("#377EB8","#969696", "#E41A1C"))+
  geom_hline(yintercept = 0.5,lty=4,lwd=0.6,alpha=0.8) + geom_hline(yintercept = -0.5,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+theme_bw(base_size = 18)+ylim(-3,3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  labs(x="Log2 average normalized read counts", y="Log2 fold change") 
dev.off()

# feature select
DE.miRNA.CP.fs <- feature.select(DE.miRNA.CP.expr, 60, 53)
DE.miRNA.CP.fs.expr2 <- limmaBatch[rownames(limmaBatch)%in%DE.miRNA.C.RP.fs$miRNA,]

# Confusion table 
svmModle5.2 <- acc.dis.singleModle.2(DE.miRNA.CP.fs.expr2, 60, 53, seed = 4)
roc.plot3(DE.miRNA.CP.fs.expr2, 60, 53, "C.P.miRNA", r=200, svmModle5.2)

acc_results.CP <- acc.dis.singleModle(DE.miRNA.CP.fs.expr2, 60, 53, seed = 4, info = TRUE)

# train
tmp <- acc_results.CP$train_table
tmp <- matrix(tmp, ncol = 2, dimnames = list(c("Healthy", "PD"), c("Healthy", "PD")))
pheatmap(tmp, cluster_rows = F, cluster_cols = F, display_numbers = T, angle_col=0,
         number_format = "%.0f", number_color = "black", fontsize_number = 20, fontsize = 15,
         border_color = NA, filename = "C.P.train_matrix.pdf", height = 4, width = 4.8,
         color = colorRampPalette(brewer.pal(n =6, name ="Reds"))(11))
# test
tmp <- acc_results.CP$test_table
tmp <- matrix(tmp, ncol = 2, dimnames = list(c("Healthy", "PD"), c("Healthy", "PD")))
pheatmap(tmp, cluster_rows = F, cluster_cols = F, display_numbers = T, angle_col=0,
         number_format = "%.0f", number_color = "black", fontsize_number = 20, fontsize = 15,
         border_color = NA, filename = "C.P.test_matrix.pdf", height = 4, width = 4.8,
         color = colorRampPalette(brewer.pal(n = 6, name ="Blues"))(11))


###################################################################################
# Healthy vs iRBD
condition <-factor(c(rep("C", 60),rep("iRBD",56)))
design <- model.matrix(~condition)
fit <- lmFit(limmaBatch[,c(1:60,61:116)], design)
fit <- eBayes(fit)
topTable(fit, coef = 2)

DE.miRNA.CR <- topTable(fit, coef = 2,sort.by = "P", number = nrow(limmaBatch))
DE.miRNA.CR$miRNA <- rownames(DE.miRNA.CR)
DE.miRNA.CR.sig <- DE.miRNA.CR　%>%　as_tibble() %>% filter(P.Value<0.05) %>% filter(abs(logFC)>0.6)
DE.miRNA.CR.expr <-  limmaBatch[rownames(limmaBatch)%in% DE.miRNA.CR.sig$miRNA, c(1:60, 61:116)]

# heatmap
df <- data.frame(Group=rep(c("Healthy", "iRBD"),c(60, 56)), Gender=PD.sex.age[c(1:60, 61:116), 2],  Age=PD.sex.age[c(1:60, 61:116), 3])
rownames(df) <- colnames(DE.miRNA.CR.expr)

nc<-brewer.pal(9,"Set1")
nc2 <-  brewer.pal(9,"Greys")
anno_color = list(Group = c("Healthy" = nc[1], "iRBD" = nc[2]), Gender=c("Female"=nc[3], "Male"=nc[4]),
  Age = c("35-50"=nc2[3], "51-60"=nc2[5], "61-70"=nc2[7], "71-80"=nc2[9]))
pheatmap(DE.miRNA.CR.expr, cluster_rows = T, cluster_cols = T, scale = "row", fontsize_row = 5, 
         display_numbers=F, number_color = "black",  fontsize_col = 4, 
         annotation_col = df, annotation_colors = anno_color,
         filename = "DE.miRNA.CR.expr.pdf", height = 6, width = 5, show_rownames = F, show_colnames = F,border_color = NA,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(21))

# PCA
PCA.point(DE.miRNA.CR.expr, "Healthy", "iRBD", 60, 56, output = "DE.miRNA.CR", angle = 60)

# MA-plot
dMA <- DE.miRNA.CR
dMA$change <- as.factor(ifelse(dMA$adj.P.Val < 0.05 & abs(dMA$logFC) > 0.6,ifelse(dMA$logFC > 0.6,'Up','Down'),'Stable'))
dMA$gene <- rownames(dMA)
pdf("DE.miRNA.CR.MA.pdf", width = 6,height = 4)
ggplot(dMA, aes(x=AveExpr, y=logFC, color=change)) + geom_point(alpha=0.8, size = 3)+
  scale_color_manual(values =c("#377EB8","#969696", "#E41A1C"))+
  geom_hline(yintercept = 0.6,lty=4,lwd=0.6,alpha=0.8) + geom_hline(yintercept = -0.6,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+theme_bw(base_size = 18)+ylim(-3,3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  labs(x="Log2 average normalized read counts", y="Log2 fold change") 
dev.off()


# feature selection 
DE.miRNA.CR.fs <- feature.select(DE.miRNA.CR.expr, 56, 53)
DE.miRNA.C.R.fs.expr4 <- limmaBatch[rownames(limmaBatch)%in%DE.miRNA.CR.fs$miRNA,]

svmModle4.2 <- acc.dis.singleModle.2(DE.miRNA.C.R.fs.expr4, 60, 56, seed = 13)
roc.plot3(DE.miRNA.C.R.fs.expr4, 60, 56, "C.R.miRNA.3", r=100, svmModle4.2)

acc_results.CR <- acc.dis.singleModle(DE.miRNA.C.R.fs.expr4, 60, 56, seed = 13, info = TRUE)
# train
tmp <- acc_results$train_table
tmp <- matrix(tmp, ncol = 2, dimnames = list(c("Healthy", "iRBD"), c("Healthy", "iRBD")))
pheatmap(tmp, cluster_rows = F, cluster_cols = F, display_numbers = T, angle_col=0,
         number_format = "%.0f", number_color = "black", fontsize_number = 20, fontsize = 15,
         border_color = NA, filename = "C.R.train_matrix.pdf", height = 4, width = 4.8,
         color = colorRampPalette(brewer.pal(n =5, name ="Reds"))(11))
# test
tmp <- acc_results$test_table
tmp <- matrix(tmp, ncol = 2, dimnames = list(c("Healthy", "iRBD"), c("Healthy", "iRBD")))
pheatmap(tmp, cluster_rows = F, cluster_cols = F, display_numbers = T, angle_col=0,
         number_format = "%.0f", number_color = "black", fontsize_number = 20, fontsize = 15,
         border_color = NA, filename = "C.R.test_matrix.pdf", height = 4, width = 4.8,
         color = colorRampPalette(brewer.pal(n =5, name ="Blues"))(11))


###################################################################################
# iRBD vs PD
condition <-factor(c(rep("iRBD", 56),rep("PD", 53)))
design <- model.matrix(~condition)
fit <- lmFit(limmaBatch[,c(61:169)], design)
fit <- eBayes(fit)
topTable(fit, coef = 2)

DE.miRNA.RP <- topTable(fit, coef = 2,sort.by = "P", number = nrow(limmaBatch))
DE.miRNA.RP$miRNA <- rownames(DE.miRNA.RP)
DE.miRNA.RP.sig <- DE.miRNA.RP　%>%　as_tibble() %>% filter(P.Value <0.05) %>% filter(abs(logFC)>0.5)
DE.miRNA.RP.expr <-  limmaBatch[rownames(limmaBatch)%in% DE.miRNA.RP.sig$miRNA, c(61:169)]

# heatmap
df <- data.frame(Group=rep(c("iRBD", "PD"),c(56, 53)), Gender=PD.sex.age[c(61:169),2] , Age=PD.sex.age[c(61:169),3])
rownames(df) <- colnames(DE.miRNA.RP.expr)

nc<-brewer.pal(9,"Set1")
nc2 <-  brewer.pal(9,"Greys")
anno_color = list(Group = c("iRBD" = nc[1], "PD" = nc[2]), Gender=c("Female"=nc[3], "Male"=nc[4]),
                  Age = c("35-50"=nc2[3], "51-60"=nc2[5], "61-70"=nc2[7], "71-80"=nc2[9]))

pheatmap(DE.miRNA.RP.expr, cluster_rows = T, cluster_cols = T, scale = "row", fontsize_row = 5, 
         display_numbers=F, number_color = "black",  fontsize_col = 4, 
         annotation_col = df, annotation_colors = anno_color, clustering_method = "complete",
         filename = "DE.miRNA.RP.expr.pdf", height = 6, width = 5, show_rownames = F, show_colnames = F,border_color = NA,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(21))

# PAC
PCA.point(DE.miRNA.RP.expr, "iRBD", "PD", 56, 53, output = "DE.miRNA.RP", angle = 60)

# MA-plot
dMA <- DE.miRNA.RP
dMA$change <- as.factor(ifelse(dMA$P.Value < 0.05 & abs(dMA$logFC) > 0.5,ifelse(dMA$logFC > 0.5,'Up','Down'),'Stable'))# 42
dMA$gene <- rownames(dMA)

pdf("DE.miRNA.RP.MA.pdf", width = 6,height = 4)
ggplot(dMA, aes(x=AveExpr, y=logFC, color=change)) + geom_point(alpha=0.8, size = 3)+
  scale_color_manual(values =c("#377EB8","#969696", "#E41A1C"))+
  geom_hline(yintercept = 0.5,lty=4,lwd=0.6,alpha=0.8) + geom_hline(yintercept = -0.5,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+theme_bw(base_size = 18)+ylim(-3,3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  labs(x="Log2 average normalized read counts", y="Log2 fold change") 
dev.off()

# feature select
DE.miRNA.RP.fs <- feature.select(DE.miRNA.RP.expr, 56, 53)
DE.miRNA.RP.fs.expr2 <- limmaBatch[rownames(limmaBatch)%in%DE.miRNA.RP.fs$miRNA,]

svmModle6.2 <- acc.dis.singleModle.2(DE.miRNA.RP.fs.expr2, 56, 53, seed = 2)
roc.plot3(DE.miRNA.RP.fs.expr2, 56, 53, "R.P.miRNA", r=200, svmModle6.2)

acc_results.RP <- acc.dis.singleModle(DE.miRNA.RP.fs.expr2, 56, 53, seed = 2, info = TRUE)

# train
tmp <- acc_results.RP$train_table
tmp <- matrix(tmp, ncol = 2, dimnames = list(c("iRBD", "PD"), c("iRBD", "PD")))
pheatmap(tmp, cluster_rows = F, cluster_cols = F, display_numbers = T, angle_col=0,
         number_format = "%.0f", number_color = "black", fontsize_number = 20, fontsize = 15,
         border_color = NA, filename = "R.P.train_matrix.pdf", height = 4, width = 4.8,
         color = colorRampPalette(brewer.pal(n =6, name ="Reds"))(11))
# test
tmp <- acc_results.RP$test_table
tmp <- matrix(tmp, ncol = 2, dimnames = list(c("iRBD", "PD"), c("iRBD", "PD")))
pheatmap(tmp, cluster_rows = F, cluster_cols = F, display_numbers = T, angle_col=0,
         number_format = "%.0f", number_color = "black", fontsize_number = 20, fontsize = 15,
         border_color = NA, filename = "R.P.test_matrix.pdf", height = 4, width = 4.8,
         color = colorRampPalette(brewer.pal(n = 6, name ="Blues"))(11))


### iRBD, PD-iRBD
PD.iRBD.id <- read.csv("PD_with_iRBD.csv", header = T)

PD <- limmaBatch[, 117:169]
iRBD <- limmaBatch[, 61:116]
PD.iRBD <- PD[, colnames(PD) %in% PD.iRBD.id$V1]
PD.no <- PD[, !colnames(PD) %in% PD.iRBD.id$V1]

d1 <- cbind(iRBD, PD.iRBD)
condition <-factor(c(rep("C", 56),rep("PD", 18)))
design <- model.matrix(~condition)
fit <- lmFit(d1, design)
fit <- eBayes(fit)
topTable(fit, coef = 2)

DE.miRNA.d1 <- topTable(fit, coef = 2,sort.by = "P", number = nrow(limmaBatch))
DE.miRNA.d1$miRNA <- rownames(DE.miRNA.d1)
DE.miRNA.d1.sig <- DE.miRNA.d1　%>%　as_tibble() %>% filter(P.Value<0.05) %>% 
  filter(abs(logFC)>0.5) 
DE.miRNA.d1.expr <-  d1[rownames(d1)%in% DE.miRNA.d1.sig$miRNA, ]

write.table(DE.miRNA.d1.sig, file = "DE.miRNA.iRBD.PD-with-iRBD.sig.txt", sep = "\t", quote = F)

# MA-plot
dMA <- DE.miRNA.d1
dMA$change <- as.factor(ifelse(dMA$P.Value < 0.05 & abs(dMA$logFC) > 0.5,ifelse(dMA$logFC > 0.5,'Up','Down'),'Stable'))

dMA$gene <- rownames(dMA)

pdf("DE.miRNA.d1.MA.pdf", width = 6,height = 4)
ggplot(dMA, aes(x=AveExpr, y=logFC, color=change)) + geom_point(alpha=0.8, size = 3)+
  scale_color_manual(values =c("#377EB8","#969696", "#E41A1C"))+
  geom_hline(yintercept = 0.6,lty=4,lwd=0.6,alpha=0.8) + geom_hline(yintercept = -0.6,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+theme_bw(base_size = 18)+ylim(-3,3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  labs(x="Log2 average normalized read counts", y="Log2 fold change") 
dev.off()

# heatmap
df <- data.frame(group=rep(c("iRBD", "PD-iRBD"),c(56, 18)))
rownames(df) <- colnames(DE.miRNA.d1.expr)

nc<-brewer.pal(9,"Set1")
anno_color = list(group = c("iRBD" = nc[1], "PD-iRBD" = nc[2]))
pheatmap(DE.miRNA.d1.expr, cluster_rows = T, cluster_cols = F, scale = "row", fontsize_row = 5, 
         display_numbers=F, number_color = "black",  fontsize_col = 4, 
         annotation_col = df, annotation_colors = anno_color,
         filename = "DE.miRNA.d1.expr.pdf", height = 3, width = 5, 
         show_rownames = F, show_colnames = F,border_color = NA,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(21))


### PD-iRBD, PD-no-iRBD
d3 <- cbind(PD.iRBD, PD.no)
condition <-factor(c(rep("C", 18),rep("PD", 35)))
design <- model.matrix(~condition)
fit <- lmFit(d3, design)
fit <- eBayes(fit)
topTable(fit, coef = 2)

DE.miRNA.d3 <- topTable(fit, coef = 2,sort.by = "P", number = nrow(limmaBatch))
DE.miRNA.d3$miRNA <- rownames(DE.miRNA.d3)
DE.miRNA.d3.sig <- DE.miRNA.d3　%>%　as_tibble() %>% filter(P.Value<0.05) %>% filter(abs(logFC)>0.5) 
DE.miRNA.d3.expr <-  d3[rownames(d3)%in% DE.miRNA.d3.sig$miRNA, ]

DE.miRNA.d3.fs <- feature.select(DE.miRNA.d3.expr, 18, 35)
DE.miRNA.d3.fs.expr <- d3[rownames(d3)%in%DE.miRNA.d3.fs$miRNA,]
DE.miRNA.d3.fs$miRNA

write.table(DE.miRNA.d3.sig, file = "DE.miRNA.PD-iRBD.PD-no-iRBD.sig.txt", sep = "\t", quote = F)
write.table(DE.miRNA.d3.fs, file = "DE.miRNA.PD-iRBD.PD-no-iRBD.biomarkers.txt", sep = "\t", quote = F)

PCA.point(DE.miRNA.d3.expr, "PD-iRBD", "PD-no-iRBD", 18, 35, output = "DE.miRNA.d3", angle = 45)

# MA-plot
dMA <- DE.miRNA.d3
dMA$change <- as.factor(ifelse(dMA$P.Value < 0.05 & abs(dMA$logFC) > 0.5,ifelse(dMA$logFC > 0.5,'Up','Down'),'Stable'))
table(dMA$change)

dMA$gene <- rownames(dMA)
brewer.pal(9,"Set1")
pdf("DE.miRNA.d3.MA.pdf", width = 6,height = 4)
ggplot(dMA, aes(x=AveExpr, y=logFC, color=change)) + geom_point(alpha=0.8, size = 3)+
  scale_color_manual(values =c("#377EB8","#969696", "#E41A1C"))+
  geom_hline(yintercept = 0.6,lty=4,lwd=0.6,alpha=0.8) + geom_hline(yintercept = -0.6,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+theme_bw(base_size = 18)+ylim(-3,3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  labs(x="Log2 average normalized read counts", y="Log2 fold change") 
dev.off()

# heatmap
df <- data.frame(group=rep(c("PD-iRBD", "PD-no-iRBD"),c(18, 35)))
rownames(df) <- colnames(DE.miRNA.d3.expr)

nc<-brewer.pal(9,"Set1")
anno_color = list(group = c("PD-iRBD" = nc[1], "PD-no-iRBD" = nc[2]))
pheatmap(DE.miRNA.d3.expr, cluster_rows = T, cluster_cols = F, scale = "row", fontsize_row = 5, 
         display_numbers=F, number_color = "black",  fontsize_col = 4, 
         annotation_col = df, annotation_colors = anno_color,
         filename = "DE.miRNA.d3.expr.pdf", height = 2, width = 5, 
         show_rownames = F, show_colnames = F,border_color = NA,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(21))

### boxplot 0f iRBD, PD-iRBD and PD-no-iRBD
d4 <- cbind(d1, PD.no)
miRNA.set <- c("miR-10a-5p", "miR-4485-3p", "miR-369-5p", "miR-625-5p", "miR-874-3p")
boxp <- d4[rownames(d4) %in% miRNA.set, ]
mirna.num <- nrow(boxp)
boxp$gene <- factor(rownames(boxp),levels = rownames(boxp))
boxp.melt <- melt(boxp, id.vars="gene")

boxp.melt$group <- factor(c(rep("iRBD", mirna.num*56), rep("PD.iRBD", mirna.num*18), rep("PD.no", mirna.num*35)),
                          levels = c("iRBD", "PD.iRBD", "PD.no"))
my_comparisons <- list( c("iRBD", "PD.iRBD"), c("PD.iRBD", "PD.no"), c("iRBD", "PD.no") )
pdf("boxplot.pdf", width = 9,height = 6)
ggplot(boxp.melt, aes(x=group, y=value, color=group)) + geom_boxplot(outlier.shape = 18,outlier.size = 1)+ 
  labs(x="", y="log2(normalized read count)") + 
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")+
  facet_wrap(.~gene, scales = "free", ncol = 3) +theme(legend.position="none")+
  stat_compare_means( comparisons = my_comparisons, size=3, method = "wilcox.test")
