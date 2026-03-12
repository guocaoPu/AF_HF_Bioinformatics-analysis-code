rm(list = ls())
options(stringsAsFactors = F)

###DEGs
if(T){
exprset <- read.table('GSE79768_series_matrix.txt',sep = '\t',header = T,row.names = 1,check.names = F)
pdata <- read.table('GSE79768_Group.txt',sep = '\t',header = T,check.names = F)
group <- pdata$Group
library(limma)
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
row.names(design) = colnames(exprset)
design
contrast.matrix <- makeContrasts(paste0(unique(group),collapse = '-'),levels = design)
contrast.matrix
fit <- lmFit(exprset,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput =topTable(fit2,coef = 1,n = Inf)#筛选条件
nrDEG = na.omit(tempOutput)
head(nrDEG)
write.table(nrDEG,'GSE79768_DEGs.txt',sep = '\t')


library(oligo)
library(affy)
celfile = list.celfiles('GSE14975/',full.name = T)
celfile
raw.cell = oligo::read.celfiles(celfile)
#image(raw.cell)
rawData <- raw.cell
exprs(rawData)[1:6,1:6]
#rna normalize
boxplot(rawData)
normData <- oligo::rma(rawData,background = T,normalize = T)
normData
boxplot(normData)
head(exprs(normData))
data <- exprs(normData)
write.table(data,'GSE14975_Normalize.txt',sep = '\t')

exprset <- read.table('GSE14975_Normalize.txt',sep = '\t',header = T,row.names = 1,check.names = F)
group <- c(rep('AF',5),rep('Ctrl',5))
library(limma)
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
row.names(design) = colnames(exprset)
design
contrast.matrix <- makeContrasts(paste0(unique(group),collapse = '-'),levels = design)
contrast.matrix
fit <- lmFit(exprset,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput =topTable(fit2,coef = 1,n = Inf)#筛选条件
nrDEG = na.omit(tempOutput)
head(nrDEG)
write.table(nrDEG,'GSE14975_DEGs.txt',sep = '\t')

exprset <- read.table('GSE141910_Count.txt',sep = '\t',header = T,row.names = 1,check.names = F)
pdata <- read.table('GSE141910_Group.txt',sep = '\t',header = T,check.names = F)
exprset <- round(2^exprset, digits = 0)
group <- pdata$Group
library(edgeR)
dge <- DGEList(counts = exprset,group = group)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
et <- exactTest(dge)
topTags(et)
order_tags <- topTags(et,n=100000)
order_tags <- data.frame(order_tags)
write.table(order_tags,file="GSE141910_DEGs.txt",sep='\t')
write.table(log2(dge$pseudo.counts + 1),'GSE141910_Normalize.txt',sep = '\t')


library(affy)
celfile = list.celfiles('GSE17800/',full.name = T)
celfile
raw.cell = oligo::read.celfiles(celfile)
#image(raw.cell)
rawData <- raw.cell
exprs(rawData)[1:6,1:6]
pData(rawData)
pData(rawData)$filename <- sampleNames(rawData)
boxplot(rawData)
normData <- oligo::rma(rawData,background = T,normalize = T)
normData
boxplot(normData)
head(exprs(normData))
data <- exprs(normData)
write.table(data,'GSE17800_Normalized.txt',sep = '\t')

exprset <- read.table('GSE17800_Normalized.txt',sep = '\t',header = T,row.names = 1,check.names = F)
pdata <- read.table('GSE17800_Group.txt',sep = '\t',header = T,check.names = F)
group <- pdata$Group
library(limma)
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
row.names(design) = colnames(exprset)
design
contrast.matrix <- makeContrasts(paste0(unique(group),collapse = '-'),levels = design)
contrast.matrix
fit <- lmFit(exprset,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput =topTable(fit2,coef = 1,n = Inf)
nrDEG = na.omit(tempOutput)
head(nrDEG)
write.table(nrDEG,'GSE17800_DEGs.txt',sep = '\t')

Count <- list.files('GSE141910/',pattern="*.gz")
setwd(dir = 'GSE141910/')
sample <-  lapply(Count, function(x) read.table(x,sep = ",",comment.char = "#", stringsAsFactors = F,header = F, fill=TRUE))
sample1 <- sample[[1]]
y <- sample[[1]]
x <- sample[[2]]
for (i in 2:366) {
  x <- sample[[i]]
  table(y$V1%in%x$V1) 
  x = x[match(y$V1,x$V1),] 
  y <- cbind(y,x$V2)
}
colnames(y) <- Count
write.table(y,'C:/Myproject/20240301_AF_HF/GSE141910_Count.txt',sep = '\t')

exprset <- read.table('GSE141910_Count.txt',sep = '\t',header = T,row.names = 1,check.names = F)
pdata <- read.table('GSE141910_Group.txt',sep = '\t',header = T,check.names = F)
group <- pdata$Group
library(limma)
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
row.names(design) = colnames(exprset)
design
contrast.matrix <- makeContrasts(paste0(unique(group),collapse = '-'),levels = design)
contrast.matrix
fit <- lmFit(exprset,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput =topTable(fit2,coef = 1,n = Inf)
nrDEG = na.omit(tempOutput)
head(nrDEG)
write.table(nrDEG,'GSE141910_DEGs.txt',sep = '\t')

exprset <- read.table('GSE115574_Normalize.txt',sep = '\t',header = T,row.names = 1,check.names = F)
pdata <- read.table('GSE115574_Group.txt',sep = '\t',header = T,check.names = F)
group <- pdata$Group
library(limma)
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
row.names(design) = colnames(exprset)
design
contrast.matrix <- makeContrasts(paste0(unique(group),collapse = '-'),levels = design)
contrast.matrix
fit <- lmFit(exprset,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput =topTable(fit2,coef = 1,n = Inf)
nrDEG = na.omit(tempOutput)
head(nrDEG)
write.table(nrDEG,'GSE115574_DEGs.txt',sep = '\t')

library(ggplot2)
data <- read.table('C:/Myproject/20240301_AF_HF/AF_HF_Result20240627/1.AF_LV_DEGs/GSE115574_DEGs.csv',sep = ',',header = T)
data$logFC
ggplot(data,aes(logFC,-log10(P.Value))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(0),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(P.Value),color = -log10(P.Value))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

library(ggplot2)
data <- read.table('C:/Myproject/20240301_AF_HF/AF_HF_Result20240627/1.AF_LV_DEGs/GSE79768_DEGs.csv',sep = ',',header = T)
data$logFC
ggplot(data,aes(logFC,-log10(P.Value))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(0),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(P.Value),color = -log10(P.Value))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

library(ggplot2)
data <- read.table('C:/Myproject/20240301_AF_HF/AF_HF_Result20240627/5.HF_LV_DEGs/GSE141910_DEGs.txt',sep = '\t',header = T)
data$logFC
ggplot(data,aes(logFC,-log10(P.Value))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(0),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(P.Value),color = -log10(P.Value))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

library(ggplot2)
data <- read.table('C:/Myproject/20240301_AF_HF/AF_HF_Result20240627/5.HF_LV_DEGs/GSE17800_DEGs.csv',sep = ',',header = T)
data$logFC
ggplot(data,aes(-logFC,-log10(P.Value))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(0),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(P.Value),color = -log10(P.Value))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())
}

df <- read.table("GSE115574_ROC.txt",head=T,sep="\t",check.names = F)
library("pROC")
mycol <- c("#e0bc58","#64abc0","#fab37f","#e98741","#8fc0dc",
           "#967568","#f2d3ca","#eebd85","#82c785","#edeaa4",
           "#cdaa9f","#794976","#bcacd3","#889b5d","#4e9592",
           "#dbad5f","#64ae79","#ac5092")
pdf("ROC.pdf",height=8,width=8)
auc.out <- c()
x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=T, 
              ci=TRUE, 
              main="",
              print.thres="best", 
              col=mycol[2],
              lwd=2,
              legacy.axes=T)
ci.lower <- round(as.numeric(x$ci[1]),3) 
ci.upper <- round(as.numeric(x$ci[3]),3) 
auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
auc.out <- rbind(auc.out,auc.ci)

for (i in 3:ncol(df)){
  x <- plot.roc(df[,1],df[,i],
                add=T, 
                smooth=T,
                ci=TRUE,
                col=mycol[i],
                lwd=2,
                legacy.axes=T)
  
  ci.lower <- round(as.numeric(x$ci[1]),3)
  ci.upper <- round(as.numeric(x$ci[3]),3)
  
  auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
}

p.out <- c()
for (i in 2:(ncol(df)-1)){
  for (j in (i+1):ncol(df)){
    p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
    p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
    p.out <- rbind(p.out,p.tmp)
  }
}

p.out <- as.data.frame(p.out)
colnames(p.out) <- c("ROC1","ROC2","p.value")
write.table(p.out,"pvalue_output.xls",sep="\t",quote=F,row.names = F,col.names = T)

auc.out <- as.data.frame(auc.out)
colnames(auc.out) <- c("Name","AUC","AUC CI")
write.table(auc.out,"auc_output.xls",sep="\t",quote = F,row.names = F,col.names = T)

legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
legend("bottomright", 
       legend=legend.name,
       col = mycol[2:length(df)],
       lwd = 2,
       bty="n")
dev.off()

df <- read.table("GSE17800_ROC.txt",head=T,sep="\t",check.names = F)
library("pROC")
mycol <- c("#e0bc58","#64abc0","#fab37f","#e98741","#8fc0dc",
           "#967568","#f2d3ca","#eebd85","#82c785","#edeaa4",
           "#cdaa9f","#794976","#bcacd3","#889b5d","#4e9592",
           "#dbad5f","#64ae79","#ac5092")
pdf("ROC.pdf",height=8,width=8)
auc.out <- c()
x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=T, 
              ci=TRUE, 
              main="",
              #print.thres="best",
              col=mycol[2],
              lwd=2, 
              legacy.axes=T)
ci.lower <- round(as.numeric(x$ci[1]),3) 
ci.upper <- round(as.numeric(x$ci[3]),3) 
auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
auc.out <- rbind(auc.out,auc.ci)

for (i in 3:ncol(df)){
  x <- plot.roc(df[,1],df[,i],
                add=T, 
                smooth=T,
                ci=TRUE,
                col=mycol[i],
                lwd=2,
                legacy.axes=T)
  
  ci.lower <- round(as.numeric(x$ci[1]),3)
  ci.upper <- round(as.numeric(x$ci[3]),3)
  
  auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
}

p.out <- c()
for (i in 2:(ncol(df)-1)){
  for (j in (i+1):ncol(df)){
    p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
    p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
    p.out <- rbind(p.out,p.tmp)
  }
}
p.out <- as.data.frame(p.out)
colnames(p.out) <- c("ROC1","ROC2","p.value")
write.table(p.out,"pvalue_output.xls",sep="\t",quote=F,row.names = F,col.names = T)

auc.out <- as.data.frame(auc.out)
colnames(auc.out) <- c("Name","AUC","AUC CI")
write.table(auc.out,"auc_output.xls",sep="\t",quote = F,row.names = F,col.names = T)
legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
legend("bottomright", 
       legend=legend.name,
       col = mycol[2:length(df)],
       lwd = 2,
       bty="n")
dev.off()


###GSE115574_WGCNA
if(T){
  exprset <- read.table('GSE115574_WGCNA.txt',sep = '\t',header = T,row.names = 1,check.names = F)
  pdata <- read.table('GSE115574_WGCNA_Group.txt',sep = '\t',header = T,row.names = 1,check.names = F)
  library(WGCNA)
  library(tidyverse)
  library(corrplot)
  library(pheatmap)
 
  enableWGCNAThreads()
  options(stringsAsFactors = FALSE)
  
  datExpr <- t(exprset)
  
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      print(paste("remove", sum(!gsg$goodGenes), "low sample"))
      datExpr <- datExpr[, gsg$goodGenes]
    }
    if (sum(!gsg$goodSamples) > 0) {
      print(paste("remove", sum(!gsg$goodSamples), "low sample"))
      datExpr <- datExpr[gsg$goodSamples, ]
    }
  }
  
  sampleTree <- hclust(dist(datExpr), method = "average")

  pdf("WGCNA_sample_clustering.pdf", width = 12, height = 8)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", 
       sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()

  powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))

  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

  pdf("WGCNA_soft_threshold_selection.pdf", width = 10, height = 6)
  par(mfrow = c(1, 2))
  cex1 <- 0.9

  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  
  dev.off()
  

  softPower <- sft$powerEstimate
  if (is.na(softPower)) {
    softPower <- 6  
  }
  print(paste("选择的软阈值功率:", softPower))
  
  net <- blockwiseModules(
    datExpr,
    power = softPower,
    TOMType = "unsigned", 
    minModuleSize = 30,          
    deepSplit = 2,               
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25,      
    numericLabels = TRUE,         
    saveTOMs = TRUE,
    saveTOMFileBase = "WGCNATOM",
    verbose = 3
  )
  
  moduleColors <- labels2colors(net$colors)

  pdf("WGCNA_module_dendrogram.pdf", width = 12, height = 9)
  plotDendroAndColors(
    net$dendrograms[[1]], 
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05
  )
  dev.off()
  
  MEs <- net$MEs
  MEs <- orderMEs(MEs) 
  
  pdata$COPD <- ifelse(pdata$Group == "COPD", 1, 0)
  pdata$Ctrl <- ifelse(pdata$Group == "Ctrl", 1, 0)
  traits <- pdata %>% 
    select(where(is.numeric))  

  moduleTraitCor <- cor(MEs, traits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  pdf("WGCNA_module_trait_relationships.pdf", width = 8, height = 6)
  
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(traits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  dev.off()

  trait <- traits[,1]
  traitName <- colnames(traits)[1]
  
  GS <- as.numeric(cor(datExpr, trait, use = "p"))
  GeneSignificance <- abs(GS)

  modNames <- substring(names(MEs), 3)
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
  
  names(geneModuleMembership) <- paste("MM", modNames, sep = "")
  names(MMPvalue) <- paste("p.MM", modNames, sep = "")
  
  geneInfo <- data.frame(
    Gene = colnames(datExpr),
    Module = moduleColors,
    GeneSignificance = GS,
    stringsAsFactors = FALSE
  )
  
  for (mod in modNames) {
    geneInfo[[paste0("MM.", mod)]] <- geneModuleMembership[[paste0("MM", mod)]]
    geneInfo[[paste0("p.MM.", mod)]] <- MMPvalue[[paste0("p.MM", mod)]]
  }
  
  geneInfo <- geneInfo[order(-abs(geneInfo$GeneSignificance)), ]

  moduleInfo <- data.frame(
    Module = unique(moduleColors),
    Size = as.numeric(table(moduleColors)),
    stringsAsFactors = FALSE
  )
  
  write.csv(geneInfo, "WGCNA_gene_info.csv", row.names = FALSE)
  write.csv(moduleInfo, "WGCNA_module_info.csv", row.names = FALSE)
  write.csv(moduleTraitCor, "WGCNA_module_trait_correlation.csv")
  write.csv(moduleTraitPvalue, "WGCNA_module_trait_pvalue.csv")

  corThreshold <- 0.5  
  pvalThreshold <- 0.05 
  
  significantModules <- which(apply(moduleTraitPvalue < pvalThreshold & 
                                      abs(moduleTraitCor) > corThreshold, 1, any))
  
  if (length(significantModules) > 0) {
    for (i in significantModules) {
      module <- modNames[i]
      moduleColor <- moduleColors[net$colors == as.numeric(module)]
      
      moduleGenes <- colnames(datExpr)[moduleColors == moduleColor]

      write.table(moduleGenes, 
                  paste0("module_", moduleColor, "_genes.txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)

      pdf(paste0("module_", moduleColor, "_connectivity.pdf"), width = 8, height = 6)
      verboseBarplot(geneModuleMembership[[paste0("MM", module)]],
                     moduleColors,
                     main = paste("Module membership in", moduleColor, "module"),
                     xlab = "", ylab = "Module Membership")
      dev.off()
    }
  }

  if (length(significantModules) > 0) {
    module <- modNames[significantModules[1]]
    moduleColor <- labels2colors(as.numeric(module))
    
    probes <- colnames(datExpr)
    inModule <- (moduleColors == moduleColor)
    modProbes <- probes[inModule]
    
    modTOM <- TOM[inModule, inModule]
    dimnames(modTOM) <- list(modProbes, modProbes)
    
    cyt <- exportNetworkToCytoscape(
      modTOM,
      edgeFile = paste("CytoscapeInput-edges-", moduleColor, ".txt", sep = ""),
      nodeFile = paste("CytoscapeInput-nodes-", moduleColor, ".txt", sep = ""),
      weighted = TRUE,
      threshold = 0.02,
      nodeNames = modProbes,
      nodeAttr = moduleColors[inModule]
    )
  }
  
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  heatmap_data <- moduleTraitCor
  head(heatmap_data)
  module_color_mapping <- data.frame(
    ME_number = paste0("ME", 1:length(unique(moduleColors))),
    Module_color = unique(moduleColors)
  )
  
  print(module_color_mapping)
  
  print(rownames(heatmap_data))
  
  me_to_color <- setNames(module_color_mapping$Module_color, module_color_mapping$ME_number)
  
  if (all(rownames(heatmap_data) %in% names(me_to_color))) {
    rownames(heatmap_data) <- me_to_color[rownames(heatmap_data)]
    rownames(moduleTraitPvalue) <- me_to_color[rownames(moduleTraitPvalue)]
    
    print(rownames(heatmap_data))
  } else {
    warning("ME nocolor")
    matching_rows <- rownames(heatmap_data) %in% names(me_to_color)
    rownames(heatmap_data)[matching_rows] <- me_to_color[rownames(heatmap_data)[matching_rows]]
    rownames(moduleTraitPvalue)[matching_rows] <- me_to_color[rownames(moduleTraitPvalue)[matching_rows]]
  }

  annotation_matrix <- matrix("", nrow = nrow(heatmap_data), ncol = 3)
  rownames(annotation_matrix) <- rownames(heatmap_data)
  colnames(annotation_matrix) <- colnames(heatmap_data)

  for (i in 1:nrow(heatmap_data)) {
    pval <- moduleTraitPvalue[i, 1]
    if (pval < 0.001) {
      annotation_matrix[i, 1] <- "***"
    } else if (pval < 0.01) {
      annotation_matrix[i, 1] <- "**"
    } else if (pval < 0.05) {
      annotation_matrix[i, 1] <- "*"
    }
  }
  
  pdf("WGCNA_module_COPD_correlation_heatmap.pdf", width = 6, height = 8)
  
  pheatmap(heatmap_data,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks = seq(-1, 1, length.out = 101),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           display_numbers = annotation_matrix,
           number_color = "white",
           number_format = "%.2f",
           fontsize_number = 9,  
           fontsize_row = 10,
           fontsize_col = 11,
           angle_col = 45, 
           main = "Module-Trait Correlation: COPD Status",
           border_color = NA,
           cellwidth = 45,  
           cellheight = 16,
           legend = TRUE,
           show_rownames = TRUE,
           gaps_col = NULL,
           treeheight_row = 30,
           treeheight_col = 0)
  
  dev.off()
}


###GSE17800_WGCNA
if(T){
  exprset <- read.table('GSE17800_WGCNA.txt',sep = '\t',header = T,row.names = 1,check.names = F)
  pdata <- read.table('GSE17800_WGCNA_Group.txt',sep = '\t',header = T,row.names = 1,check.names = F)
  library(WGCNA)
  library(tidyverse)
  library(corrplot)
  library(pheatmap)
 
  enableWGCNAThreads()
  options(stringsAsFactors = FALSE)
  
  datExpr <- t(exprset)
  
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      print(paste("remove", sum(!gsg$goodGenes), "low sample"))
      datExpr <- datExpr[, gsg$goodGenes]
    }
    if (sum(!gsg$goodSamples) > 0) {
      print(paste("remove", sum(!gsg$goodSamples), "low sample"))
      datExpr <- datExpr[gsg$goodSamples, ]
    }
  }
  
  sampleTree <- hclust(dist(datExpr), method = "average")

  pdf("WGCNA_sample_clustering.pdf", width = 12, height = 8)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", 
       sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()

  powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))

  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

  pdf("WGCNA_soft_threshold_selection.pdf", width = 10, height = 6)
  par(mfrow = c(1, 2))
  cex1 <- 0.9

  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  
  dev.off()
  

  softPower <- sft$powerEstimate
  if (is.na(softPower)) {
    softPower <- 6  
  }
  print(paste("选择的软阈值功率:", softPower))
  
  net <- blockwiseModules(
    datExpr,
    power = softPower,
    TOMType = "unsigned", 
    minModuleSize = 30,          
    deepSplit = 2,               
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25,      
    numericLabels = TRUE,         
    saveTOMs = TRUE,
    saveTOMFileBase = "WGCNATOM",
    verbose = 3
  )
  
  moduleColors <- labels2colors(net$colors)

  pdf("WGCNA_module_dendrogram.pdf", width = 12, height = 9)
  plotDendroAndColors(
    net$dendrograms[[1]], 
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05
  )
  dev.off()
  
  MEs <- net$MEs
  MEs <- orderMEs(MEs) 
  
  pdata$COPD <- ifelse(pdata$Group == "COPD", 1, 0)
  pdata$Ctrl <- ifelse(pdata$Group == "Ctrl", 1, 0)
  traits <- pdata %>% 
    select(where(is.numeric))  

  moduleTraitCor <- cor(MEs, traits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  pdf("WGCNA_module_trait_relationships.pdf", width = 8, height = 6)
  
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(traits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  dev.off()

  trait <- traits[,1]
  traitName <- colnames(traits)[1]
  
  GS <- as.numeric(cor(datExpr, trait, use = "p"))
  GeneSignificance <- abs(GS)

  modNames <- substring(names(MEs), 3)
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
  
  names(geneModuleMembership) <- paste("MM", modNames, sep = "")
  names(MMPvalue) <- paste("p.MM", modNames, sep = "")
  
  geneInfo <- data.frame(
    Gene = colnames(datExpr),
    Module = moduleColors,
    GeneSignificance = GS,
    stringsAsFactors = FALSE
  )
  
  for (mod in modNames) {
    geneInfo[[paste0("MM.", mod)]] <- geneModuleMembership[[paste0("MM", mod)]]
    geneInfo[[paste0("p.MM.", mod)]] <- MMPvalue[[paste0("p.MM", mod)]]
  }
  
  geneInfo <- geneInfo[order(-abs(geneInfo$GeneSignificance)), ]

  moduleInfo <- data.frame(
    Module = unique(moduleColors),
    Size = as.numeric(table(moduleColors)),
    stringsAsFactors = FALSE
  )
  
  write.csv(geneInfo, "WGCNA_gene_info.csv", row.names = FALSE)
  write.csv(moduleInfo, "WGCNA_module_info.csv", row.names = FALSE)
  write.csv(moduleTraitCor, "WGCNA_module_trait_correlation.csv")
  write.csv(moduleTraitPvalue, "WGCNA_module_trait_pvalue.csv")

  corThreshold <- 0.5  
  pvalThreshold <- 0.05 
  
  significantModules <- which(apply(moduleTraitPvalue < pvalThreshold & 
                                      abs(moduleTraitCor) > corThreshold, 1, any))
  
  if (length(significantModules) > 0) {
    for (i in significantModules) {
      module <- modNames[i]
      moduleColor <- moduleColors[net$colors == as.numeric(module)]
      
      moduleGenes <- colnames(datExpr)[moduleColors == moduleColor]

      write.table(moduleGenes, 
                  paste0("module_", moduleColor, "_genes.txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)

      pdf(paste0("module_", moduleColor, "_connectivity.pdf"), width = 8, height = 6)
      verboseBarplot(geneModuleMembership[[paste0("MM", module)]],
                     moduleColors,
                     main = paste("Module membership in", moduleColor, "module"),
                     xlab = "", ylab = "Module Membership")
      dev.off()
    }
  }

  if (length(significantModules) > 0) {
    module <- modNames[significantModules[1]]
    moduleColor <- labels2colors(as.numeric(module))
    
    probes <- colnames(datExpr)
    inModule <- (moduleColors == moduleColor)
    modProbes <- probes[inModule]
    
    modTOM <- TOM[inModule, inModule]
    dimnames(modTOM) <- list(modProbes, modProbes)
    
    cyt <- exportNetworkToCytoscape(
      modTOM,
      edgeFile = paste("CytoscapeInput-edges-", moduleColor, ".txt", sep = ""),
      nodeFile = paste("CytoscapeInput-nodes-", moduleColor, ".txt", sep = ""),
      weighted = TRUE,
      threshold = 0.02,
      nodeNames = modProbes,
      nodeAttr = moduleColors[inModule]
    )
  }
  
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  heatmap_data <- moduleTraitCor
  head(heatmap_data)
  module_color_mapping <- data.frame(
    ME_number = paste0("ME", 1:length(unique(moduleColors))),
    Module_color = unique(moduleColors)
  )
  
  print(module_color_mapping)
  
  print(rownames(heatmap_data))
  
  me_to_color <- setNames(module_color_mapping$Module_color, module_color_mapping$ME_number)
  
  if (all(rownames(heatmap_data) %in% names(me_to_color))) {
    rownames(heatmap_data) <- me_to_color[rownames(heatmap_data)]
    rownames(moduleTraitPvalue) <- me_to_color[rownames(moduleTraitPvalue)]
    
    print(rownames(heatmap_data))
  } else {
    warning("ME nocolor")
    matching_rows <- rownames(heatmap_data) %in% names(me_to_color)
    rownames(heatmap_data)[matching_rows] <- me_to_color[rownames(heatmap_data)[matching_rows]]
    rownames(moduleTraitPvalue)[matching_rows] <- me_to_color[rownames(moduleTraitPvalue)[matching_rows]]
  }

  annotation_matrix <- matrix("", nrow = nrow(heatmap_data), ncol = 3)
  rownames(annotation_matrix) <- rownames(heatmap_data)
  colnames(annotation_matrix) <- colnames(heatmap_data)

  for (i in 1:nrow(heatmap_data)) {
    pval <- moduleTraitPvalue[i, 1]
    if (pval < 0.001) {
      annotation_matrix[i, 1] <- "***"
    } else if (pval < 0.01) {
      annotation_matrix[i, 1] <- "**"
    } else if (pval < 0.05) {
      annotation_matrix[i, 1] <- "*"
    }
  }
  
  pdf("WGCNA_module_COPD_correlation_heatmap.pdf", width = 6, height = 8)
  
  pheatmap(heatmap_data,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks = seq(-1, 1, length.out = 101),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           display_numbers = annotation_matrix,
           number_color = "white",
           number_format = "%.2f",
           fontsize_number = 9,  
           fontsize_row = 10,
           fontsize_col = 11,
           angle_col = 45, 
           main = "Module-Trait Correlation: COPD Status",
           border_color = NA,
           cellwidth = 45,  
           cellheight = 16,
           legend = TRUE,
           show_rownames = TRUE,
           gaps_col = NULL,
           treeheight_row = 30,
           treeheight_col = 0)
  
  dev.off()
}

###scRNA data
if(T){
  library(Seurat)
  library(tidyverse)
  library(ggsci)
  library(clustree)
  library(harmony)
  library(patchwork)
  #library(MySeuratWrappers)
  library(VISION)
  library(monocle)
  library(scMetabolism)
  library(glmGamPoi)
  Donor1 <- Read10X(data.dir = 'GSE226314/H_ZC-11-292/filtered_feature_bc_matrix/')
  Donor2 <- Read10X(data.dir = 'GSE226314/TWCM-11-103/filtered_feature_bc_matrix/')
  Donor3 <- Read10X(data.dir = 'GSE226314/TWCM-11-192/filtered_feature_bc_matrix/')
  Donor4 <- Read10X(data.dir = 'GSE226314/TWCM-11-41/filtered_feature_bc_matrix/')
  Donor5 <- Read10X(data.dir = 'GSE226314/TWCM-11-74/filtered_feature_bc_matrix/')
  Donor6 <- Read10X(data.dir = 'GSE226314/TWCM-11-78/filtered_feature_bc_matrix/')
  Donor7 <- Read10X(data.dir = 'GSE226314/TWCM-11-82/filtered_feature_bc_matrix/')
  Donor8 <- Read10X(data.dir = 'GSE226314/TWCM-13-1/filtered_feature_bc_matrix/')
  Donor9 <- Read10X(data.dir = 'GSE226314/TWCM-13-104/filtered_feature_bc_matrix/')
  Donor10 <- Read10X(data.dir = 'GSE226314/TWCM-13-152/filtered_feature_bc_matrix/')
  Donor11 <- Read10X(data.dir = 'GSE226314/TWCM-13-168/filtered_feature_bc_matrix/')
  Donor12 <- Read10X(data.dir = 'GSE226314/TWCM-13-192/filtered_feature_bc_matrix/')
  Donor13 <- Read10X(data.dir = 'GSE226314/TWCM-13-80/filtered_feature_bc_matrix/')
  HF1 <- Read10X(data.dir = 'GSE226314/TWCM-190-R-pre/filtered_feature_bc_matrix/')
  HF2 <- Read10X(data.dir = 'GSE226314/TWCM-229-R-pre/filtered_feature_bc_matrix/')
  HF3 <- Read10X(data.dir = 'GSE226314/TWCM-239-R-pre/filtered_feature_bc_matrix/')
  HF4 <- Read10X(data.dir = 'GSE226314/TWCM-296-R-pre/filtered_feature_bc_matrix/')
  HF5 <- Read10X(data.dir = 'GSE226314/TWCM-359-NR-pre/filtered_feature_bc_matrix/')
  HF6 <- Read10X(data.dir = 'GSE226314/TWCM-363-NR-pre/filtered_feature_bc_matrix/')
  HF7 <- Read10X(data.dir = 'GSE226314/TWCM-373-NR-pre/filtered_feature_bc_matrix/')
  HF8 <- Read10X(data.dir = 'GSE226314/TWCM-376-NR-pre/filtered_feature_bc_matrix/')
  HF9 <- Read10X(data.dir = 'GSE226314/TWCM-378-NR-pre/filtered_feature_bc_matrix/')
  HF10 <- Read10X(data.dir = 'GSE226314/TWCM-388-NR-pre/filtered_feature_bc_matrix/')
  HF11 <- Read10X(data.dir = 'GSE226314/TWCM-397-NR-pre/filtered_feature_bc_matrix/')
  HF12 <- Read10X(data.dir = 'GSE226314/TWCM-410-NR-pre/filtered_feature_bc_matrix/')
  HF13 <- Read10X(data.dir = 'GSE226314/TWCM-463-R-pre/filtered_feature_bc_matrix/')
  AF1 <- Read10X(data.dir = 'GSE226314/77/')
  AF2 <- Read10X(data.dir = 'GSE226314/81/')
  AF3 <- Read10X(data.dir = 'GSE226314/85/')
  
  Donor1 <- CreateSeuratObject(counts = Donor1,project = 'Donor',min.cells = 3,min.features = 200)
  Donor2 <- CreateSeuratObject(counts = Donor2,project = 'Donor',min.cells = 3,min.features = 200)
  Donor3 <- CreateSeuratObject(counts = Donor3,project = 'Donor',min.cells = 3,min.features = 200)
  Donor4 <- CreateSeuratObject(counts = Donor4,project = 'Donor',min.cells = 3,min.features = 200)
  Donor5 <- CreateSeuratObject(counts = Donor5,project = 'Donor',min.cells = 3,min.features = 200)
  Donor6 <- CreateSeuratObject(counts = Donor6,project = 'Donor',min.cells = 3,min.features = 200)
  Donor7 <- CreateSeuratObject(counts = Donor7,project = 'Donor',min.cells = 3,min.features = 200)
  Donor8 <- CreateSeuratObject(counts = Donor8,project = 'Donor',min.cells = 3,min.features = 200)
  Donor9 <- CreateSeuratObject(counts = Donor9,project = 'Donor',min.cells = 3,min.features = 200)
  Donor10 <- CreateSeuratObject(counts = Donor10,project = 'Donor',min.cells = 3,min.features = 200)
  Donor11 <- CreateSeuratObject(counts = Donor11,project = 'Donor',min.cells = 3,min.features = 200)
  Donor12 <- CreateSeuratObject(counts = Donor12,project = 'Donor',min.cells = 3,min.features = 200)
  Donor13 <- CreateSeuratObject(counts = Donor13,project = 'Donor',min.cells = 3,min.features = 200)
  HF1 <- CreateSeuratObject(counts = HF1,project = 'HF',min.cells = 3,min.features = 200)
  HF2 <- CreateSeuratObject(counts = HF2,project = 'HF',min.cells = 3,min.features = 200)
  HF3 <- CreateSeuratObject(counts = HF3,project = 'HF',min.cells = 3,min.features = 200)
  HF4 <- CreateSeuratObject(counts = HF4,project = 'HF',min.cells = 3,min.features = 200)
  HF5 <- CreateSeuratObject(counts = HF5,project = 'HF',min.cells = 3,min.features = 200)
  HF6 <- CreateSeuratObject(counts = HF6,project = 'HF',min.cells = 3,min.features = 200)
  HF7 <- CreateSeuratObject(counts = HF7,project = 'HF',min.cells = 3,min.features = 200)
  HF8 <- CreateSeuratObject(counts = HF8,project = 'HF',min.cells = 3,min.features = 200)
  HF9 <- CreateSeuratObject(counts = HF9,project = 'HF',min.cells = 3,min.features = 200)
  HF10 <- CreateSeuratObject(counts = HF10,project = 'HF',min.cells = 3,min.features = 200)
  HF11 <- CreateSeuratObject(counts = HF11,project = 'HF',min.cells = 3,min.features = 200)
  HF12 <- CreateSeuratObject(counts = HF12,project = 'HF',min.cells = 3,min.features = 200)
  HF13 <- CreateSeuratObject(counts = HF13,project = 'HF',min.cells = 3,min.features = 200)
  AF1 <- CreateSeuratObject(counts = AF1,project = 'AF',min.cells = 3, min.features = 200)
  AF2 <- CreateSeuratObject(counts = AF2,project = 'AF',min.cells = 3, min.features = 200)
  AF3 <- CreateSeuratObject(counts = AF3,project = 'AF',min.cells = 3, min.features = 200)
  
  sce <- merge(Donor1,y=c(Donor2,Donor3,Donor4,Donor5,Donor6,Donor7,Donor8,Donor9,Donor10,Donor11,Donor12,Donor13,
                          HF1,HF2,HF3,HF4,HF5,HF6,HF7,HF8,HF9,HF10,HF11,HF12,HF13,
                          AF1,AF2,AF3))
  sce <- JoinLayers(sce)
  sce[["RNA"]] <- split(sce[["RNA"]], f = sce$orig.ident)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce,pattern = "^MT-")
  sce <- subset(sce,  subset =  
                  nFeature_RNA > 500 &
                  nFeature_RNA < 4000 & 
                  nCount_RNA > 1000 & 
                  nCount_RNA < 25000 &
                  percent.mt < 5 
  )
  
  sce <- sce%>%NormalizeData(normalization.method = "LogNormalize",scale.factor = 1e4)
  sce <- sce%>%FindVariableFeatures(selection.method = "vst",nfeatures = 2000) 
  
  sce <- sce%>%ScaleData(verbose = FALSE)
  sce <- sce%>%RunPCA(npcs = 50, verbose = FALSE)
  
  sce <- IntegrateLayers(object = sce, method = CCAIntegration, orig.reduction = "pca",
                         new.reduction = 'integrated.cca', verbose = FALSE)
  saveRDS(sce,"sce.cca.RDS")
  
  sce <- readRDS('sce.cca.RDS')
  dim.use <- 1:30
  sce <- sce%>%FindNeighbors(reduction = "integrated.cca", dims = dim.use)
  sce <- Seurat::FindClusters(sce,resolution = 0.6)
  sce <- Seurat::RunTSNE(sce,dims = dim.use,reduction = "integrated.cca")
  
  DimPlot(sce,label = T,reduction = 'tsne') + NoLegend() + labs(x = 'TSNE1', y = 'TSNE2') +
    theme(panel.border = element_rect(fill = NA,color = 'black',size = 1,linetype = 'solid'),
          axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank())
  
  DimPlot(sce,label = T,reduction = 'tsne',group.by = 'orig.ident',cols = c('#F7903D','#4D85BD','#FA7F6F')) + labs(x = 'TSNE1', y = 'TSNE2') +
    theme(panel.border = element_rect(fill = NA,color = 'black',size = 1,linetype = 'solid'),
          axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank())
  
  sce <- JoinLayers(sce)
  DefaultAssay(sce) <- "RNA"
  all.markers <- FindAllMarkers(sce, 
                                only.pos = TRUE,# 同样，默认情况下，此函数将返回正向和负向表达变化的基因。我们可以使用only.pos观察阳性变化
                                logfc.threshold = 0.25,
                                # min.diff.pct = 0.3, #某个cluster中表达基因的细胞百分比与其他所有clusters中该表达基因的细胞百分比之间的最小百分比差异。
                                min.pct = 0.1 #在两个clusters中任一clusters的最少比例细胞中检测到的基因。默认值为0.1
  )
  write.table(all.markers,file=paste0("total_marker_genes_",max(dim.use),"PC.txt"),sep="\t",quote = F,row.names = F)
  
  saveRDS(sce,'scRNA_Result.rds')
  
  sce <- readRDS('scRNA_Result.rds')
  scemetadata <- read.csv('scemetadata.csv',header = T,check.names = F)
  sce@meta.data$Celltype <- scemetadata$Celltype
  
  library(scRNAtoolVis)
  clusterCornerAxes(object = sce,
                    reduction = 'tsne',
                    noSplit = T,
                    cornerTextSize = 3.5,
                    themebg = 'bwCorner',
                    cicAlpha = 0.2,
                    nbin = 200,
                    pSize = 0.05)
  
  clusterCornerAxes(object = sce,
                    reduction = 'tsne',
                    noSplit = T,
                    cornerTextSize = 3.5,
                    themebg = 'bwCorner',
                    cicAlpha = 0.2,
                    nbin = 200,
                    clusterCol='Celltype',
                    pSize = 0.05)
  
  clusterCornerAxes(object=sce,
                    reduction='tsne',
                    noSplit=F,
                    groupFacet='orig.ident',
                    aspect.ratio=1,
                    relLength=0.5,
                    lineTextcol='grey50',
                    clusterCol='Celltype',
                    pSize = 0.05)
  
  #CUX1	CACNA1E	SEC62	RUSC1.AS1	KDELR3	NRP2	CSNK2A1	NES	LAMA4	COTL1	TLR4	SERPINB8	WASL	POLR1C	EXT1	PLXDC2
  selectgene <-c('CUX1','CACNA1E','SEC62','RUSC1-AS1','KDELR3','NRP2','CSNK2A1','NES','LAMA4','COTL1','TLR4','SERPINB8','WASL','POLR1C','EXT1','PLXDC2') 
  
  AverageHeatmap(objec = sce,
                 markerGene = selectgene,
                 clusterAnnoName = F,
                 showRowNames = F,#不展示所有基因
                 markGenes = selectgene,#标记基因
                 cluster_columns = F,
                 cluster_rows = T,
                 width=8,height=10
  )
  
  Idents(sce) = sce@meta.data$Celltype
  DotPlot(sce, features = selectgene)+coord_flip()+theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
  }

###Cellchat
if(T){
  library(tidyverse)
  library(Seurat)
  library(CellChat)
  sce <- readRDS('scRNA_Result.rds')
  scemetadata <- read.csv('scemetadata.csv',header = T,row.names = 1,check.names = F)
  colnames(scemetadata) <- c('Group1','Count1','Feature1','mt','snn','cluster1','Celltype')
  sce <- AddMetaData(sce,scemetadata)
  table(sce@meta.data$Group1)
  
  sceAF <- subset(x = sce,subset = Group1 == "AF")
  data.input <- sceAF@assays$RNA$data
  head(data.input[1:40,1:40])
  meta <- sceAF@meta.data
  cellchatAF <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")
  cellchat <- cellchatAF
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  
  write.csv(df.net, "AF_cell-cell_communications.all.csv")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  cellchatAF <- cellchat
  saveRDS(cellchatAF,'AF_cellchat.RDS')
  
  sceDonor <- subset(x = sce,subset = Group1 == "Donor")
  data.input <- sceDonor@assays$RNA$data
  meta <- sceDonor@meta.data
  cellchatDonor <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")
  cellchat <- cellchatDonor
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  
  write.csv(df.net, "Donor_cell-cell_communications.all.csv")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  cellchatDonor <- cellchat
  saveRDS(cellchatDonor,'Donor_cellchat.RDS')
  
  sceHF <- subset(x = sce,subset = Group1 == "HF")
  data.input <- sceHF@assays$RNA$data
  meta <- sceHF@meta.data
  cellchatHF <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")
  cellchat <- cellchatHF
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  
  write.csv(df.net, "HF_cell-cell_communications.all.csv")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  cellchatHF <- cellchat
  saveRDS(cellchatHF,'HF_cellchat.RDS')
  
  cellchat.list <- list(HF = cellchatHF,AF = cellchatAF,Donor = cellchatDonor)
  cellchat <- mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = TRUE)
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "count")
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
  p <- gg1 + gg2
  p
  ggsave("Overview number strength.pdf",p,width = 6, height = 4)
  
  par(mfrow = c(1,3))
  netVisual_diffInteraction(cellchat,weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  
  par(mfrow = c(1,1))
  h1 <- netVisual_heatmap(cellchat)
  h2 <- netVisual_heatmap(cellchat, measure = "weight")
  h1 + h2
  
  
  # 获取最大权重（用于统一比例）
  weight.max <- getMaxWeight(cellchat.list, attribute = c("idents", "count"))
  # 设置输出目录（可选）
  output_dir <- "CellChat_plots"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # 循环导出每个对象的 PDF
  for (i in 1:length(cellchat.list)) {
    # 设置 PDF 文件名（按组名命名）
    pdf_file <- file.path(output_dir, paste0("CellChat_Network_", names(cellchat.list)[i], ".pdf"))
    
    # 开启 PDF 设备
    pdf(pdf_file, width = 8, height = 6)  # 可调整宽高
    
    # 绘制网络图
    netVisual_circle(
      cellchat.list[[i]]@net$count,
      weight.scale = TRUE,
      label.edge = FALSE,
      edge.weight.max = weight.max[2],
      edge.width.max = 6,
      title.name = paste0("Number of interactions ~ ", names(cellchat.list)[i])
    )
    
    # 关闭 PDF 设备
    dev.off()
    
    # 提示保存成功
    message("Saved: ", pdf_file)
  }
  
  pathway.union <- union(cellchat.list[[1]]@netP$pathways,cellchat.list[[2]]@netP$pathways)
  pathway.union <- union(pathway.union,cellchat.list[[3]]@netP$pathways)
  
  pathway_scores_list <- list()
  
  for(i in 1:length(cellchat.list)) {
    # 提取通路活性得分矩阵
    score_mat <- cellchat.list[[i]]@netP$prob  # 或者使用其他适当的矩阵
    # 计算每个通路的平均得分
    pathway_scores_list[[i]] <- rowMeans(score_mat)
  }
  
  # 创建数据框计算变化
  pathway_changes <- data.frame(
    pathway = pathway.union
  )
  
  # 填充每个样本的得分
  for(i in 1:length(cellchat.list)) {
    pathway_changes[[paste0("sample", i)]] <- 
      pathway_scores_list[[i]][match(pathway.union, names(pathway_scores_list[[i]]))]
  }
  
  # 计算方差（变化程度）
  pathway_changes$variance <- apply(pathway_changes[, -1], 1, var, na.rm = TRUE)
  
  # 按方差排序，选择前40个
  top40_pathways <- pathway_changes %>%
    arrange(desc(variance)) %>%
    head(40) %>%
    pull(pathway)
  
  # 使用筛选后的通路绘制热图
  ht1 = netAnalysis_signalingRole_heatmap(cellchat.list[[1]], pattern = "all", signaling = top40_pathways,
                                          title = names(cellchat.list)[1], width = 6, height = 16)
  ht2 = netAnalysis_signalingRole_heatmap(cellchat.list[[2]], pattern = "all", signaling = top40_pathways,
                                          title = names(cellchat.list)[2], width = 6, height = 16)
  ht3 = netAnalysis_signalingRole_heatmap(cellchat.list[[3]], pattern = "all", signaling = top40_pathways,
                                          title = names(cellchat.list)[3], width = 6, height = 16)
  
  ht1 + ht2 + ht3
  
  pdf("CellChat_plots/Cellchat_HF_Overall signaling patterns.pdf", width = 6, height = 8)
  draw(ht1)  # 绘制 Heatmap
  dev.off()  # 关闭图形设备
  
  pdf("CellChat_plots/Cellchat_AF_Overall signaling patterns.pdf", width = 6, height = 8)
  draw(ht1)  # 绘制 Heatmap
  dev.off()  # 关闭图形设备
  
  pdf("CellChat_plots/Cellchat_Donor_Overall signaling patterns.pdf", width = 6, height = 8)
  draw(ht1)  # 绘制 Heatmap
  dev.off()  # 关闭图形设备
  
  if(T){
    # 2. 设置图形参数以获得更好的视觉效果
    par(mfrow = c(1, 3), mar = c(2, 2, 4, 2))  # 调整边距
    
    # 3. 计算所有组的最大交互数以统一比例
    max_count <- max(sapply(cellchat.list, function(x) max(x@net$count)))
    
    # 4. 绘制三个组的网络图
    for (i in 1:length(cellchat.list)) {
      # 获取当前组的细胞类型数量
      groupSize <- as.numeric(table(cellchat.list[[i]]@idents))
      
      # 绘制网络图
      netVisual_circle(cellchat.list[[i]]@net$count, 
                       vertex.weight = groupSize,
                       weight.scale = TRUE, 
                       label.edge = FALSE,
                       edge.weight.max = max_count,  # 使用统一的最大值
                       title.name = paste0(names(cellchat.list)[i], 
                                           "\n(Total: ", sum(cellchat.list[[i]]@net$count), ")"))
    }
    # 5. 重置图形参数
    par(mfrow = c(1, 1))
  }
  
  if(T){
    # 2. 设置图形参数以获得更好的视觉效果
    par(mfrow = c(1, 3), mar = c(2, 2, 4, 2))  # 调整边距
    
    # 3. 计算所有组的最大交互数以统一比例
    max_count <- max(sapply(cellchat.list, function(x) max(x@net$count)))
    
    # 4. 绘制三个组的网络图
    for (i in 1:length(cellchat.list)) {
      # 获取当前组的细胞类型数量
      groupSize <- as.numeric(table(cellchat.list[[i]]@idents))
      
      # 绘制网络图
      netVisual_circle(cellchat.list[[i]]@net$weight, 
                       vertex.weight = groupSize,
                       weight.scale = TRUE, 
                       label.edge = FALSE,
                       edge.weight.max = max_count,  # 使用统一的最大值
                       title.name = paste0(names(cellchat.list)[i], 
                                           "\n(Total: ", sum(cellchat.list[[i]]@net$weight), ")"))
    }
    # 5. 重置图形参数
    par(mfrow = c(1, 1))
  }
  
  ###hf.vs.af
  if(T){
    cellchat.list <- list(HF = cellchatHF,AF = cellchatAF)
    cellchat <- mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = TRUE)   
    par(mfrow = c(1,2))
    h1 <- netVisual_heatmap(cellchat)
    h2 <- netVisual_heatmap(cellchat, measure = "weight")
    p1 <- h1 + h2
    pdf("CellChat_plots/Cellchat_hf.vs.af_netVisual_heatmap.pdf", width = 10, height = 8)
    draw(p1)  # 绘制 Heatmap
    dev.off()  # 关闭图形设备
    
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T,do.stat = TRUE)
    ggsave("CellChat_plots/Cellchat_hf.vs.af.stacked.Compare_pathway_strengh.pdf",gg1,width = 4, height = 8)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F,do.stat = TRUE)
    ggsave("CellChat_plots/Cellchat_hf.vs.af.Barplot.Compare_pathway_strengh.pdf",gg2,width = 4, height = 8)
  }
}

library(ClusterGVis)
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)
sce@active.ident <- as.factor(sce@meta.data$Celltype)
Idents(sce) <- sce@meta.data[["Celltype"]]
pbmc.markers.all <- Seurat::FindAllMarkers(sce, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
pbmc.markers.all <- read.table('pbmc.markers.all.txt',sep = '\t',header = T,row.names = 1,check.names = F)
pbmc.markers <- pbmc.markers.all %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 40, wt = avg_log2FC)
st.data <- prepareDataFromscRNA(object = sce,diffData = pbmc.markers,showAverage = TRUE)
enrich <- enrichCluster(object = st.data,OrgDb = org.Hs.eg.db,type = "BP",organism = "hsa",pvalueCutoff = 1,topn = 3, seed = 0627)
markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,replace = F)]
visCluster(object = st.data, plot.type = "line")

pdf('sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 70,
           markGenes = markGenes,
           cluster.order = c(1:13))
dev.off()

pdf('sc2.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 70,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:12),
           go.col = rep(jjAnno::useMyCol("stallion",n = 12),each = 3),
           add.bar = T)
dev.off()

saveRDS(sce,"heatmap.rds")
write.table(pbmc.markers.all,'pbmc.markers.all.txt',sep = '\t')


