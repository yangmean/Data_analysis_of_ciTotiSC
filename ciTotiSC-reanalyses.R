library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(scales)
library(xlsx)
library(limma)
library(org.Mm.eg.db)
library(factoextra)
library(FactoMineR)
library(sva)

##################Figure 1 Related to Fig.1h###############################
gene.fpkm <- read.table("genes.fpkm_table.csv", sep = ",", header = T, row.names = 1, check.names = F)
gene.fpkm.log <- log2(gene.fpkm + 1)
gene.fpkm.log <- gene.fpkm.log[, c(1,3,5,7,9,10,12,14:20,25)]
Maternal.gene <- c("Ccdc122", "Cpb2", "Echdc3", "Gfra1", "Obox2", "Rhobtb1", "Slco3a1", "Thap11", "Ust", "Hoxa1")
ZGA.gene <- c("Gm13078", "Gm2016", "Gm5662", "Gm8300", "Obox3", "Pramef6", "Spz1", "Tdpoz1", "Cd74", "Lysmd1", "Parp10", 
              "Plk3", "Tacr3", "Tdpoz3", "Tdpoz4", "Iqcf1")
Toti.gene <- c("BC080695", "Chrm3", "Duoxa2", "Fgf7", "Gm11487", "Gm12794", "Gm13083", "Olfr376", "Tmem132c", "Zfp352", 
               "Zfp760", "Zscan4c", "Zscan4d", "Zscan4f")
gene.fpkm.select <- gene.fpkm.log[c(Maternal.gene, ZGA.gene, Toti.gene), ]
gene.fpkm.select.s <- data.frame(t(apply(gene.fpkm.select, 1, scale)), check.names = F)
colnames(gene.fpkm.select.s) <- colnames(gene.fpkm.select)
color <- brewer.pal(11, "RdBu")
color <- color[c(1:5, 7:10)]
gene.fpkm.select.s$gene <- rownames(gene.fpkm.select.s)
gene.fpkm.select.m <- melt(gene.fpkm.select.s)
gene.fpkm.select.m$gene <- factor(gene.fpkm.select.m$gene, levels = c(Maternal.gene, ZGA.gene, Toti.gene))
gene.fpkm.select.m$variable <- factor(gene.fpkm.select.m$variable, levels = rev(colnames(gene.fpkm.select)))
ggplot() + 
  geom_point(gene.fpkm.select.m, mapping = aes(x = gene, y = variable, color = value), size = 3) + 
  scale_color_gradientn(colors = colorRampPalette(colors = rev(color))(100)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        axis.text = element_text(size = 15))

##################Figure 2 Related to Fig.2a###############################
#Read datasets
TBLC.nor <- read.table("TBLC-passage-CPM-normalized.csv", sep = ",", header = T)
TBLC.removal <- read.table("TBLC-converted-PSC-CPM-normalized.csv", sep = ",", header = T)
colnames(TBLC.removal)[1] <- "gene"
TBLC.nor <- merge(TBLC.nor, TBLC.removal, by = "gene", all = F)
TBLC.nor <- TBLC.nor[, c(1:10,13,14)]
mouse.nor <- read.table("Mouse-embryo-CPM-normalized.csv", sep = ",", header = T)
ciTotiSC.nor <- read.table("ciTotiSC-CPM-normalized.csv", sep = ",", header = T)
#merge all datasets
TBLC.ciTotiSC.embryo <- Reduce(function(x, y){merge(x, y, by = "gene", all = F)}, 
                                    list(TBLC.nor, ciTotiSC.nor, mouse.nor))
rownames(TBLC.ciTotiSC.embryo) <- TBLC.ciTotiSC.embryo$gene
TBLC.ciTotiSC.embryo$gene <- NULL
TBLC.ciTotiSC.embryo <- TBLC.ciTotiSC.embryo[rowSums(TBLC.ciTotiSC.embryo)>0, ]
#using ComBat to remove batch effects
batch <- c(rep("1", 11), rep("2", 10), rep("3", 10))
TBLC.ciTotiSC.embryo.batch <- data.frame(ComBat(as.matrix(TBLC.ciTotiSC.embryo), batch = batch), check.names = F)
TBLC.ciTotiSC.embryo.batch.s <- TBLC.ciTotiSC.embryo.batch[, c(1,5:11,12,14,16,18,20,22:31)]
TBLC.ciTotiSC.embryo.pca <- PCA(t(TBLC.ciTotiSC.embryo.batch.s), graph = F)
fviz_pca_ind(TBLC.ciTotiSC.embryo.pca)

##################Figure 3 Related to Fig.2b###############################
#correlation calculation
TBLC.ciTotiSC.embryo.sub <- TBLC.TLSC.ciTotiSC.embryo[, c(1:11,12,14,16,18,20,22:31)]
TBLC.ciTotiSC.embryo.fc <- data.frame(TBLC.ciTotiSC.embryo.sub[, 2:11]/(TBLC.ciTotiSC.embryo.sub$P0.x+1), 
                                      TBLC.ciTotiSC.embryo.sub[, 13:16]/(TBLC.ciTotiSC.embryo.sub$OG2_mESC_rep1+1), 
                                      TBLC.ciTotiSC.embryo.sub[, 17:25]/(TBLC.ciTotiSC.embryo.sub$lateblast+1))
TBLC.ciTotiSC.embryo.cor <- cor(TBLC.ciTotiSC.embryo.fc)
TBLC.ciTotiSC.embryo.cor <- TBLC.ciTotiSC.embryo.cor[1:14, 15:23]
TBLC.ciTotiSC.embryo.cor.m <- melt(TBLC.ciTotiSC.embryo.cor)
TBLC.ciTotiSC.embryo.cor.m$Var1 <- factor(TBLC.ciTotiSC.embryo.cor.m$Var1, levels = rownames(TBLC.ciTotiSC.embryo.cor))
TBLC.ciTotiSC.embryo.cor.m$Var2 <- factor(TBLC.ciTotiSC.embryo.cor.m$Var2, levels = colnames(TBLC.ciTotiSC.embryo.cor))
color <- brewer.pal(11, "RdBu")
ggplot() + 
  geom_point(TBLC.ciTotiSC.embryo.cor.m, mapping = aes(x = Var1, y = Var2, color = value, size = value)) + 
  scale_color_gradientn(colors = colorRampPalette(colors = rev(color))(100)) + 
  scale_size(range = c(0, 6)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        axis.text = element_text(size = 15))





















