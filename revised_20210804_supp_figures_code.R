### packages
library(ggplot2)
library(gplots)
library(ggpubr)
library(patchwork)
library(dplyr)
library(PharmacoGx)
library(Biobase)
library(GEOquery)
library(annotate)
library(org.Hs.eg.db)
library(EnhancedVolcano)


#-------------------------------------------------------------
### S1a
plot <- read.csv("S1a pearson coefficient.csv", header = TRUE, row.names = 1)
plot <- as.matrix(plot)
tiff("1e pearson coefficient.tiff", units = "in", width = 4, height = 4, res = 300)
heatmap.2(plot,
          scale = "none",
          trace = "none",
          density.info = "none",
          dendrogram = "none",
          col = colorpanel(100, low = "blue", mid = "white", high = "red"),
          #breaks = seq(-4,4,length.out = 101),
          margins = c(8,12),
          keysize = 0.5,
          cexRow = 0.8,
          cexCol = 0.8,
          srtCol = 45,
          Colv = FALSE,
          Rowv = FALSE)
dev.off()

#-------------------------------------------------------------
### S1b
plot <- read.csv("S1b_pathways_heatmap.csv", header = TRUE, row.names = 1)
plot <- as.matrix(plot)
tiff("S1b_pathways.heatmap.tiff", units = "in", width = 7, height = 9, res = 300)
heatmap.2(plot,
          scale = "none",
          trace = "none",
          density.info = "none",
          dendrogram = "none",
          col = colorpanel(100, low = "blue", mid = "white", high = "red"),
          margins = c(12,18.5),
          keysize = 1,
          cexRow = 1.2,
          cexCol = 1,
          srtCol = 45,
          Colv = FALSE,
          Rowv = FALSE,
          key.xlab = "z-score",
          na.color = "grey",
          lhei = c(1,7))
dev.off()

#-------------------------------------------------------------
### S1c
tiff("pal.tiff", units = "in", width = 7, height = 6, res = 300)
plot <- read.csv('plot_diff.csv', header = TRUE)
p <- ggdotchart(plot, x = "pathway", y = "logp",
                color = "group",
                #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                sorting = "ascending",
                add = "segments",
                add.params = list(color = "group"),
                rotate = TRUE,
                group = "group",
                dot.size = 5,
                label = round(plot$average.difference),
                font.label = list(color = "white", size = 8, 
                                  vjust = 0.5, style = "bold"),
                ggtheme = theme_pubr(base_size = 2),
                ylab = "-log(pvalue)",
                xlab = "Pathway") +
  font("y.text", size = 9, face = "bold")
ggpar(p,
      legend.title = "",
      font.legend = c(7, "bold", "black"),
      font.x = c(12, "bold"),
      font.y = c(12,"bold"))
dev.off()




#-------------------------------------------------------------
### S1d
# GO Terms Analysis
plot <- read.csv("plot_GO.csv", header = TRUE)
plot <- subset(plot, gene == "up")
plot$Term <- factor(plot$Term, levels = plot$Term, ordered = TRUE)
up <- ggplot(plot, aes(x = Term, y = log)) +
  geom_bar(stat = "identity",
           fill = "coral2") +
  coord_flip() +
  ggtitle("Age Upregulated genes") +
  xlab("GO Biological Process") +
  ylab("-log(FDR)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", face = "bold", size = 16),
        plot.title = element_text(face = "bold", size = 16))

plot <- read.csv("plot_GO.csv", header = TRUE)
plot <- subset(plot, gene == "down")
plot$Term <- factor(plot$Term, levels = plot$Term, ordered = TRUE)
down <- ggplot(plot, aes(x = Term, y = log)) +
  geom_bar(stat = "identity",
           fill = "cornflowerblue") +
  coord_flip() +
  ggtitle("Age Downregulated genes") +
  xlab("GO Biological Process") +
  ylab("-log(FDR)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", face = "bold", size = 16),
        plot.title = element_text(face = "bold", size = 16))

tiff("GO.tiff", height = 6, width = 8, units = "in", res = 300)
up/down
dev.off()

#-------------------------------------------------------------
### S1e
# gh2ax cryo
tiff("gh2ax cryo rep.tiff", height = 3, width = 3, units = "in", res = 300)
plot <- read.csv("cryo_plot_rep.csv", header = TRUE)
plot$time <- factor(plot$time, levels = c("short", "long"), ordered = TRUE)
plot$ID <- factor(plot$ID, levels = c(6,7,9,32,36,42), ordered = TRUE)
ggplot(plot, aes(x = time, y = count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),
             aes(color = ID)) +
  ylim(c(-0.5, 25)) +
  stat_compare_means(method = "t.test",
                     label.y = 20,
                     label.x.npc = "left",
                     label = "..p.format..",
                     size = 5) +
  ylab("γh2ax + cells per HPF") +
  xlab("Time to Cryopreservation (hr)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 10, family = "Arial"),
        axis.title = element_text(size = 13, family = "Arial", face = "bold"))
dev.off()

#-------------------------------------------------------------
### S1f
# gh2ax old vs young 3 dots
tiff("gh2ax_3dots.tiff", height = 3, width = 3, units = "in", res = 300)
plot <- read.csv("S1f_3dots.csv", header = TRUE)
plot$donor <- factor(plot$donor, levels = c("short1", "short2", "short3", "long1", "long2", "long3"), ordered = TRUE)
ggplot(plot, aes(x = TL, y = count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),
             aes(color = donor)) +
  scale_color_brewer(palette = "Spectral") +
  ylim(c(0, 4)) +
  stat_compare_means(method = "t.test",
                     label.y = 3.5,
                     label.x.npc = "left",
                     label = "..p.format..",
                     size = 5) +
  ylab("Avg # γh2ax + cells per HPF") +
  xlab("Telomere Length") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 10, family = "Arial"),
        axis.title = element_text(size = 13, family = "Arial", face = "bold"))
dev.off()

#-------------------------------------------------------------
### S1g
# time vs. age
plot <- read.csv("time.vs.age.csv", header = TRUE, row.names = 1)
formula <- y ~ x
cor.test(plot$time, plot$age, method = "pearson")
tiff("time.vs.age.tiff", units="in", width=3.5, height=3.5, res=300)
ggplot(plot, aes(x=age, y=time)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  xlab("Age (yr)") +
  ylab("Time to Cryopreservation (hr)") +
  ggtitle("Age vs. Time") +
  ylim(c(0,60)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(size = 14, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 21, y = 56,
           label = "r = 0.012\np = 0.94", hjust = 0,
           size = 6.5)
dev.off()
# time vs. CDKN2A
plot <- read.csv("time.ncounts.plot.csv", header = TRUE, row.names = 1)
formula <- y ~ x
cor.test(log10(plot$CDKN2A), plot$time, method = "pearson")$estimate
tiff("time_cdkn2a.tiff", units="in", width=3.5, height=3.5, res=300)
ggplot(plot, aes(x=time, y=CDKN2A)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  scale_y_log10(name = "Normalized Count", breaks=c(25,100,400), limits = c(6, 257)) +
  scale_x_continuous(name = "Time to Cryopreservation (hr)") +
  ggtitle("CDKN2A") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 14, family = "Arial"),
        axis.title = element_text(size = 14, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 25, y = 190,
           label = "r = -0.009\np = 0.95", hjust = 0,
           size = 6.5)
dev.off()
# time vs. CDKN1A
plot <- read.csv("time.ncounts.plot.csv", header = TRUE, row.names = 1)
formula <- y ~ x
cor.test(log10(plot$CDKN1A), plot$time, method = "pearson")
tiff("time_cdkn1a.tiff", units="in", width=3.5, height=3.5, res=300)
ggplot(plot, aes(x=time, y=CDKN1A)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  scale_y_log10(name = "Normalized Count", breaks = c(3000,5000,10000), limits = c(1000, 21000)) +
  scale_x_continuous(name = "Time to Cryopreservation (hr)") +
  ggtitle("CDKN1A") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(size = 14, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 25, y = 17000,
           label = "r = -0.02\np = 0.90", hjust = 0,
           size = 6.5)
dev.off()
# time vs apoptosis
plot <- read.csv("score.apoptosis.csv", header = TRUE)
tiff("apoptosis.tiff", units="in", width=3.5, height=3.5, res=300)
ggplot(plot, aes(x=time, y=TotalScore)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  xlab("Time to Cryopreservation (hr)") +
  ggtitle("Time vs. Apoptosis Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 16, vjust = 2),
        axis.text = element_text(size = 16, family = "Arial"),
        axis.title = element_text(size = 14, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 32, y = 0.035,
           label = "r = -0.07\np = 0.64", hjust = 0,
           size = 7)
dev.off()

cor.test(plot$time, plot$TotalScore, method = "pearson")

# time vs necrosis
plot <- read.csv("score.necrosis.csv", header = TRUE)
tiff("necrosis.tiff", units="in", width=3.5, height=3.5, res=300)
ggplot(plot, aes(x=time, y=TotalScore)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  xlab("Time to Cryopreservation (hr)") +
  ggtitle("Time vs. Necrosis Score") +
  ylim(c(-0.04,0.04)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 16, vjust = 2),
        axis.text = element_text(size = 16, family = "Arial"),
        axis.title = element_text(size = 14, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 32, y = 0.033,
           label = "r = -0.26\np = 0.09", hjust = 0,
           size = 7)
dev.off()

cor.test(plot$time, plot$TotalScore, method = "pearson")


#-------------------------------------------------------------
## Figure S2a
plot <- read.csv("LAC Cell Type Proportions.csv", header = TRUE, row.names = 1)

cor.test(plot$Macrophages, plot$Age, method = "pearson")
lacmac <- ggplot(plot, aes(x=Age, y=Macrophages)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  ggtitle("LAC Macrophages") +
  ylab("Proportion") +
  ylim(0,0.45) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 16, y = 0.4,
           label = "r = 0.02\np = 0.85", hjust = 0,
           size = 7)
cor.test(plot$Bcells, plot$Age, method = "pearson")
lacb <- ggplot(plot, aes(x=Age, y=Bcells)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  ggtitle("LAC B Cells") +
  ylab("Proportion") +
  ylim(c(0,0.055)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 16, y = 0.049,
           label = "r = 0.17\np = 0.13", hjust = 0,
           size = 7)
cor.test(plot$Endothelial.Cells, plot$Age, method = "pearson")
lacendo <- ggplot(plot, aes(x=Age, y=Endothelial.Cells)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  ggtitle("LAC Endothelial Cells") +
  ylab("Proportion") +
  ylim(c(0,0.28)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 16, y = 0.255,
           label = "r = -0.08\np = 0.45", hjust = 0,
           size = 7)
cor.test(plot$Monocytes, plot$Age, method = "pearson")
lacmono <- ggplot(plot, aes(x=Age, y=Monocytes)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  ggtitle("LAC Monocytes") +
  ylab("Proportion") +
  ylim(c(0,0.06)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           label = "r = 0.05\np = 0.65",
           hjust = 0,
           x = 16,
           y = 0.05,
           size = 7)
cor.test(plot$NK.cells, plot$Age, method = "pearson")
lacnk <- ggplot(plot, aes(x=Age, y=NK.cells)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  ggtitle("LAC NK Cells") +
  ylab("Proportion") +
  ylim(c(0,0.15)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           label = "r = -0.04\np = 0.7207",
           hjust = 0,
           x = 16,
           y = 0.13,
           size = 7)

#gtex
plot <- read.csv("gtex Cell Type Proportions.csv", header = TRUE, row.names = 1)
plot$Age <- factor(plot$Age, levels = c("20","30","40","50","60","70"), labels = c("20","30","40","50","60","70"))
gtexmac <- ggplot(plot, aes(x = Age, y = Macrophages)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_classic() +
  ggtitle("GTEx Macrophages") +
  ylab("Proportion") +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  #geom_jitter(width = 0.1, size = 1) +
  annotate(geom = "text",
           label = "p = 0.009",
           x = "20",
           y = 0.35,
           hjust = 0,
           size = 7)
gtexb <- ggplot(plot, aes(x = Age, y = Bcells)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_classic() +
  ggtitle("GTEx B Cells") +
  ylab("Proportion") +
  ylim(c(0,0.07)) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  #geom_jitter(width = 0.1, size = 1) +
  annotate(geom = "text",
           label = "p = 0.48",
           x = "20",
           y = 0.065,
           hjust = 0,
           size = 7)
gtexmono <- ggplot(plot, aes(x = Age, y = Monocytes)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_classic() +
  ggtitle("GTEx Monocytes") +
  ylab("Proportion") +
  ylim(c(0,0.06)) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  #geom_jitter(width = 0.1, size = 1) +
  annotate(geom = "text",
           label = "p = 0.66",
           x = "20",
           y = 0.058,
           hjust = 0,
           size = 7)
gtexendo <- ggplot(plot, aes(x = Age, y = Endothelial.cells)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_classic() +
  ggtitle("GTEx Endothelial Cells") +
  ylab("Proportion") +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  #geom_jitter(width = 0.1, size = 1) +
  annotate(geom = "text",
           label = "p = 0.13",
           x = "20",
           y = 0.7,
           hjust = 0,
           size = 7)
gtexnk <- ggplot(plot, aes(x = Age, y = NK.cells)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_classic() +
  ggtitle("GTEx NK Cells") +
  ylab("Proportion") +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  #geom_jitter(width = 0.1, size = 1) +
  annotate(geom = "text",
           label = "p = 0.28",
           x = "20",
           y = 0.1,
           hjust = 0,
           size = 7)
tiff("S2a_b", units = "in", width = 3.5, height = 7, res = 300)
lacb/gtexb
dev.off()
tiff("S2a_endo", units = "in", width = 3.5, height = 7, res = 300)
lacendo/gtexendo
dev.off()
tiff("S2a_mono", units = "in", width = 3.5, height = 7, res = 300)
lacmono/gtexmono
dev.off()
tiff("S2a_nk", units = "in", width = 3.5, height = 7, res = 300)
lacnk/gtexnk
dev.off()
tiff("S2a_mac", units = "in", width = 3.5, height = 7, res = 300)
lacmac/gtexmac
dev.off()

#-------------------------------------------------------------
### Figure S2b
# lac
plot <- read.csv("trav.deconv.lac.weighted.prop.csv", header = TRUE)
cor.test(plot$Age, plot$Alveolar.Epithelial.Type.2, method = "pearson")
lact2 <- ggplot(plot, aes(x = Age, y = Alveolar.Epithelial.Type.2)) +
  geom_point(position = position_jitter(w = 0.1, h = 0),
             color = "black", 
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  theme_classic() +
  ggtitle("LAC Alveolar Epithelial T2") +
  xlab("Age") +
  ylab("Proportion") +
  ylim(c(0,0.55)) +
  geom_jitter(width = 0.1) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           label = "r = -0.17\np = 0.13",
           hjust = 0,
           x = 16,
           y = 0.47,
           size = 7)
cor.test(plot$Age, plot$Alveolar.Fibroblast, method = "pearson")
lacfib <- ggplot(plot, aes(x = Age, y = Alveolar.Fibroblast)) +
  geom_point(position = position_jitter(w = 0.1, h = 0),
             color = "black", 
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  theme_classic() +
  ggtitle("LAC Alveolar Fibroblast") +
  xlab("Age") +
  ylab("Proportion") +
  ylim(c(-0.01,0.32)) +
  geom_jitter(width = 0.1) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           label = "r = 0.12\np = 0.28",
           hjust = 0,
           x = 16,
           y = 0.28,
           size = 7)

# gtex
plot <- read.csv("trav.deconv.gtex.weighted.prop.csv", header = TRUE)
plot$Age <- factor(plot$Age)
gtext2 <- ggplot(plot, aes(x = Age, y = Alveolar.Epithelial.Type.2)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("GTEx Alveolar Epithelial T2") +
  ylab("Proportion") +
  ylim(c(-0.005,0.24)) +
  geom_jitter(width = 0.1) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           label = "p = 0.0061",
           x = "20",
           y = 0.22,
           hjust = 0,
           size = 7)
gtexfib <- ggplot(plot, aes(x = Age, y = Alveolar.Fibroblast)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("GTEx Alveolar Fibroblast") +
  ylab("Proportion") +
  ylim(c(0,0.6)) +
  geom_jitter(width = 0.1) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           label = "p = 0.0067",
           x = "20",
           y = 0.58,
           hjust = 0,
           size = 7)

tiff("S2b_t2.tiff", width = 9, height = 4.5, units = "in", res = 300)
lact2+gtext2
dev.off()
tiff("S2b_fib.tiff", width = 9, height = 4.5, units = "in", res = 300)
lacfib+gtexfib
dev.off()

#-------------------------------------------------------------
### Figure S2c
plot <- read.csv("skin.deconv.weighted.prop.csv", header = TRUE)
plot$Age <- factor(plot$Age)
skinepi <- ggplot(plot, aes(x = Age, y = Epithelial.cells)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("GTEx Sun-exposed Skin\nEpithelial Cells") +
  ylab("Proportion") +
  ylim(c(-0.01,0.5)) +
  geom_jitter(width = 0.1) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 10, family = "Arial", color = "black"),
        axis.title = element_text(size = 14, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           label = "p = 0.024",
           x = "20-29",
           y = 0.45,
           hjust = 0,
           size = 7)
skinfib <- ggplot(plot, aes(x = Age, y = Fibroblasts)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("GTEx Sun-exposed Skin\nFibroblasts") +
  ylab("Proportion") +
  ylim(c(-0.01,0.4)) +
  geom_jitter(width = 0.1) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 10, family = "Arial", color = "black"),
        axis.title = element_text(size = 14, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           label = "p = 0.0072",
           x = "20-29",
           y = 0.37,
           hjust = 0,
           size = 7)
tiff("S2c_combined.tiff", width = 3.5, height = 7, units = "in", res = 300)
skinepi/skinfib
dev.off()

#-------------------------------------------------------------
### Figure S3b
#sun exposed skin
plot <- read.csv("S3b_gtexsunskin_enriched.collagen.csv", header = TRUE)
plot$Age <- factor(plot$Age, levels = c("young", "old"), ordered = TRUE)
sskin <- ggplot(plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("GTEx Skin: Sun Exposed") +
  ylim(0.3,0.41) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold"),
        axis.text = element_text(size = 16,family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold"),
        legend.title = element_text(size = 15, family = "Arial"),
        legend.text = element_text(size = 15, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label.y = 0.405,
                     label = "..p.format..",
                     size = 6)
#not sun exposed skin
plot <- read.csv("S3b_gtexNOsunskin_enriched.collagen.csv", header = TRUE)
plot$Age <- factor(plot$Age, levels = c("young", "old"), ordered = TRUE)
nsskin <- ggplot(plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("GTEx Skin: Not Sun Exposed") +
  ylim(0.29,0.415) +
  theme(plot.title = element_text(hjust = 0.5, size = 15.8, family = "Arial", face = "bold"),
        axis.text = element_text(size = 16,family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold"),
        legend.title = element_text(size = 15, family = "Arial"),
        legend.text = element_text(size = 15, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label.y = 0.405,
                     label = "..p.format..",
                     size = 6)
#heart
plot <- read.csv("S3b_gtexheart_enriched.collagen.csv", header = TRUE)
plot$Age <- factor(plot$Age, levels = c("young", "old"), ordered = TRUE)
heart <- ggplot(plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("GTEx Heart") +
  ylim(0.24,0.415) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold"),
        axis.text = element_text(size = 16,family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold"),
        legend.title = element_text(size = 15, family = "Arial"),
        legend.text = element_text(size = 15, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label.y = 0.405,
                     label = "..p.format..",
                     size = 6)
#liver
plot <- read.csv("S3b_gtexliver_enriched.collagen.csv", header = TRUE)
plot$Age <- factor(plot$Age, levels = c("young", "old"), ordered = TRUE)
liver <- ggplot(plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("GTEx Liver") +
  ylim(0.23,0.415) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold"),
        axis.text = element_text(size = 16,family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold"),
        legend.title = element_text(size = 15, family = "Arial"),
        legend.text = element_text(size = 15, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label.y = 0.405,
                     label = "..p.format..",
                     size = 6)
#kidney
plot <- read.csv("S3b_gtexkidney_enriched.collagen.csv", header = TRUE)
plot$Age <- factor(plot$Age, levels = c("young", "old"), ordered = TRUE)
kidney <- ggplot(plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("GTEx Kidney") +
  ylim(0.24,0.415) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold"),
        axis.text = element_text(size = 16,family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold"),
        legend.title = element_text(size = 15, family = "Arial"),
        legend.text = element_text(size = 15, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label.y = 0.405,
                     label = "..p.format..",
                     size = 6)
tiff("S3b_skins.tiff", width = 7, height = 3.5, units = "in", res = 300)
sskin + nsskin
dev.off()

tiff("S3b_others.tiff", width = 10.5, height = 3.5, units = "in", res = 300)
kidney + liver + heart
dev.off()

#-------------------------------------------------------------
### Figure S3c
genes = read.table('cmap_genes.txt',row.names = NULL,fill = T,header=T)
Cmap.pertub <- downloadPertSig("CMAP_2016")
rownames(Cmap.pertub) <- gsub('_at', '', rownames(Cmap.pertub))
# Aging Signature CMAP
g1 = genes$Aging_up
g2 = genes$Aging_dn
g1 = g1[g1!='']
g2 = g2[g2!='']
ens1 <- mapIds(org.Hs.eg.db, keys = g1, keytype="SYMBOL",column = "ENSEMBL" )
ens2 <- mapIds(org.Hs.eg.db, keys = g2, keytype="SYMBOL",column = "ENSEMBL" )
ens1 = ens1[!is.na(ens1)]
ens2 = ens2[!is.na(ens2)]
ens = data.frame(row.names=as.character(c(ens1,ens2)),
                 feature=as.character(c(ens1,ens2)),
                 direction=c(rep(1,length(ens1)),rep(-1,length(ens2))))
ens = ens[ens$feature %in% rownames(Cmap.pertub),]
res_aging <- apply(Cmap.pertub[,,c("tstat", "fdr")],
                   2, function(x, genes){
                     return(connectivityScore(x=x,
                                              y=ens[,2,drop=F],
                                              method="fgsea", nperm=100))
                   }, genes=ens)
res_aging = as.data.frame(t(res_aging))
write.csv(res_aging, "res_aging.csv")

tiff("aging_cmap.tiff", width = 6, height = 8, units = "in", res = 300)
ageplot <- EnhancedVolcano(res_aging,
                           lab = rownames(res_aging),
                           x = 'V1',
                           y = 'V2',
                           title = 'Aging Signature',
                           pCutoff = 0.05,
                           FCcutoff = 0.5,
                           pointSize = 2.0,
                           labSize = 4.0,ylim = c(0,2))
ageplot
dev.off()

# Collagen CMAP
g1 = genes$Collagen
g1 = g1[g1!='']
ens1 <- mapIds(org.Hs.eg.db, keys = g1, keytype="SYMBOL",column = "ENSEMBL" )
ens1 = ens1[!is.na(ens1)]
ens = data.frame(row.names=as.character(c(ens1)),
                 feature=as.character(c(ens1)),
                 direction=c(rep(1,length(ens1))))
ens = ens[ens$feature %in% rownames(Cmap.pertub),]
res_collagen <- apply(Cmap.pertub[,,c("tstat", "fdr")],
                      2, function(x, genes){
                        return(connectivityScore(x=x,
                                                 y=ens[,2,drop=F],
                                                 method="fgsea", nperm=100))
                      }, genes=ens)
res_collagen = as.data.frame(t(res_collagen))
write.csv(res_collagen, "res_collagen.csv")
tiff("collagen_cmap.tiff", width = 6, height = 8, units = "in", res = 300)
colplot <- EnhancedVolcano(res_collagen,
                           lab = rownames(res_collagen),
                           x = 'es',
                           y = 'p',
                           title = 'Collagen Signature',
                           pCutoff = 0.05,
                           FCcutoff = 0.5,
                           pointSize = 2.0,
                           labSize = 4.0,ylim = c(0,2))
colplot
dev.off()


########SUPP FIGURES END#####################################










