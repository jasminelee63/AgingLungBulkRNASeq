### packages
library(ggplot2)
library(gplots)
library(ggpubr)
library(patchwork)
library(dplyr)
library(DOSE)
library(clusterProfiler)
library(enrichplot)

#-------------------------------------------------------------
### Figure 1c
# smoking
smoking <- data.frame(
  group = c("Yes", "No"),
  value = c(47.7, 52.3),
  label = c("Yes\n47.7%", "No\n52.3%")
)
smoking_fig <- ggplot(smoking, aes(x = "", y = value, fill = group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  ggtitle("Ever Smoked") +
  theme(axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5)) +
  geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                label = label), size=10)
# gender
gender <- data.frame(
  group = c("Male", "Female"),
  value = c(57, 43),
  label = c("Male\n57%", "Female\n43%")
)
gender_fig <- ggplot(gender, aes(x = "", y = value, fill = group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal() +
  ggtitle("Sex") +
  theme(axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5)) +
  geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                label = label), size=10)
# ethnicity
ethnicity <- data.frame(
  group = c("White: 61.6%", "Hispanic/Latino : 25.6%", "Asian : 7%", "Black/African American : 3.5%", "Other : 2.3%"),
  value = c(61.6, 25.6, 7, 3.5, 2.3),
  label = c("W", "H", "A", "B", "O")
)
ethnicity <- ethnicity %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(ethnicity$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
ethnicity_fig <- ggplot(ethnicity, aes(x = "", y = prop, fill = group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  ggtitle("Ethnicity") +
  theme(axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 25),
        legend.title = element_blank(),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5))
tiff("combined_demo.tiff", width = 20, height = 6, units = "in", res = 300, bg = "transparent")
gender_fig+smoking_fig+ethnicity_fig
dev.off()

#-------------------------------------------------------------
### Figure 1d
significant <- read.csv("1d significant.ncounts.raw.csv", header = TRUE, row.names = 1)
significant <- as.matrix(significant)
tiff("1d_ncounts_heatmap.tiff", units="in", width=8, height=6, res=300, bg = "transparent")
heatmap.2(significant,
          scale = "row",
          trace = "none",
          density.info = "none",
          dendrogram = "none",
          col = colorpanel(100, low = "blue", mid = "white", high = "red"),
          breaks = seq(-4,4,length.out = 101),
          margins = c(8,8),
          keysize = 1.1,
          cexRow = 1.2,
          Colv = FALSE,
          Rowv = FALSE,
          labCol = FALSE)
          #lmat=rbind(c(4, 2), c(1, 3)),
          #lhei=c(1, 3),
          #lwid=c(4, 1),
          #key.par=list(mar=c(4,4,4,8)))
dev.off()

#-------------------------------------------------------------
### Figure 2a
plot <- read.csv("2a_cdkn2a.csv", header = TRUE, row.names = 1)
formula <- y ~ x
tiff("2a_cdkn2a.tiff", units="in", width=3.5, height=3.5, res=300)
ggplot(plot, aes(x=Age, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "coral2",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  scale_y_log10(name = "log(Normalized Count)", breaks=c(25,100,400), limits = c(6, 400)) +
  scale_x_continuous(name = "Age (yr)", breaks=c(20,30,40,50,60,70)) +
  ggtitle("CDKN2A") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 16, y = 280,
           label = "r = 0.54\np = 6.61e-08", hjust = 0,
           size = 7)
dev.off()
cor.test(log10(plot$count), plot$Age, method = "pearson")

#-------------------------------------------------------------
### Figure 2b
## GSEA
# ranked gene list, delete duplicated gene names
d <- read.csv("2b_corr.csv", header = TRUE)
d <- d[!duplicated(d$ID),]
# rank by pearson correlation
dcor <- d[,2]
names(dcor) <- as.character(d[,1])
dcor <- sort(dcor, decreasing = TRUE)
# gene sets
gs <- read.csv("2b_geneset.csv", header = TRUE)
collagen <- subset(gs, gs_name == "Collagen_Processing_Genes")
senescence <- subset(gs, gs_name == "Senescence_Consensus_Genes")
# run GSEA
em <- GSEA(dcor, TERM2GENE = gs)
p <- lapply(1:1, function(i) {
  anno <- em[i, c("NES", "pvalue", "p.adjust")]
  anno$NES <- round(anno$NES, 3)
  anno$pvalue <- formatC(anno$pvalue, format = "e", digits = 2)
  anno$p.adjust <- formatC(anno$p.adjust, format = "e", digits = 2)
  lab <- paste0(names(anno), "=", anno, collapse="\n")
  
  gseaplot(em, geneSetID = i, by = "runningScore", title = em$Description[i], color.line = "#C50000") +
    annotate("text", 21000, em[i, "enrichmentScore"] * .7, label = lab, hjust=0, vjust=0)
})
# plot
tiff("collagen_processing.tiff", res = 300, width = 3.5, height = 3.5, units = "in")
p[1]
dev.off()
tiff("senescence.tiff", res = 300, width = 3.5, height = 3.5, units = "in")
p[2]
dev.off()

## plot singscore
age.plot <- read.csv("2b_LAC_consensus_score.csv", header = TRUE)
age.plot$Age <- factor(age.plot$Age, levels = c("young", "old"), ordered = TRUE)
lac <- ggplot(age.plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("Lung Aging Cohort") +
  ylim(0.325,0.38) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold"),
        axis.text = element_text(size = 10, family = "Arial", color = "black"),
        axis.title = element_text(size = 12, family = "Arial", face = "bold"),
        legend.title = element_text(size = 10, family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label = "..p.format..",
                     label.y = 0.377,
                     size = 3.4)
gtexlung.plot <- read.csv("2b_gtexlung_consensus_score.csv", header = TRUE)
gtexlung.plot$Age <- factor(gtexlung.plot$Age, levels = c("young", "old"), ordered = TRUE)
gtex <- ggplot(gtexlung.plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("GTEx Lung") +
  ylim(0.32,0.405) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold"),
        axis.text = element_text(size = 10, family = "Arial", color = "black"),
        axis.title = element_text(size = 12, family = "Arial", face = "bold"),
        legend.title = element_text(size = 10, family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label = "..p.format..",
                     label.y = 0.40,
                     size = 3.4)
tiff("2b.tiff", units="in", width=5, height=2.5, res=300)
lac+gtex
dev.off()

#-------------------------------------------------------------
### Figure 2c
plot <- read.csv("2c_LAC_IPA_selected_pathways.csv", header = TRUE)
plot$name<- factor(plot$name, levels = plot$name)
tiff("2c_pathways.tiff", units = "in", width = 10, height = 12, res = 300)
ggplot(plot, aes(x = name, y = zscore, fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"YlGnBu"),
                       trans = "reverse") +
  theme_classic() +
  ggtitle("LAC Significant Pathways") +
  labs(x = "Pathway/Function/Upstream Regulator",
       y = "z-score") +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, size = 25, family = "Arial", face = "bold"),,
        axis.text = element_text(size = 21, color = "black"),
        axis.title = element_text(size = 21, face = "bold"),
        legend.title = element_text(size = 21, family = "Arial"),
        legend.text = element_text(size = 21, family = "Arial"))
dev.off()

#-------------------------------------------------------------
### Figure 2d
# Telomere qPCR
plot <- read.csv("2d_TL.csv", header = TRUE, row.names = 1)
formula <- y ~ x
tiff("2d_TL.tiff", units="in", width=3.5, height=3.5, res=300)
ggplot(plot, aes(x=age, y=TL)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "cornflowerblue",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  ylab("Average Length (bp)")+
  ylim(6500,9800) +
  scale_x_continuous(name = "Age (yr)", breaks=c(20,30,40,50,60,70)) +
  ggtitle("Telomere Length") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 18, family = "Arial"),
        axis.title = element_text(size = 18, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x =35, y = 9590,
           label = "r = -0.41\np = 1.034e-04", hjust = 0,
           size = 6)
dev.off()
cor.test(plot$TL, plot$age, method = "pearson")
# gh2ax IHC count
# order : N1603, N1610, N1731, N1861, N1863, N1869
library(viridisLite)
library(grDevices)
plot <- read.csv("2d_teloIHCcount.csv", header = TRUE)
plot$donor <- factor(plot$donor, levels = c("short1", "short2", "short3", "long1", "long2", "long3"), ordered = TRUE)
plot$TL <- factor(plot$TL, levels = c("long", "short"), ordered = TRUE)
hcl <- hcl.colors(n = 6, palette = "Dynamic")
hcl <- as.character(c(hcl[c(8:10)], hcl[c(1:3)]))
tiff("2d_teloIHC.tiff", units="in", width=3.5, height=3.5, res=500)
ggbarplot(plot, x = "TL", y = "count",
          ylim = c(0,8.5),
          add = c("mean_se"),
          alpha = 0.5,
          legend = "right",
          xlab = "Telomere Length",
          ylab = "γH2AX Positive Cells per HPF",
          order = c("long", "short")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold"),
        axis.text = element_text(size = 18,family = "Arial", color = "black"),
        axis.title = element_text(size = 14, family = "Arial", face = "bold"),
        legend.title = element_text(size = 15, family = "Arial"),
        legend.text = element_text(size = 15, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = donor), size = 2) +
  #scale_color_manual(values = hcl) +
  scale_color_brewer(palette = "Spectral") +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label = "..p.format..",
                     label.y = 8,
                     size = 7)
dev.off()

#-------------------------------------------------------------
### Figure 3a
#LAC epi
plot <- read.csv("3a_Reyfman.young.LAC.weighted.csv", header = TRUE, row.names = 1)
cor.test(plot$Epithelial.cells, plot$Age, method = "pearson")
lacepi <- ggplot(plot, aes(x=Age, y=Epithelial.cells)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  ggtitle("LAC Epithelial Cells") +
  ylim(0,0.8) +
  ylab("Proportion") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 13, family = "Arial", color = "black"),
        axis.title = element_text(size = 13, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 16, y = 0.75,
           label = "r = -0.14\np = 0.19", hjust = 0,
           size = 7)
#LAC fib
cor.test(plot$Fibroblasts, plot$Age, method = "pearson")
lacfib <- ggplot(plot, aes(x=Age, y=Fibroblasts)) +
  geom_point(position=position_jitter(w=0.1,h=0),
             color = "black",
             size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  ggtitle("LAC Fibroblasts") +
  ylim(0.2,0.85) +
  ylab("Proportion") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 13, family = "Arial", color = "black"),
        axis.title = element_text(size = 13, family = "Arial", face = "bold")) +
  annotate(geom = "text",
           x = 16, y = 0.79,
           label = "r = 0.20\np = 0.06", hjust = 0,
           size = 7)
#gtexlung epi
plot <- read.csv("3a_Reyfman.young.gtex.weighted.csv", header = TRUE, row.names = 1)
plot$Age <- as.factor(plot$Age)
gtexepi <- ggplot(plot, aes(x = Age, y = Epithelial.cells)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("GTEx Epithelial Cells") +
  ylim(c(0,0.77)) +
  ylab("Proportion") +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 13, family = "Arial", color = "black"),
        axis.title = element_text(size = 13, family = "Arial", face = "bold")) +
  #geom_jitter(width = 0.1, size = 1) +
  annotate(geom = "text",
           x = "20-29", y = 0.72,
           label = "p = 1.8e-04", hjust = 0,
           size = 7)
#gtexlung fib
gtexfib <- ggplot(plot, aes(x = Age, y = Fibroblasts)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("GTEx Fibroblasts") +
  ylab("Proportion") +
  ylim(c(0.25,1.1)) +
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 18, vjust = 2),
        axis.text = element_text(size = 13, family = "Arial", color = "black"),
        axis.title = element_text(size = 13, family = "Arial", face = "bold")) +
  #geom_jitter(width = 0.1, size = 1) +
  annotate(geom = "text",
           x = "20-29", y = 1.1,
           label = "p = 8.2e-07", hjust = 0,
           size = 7)

tiff("3a_combined.tiff", units="in", width=8, height=8, res=300)
(lacepi + gtexepi) / (lacfib + gtexfib)
dev.off()

#-------------------------------------------------------------
### Figure 3b
# order : N1854(18), N1724(19), N1744(20), N1216(22), N1865(68), N1329(72), N1609(76), N1728(76)
donor <- factor(c(rep("18M",6), rep("19M",6), rep("20F",6), rep("22F",6), rep("68M",6), rep("72F",6), rep("76F",6), rep("76M",6)), ordered = TRUE)
agegroup <- factor(c(rep("young",24), rep("old",24)))
count <- c(57,51,16,18,15,17,38,47,39,18,26,23,32,21,33,20,56,38,41,52,31,39,41,58,23,19,20,20,28,21,5,4,8,5,6,4,7,7,6,9,8,6,1,0,2,2,0,2)
ratio <- c(0.126666667,0.142458101,0.089385475,0.125,0.130434783,0.129770992,0.087962963,0.11407767,0.065878378,0.029363785,0.050290135,0.058823529,0.118081181,0.094170404,0.099697885,0.052493438,0.092868988,0.092457421,0.109625668,0.096118299,0.065539112,0.079107505,0.138513514,0.142857143,0.074918567,0.053824363,0.051813472,0.08097166,0.09929078,0.07,0.026041667,0.020618557,0.037914692,0.034482759,0.037974684,0.026143791,0.046052632,0.094594595,0.092307692,0.069767442,0.06557377,0.133333333,0.006666667,0,0.008658009,0.00913242,0,0.00896861)
percent <- ratio*100
permicron <- c(0.022678954,
             0.022889433,
             0.013544562,
             0.012132643,
             0.012593645,
             0.014625176,
             0.021697148,
             0.023943264,
             0.015447515,
             0.006741696,
             0.01639348,
             0.008331297,
             0.016162801,
             0.017024442,
             0.022291733,
             0.011180317,
             0.025990933,
             0.01946071,
             0.015255498,
             0.017519802,
             0.013160251,
             0.01829702,
             0.02144251,
             0.02658708,
             0.011325388,
             0.010768799,
             0.011628018,
             0.011729987,
             0.011989786,
             0.010438246,
             0.002727138,
             0.002759754,
             0.005471801,
             0.005080953,
             0.004877476,
             0.003310052,
             0.004933321,
             0.010068917,
             0.006663811,
             0.006181675,
             0.007313347,
             0.010605288,
             0.001093588,
             0,
             0.001044211,
             0.001251597,
             0,
             0.001426001)
permm <- permicron*1000
df <- data.frame("donor" = donor,
                 "agegroup" = agegroup,
                 "count" = count,
                 "permm" = permm)
df$agegroup <- factor(df$agegroup, levels = c("young", "old"), ordered = TRUE)
df$donor <- factor(df$donor, levels = c("18M", "19M", "20F", "22F", "68M", "72F", "76F", "76M"), ordered = TRUE)
dfyoung <- df[1:24, ]
dfold <- df[25:48, ]

tiff("3b_spc_perunitlength.tiff", units="in", width=5, height=3.5, res=300)
ggplot(df, aes(x = agegroup, y = permm)) +
  geom_boxplot() +
  theme_classic() +
  ylab("SPC + cells/\nmm alveolar surface") +
  xlab("Age") +
  #ylim(0,1.6) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold"),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 16, family = "Arial", face = "bold"),
        legend.title = element_text(size = 18, family = "Arial"),
        legend.text = element_text(size = 15, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = donor), size = 1.5) +
  scale_color_brewer(palette = "Paired") +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label.y = 30,
                     label = "..p.format..",
                     size = 7,
                     vjust = 0.8)
dev.off()

#-------------------------------------------------------------
### Figure 3c
#LAC
age.plot <- read.csv("3c_LAC_enriched.collagen_score.csv", header = TRUE)
age.plot$Age <- factor(age.plot$Age, levels = c("young", "old"), ordered = TRUE)
lac <- ggplot(age.plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("Lung Aging Cohort") +
  ylim(0.26,0.362) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold"),
        axis.text = element_text(size = 10, family = "Arial", color = "black"),
        axis.title = element_text(size = 12, family = "Arial", face = "bold"),
        legend.title = element_text(size = 10, family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label.y = 0.36,
                     label = "..p.format..",
                     size = 3.2)
#gtexlung
gtexlung.plot <- read.csv("3c_gtexlung_enriched.collagen_score.csv", header = TRUE)
gtexlung.plot$Age <- factor(gtexlung.plot$Age, levels = c("young", "old"), ordered = TRUE)
gtex <- ggplot(gtexlung.plot, aes(x = Age, y = TotalScore)) +
  geom_boxplot(aes(color = Age)) +
  theme_classic() +
  ggtitle("GTEx Lung") +
  ylim(0.25,0.40) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold"),
        axis.text = element_text(size = 10, family = "Arial", color = "black"),
        axis.title = element_text(size = 12, family = "Arial", face = "bold"),
        legend.title = element_text(size = 10, family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial")) +
  geom_jitter(width = 0.1, aes(color = Age), size = 1) +
  scale_color_manual(values = c("cornflowerblue", "coral2")) +
  stat_compare_means(method = "t.test",
                     label.x.npc = "left",
                     label.y = 0.395,
                     label = "..p.format..",
                     size = 3.2)
tiff("3c_combined.tiff", units="in", width=5, height=2.5, res=300)
lac+gtex
dev.off()

#####MAIN FIGURES END#####################################