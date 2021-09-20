library(singscore)
library(GSEABase)
library(plotly)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(dplyr)
library(ggpubr)

#For Singscore: start with RAW COUNTS table. You will need sample information and will have to
#transform the data so that the row names are in HGNC symbols, and no duplicate gene names.
#Will transform into tpm. Will need gene sets file containing each gene set in each column.


## GTEx Data Transformation for singscore
#Extracted quintiles samples from sample.info
#Extract corresponding counts
df <- read.csv("fgtex.sunNexposed.skin.counts.csv", header = TRUE, row.names = 1)
names <- read.csv("gtexsunNexposed.skin.pent.sample.info.csv", header = TRUE)
names <- names$ID
names <- as.character(names)
df <- as.data.frame(df)
sing <- df[colnames(df) %in% names == TRUE]
write.csv(as.data.frame(sing), "gtexsunNexposed.skin.pent.counts.csv")

#Converted to HGNC symbols in excel
#Added Gene lengths
#Remove duplicate genes
sing <- read.csv("gtexsunNexposed.skin.pent.counts.csv", header = TRUE)
sing <- sing[!duplicated(sing$X),]
row.names(sing) <- sing$X
sing <- sing[,-1]


#Generate tpm
#--------------------------------------------------------------------------
#counts <- read.csv("age.pent.counts.nodup.csv", row.names = 1)
onlycounts <- sing[,-1]
tpm <- calculateTPM(onlycounts, lengths = as.numeric(sing$lengths))
write.csv(as.data.frame(tpm), file = "age.pent.tpm.csv")
#--------------------------------------------------------------------------

#Gene Expression Data
#--------------------------------------------------------------------------
#exptpm <- read.csv("age.pent.tpm.csv", header = TRUE, row.names = 1)

##Sample Scoring
#Rank genes: genes are ranked by increasing mRNA abundance
rankData <- rankGenes(tpm)
#--------------------------------------------------------------------------

#Load Gene Set Pair
#--------------------------------------------------------------------------
GeneSets <- read.csv("genesets.csv", header = TRUE, row.names = NULL)

NegSet <- as.character(GeneSets$aging_dn)
NegSet <- NegSet[NegSet != ""]

PosSet <- as.character(GeneSets$aging_up)
PosSet <- PosSet[PosSet != ""]
#--------------------------------------------------------------------------


#Score by signature
#--------------------------------------------------------------------------
scoredf <- simpleScore(rankData, upSet = PosSet, knownDirection = TRUE)

#Remove missing genes if any
NegRemove <- c("")
NegSet <- NegSet[!NegSet %in% NegRemove]
PosRemove <- c("")
PosSet <- PosSet[!PosSet %in% PosRemove]

#Add column 'Age'
age <- read.csv("sample.info.pent.csv", header=TRUE, row.names = 1)
scoredf <- mutate(scoredf, age$Age)
row.names(scoredf) <- rownames(age)

#Export scores as pdf
write.csv(as.data.frame(scoredf), file = "score.aging.csv")
#--------------------------------------------------------------------------
