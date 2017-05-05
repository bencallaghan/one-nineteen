
library("org.Hs.eg.db")
library(Rgraphviz)

mssng.hitlist$ids <- mapIds(org.Hs.eg.db, levels(mssng.hitlist$GeneSymbol)[mssng.hitlist$GeneSymbol], "ENTREZID", "ALIAS", multiVals = "first")

geneNames <- mssng.hitlist$ids
mssng.hitlist %>% filter(Freq.corr > 0) %>% dplyr::select(ids) -> myInterestingGenes
myInterestingGenes <- as.character(myInterestingGenes$ids)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)
geneList


# GEA with topGO

data(ALL) # Using Acute Lymphoblastic Leukemia GE data 
data(geneList)
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)
sum(topDiffGenes(geneList))

sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)

sampleGOdata

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher

resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
        elim = pValue.elim[sel.go],
        classic = pValue.classic[sel.go])
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all',)













# Creating tables for results section
library(xtable)
data(tli)
fm2 <- lm(tlimth ~ sex * ethnicty, data = tli)
print(xtable(anova(fm2)), type="html")


# De novo 
# 5 of 69 prioritized variants were confirmed de novo by Sanger sequencing
# SCN2A (LOF + MS)x
# NIPBL
# WDR45x
# ARID2x
toplof31 %>% select(SAMPLE_ID,Gene,X.CHROM,POS, REF, ALT, EFF.0..EFFECT) %>% filter(row.names(toplof31) %in% c('2','8','13'))-> denovo.LOF
names(denovo.LOF) <- c("Subject_ID","Gene", "Chrom", "Pos", "Ref", "Alt", "Effect")
denovo.LOF$Effect <- gsub("_", " ", denovo.LOF$Effect)
denovo.LOF$Effect <- gsub("&.+", "", denovo.LOF$Effect)
write.csv(x = denovo.LOF, "outputs/lof.dn.variants.csv", row.names = F)

denovo.LOF
print(xtable(denovo.LOF))

print(xtable(denovo.LOF), type = "html")

topms66 %>% select(SAMPLE_ID,GENE,CHROM,POS, REF, ALT, HGVS_P, CADD_PHRED) %>% filter(row.names(topms66) %in% c('24','42'))-> denovo.MS
names(denovo.MS) <- c("Subject_ID", "Gene", "Chrom", "Pos", "Ref", "Alt", "Effect", "Cadd")
denovo.MS$Effect <- gsub("_", " ", denovo.MS$Effect)
write.csv(x = denovo.MS, "outputs/ms.dn.variants.csv", row.names = F)
xtable(denovo.MS)
print(xtable(denovo.MS), type = "html")


library(stargazer)
stargazer(denovo.LOF, summary = FALSE)


# Fisher's exact test - testing mutation burden analysis

TeaTasting <-
  matrix(c(3, 1, 1, 3),
         nrow = 2,
         dimnames = list(Guess = c("Milk", "Tea"),
                         Truth = c("Milk", "Tea")))
fisher.test(TeaTasting, alternative = "greater")


library(ggplot2)
ggplot(diamonds, aes(x = depth, colour = cut)) +
  geom_density() +
  xlim(55, 70)
diamonds

############### 22-02 Analysis ----------------------

rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}

adf <- rbind.all.columns(lof.df, miss.df)


a <- c("a", "ab", "bafds", "afdsb", "adfs", "afdsbfds", "bfjkdls", "fdjsl")
a[grep("^a[[:alnum:]]*[^b]+[[:alnum:]]*",a, perl = TRUE)]
b <- a[grep("a.*",a)]
c <- b[!grepl("b",b)]
d <- grep(".*b.*a.*")
c


