# Analysis script for the 119 project
# Run load / clean / func scripts first
## Summary plots of mutations ---------------------------------------------

# LOF counts per subject 
# TODO: rotate sample ID names
p <- ggplot(lof.df.filt.sort, aes(factor(SAMPLE_ID))) + geom_point(stat = "count")
p

# Exac frequency of LOF variants 
p <- ggplot(lof.df.filt, aes(exac03)) + geom_histogram()
print(p)
ggsave("plots/ExacDistributionLOF",plot=p,device="png")

# Exac Frequency of LOF mutation / LOF containing gene
# This doesn't mean much, really need all genes to get the full picture
p <- ggplot(lof.df.filt, aes(ExACgeneLOF)) + geom_histogram(bins=100)
print(p)
# ggsave("plots/ExacLOFCountsperGene",plot=p,device="png")

# GDI of LOF containing genes (phred scaled)
# Again, this doesn't mean much as we need all the genes for full picture
p <- ggplot(lof.df.filt, aes(GDI.Phred)) + geom_histogram() + geom_vline(xintercept = 13.84, colour = "red")
print(p)
# ggsave("plots/GDIscores",plot=p,device="png") 

#Exac Constraint scores / LOF containing gene
#See above
p <- ggplot(lof.df.filt, aes(Exac.pLI)) + geom_histogram()
print(p)
ggsave("plots/ExacPLIdistribution",plot=p,device="png") 

#Exac n for LOF / gene
#
p <- ggplot(lof.df.filt, aes(Exac.n.LOF)) + geom_histogram()
print(p)
ggsave("plots/ExacNLOF",plot=p,device="png") 

# Exac pLI plots
ggplot(lof.df.filt, aes(Exac.pLI)) + geom_histogram() + geom_vline(xintercept = .90, colour = "red")
ggplot(lof.df.filt, aes(x = rank(Exac.pLI), y = GDI.score)) + geom_point()
ggplot(lof.df.filt, aes(y = Exac.pLI, x = GDI.Phred)) + geom_point()
ggplot(lof.df.filt, aes(x = -10*log10(Exac.pLI), y = GDI.Phred)) + geom_point()
ggplot(lof.df.filt, aes(x = -10*log10(Exac.pLI), y = log10(GDI.score))) + geom_point() 
p

## EXAC VS GDI

ggplot(gdi.exac, aes(x = pLI, y = GDI.Phred)) + geom_point(alpha = 0.05) 
ggplot(gdi.exac, aes(y = pLI, x = GDI.Phred)) + geom_point(alpha = 0.05)
ggplot(gdi.exac, aes(x = pLI, y = GDI)) + geom_point(alpha = 0.05) + scale_y_log10()



# By subject Variant prioritisation ---------------------------------------
# MS variants for by-subject prioritisation steps
miss.df %>% dplyr::select(c(5,1,2,3,7,4,6,12,13,15,18,24,31,47,40,42)) %>% 
  filter(dbNSFP_CADD_phred >= 25) %>% 
  mutate(Priority = rank(
    rank((exac03), ties = "average", na.last="keep" ) +
      rank(desc(Exac.MSZ), na.last="keep") +
      rank(desc(dbNSFP_CADD_phred), na.last = "keep") +
      rank(GDI.Phred, na.last = "keep"), na.last = "keep")) -> miss.df.iter
names(miss.df.iter) <- c("SAMPLE_ID", "GENE", "CHROM", "POS", "REF","ALT","ID","HGVS_C","HGVS_P","GENOTYPE","EFFECT","CADD_PHRED","NATEDB_LOF_SCORE","EXAC03_FREQ","GDI_PHRED","EXAC_MSZ","PRIORITY")

# Write CSV of top MS vars
miss.df.iter %>% filter(CADD_PHRED > 25 & EXAC_MSZ >3 & GDI_PHRED < 12 & EXAC03_FREQ == 0) %>% arrange(PRIORITY) -> topms66
write.table(topms66,file='outputs/topms66.tsv', sep = '\t', quote=FALSE)

# LOF variants for by-subject prioritisation steps
lof.df.filt %>% dplyr::select(c(2,1,3,4,5,6,7,10,11,13,15,16,17,19,21,24,25)) %>%   
  mutate(Priority = rank(
    rank((exac03), ties = "average", na.last="keep" ) +
      rank(desc(Exac.pLI), na.last="keep") +
      rank(desc(DnLOF), na.last = "keep") +
      rank(GDI.Phred, na.last = "keep"), na.last = "keep")) -> lof.df.iter

# Write CSV of top LOF vars
lof.df.iter %>% filter(Exac.pLI > 0.9) %>% arrange(Priority) -> toplof31
write.table(toplof31,file='outputs/toplof31.tsv', sep = '\t', quote=FALSE)

# View variants per subject
# view.vars("18436-29827")
# 
# lof.df.filt %>% filter(SAMPLE_ID == "15210-25677")
# miss.df.iter %>% filter(SAMPLE_ID == "15210-25677") %>% arrange(desc(CADD_PHRED * EXAC_MSZ))


# toplof31 %>% filter(SAMPLE_ID %in% c("21705-34281" ))

# MARV counts ---------------------------------#

exac.constraints <- read.table("inputs/exac_constraints.txt", header = T)
exac.constraints$gene <- levels(exac.constraints$gene)[exac.constraints$gene]

marvdb <- read.csv("/home/bcallaghan/Projects/Onenineteen/temp/marv.anno.hg19_multianno.csv")
marvdb %>% dplyr::select(everything(), inheritance = Otherinfo) -> marv.wdn


# marv.wdn <- read.table("inputs/marv_all_affected_variants_annotated.tsv", sep = "\t", header = T)
marv.wdn$exonicFunc <- levels(marv.wdn$exonicFunc)[marv.wdn$exonicFunc]
marv.wdn$exonicFunc[grepl(".+\\+[1-2][ACGT]>.+", marv.wdn$GeneDetail.refGene)] <- "splicing"
marv.wdn$exonicFunc[grepl(".+\\-[1-2][ACGT]>.+", marv.wdn$GeneDetail.refGene)] <- "splicing"
marv.wdn$CADD_phred <- as.numeric(levels(marv.wdn$CADD_phred)[marv.wdn$CADD_phred])
marv.wdn$ExAC_ALL[marv.wdn$ExAC_ALL == "."] <- 0

marv.wdn %>% filter(CADD_phred > 0) -> marv.c0
marv.c0 <- merge(marv.c0, exac.constraints, by.x = 'Gene.refGene', by.y = 'gene')
marv.c0 %>% filter(inheritance == "yes") -> marv.c0.dn
marv.c0.dn %>% filter(CADD_phred > 25) -> marv.c25.dn
marv.c0.dn %>% filter(CADD_phred > 25 & ExAC_ALL == 0 & (mis_z > 3 | pLI > 0.9)) -> marv.c25.ex.con

marvcounts.c25.ex.con <- as.data.frame(table(marv.c25.ex.con$Gene.refGene))
marvcounts.c25.ex.con.mrg <- merge(marvcounts.c25.ex.con, exac.constraints, by.x = 'Var1', by.y = 'gene', all = T)
marvcounts.c25.ex.con.mrg %>% dplyr::rename(MARV.Freq = Freq) %>% mutate(MARV.Freq.corr = MARV.Freq / bp) %>% arrange(desc(MARV.Freq.corr)) -> marvcounts.c25.ex.con.mrg.cor
head(marvcounts.c25.ex.con.mrg.cor)


# MSSNG counts --------------------------------#

mssng.wdn <- read.table("inputs/MSSNG_all_affected_variants_annotated.tsv", sep = "\t", header = T)
mssng.wdn$exonicFunc <- levels(mssng.wdn$exonicFunc)[mssng.wdn$exonicFunc]
mssng.wdn$exonicFunc[grep(".+\\+[1-2][ACGT]>.+", mssng.wdn$GeneDetail.refGene)] <- "splicing"
mssng.wdn$exonicFunc[grep(".+\\-[1-2][ACGT]>.+", mssng.wdn$GeneDetail.refGene)] <- "splicing"

# mssng.wdn %>% filter(!is.na(exonicFunc)) %>% filter(exonicFunc != "synonymous SNV") -> mssng.lgd

## filter(CADD_PHRED > 25 & EXAC_MSZ > 3 & GDI_PHRED < 12 & EXAC03_FREQ == 0) %>% arrange(PRIORITY) 
## Filtering criteria for variants in the 119, should try and keep criteria the same
mssng.wdn %>% filter(CADD_phred > 0) -> mssng.c0

mssng.c0 <- merge(mssng.c0, exac.constraints, by.x = 'Gene.refGene', by.y = 'gene')

mssng.c0 %>% filter(inheritance == "0/0") -> mssng.c0.dn
mssng.c0.dn %>% filter(CADD_phred > 25) -> mssng.c25.dn
mssng.c0.dn %>% filter(CADD_phred > 25 & ExAC_ALL == 0 & (mis_z > 3 | pLI > 0.9)) -> mssng.c25.ex.con

# Gene-wise counts of higher priority mssng variants
# mssngcounts.c25 <- as.data.frame(table(mssng.c25.dn$Gene.refGene))
# mssngcounts.c25.mrg <- merge(mssngcounts.c25, exac.constraints, by.x = 'Var1', by.y = 'gene', all = T)
# mssngcounts.c25.mrg %>% mutate(Freq.corr = Freq / bp) %>% arrange(desc(Freq.corr)) -> mssngcounts.c25.mrg.cor

mssngcounts.c25.ex.con <- as.data.frame(table(mssng.c25.ex.con$Gene.refGene))
mssngcounts.c25.ex.con.mrg <- merge(mssngcounts.c25.ex.con, exac.constraints, by.x = 'Var1', by.y = 'gene', all = T)
mssngcounts.c25.ex.con.mrg %>% dplyr::rename(MSSNG.Freq = Freq) %>%  mutate(MSSNG.Freq.corr = MSSNG.Freq / bp) %>% arrange(desc(MSSNG.Freq.corr)) -> mssngcounts.c25.ex.con.mrg.cor
head(mssngcounts.c25.ex.con.mrg.cor)

# Merge 119 variants with MARV and MSSNG counts for other sources of "interestingness"
miss.df.mssng <- merge(miss.df, mssngcounts.c25.ex.con.mrg.cor, by.x = 'Gene', by.y = 'Var1')
lof.df.mssng <- merge(lof.df, mssngcounts.c25.ex.con.mrg.cor, by.x = 'Gene', by.y = 'Var1')

miss.df.mssng.marv <- merge(miss.df.mssng, marvcounts.c25.ex.con.mrg.cor, by.x = 'Gene', by.y = 'Var1')
lof.df.mssng.marv <- merge(lof.df.mssng, marvcounts.c25.ex.con.mrg.cor, by.x = 'Gene', by.y = 'Var1')


# For Sanja's Poster / Evaluating the impact of MSSNG variants on the 119 project
# How do MSSNG variants change our conceptions of what the "interesting" variants of the 119 project are?
# Do they cause us to re-evaluate an interesting gene as not interesting (no, I doubt it)
# Could they cause us to re-evaluate a less interesting gene as more interesting, maybe even to the point
# whereby why should have resequenced it? ( << yeah, more this one )

# Check Loss of function variants for interestingness according to MSSNG recurrence (Freq.corr)
lof.df.mssng.marv %>% 
  mutate(gchange = paste0(X.CHROM, ":", POS, "-", REF, "/", ALT)) %>% 
  mutate(Family = ifelse(SAMPLE_ID %in% samples.spx,"spx","mpx")) %>%
  mutate(Inheritance = ifelse(gchange %in% sanger.denovo, yes = "de novo", # Why do you need to search within factors? **Find out
                              no = ifelse(gchange %in% sanger.inherited, yes = "inherited", no = NA))) %>% 
  dplyr::select(Gene, Sample_ID = SAMPLE_ID, Family, gchange, Inheritance, Effect = EFF.0..EFFECT, 
                MSSNG.Freq, MARV.Freq, ExAC_Freq = exac03, MSSNG.Freq.corr, MARV.Freq.corr, pLI.y, mis_z.x) %>% 
  arrange(desc(MARV.Freq.corr)) 


# Check missense variants for interestingness according to MSSNG recurrence (Freq.corr)
miss.df.mssng.marv %>%
  filter(dbNSFP_CADD_phred > 25) %>% filter(exac03 == 0) %>% filter(mis_z.x > 3) %>% 
  mutate(gchange = paste0(X.CHROM, ":", POS, "-", REF, "/", ALT)) %>% 
  mutate(Family = ifelse(SAMPLE_ID %in% samples.spx,"spx","mpx")) %>% 
  mutate(Inheritance = ifelse(gchange %in% sanger.denovo, yes = "de novo", # Why do you need to search within factors? **Find out
                              no = ifelse(gchange %in% sanger.inherited, yes = "inherited", no = NA))) %>% 
  dplyr::select(Gene, Sample_ID = SAMPLE_ID, Family, gchange, Inheritance, Effect = EFF.0..EFFECT, 
                MSSNG.Freq, MARV.Freq, CADD.phred = dbNSFP_CADD_phred, ExAC_Freq = exac03, MSSNG.Freq.corr, MARV.Freq.corr, pLI.x, mis_z.x) %>% 
  arrange(desc(MARV.Freq.corr), desc(MSSNG.Freq.corr))








# GEA with topGO ----------------

# Data Prep
# Running the enrichment tests
# Analysis of the results





