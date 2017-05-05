

# LOF variants dataframe cleaning
lof.df %>% dplyr::select(c(1:9,15:21,26,46,53:61)) -> lof.df.filt
lof.df.filt$DnLOF[is.na(lof.df.filt$DnLOF)] <- 0

# MS variants dataframe cleaning
miss.df$dbNSFP_CADD_phred <- as.numeric(as.character(miss.df$dbNSFP_CADD_phred))
miss.df$SAMPLE_ID <- as.character(miss.df$SAMPLE_ID)
miss.df.filtered <- miss.df[which(miss.df$dbNSFP_CADD_phred >= 25),]

missense4anno <- miss.df[,c(2,3,3,7,4)]
names(missense4anno) <- c("CHR","start","stop","ref","alt")
write.table(missense4anno, file="temp/missense4anno", quote=FALSE, sep = "\t", row.names= FALSE)
# TODO: Run annovar here (run.remote())...for now just use old output (missense.exac.hg19_multianno.csv)

miss.exac <- miss.exac[-1,]
miss.exac$exac03[which(miss.exac$exac03 == ".")] <- 0
miss.df$exac03 <- as.numeric(as.character(miss.exac$exac03))
miss.df$natedb.lof.score[which(is.na(miss.df$natedb.lof.score))] <- 0


# Clean for mutation summary plotting
lof.df.filt.sort <- within(lof.df.filt, 
                           SAMPLE_ID <- factor(SAMPLE_ID, 
                                               levels=names(sort(table(SAMPLE_ID), 
                                                                 decreasing=TRUE))))

# Clean gene scores dfs
exac.all.genes %>% dplyr::select(c(2,20)) -> exac.all.genes
gdi.all.genes %>% dplyr::select(c(1,2,3)) -> gdi.all.genes

exac.all.genes$gene <- as.character(exac.all.genes$gene)
gdi.all.genes$Gene <- as.character(gdi.all.genes$Gene)

# I should merge("complete.obs")
# TODO: use merge
gdi.all.genes <- gdi.all.genes[which(gdi.all.genes$Gene %in% exac.all.genes$gene),]
exac.all.genes <- exac.all.genes[which(exac.all.genes$gene %in% gdi.all.genes$Gene),]
gdi.all.genes <- gdi.all.genes[order(gdi.all.genes$Gene),]
exac.all.genes <- exac.all.genes[order(exac.all.genes$gene),]
gdi.all.genes$Gene == exac.all.genes$gene
gdi.exac <- cbind(gdi.all.genes, exac.all.genes)


# sanger.results
# sanger.results$Gene.1 %>% 
# gsub(pattern = "-F" , replacement = "" , x = . ) %>% 
# gsub(pattern = "[0-9]+-", replacement = "", x = ., perl = TRUE) -> sanger.results$Gene.1
sanger.results %>% 
  mutate(gchange = paste0(Location,"-", BaseChange)) %>% 
  dplyr::select(-Location, -BaseChange) -> sanger.results

# Sanger results, dn / inherited
sanger.denovo <- as.factor(sanger.results[sanger.results$inheritance == "de novo",]$gchange)
sanger.inherited <- as.factor(sanger.results[sanger.results$inheritance == "inherited",]$gchange)

# DNA availability, mpx/spx
dna.info$Sample.Name[grep("Y", dna.info$ASD.affected.Sib)] -> samples.mpx
dna.info$Sample.Name[grep("N", dna.info$ASD.affected.Sib)] -> samples.spx
dna.info %>% filter(Both.parents != "Y") %>% dplyr::select(Sample.Name) -> samples.notrio # Samples for which we do not have full trios
dna.info %>% filter(Both.parents == "Y") %>% dplyr::select(Sample.Name) -> samples.trio # Samples for which we have full trios
