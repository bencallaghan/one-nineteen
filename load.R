#Load the loss of function variants (annotated with annovar)## What was the original file?
lof.df <- read.table("inputs/lofAllVariantsMergedFiltered_GDI_Exac.csv", header = T, row.names = NULL, sep = "\t")

# Load the missense variants (annotated with annovar) ## What was the original file?
miss.df <- read.csv("inputs/MS_GDI_Exac.csv", header = T, sep = "\t")
miss.exac <- read.csv("inputs/missense.exac.hg19_multianno.csv", header = T) #rerun.. see note in clean.R

# Load MARV LOF counts / gene ## Is this still up-to-date?
# nateByGene <- read.table("inputs/nate_db_LOF_by_gene", sep = ",", header = T)
marvdb <- read.csv("inputs/marvdbdump-03-07-17.csv", header = T, sep = "\t", fill = T)
marvdb %>% dplyr::select(chromosome,start_hg19, stop_hg19, ref, alt, denovo) -> marvdb.anno.in
colnames(marvdb.anno.in) <- c("chr", "start", "stop", "ref", "alt", "denovo")
write.table(x = marvdb.anno.in, file = "temp/marv.anno.in" , quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
marvdb <- read.csv("temp/marv.anno.in.hg19_multianno.csv")

if(FALSE){
  annocmd <- paste0("perl /space/bin/annovar/table_annovar.pl ", "temp/marv.anno.in", " /space/bin/annovar/humandb/ -buildver hg19 -out ", "outputs/marv.anno.out", " -remove -protocol refGene,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,ljb_all,exac03,cadd -otherinfo -operation g,r,f,f,f,f,f,f -nastring . -csvout")
  cmd.out <- run.remote(cmd=annocmd , remote= "apu")
  annocmd2 <- paste0("perl /space/bin/annovar/annotate_variation.pl ", "temp/marv.anno.in", " /space/bin/annovar/humandb -filter -dbtype cadd -buildver hg19 -out outputs/marv.anno.out2 -otherinfo")
  cmd.out <- run.remote(cmd=annocmd2 , remote= "-q apu",stderr.redirect=F)
}

# Demographics and DNA availability
dna.info <- read.table("inputs/DNA_availability.csv", sep = ",", header = T) # DNA availability info (from Kristina)



# annocmd <- paste0("perl /space/bin/annovar/table_annovar.pl ", "temp/marv.anno.in",
#                   " /space/bin/annovar/humandb/ -buildver hg19 -out ", "outputs/marv.anno.out",
#                   " -remove -protocol refGene,genomicSuperDups,esp6500si_all,",
#                   "1000g2012apr_all,snp135,ljb_all,exac03,cadd -otherinfo -operation g,r,f,f,f,f,f,f -nastring . -csvout")
# cmd.out <- run.remote(cmd=annocmd , remote= "-q troy", stderr.redirect = F)

# run_annovar_on_vcf_df(marvdb.anno.in, tempPath = "temp/marv.anno.in", outputPath = "outputs/marv.anno.out")

# Sander's genes list
sandy.genes <- c("ADNP", "ANK2", "ARID1B", "ASH1L", "CHD2", "CHD8", "CUL3", "DSCAM", "DYRK1A", "GRIN2B", "KATNAL2", "KDM5B", "KMT2C", 
                 "NCKAP1", "POGZ", "SCN2A", "SUV420H1", "SYNGAP1", "TBR1", "TCF7L2", "TNRC6B", "WAC", "NRXN1", "PTEN", "SETD5", 
                 "SHANK2", "SHANK3", "TRIP12", "BCL11A", "FOXP1", "GIGYF1", "ILF2", "KDM6B", "PHF2", "RANBP17", "SPAST", "WDFY3", 
                 "DNMT3A", "GABRB3", "KAT2B", "MFRP", "MYT1L", "P2RX5", "MIB1", "SLC6A1", "ZNF559")


# Exac constraints list (all genes)
exac.all.genes <- read.table("inputs/exac_constraints.txt", header = T)
gdi.all.genes <- read.table("inputs/GDI_full.txt", header = T, sep = "\t")


# DNA availability info (from Kristina)
dna.info <- read.table("inputs/DNA_availability.csv", sep = ",", header = T)
maybes <- c("2047-24322","48-14470","1989-24099","21758-34760","1588-22326","1937-23920","19991-31601","779-18300")
dna.info %>% filter(Both.parents != "Y") %>% dplyr::select(Sample.Name)

# Sanger validation results
sanger.results <- read.csv("inputs/sanger_results2.csv", header = T)


# MSSNG
# mssng.wdn <- read.table("inputs/MSSNG_all_affected_variants_annotated.tsv", sep = "\t", header = T)
# head(mssng.wdn$GeneDetail.refGene[ grep("+", mssng.wdn$GeneDetail.refGene, perl = TRUE) ] )
# mssng.wdn %>% head

