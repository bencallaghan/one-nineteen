view.vars <- function(sample,pngout = TRUE){
  
  lof.df.iter %>%
    filter(SAMPLE_ID == sample) %>%
    arrange(PRIORITY) -> p
  miss.df.iter %>%
    filter(SAMPLE_ID == sample) %>%
    arrange(PRIORITY) -> q
  
  print("### LOF MUTATIONS ###")
  print(p)
  print("")
  print("### MS MUTATIONS ###")
  print(q)
  
  if(pngout == TRUE){
    if(nrow(p) > 0){
      png(paste0("out/",sample,".png"), width = 1500, height= (25 * (nrow(p)+1)))
      grid.arrange(tableGrob(p))
      dev.off()
    }else{
      print("No LOF Mutations!")
    }
    
    if(nrow(q) > 0 ){
      png(paste0("out/",sample,"ms.png"), width = 1500, height= (25 * nrow(q)+1))
      grid.arrange(tableGrob(q))
      dev.off()
    }else{
      print("No MS Mutations!")
    }
  }
}

BM_gene_2_tx <- function(hgnc_symbols,mode = "ensembl"){
  
  BMattributes <- c("hgnc_symbol","uniprot_swissprot_accession","refseq_mrna","ensembl_transcript_id","transcript_start","transcript_end" ,
                    "transcript_status","transcript_count", "transcript_biotype", "transcript_source",
                    "transcript_version","transcript_length")
  BMtranscript <- getBM(attributes = BMattributes, 
                        filters = c("hgnc_symbol"),values = list(hgnc_symbols), mart = mart, verbose = FALSE)
  BMtranscript %>% filter(transcript_status == "KNOWN", transcript_biotype == "protein_coding") %>% 
    arrange(desc(transcript_length)) %>% arrange(hgnc_symbol) %>%
    filter(refseq_mrna != "") %>% filter(duplicated(hgnc_symbol) == FALSE) -> BMtranscript.sort
  #   TRANSCRIPT <- BMtranscript.sort$ensembl_transcript_id[1]
  #   return(TRANSCRIPT)
  
  if(mode == 'ensembl'){
    return(BMtranscript.sort)
  }
  
}