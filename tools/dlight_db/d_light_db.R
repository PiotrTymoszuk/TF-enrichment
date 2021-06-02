# this script clears the raw output of the D-Light TF scanning software:
# The table with all transcripts and binding sites in their promoters

## libraries -----

  library(plyr)
  library(tidyverse)
  library(stringi)

# data container ----

  dlight_db <- list()
  
# functions ----

  
  dlight_db$count_tf_sites <- function(transcript, tf, db_tbl = dlight_db$raw_db) {
    
    ## counts the binding sites of a given TF or TFs in the transcript promoter

    tf_count <- db_tbl %>% 
      filter(transcript_id == transcript, 
             feature_acc %in% tf) %>% 
      count(feature_acc)
    
    return(tf_count)
    
  }
  
  dlight_db$count_tf_list <- function(transcript_vec, tf, db_tbl_lst) {
    
    ## a wrapper for handling vectors of transcripts and tables
    
    start_time <- Sys.time()
    
    message(paste('Counting TFs for', length(transcript_vec), 'transcripts'))

    res_lst <- list(transcript = transcript_vec, 
                    db_tbl = db_tbl_lst) %>% 
      pmap(dlight_db$count_tf_sites, 
           tf = tf)
    
    res_lst <- map2(res_lst, 
                    transcript_vec, 
                    function(x, y) set_names(x, c('feature_acc', y)))
 
    res_lst <- res_lst %>% 
      map(function(x) left_join(tibble(feature_acc = tf), 
                                x, 
                                by = 'feature_acc')[, 2]) %>% 
      do.call(cbind, .)

    res_lst <- res_lst %>% 
      map_dfc(function(x) ifelse(is.na(x), 0, x)) %>% 
      cbind(feature_acc = tf, .) %>% 
      as_tibble
    
    message(paste('Elapsed:', Sys.time() - start_time))
    
    return(res_lst)
    
  }
  
  
# reading the program output, identifying unique transcripts and TFs ----
  
  dlight_db$raw_db <- read_tsv('./tools/dlight_db/all_gene_d_light.csv') %>% 
    mutate(transcript_id = stri_replace(gene_acc, 
                                        regex = '\\.\\d+$', 
                                        replacement = '')) ## transcript ids without entry version
  
  dlight_db$feature_vec <- unique(dlight_db$raw_db$feature_acc) ## unique TFs identified by the software
  
  dlight_db$feature_vec <- dlight_db$feature_vec[!dlight_db$feature_vec %in% c('TSS', 'CSS')] ## we don't need the transcriptiion starting sites
  
  dlight_db$trans_vec <- unique(dlight_db$raw_db$transcript_id) ## unique transcripts identified by the software
  
  
# making a DB table split into unique transcripts to accelerate counting, restricting to unique transcripts -----
  
  dlight_db$split_db <- dlight_db$raw_db %>% 
    dlply(.(transcript_id), as_tibble)

# for each transcript, the the number of TF binding sites for each TF in the database is counted -----

  dlight_db$count_tbl <- dlight_db$count_tf_list(transcript_vec = dlight_db$trans_vec, 
                                                 tf = dlight_db$feature_vec, 
                                                 db_tbl_lst = dlight_db$split_db[dlight_db$trans_vec], 
                                                 .parallel = T)

  dlight_db$count_tbl <- dlight_db$count_tbl %>% 
    column_to_rownames('feature_acc') %>% 
    t %>% 
    as.data.frame 
    
  
# annotation table, re-naming the count table with entrez IDs -----
  
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  
  dlight_db$annotation_vec <- mapIds(org.Mm.eg.db, 
                                     keys = rownames(dlight_db$count_tbl), 
                                     keytype = 'REFSEQ', 
                                     column = 'ENTREZID', 
                                     multiVals = 'first')
  
  detach('package:org.Mm.eg.db', unload = T, force = T)
  detach('package:AnnotationDbi', unload = T, force = T)
  detach('package:IRanges', unload = T, force = T)
  detach('package:S4Vectors', unload = T, force = T)
  detach('package:Biobase', unload = T, force = T)
  detach('package:BiocGenerics', unload = T, force = T)
  
  ## restricting the transcripts to those with valid, non-duplicated entrez ids
  
  dlight_db$annotation_vec <- dlight_db$annotation_vec[!is.na(dlight_db$annotation_vec)]
  
  dlight_db$annotation_vec <- dlight_db$annotation_vec[!duplicated(dlight_db$annotation_vec)]

  ## renaming the count table
  
  dlight_db$count_tbl <- dlight_db$count_tbl[names(dlight_db$annotation_vec), ]
  
  rownames(dlight_db$count_tbl) <- dlight_db$annotation_vec[rownames(dlight_db$count_tbl)]
  
# Annotation of the TFs -----
  
  dlight_db$annotation_tf <- read_tsv('./tools/dlight_db/tf_annotation.csv') %>% 
    set_names(c('feature_acc', 
                'tf_symbol'))
  
# saving the database on the disc ----
  
  dlight_db_lite <- dlight_db[c('count_tbl', 
                                'annotation_vec',
                                'annotation_tf')]
  
  save(dlight_db_lite, 
       file = './tools/dlight_db/dlight_db.RData')

# END -----
  