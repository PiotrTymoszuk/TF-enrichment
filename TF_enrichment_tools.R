# A medley of functions to analyze TF binding site enrichment in the given gene set
# The database of whole-genome promoter boinding sites was done with D-Light
# https://pbwww.services.came.sbg.ac.at/?page_id=40

# data and tools ------

  load('./tools/dlight_db/dlight_db.RData')


  require(plyr)
  require(tidyverse)
  require(furrr)

# functions ------
  
  count_tf_gene_set <- function(entrez_id_vec, dlight_db = dlight_db_lite) {
    
    ## gives back the number of binding sites in the gene set defined by entrez_id_vec
    ## for each TF present in the given database
    ## entrez IDs absent from the TF database are skipped silently
    
    gen_set_sites <- dlight_db$count_tbl[entrez_id_vec, ] %>% 
      map_dfc(sum, na.rm = T)

    return(gen_set_sites)
    
  }
  
  compare_count_sets <- function(count_tbl1, count_tbl2, dlight_db = dlight_db_lite, .parallel = F) {
    
    ## checks if the counts of binding sites in the count table 2
    ## are larger than in the count table 1
    ## it a list of count tables 2 is provided, the sum of cases when
    ## the number of TF binding sites is larger than in the cont table 1
    ## is returned
    
    if(any(class(count_tbl2) == 'list')) {
      
      ## comparison results
      
      if(.parallel) {
        
        plan('multisession')
        
        com_lst <- count_tbl2 %>% 
          future_map(compare_count_sets, 
                     count_tbl1 = count_tbl1)
        
        plan('sequential')
        
      } else {
        
        com_lst <- count_tbl2 %>% 
          map(compare_count_sets, 
              count_tbl1 = count_tbl1)
        
      }
      
      comparison <- com_lst %>% 
        map_dfr(function(x) x$comparison) %>% 
        map_dfc(sum)
      
      ## bootstrap p value: the wors-case scenario is assumed (1 if count of cases is 0)
      
      p_values <- comparison %>% 
        map(function(x) if(x[1] == 0) 1/length(com_lst) else x/length(com_lst)) %>% 
        tibble(feature_acc = names(.), 
               p_raw = unlist(.)) %>% 
        select(feature_acc, 
               p_raw)
      
      ## enrichment over the whole genome occurence
      
      enrichment <- com_lst %>% 
        map_dfr(function(x) x$enrichment)
      
      enrichment_stats <- enrichment %>% 
        map(function(x) tibble(median = median(x, na.rm = T), 
                               lower_ci = quantile(x, 0.025, na.rm = T), 
                               upper_ci = quantile(x, 0.975, na.rm = T))) %>% 
        map2_dfr(., names(.), 
                 function(x, y) mutate(x, feature_acc = y))
      
      enrichment_stats <- left_join(enrichment_stats, 
                                    p_values, 
                                    by = 'feature_acc') %>% 
        left_join(., 
                  dlight_db$annotation_tf, 
                  by = 'feature_acc')
      
      return(enrichment_stats)
      
    }
    
    comp_res <- map2_dfc(count_tbl1, 
                         count_tbl2, 
                         function(x, y) as.numeric(x < y)) %>% 
      set_names(names(count_tbl1)) ## number of cases
    
    enrich_res <- map2_dfc(count_tbl1, 
                           count_tbl2, 
                           function(x, y) if(y > 0) x/y else x) %>% 
      set_names(names(count_tbl1)) ## enrichment
    
    return(list(comparison = comp_res, 
                enrichment = enrich_res))
    
  }
  
  create_db_bootstraps <- function(boot_size = 100, boot_n = 100, dlight_db = dlight_db_lite) {
    
    ## creates a random subsets of the TF database
    
    boot_ids <- paste('boot', 1:boot_n, sep = '_')
    
    boots <- boot_ids %>% 
      map(function(x) sample(1:nrow(dlight_db$count_tbl), 
                             size = boot_size)) %>% 
      map(function(x) dlight_db$count_tbl[x, ]) %>% 
      set_names(boot_ids)
    
    return(boots)
    
  }
  
  test_tf_enrichment <- function(entrez_id_vec, boot_size = length(entrez_id_vec), boot_n = 100, 
                                 dlight_db = dlight_db_lite, .parallel = F) {
    
    ## the heartpice of the toolbox: checks how many random bootstraps of the database have more
    ## binding sites for the particular TF than the gene set of interest defined by the entrez ID vector
    
    ## random subsets of the whole-genome TF binding database
    
    start_time <- Sys.time()
    
    boots <- create_db_bootstraps(boot_size = boot_size, 
                                  boot_n = boot_n, 
                                  dlight_db = dlight_db)
    
    count_boots <- boots %>% 
      map(function(x) rownames(x)) %>% 
      map(count_tf_gene_set, 
          dlight_db = dlight_db)
    
    genes_interest_counts <- count_tf_gene_set(entrez_id_vec = entrez_id_vec, 
                                               dlight_db = dlight_db)
    
    com_table <- compare_count_sets(count_tbl1 = genes_interest_counts, 
                                    count_tbl2 = count_boots, 
                                    dlight_db = dlight_db, 
                                    .parallel = .parallel)
    
    message(paste('Elapsed:', Sys.time() - start_time))
    
    return(com_table)
    
  }
  
# END ----