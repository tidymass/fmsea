#' ####do_msea
#' #' @title do_msea
#' #' @description Annotate isotopes for known formula compounds.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param annotation_table annotation_table
#' #' @param metabolite_set metabolite_set
#' #' @param metabolite_set_name metabolite_set_name
#' #' @param exponent exponent
#' #' @param perm_num perm_num
#' #' @param min_size min_size
#' #' @param max_size max_size
#' #' @param pvalue_cutoff pvalue_cutoff
#' #' @param p_adjust_method p_adjust_method
#' #' @param seed seed
#' #' @param verbose verbose
#' #' @param ... other parameters
#' #' @export
#' 
#' 
#' no_source()
#' sxtTools::setwd_project()
#' rm(list = ls())
#' source("R/annotate_isotope.R")
#' source("R/utils.R")
#' source("R/annotate_annotation_table.R")
#' 
#' setwd("demo_data/denmark_project/metabolome/peaks/")
#' 
#' library(tidyverse)
#' 
#' load("annotation_table")
#' library(metPath)
#' load("kegg_hsa_pathway.rda")
#' 
#' metabolite_set = kegg_hsa_pathway@compound_list[[1]]$KEGG.ID
#' metabolite_set_name = kegg_hsa_pathway@pathway_name[1]
#' metabolite_set_id = kegg_hsa_pathway@pathway_id[1]
#' exponent = 1
#' perm_num = 1000
#' min_size = 5
#' max_size = 1000
#' pvalue_cutoff = 0.2
#' p_adjust_method = "BH"
#' seed = FALSE
#' verbose = TRUE
#' threads = 3
#' 
#' fmse_single_set <-
#'   function(annotation_table,
#'            metabolite_set,
#'            metabolite_set_name = "metabolite_set",
#'            metabolite_set_id = "set1",
#'            exponent = 1,
#'            perm_num = 1000,
#'            min_size = 5,
#'            max_size = 1000,
#'            pvalue_cutoff = 0.2,
#'            p_adjust_method = c("fdr",
#'                                "holm",
#'                                "hochberg",
#'                                "hommel",
#'                                "bonferroni",
#'                                "BH",
#'                                "BY",
#'                                "none"),
#'            seed = FALSE,
#'            verbose = TRUE,
#'            threads = 3,
#'            ...) {
#'     p_adjust_method = match.arg(p_adjust_method)
#'     ##filter pathways according to pathway size
#'     if (length(metabolite_set) < min_size |
#'         length(metabolite_set) > max_size) {
#'       return(NULL)
#'     }
#'     
#'     ##find where are the metabolites in the feature table
#'     ##for each metabolite in metabolite set, its position in feature table
#'     match_idx <-
#'       purrr::map(metabolite_set, function(x) {
#'         which(x == annotation_table$Lab.ID)
#'       })
#'     
#'     purrr::map2(
#'       .x = match_idx ,
#'       .y = 1:length(match_idx),
#'       .f = function(x, y) {
#'         if (length(x) == 0) {
#'           return(NULL)
#'         }
#'         data.frame(x, y)
#'       }
#'     ) %>%
#'       dplyr::bind_rows() %>%
#'       ggplot(aes(y, x)) +
#'       geom_point()
#'     
#'     ##----------------------------------------------------------------------------
#'     ##for the metabolites in the feature set,
#'     ##all the matched peaks should be remained, other metabolites
#'     ##only remain one matched peak,
#'     ##either one is fine
#'     ###idx is the index of all metabolites in metabolite set in annotation_table
#'     idx <-
#'       sort(unique(unlist(match_idx)))
#'     
#'     if (length(idx) == 0) {
#'       return(NULL)
#'     }
#'     
#'     ##for the metabolites in metabolite sets.
#'     ##annotation_table1 is the features that have the annotation in metabolite_set
#'     annotation_table1 <-
#'       annotation_table[idx, ]
#'     
#'     #for each metabolite_feature_cluster, if the peaks are in different condition_class
#'     #(for example, if the condition is fold change, some peaks fold change > 1 and some
#'     #peaks fold change < 1), we should think some method to handle this condition
#'     
#'     ##for the metabolites are not in the feature set, they only affect the rank of
#'     ##the metabolites, so only remain one time for each peak.
#'     annotation_table2 <-
#'       annotation_table[-idx, ] %>%
#'       dplyr::distinct(name, .keep_all = TRUE)
#'     
#'     ###annotation_table1 and annotation_table2 have a lot of overlapped features
#'     ##we should remainthe overlapped features in feature table2, because the annotation is not confirmed
#'     # intersect(annotation_table1$name, annotation_table2$name)
#'     
#'     annotation_table <-
#'       rbind(annotation_table1,
#'             annotation_table2) %>%
#'       dplyr::arrange(dplyr::desc(condition))
#'     
#'     if (verbose) {
#'       message("calculating observed enrichment scores...")
#'     }
#'     
#'     ###observed_info is the running enrichment score for each metabolite (peak)
#'     ###in annotation_table. And ES is the enrichment score
#'     observed_info =
#'       get_msea_score(
#'         annotation_table = annotation_table,
#'         metabolite_set = metabolite_set,
#'         exponent = exponent
#'       )
#'     
#'     observed_score <- observed_info$enrichment_score
#'     
#'     if (verbose) {
#'       message("calculating permutation scores...")
#'     }
#'     
#'     ###here, we need to do the permutation test
#'     if (seed) {
#'       seeds <- sample.int(perm_num)
#'     }
#'     
#'     if (masstools::get_os() == "windows") {
#'       bpparam =
#'         BiocParallel::SnowParam(workers = threads,
#'                                 progressbar = TRUE)
#'     } else{
#'       bpparam = BiocParallel::MulticoreParam(workers = threads,
#'                                              progressbar = TRUE)
#'     }
#'     
#'     permute_scores <-
#'       BiocParallel::bplapply(1:perm_num, function(i) {
#'         if (seed) {
#'           set.seed(seeds[i])
#'         }
#'         
#'         permute_msea_enrichment_score(
#'           annotation_table = annotation_table,
#'           metabolite_set = metabolite_set,
#'           exponent = exponent
#'         )
#'       }, BPPARAM = bpparam) %>%
#'       unlist()
#'     
#'     ####mean_pos_null_es is the mean value of all the positive null enrichment score
#'     mean_pos_null_es <- mean(permute_scores[permute_scores >= 0])
#'     ####mean_neg_null_es is the mean value of all the negative null enrichment score
#'     mean_neg_null_es <- mean(permute_scores[permute_scores < 0])
#'     
#'     ###normalize_es is the function used to normalize enrichment_score
#'     ##enrichment_score is the enrichment score
#'     normalize_es <-
#'       function(enrichment_score,
#'                mean_pos_null_es,
#'                mean_neg_null_es) {
#'         s <- sign(enrichment_score)
#'         m <- numeric(length(enrichment_score))
#'         m[s == 1] <- mean_pos_null_es[s == 1]
#'         m[s == -1] <- mean_neg_null_es[s == -1]
#'         enrichment_score / m
#'       }
#'     
#'     ###normalized_es is the normalized enrichment_score
#'     normalized_es <-
#'       normalize_es(observed_score, mean_pos_null_es, mean_neg_null_es)
#'     
#'     permute_scores = normalize_es(
#'       enrichment_score = observed_score,
#'       mean_pos_null_es = mean_pos_null_es,
#'       mean_neg_null_es = mean_neg_null_es
#'     )
#'     
#'     if (verbose) {
#'       message("calculating p values...")
#'     }
#'     
#'     ###get the p value of enrichment_score
#'     if (is.na(normalized_es)) {
#'       p_value = NA
#'     }
#'     
#'     if (normalized_es >= 0) {
#'       (sum(permute_scores[i,] >= normalized_es[i]) + 1) / (sum(permute_scores[i, ] >= 0) + 1)
#'     }
#'     
#'     p_values <- sapply(seq_along(observed_score), function(i) {
#'       if (is.na(normalized_es[i])) {
#'         NA
#'       } else if (normalized_es[i] >= 0) {
#'         (sum(permute_scores[i,] >= normalized_es[i]) + 1) / (sum(permute_scores[i, ] >= 0) + 1)
#'       } else {
#'         # normalized_es[i] < 0
#'         (sum(permute_scores[i,] <= normalized_es[i]) + 1) / (sum(permute_scores[i, ] < 0) +
#'                                                             1)
#'       }
#'       
#'     })
#'     
#'     
#'     p_values <- sapply(seq_along(observed_score), function(i) {
#'       if (is.na(normalized_es[i])) {
#'         NA
#'       } else if (normalized_es[i] >= 0) {
#'         (sum(permute_scores[i,] >= normalized_es[i]) + 1) / (sum(permute_scores[i, ] >= 0) + 1)
#'       } else {
#'         # normalized_es[i] < 0
#'         (sum(permute_scores[i,] <= normalized_es[i]) + 1) / (sum(permute_scores[i, ] < 0) +
#'                                                             1)
#'       }
#'       
#'     })
#'     
#'     p_adjust <- p.adjust(p_values, method = p_adjust_method)
#'     
#'     qvalues <- calculate_qvalue(pvals = p_values)
#'     
#'     gs.name <- names(metabolite_set)
#'     
#'     description <- stringr::str_split(gs.name, pattern = ";") %>%
#'       lapply(function(x) {
#'         x[2]
#'       }) %>%
#'       unlist()
#'     
#'     id <- stringr::str_split(gs.name, pattern = ";") %>%
#'       lapply(function(x) {
#'         x[1]
#'       }) %>%
#'       unlist()
#'     
#'     params <- list(
#'       pvalue_cutoff = pvalue_cutoff,
#'       perm_num = perm_num,
#'       p_adjust_method = p_adjust_method,
#'       exponent = exponent,
#'       min_size = min_size,
#'       max_size = max_size
#'     )
#'     
#'     
#'     res <- data.frame(
#'       ID = id,
#'       Description = description,
#'       setSize = sapply(metabolite_set, length),
#'       enrichmentScore = observed_score,
#'       normalized_es = normalized_es,
#'       pvalue = p_values,
#'       p_adjust = p_adjust,
#'       qvalues = qvalues,
#'       stringsAsFactors = FALSE
#'     )
#'     
#'     
#'     res <-
#'       res %>%
#'       dplyr::filter(!is.na(pvalue)) %>%
#'       dplyr::filter(pvalue <= pvalue_cutoff) %>%
#'       dplyr::filter(p_adjust <= pvalue_cutoff) %>%
#'       dplyr::arrange(pvalue)
#'     
#'     if (nrow(res) == 0) {
#'       message("no term enriched under specific pvalue_cutoff...")
#'       return(
#'         new(
#'           "metPathMSEA",
#'           result     = res,
#'           metabolite_set   = metabolite_set,
#'           feature_list   = feature_list,
#'           params     = params
#'         )
#'       )
#'     }
#'     
#'     row.names(res) <- res$ID
#'     observed_info <-
#'       observed_info[sapply(res$ID, function(x)
#'         grep(x, names(observed_info)))]
#'     
#'     if (verbose)
#'       message("leading edge analysis...")
#'     
#'     ledge <- get_leading_edge(observed_info)
#'     
#'     res$rank <- ledge$rank
#'     res$leading_edge <- ledge$leading_edge
#'     res$core_enrichment <-
#'       sapply(ledge$core_enrichment, paste0, collapse = '/')
#'     
#'     
#'     if (verbose)
#'       message("done...")
#'     
#'     feature_list <- annotation_table$condition
#'     names(feature_list) <- annotation_table$Lab.ID
#'     new(
#'       "metPathMSEA",
#'       result     = res,
#'       metabolite_set   = metabolite_set,
#'       feature_list   = feature_list,
#'       permute_scores = permute_scores,
#'       params     = params
#'     )
#'   }
