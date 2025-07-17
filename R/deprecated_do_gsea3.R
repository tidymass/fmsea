#' do_msea3 <- function(feature_table,
#'                      metabolite_set,
#'                      exponent = 1,
#'                      perm_num = 1000,
#'                      min_size = 5,
#'                      max_size = 1000,
#'                      pvalue_cutoff = 0.2,
#'                      p_adjust_method = c("fdr",
#'                                          "holm",
#'                                          "hochberg",
#'                                          "hommel",
#'                                          "bonferroni",
#'                                          "BH",
#'                                          "BY",
#'                                          "none"),
#'                      seed = FALSE,
#'                      verbose = TRUE,
#'                      threads = 3,
#'                      ...) {
#'   browser()
#'   p_adjust_method = match.arg(p_adjust_method)
#'   ##filter pathways according to pathway size
#'   if (length(metabolite_set) < min_size |
#'       length(metabolite_set) > max_size) {
#'     return(NULL)
#'   }
#'   
#'   ##find where are the metabolites in the feature table
#'   ##for each metabolite in metabolite set, its position in feature table
#'   match_idx <-
#'     purrr::map(metabolite_set, function(x) {
#'       which(x == feature_table$Lab.ID)
#'     })
#'   
#'   ##for the metabolites in the feature set,
#'   ##all the matched peaks should be remained, other metabolites
#'   ##only remain one matched peak,
#'   ##either one is fine
#'   ###idx is the index of all metabolites in metabolite set in feature_table
#'   idx <-
#'     sort(unique(unlist(match_idx)))
#'   
#'   ##for the metabolites in feature sets
#'   feature_table1 <-
#'     feature_table[idx,]
#'   
#'   #for each metabolite_feature_cluster, if the peaks are in different condition_class
#'   #(for example, if the condition is fold change, some peaks fold change > 1 and some
#'   #peaks fold change < 1), we should think some method to handle this condition
#'   
#'   ##for the metabolites are not in the feature set, they only affect the rank of
#'   ##the metabolites, so only remain one time for each peak.
#'   feature_table2 <-
#'     feature_table[-idx,] %>%
#'     dplyr::distinct(name, .keep_all = TRUE)
#'   
#'   
#'   ###faeture_table1 and feature_table2 have a lot of overlapped features
#'   ##we should remainthe overlapped features in feature table2, because the annotaiton is not confirmed
#'   # intersect(feature_table1$name, feature_table2$name)
#'   
#'   feature_table <-
#'     rbind(feature_table1,
#'           feature_table2) %>%
#'     dplyr::arrange(dplyr::desc(condition))
#'   
#'   if (verbose) {
#'     message("calculating observed enrichment scores...")
#'   }
#'   
#'   ###observed_info is the running enrichment score for each metabolite (peak)
#'   ###in feature_table. And ES is the enrichment score
#'   observed_info = get_msea_score2(
#'     feature_table = feature_table,
#'     metabolite_set = metabolite_set,
#'     exponent = exponent
#'   )
#'   
#'   observed_score <- observed_info$ES
#'   
#'   if (verbose) {
#'     message("calculating permutation scores...")
#'   }
#'   
#'   ###here, we need to do the permutation test
#'   if (seed) {
#'     seeds <- sample.int(perm_num)
#'   }
#'   
#'   if (masstools::get_os() == "windows") {
#'     bpparam =
#'       BiocParallel::SnowParam(workers = threads,
#'                               progressbar = TRUE)
#'   } else{
#'     bpparam = BiocParallel::MulticoreParam(workers = threads,
#'                                            progressbar = TRUE)
#'   }
#'   
#'   perm_scores <-
#'     BiocParallel::bplapply(1:perm_num, function(i) {
#'       if (seed) {
#'         set.seed(seeds[i])
#'       }
#'       
#'       perm_msea_es2(
#'         feature_table = feature_table,
#'         metabolite_set = metabolite_set,
#'         exponent = exponent
#'       )
#'     }, BPPARAM = bpparam) %>%
#'     unlist()
#'   
#'   ####mean_pos_null_es is the mean value of all the positive null enrichment score
#'   mean_pos_null_es <- mean(perm_scores[perm_scores >= 0])
#'   ####mean_neg_null_es is the mean value of all the negative null enrichment score
#'   mean_neg_null_es <- mean(perm_scores[perm_scores < 0])
#'   
#'   ###normalize_es is the function used to normalize ES
#'   ##ES is the enrichment score
#'   normalize_es <- function(ES, mean_pos_null_es, mean_neg_null_es) {
#'     s <- sign(ES)
#'     m <- numeric(length(ES))
#'     m[s == 1] <- mean_pos_null_es[s == 1]
#'     m[s == -1] <- mean_neg_null_es[s == -1]
#'     ES / m
#'   }
#'   
#'   ###normalized_es is the normalized ES
#'   normalized_es <-
#'     normalize_es(observed_score, mean_pos_null_es, mean_neg_null_es)
#'   
#'   perm_scores = normalize_es(
#'     ES = observed_score,
#'     mean_pos_null_es = mean_pos_null_es,
#'     mean_neg_null_es = mean_neg_null_es
#'   )
#'   
#'   if (verbose) {
#'     message("calculating p values...")
#'   }
#'   
#'   ###get the p value of ES
#'   if (is.na(normalized_es)) {
#'     p_value = NA
#'   }
#'   
#'   if (normalized_es >= 0) {
#'     (sum(perm_scores[i, ] >= normalized_es[i]) + 1) / (sum(perm_scores[i,] >= 0) + 1)
#'   }
#'   
#'   p_values <- sapply(seq_along(observed_score), function(i) {
#'     if (is.na(normalized_es[i])) {
#'       NA
#'     } else if (normalized_es[i] >= 0) {
#'       (sum(perm_scores[i, ] >= normalized_es[i]) + 1) / (sum(perm_scores[i,] >= 0) + 1)
#'     } else {
#'       # normalized_es[i] < 0
#'       (sum(perm_scores[i, ] <= normalized_es[i]) + 1) / (sum(perm_scores[i,] < 0) +
#'                                                            1)
#'     }
#'     
#'   })
#'   
#'   
#'   p_values <- sapply(seq_along(observed_score), function(i) {
#'     if (is.na(normalized_es[i])) {
#'       NA
#'     } else if (normalized_es[i] >= 0) {
#'       (sum(perm_scores[i, ] >= normalized_es[i]) + 1) / (sum(perm_scores[i,] >= 0) + 1)
#'     } else {
#'       # normalized_es[i] < 0
#'       (sum(perm_scores[i, ] <= normalized_es[i]) + 1) / (sum(perm_scores[i,] < 0) +
#'                                                            1)
#'     }
#'     
#'   })
#'   
#'   p_adjust <- p.adjust(p_values, method = p_adjust_method)
#'   
#'   qvalues <- calculate_qvalue(pvals = p_values)
#'   
#'   gs.name <- names(metabolite_set)
#'   
#'   description <- stringr::str_split(gs.name, pattern = ";") %>%
#'     lapply(function(x) {
#'       x[2]
#'     }) %>%
#'     unlist()
#'   
#'   id <- stringr::str_split(gs.name, pattern = ";") %>%
#'     lapply(function(x) {
#'       x[1]
#'     }) %>%
#'     unlist()
#'   
#'   params <- list(
#'     pvalue_cutoff = pvalue_cutoff,
#'     perm_num = perm_num,
#'     p_adjust_method = p_adjust_method,
#'     exponent = exponent,
#'     min_size = min_size,
#'     max_size = max_size
#'   )
#'   
#'   
#'   res <- data.frame(
#'     ID = id,
#'     Description = description,
#'     setSize = sapply(metabolite_set, length),
#'     enrichmentScore = observed_score,
#'     normalized_es = normalized_es,
#'     pvalue = p_values,
#'     p_adjust = p_adjust,
#'     qvalues = qvalues,
#'     stringsAsFactors = FALSE
#'   )
#'   
#'   
#'   res <-
#'     res %>%
#'     dplyr::filter(!is.na(pvalue)) %>%
#'     dplyr::filter(pvalue <= pvalue_cutoff) %>%
#'     dplyr::filter(p_adjust <= pvalue_cutoff) %>%
#'     dplyr::arrange(pvalue)
#'   
#'   if (nrow(res) == 0) {
#'     message("no term enriched under specific pvalue_cutoff...")
#'     return(
#'       new(
#'         "metPathMSEA",
#'         result     = res,
#'         metabolite_set   = metabolite_set,
#'         feature_list   = feature_list,
#'         params     = params
#'       )
#'     )
#'   }
#'   
#'   row.names(res) <- res$ID
#'   observed_info <-
#'     observed_info[sapply(res$ID, function(x)
#'       grep(x, names(observed_info)))]
#'   
#'   if (verbose)
#'     message("leading edge analysis...")
#'   
#'   ledge <- get_leading_edge(observed_info)
#'   
#'   res$rank <- ledge$rank
#'   res$leading_edge <- ledge$leading_edge
#'   res$core_enrichment <-
#'     sapply(ledge$core_enrichment, paste0, collapse = '/')
#'   
#'   
#'   if (verbose)
#'     message("done...")
#'   
#'   feature_list <- feature_table$condition
#'   names(feature_list) <- feature_table$Lab.ID
#'   new(
#'     "metPathMSEA",
#'     result     = res,
#'     metabolite_set   = metabolite_set,
#'     feature_list   = feature_list,
#'     perm_scores = perm_scores,
#'     params     = params
#'   )
#' }
#' 
#' 
#' 
#' ##' @name metPathMSEA-class
#' ##' @aliases metPathMSEA-class
#' ##'   show,metPathMSEA-method summary,metPathMSEA-method
#' ##'
#' ##' @docType class
#' ##' @slot result MSEA anaysis result.
#' ##' @slot metaboliteSets metaboliteSets
#' ##' @slot feature_list order rank feature_list
#' ##' @slot permScores permutation scores
#' ##' @slot params parameters
#' ##' @exportClass metPathMSEA
#' ##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
#' ##' @keywords classes
#' setClass(
#'   "metPathMSEA",
#'   representation   = representation(
#'     result          = "data.frame",
#'     metabolite_set        = "list",
#'     feature_list        = "numeric",
#'     perm_scores      = "matrix",
#'     params          = "list"
#'   )
#' )
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ###----------------------------------------------------------------------------
#' # metabolite_set <- kegg_hsa_compound_pathway
#' # feature_list <- feature_list
#' 
#' filter_metabolite_set2 <- function(metabolite_set,
#'                                    feature_table,
#'                                    min_size = 5,
#'                                    max_size = 1000) {
#'   len <- lapply(metabolite_set, function(x) {
#'     length(intersect(x, feature_table$Lab.ID))
#'   }) %>%
#'     unlist() %>%
#'     unname()
#'   
#'   remain_idx <-
#'     which(len > min_size & len < max_size)
#'   
#'   metabolite_set <- metabolite_set[remain_idx]
#'   return(metabolite_set)
#' }
#' 
#' 
#' 
#' 
#' 
#' ###get get the msea score for one metabolite set
#' get_msea_score2 <- function(feature_table,
#'                             metabolite_set,
#'                             exponent = 1,
#'                             fortify = FALSE) {
#'   ###################################################################
#'   ##    feature_table                                              ##
#'   ##                                                               ##
#'   ## 1. Rank order the N metabolites (peaks) in D                  ##
#'   ## to form L = { m_1, ... , m_N}                                 ##
#'   ##    according to the condition (fold change, correlation       ##
#'   ##    and son on),
#'   ##    r(g_j)=r_j,                  ##
#'   ##    of their expression profiles with C.                       ##
#'   ##                                                               ##
#'   ###################################################################
#'   
#'   ###################################################################
#'   ##    exponent                                                   ##
#'   ##                                                               ##
#'   ## An exponent p to control the weight of the step.              ##
#'   ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
#'   ##   the standard Kolmogorov-Smirnov statistic.                  ##
#'   ##   When p = 1, we are weighting the metabolites in S           ##
#'   ##   by their correlation with C normalized                      ##
#'   ##   by the sum of the correlations over all of the              ##
#'   ##   metabolites in S.                                           ##
#'   ##                                                               ##
#'   ###################################################################
#'   
#'   ## metabolites defined in metabolite_set should appear in feature_table.
#'   metabolite_set <- intersect(metabolite_set, feature_table$Lab.ID)
#'   
#'   ##num_feature_list is the length of feature_table
#'   num_feature_list <- nrow(feature_table)
#'   
#'   ##num_metabolite_list is the length of unique metabolites in feature_table
#'   num_metabolite_list <- length(unique(feature_table$Lab.ID))
#'   
#'   ##num_metabolite_set is the length of feature set
#'   num_metabolite_set <- length(metabolite_set)
#'   
#'   ###p_hit and p_miss are the vector to record if the metabolites in feature_table
#'   ###are in metabolite_set or not
#'   p_hit <- p_miss <- numeric(num_feature_list)
#'   
#'   ###hits is a logical vector indicates if the metabolites ID
#'   ### in feature set or not
#'   hits <- feature_table$Lab.ID %in% metabolite_set ## logical
#'   
#'   p_hit[hits] <-
#'     (abs(feature_table$condition[hits])) ^ exponent
#'   #   score <- feature_table$score[hits]
#'   #   metabolite_feature_cluster <- feature_table$metabolite_feature_cluster[hits]
#'   #
#'   #   temp_data <- data.frame(
#'   #     index = 1:sum(hits),
#'   #     p_hit = p_hit[hits],
#'   #     feature_table[hits, ],
#'   #     stringsAsFactors = FALSE
#'   #   )
#'   #
#'   # idx <- 17
#'   #   temp_data %>%
#'   #     ggplot(aes(index, p_hit)) +
#'   #     geom_point() +
#'   #     geom_point(
#'   #       aes(x = index, y = p_hit, color = metabolite_feature_cluster),
#'   #       show.legend = FALSE,
#'   #       data = temp_data %>% filter(Lab.ID == unique(Lab.ID)[idx])
#'   #     ) +
#'   #     ggrepel::geom_label_repel(
#'   #       aes(
#'   #         x = index,
#'   #         y = p_hit,
#'   #         label = paste(metabolite_feature_cluster, score, Adduct, isotope, sep = ":"),
#'   #         color = metabolite_feature_cluster
#'   #       ),
#'   #       show.legend = FALSE,
#'   #       data = temp_data %>% filter(Lab.ID == unique(Lab.ID)[idx])
#'   #     ) +
#'   #     ggsci::scale_color_aaas() +
#'   #     theme_bw()
#'   
#'   # write.csv(feature_table[hits,], "test.csv", row.names = FALSE)
#'   
#'   ###NR is the normalized runing enrichment score
#'   NR <- sum(p_hit)
#'   p_hit <- cumsum(p_hit / NR)
#'   
#'   p_miss[!hits] <- 1 / (num_feature_list - num_metabolite_set)
#'   p_miss <- cumsum(p_miss)
#'   
#'   runing_es <- p_hit - p_miss
#'   
#'   # data.frame(index = 1:length(p_hit),
#'   #            p_hit, p_miss, runing_es,
#'   #            stringsAsFactors = FALSE) %>%
#'   #   tidyr::pivot_longer(cols = -index, names_to = "class", values_to = "value") %>%
#'   #   ggplot(aes(x = index, y = value, color = class)) +
#'   #   geom_point()
#'   
#'   ## enrichment (ES) is the maximum deviation from zero of p_hit-p_miss
#'   max_es <- max(runing_es)
#'   min_es <- min(runing_es)
#'   
#'   if (abs(max_es) > abs(min_es)) {
#'     ES <- max_es
#'   } else {
#'     ES <- min_es
#'   }
#'   
#'   ####in the df data frame,x is the index of all the metabolites (peaks) in the
#'   ####feature table. running_score is the running score for each metabolites (peaks)
#'   ####in feature table. position indicates if the metabolites or peaks are in
#'   ####metabolite set or not.
#'   df <- data.frame(
#'     x = seq_along(runing_es),
#'     running_score = runing_es,
#'     position = as.integer(hits)
#'   )
#'   
#'   # df %>%
#'   #   ggplot(aes(x, running_score)) +
#'   #   geom_point(aes(color = as.character(position)))
#'   
#'   ###fortify (加强,增强)
#'   if (fortify == TRUE) {
#'     return(df)
#'   }
#'   
#'   ###if set fortify as TRUE, add the metabolites ID to each metabolite (peak)
#'   df$metabolite = feature_table$Lab.ID
#'   res <- list(ES = ES, runing_es = df)
#'   return(res)
#' }
#' 
#' 
#' 
#' 
#' #' @title perm_msea_es2
#' #' @description Permutation to get the null distribution of enrichment score
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param feature_table Feature table ordered by condition.
#' #' @param metabolite_set Metabolite set.
#' 
#' perm_msea_es2 <- function(feature_table,
#'                           metabolite_set,
#'                           exponent = 1) {
#'   permutated_feature_table <- perm_feature_table(feature_table)
#'   res =
#'     get_msea_score2(
#'       feature_table = permutated_feature_table,
#'       metabolite_set = metabolite_set,
#'       exponent = exponent,
#'       fortify = FALSE
#'     )
#'   return(res$ES)
#' }
#' 
#' #' @title perm_feature_table
#' #' @description This function is used randomly order the feature table. We just randomly order the
#' #' Lab.ID for feature_table.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param feature_table Feature table ordered by condition.
#' 
#' perm_feature_table <- function(feature_table) {
#'   perm.idx <- sample.int(nrow(feature_table))
#'   permutated_feature_table <- feature_table
#'   permutated_feature_table$condition <-
#'     permutated_feature_table$condition[perm.idx]
#'   permutated_feature_table =
#'     permutated_feature_table %>%
#'     dplyr::arrange(dplyr::desc(condition))
#'   return(permutated_feature_table)
#' }
#' 
#' 
#' 
#' 
#' 
#' calculate_qvalue <- function(pvals) {
#'   if (length(pvals) == 0)
#'     return(numeric(0))
#'   
#'   qobj <- tryCatch(
#'     qvalue::qvalue(pvals, lambda = 0.05,
#'                    pi0.method = "bootstrap"),
#'     error = function(e)
#'       NULL
#'   )
#'   
#'   if (class(qobj) == "qvalue") {
#'     qvalues <- qobj$qvalues
#'   } else {
#'     qvalues <- NA
#'   }
#'   return(qvalues)
#' }
#' 
#' 
#' ###----------------------------------------------------------------------------
#' get_leading_edge <- function(observed_info) {
#'   core_enrichment <- lapply(observed_info, function(x) {
#'     runing_es <- x$runing_es
#'     runing_es <- runing_es[runing_es$position == 1, ]
#'     ES <- x$ES
#'     if (ES >= 0) {
#'       i <- which.max(runing_es$running_score)
#'       leading_metabolite <- runing_es$metabolite[1:i]
#'     } else {
#'       i <- which.min(runing_es$running_score)
#'       leading_metabolite <- runing_es$metabolite[-c(1:(i - 1))]
#'     }
#'     return(leading_metabolite)
#'   })
#'   
#'   rank <- sapply(observed_info, function(x) {
#'     runing_es <- x$runing_es
#'     ES <- x$ES
#'     if (ES >= 0) {
#'       rr <- which.max(runing_es$running_score)
#'     } else {
#'       i <- which.min(runing_es$running_score)
#'       rr <- nrow(runing_es) - i + 1
#'     }
#'     return(rr)
#'   })
#'   
#'   tags <- sapply(observed_info, function(x) {
#'     runing_es <- x$runing_es
#'     runing_es <- runing_es[runing_es$position == 1, ]
#'     ES <- x$ES
#'     if (ES >= 0) {
#'       i <- which.max(runing_es$running_score)
#'       res <- i / nrow(runing_es)
#'     } else {
#'       i <- which.min(runing_es$running_score)
#'       res <- (nrow(runing_es) - i + 1) / nrow(runing_es)
#'     }
#'     return(res)
#'   })
#'   
#'   ll <- sapply(observed_info, function(x) {
#'     runing_es <- x$runing_es
#'     ES <- x$ES
#'     if (ES >= 0) {
#'       i <- which.max(runing_es$running_score)
#'       res <- i / nrow(runing_es)
#'     } else {
#'       i <- which.min(runing_es$running_score)
#'       res <- (nrow(runing_es) - i + 1) / nrow(runing_es)
#'     }
#'     return(res)
#'   })
#'   
#'   N <- nrow(observed_info[[1]]$runing_es)
#'   setSize <-
#'     sapply(observed_info, function(x)
#'       sum(x$runing_es$position))
#'   signal <- tags * (1 - ll) * (N / (N - setSize))
#'   
#'   tags <- paste0(round(tags * 100), "%")
#'   ll <- paste0(round(ll * 100), "%")
#'   signal <- paste0(round(signal * 100), "%")
#'   leading_edge <-
#'     paste0('tags=', tags, ", list=", ll, ", signal=", signal)
#'   
#'   res <- list(
#'     rank = rank,
#'     tags = tags,
#'     list = ll,
#'     signal = signal,
#'     leading_edge = leading_edge,
#'     core_enrichment = core_enrichment
#'   )
#'   return(res)
#' }
#' 
#' 
#' 
#' ##' show method for \code{metPathMSEA} instance
#' ##'
#' ##' @name show
#' ##' @docType methods
#' ##' @rdname show-methods
#' ##'
#' ##' @title show method
#' ##' @return message
#' ##' @importFrom methods show
#' ##' @exportMethod show
#' ##' @usage show(object)
#' ##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
#' setMethod("show", signature(object = "metPathMSEA"),
#'           function (object) {
#'             params <- object@params
#'             cat("#\n# Gene Set Enrichment Analysis\n#\n")
#'             cat("#...@feature_list", "\t")
#'             str(object@feature_list)
#'             cat("#...nPerm", "\t", params$perm_num, "\n")
#'             cat(
#'               "#...pvalues adjusted by",
#'               paste0("'", params$p_adjust_method, "'"),
#'               paste0("with cutoff <", params$pvalueCutoff),
#'               "\n"
#'             )
#'             cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
#'             str(object@result)
#'           })
#' 
#' 
#' 
#' ##' visualize analyzing result of MSEA
#' ##'
#' ##' plotting function for mseaResult
#' ##' @title msea_plot
#' ##' @rdname mseaplot
#' ##' @param x object of msea result
#' ##' @param metaboliteSetID metaboliteSet ID
#' ##' @param by one of "running_score" or "position"
#' ##' @param title plot title
#' ##' @param ... additional parameters
#' ##' @return ggplot2 object
#' ##' @export
#' ##' @examples
#' ##' library(DOSE)
#' ##' data(metaboliteList)
#' ##' x <- gseDO(metaboliteList)
#' ##' mseaplot(x, metaboliteSetID=1)
#' setGeneric(name = "msea_plot",
#'            function(x,
#'                     metabolite_set_idx = 1,
#'                     by = "all",
#'                     title = "",
#'                     color = 'black',
#'                     color.line = "#8DD3C7",
#'                     color.vline = "#FB8072",
#'                     ...) {
#'              standardGeneric("msea_plot")
#'            })
#' 
#' 
#' ##' @rdname msea_plot
#' ##' @exportMethod msea_plot
#' setMethod(f = "msea_plot", signature(x = "metPathMSEA"),
#'           function (x,
#'                     metabolite_set_idx = 1,
#'                     by = "all",
#'                     title = "",
#'                     color = 'black',
#'                     color.line = "#8DD3C7",
#'                     color.vline = "#FB8072",
#'                     ...) {
#'             msea_plot.metPathMSEA(
#'               x,
#'               metabolite_set_idx = metabolite_set_idx,
#'               by = by,
#'               title = title,
#'               color = color,
#'               color.line = color.line,
#'               color.vline = color.vline,
#'               ...
#'             )
#'           })
#' 
#' ##' @rdname msea_plot
#' ##' @param color color of line segments
#' ##' @param color.line color of running enrichment score line
#' ##' @param color.vline color of vertical line which indicating the maximum/minimal running enrichment score
#' ##' @return ggplot2 object
#' ##' @importFrom ggplot2 ggplot
#' ##' @importFrom ggplot2 geom_linerange
#' ##' @importFrom ggplot2 geom_line
#' ##' @importFrom ggplot2 geom_vline
#' ##' @importFrom ggplot2 geom_hline
#' ##' @importFrom ggplot2 xlab
#' ##' @importFrom ggplot2 ylab
#' ##' @importFrom ggplot2 xlim
#' ##' @importFrom ggplot2 aes
#' ##' @importFrom ggplot2 ggplotGrob
#' ##' @importFrom ggplot2 geom_segment
#' ##' @importFrom ggplot2 ggplot_gtable
#' ##' @importFrom ggplot2 ggplot_build
#' ##' @importFrom ggplot2 ggtitle
#' ##' @importFrom ggplot2 element_text
#' ##' @importFrom ggplot2 rel
#' ##' @importFrom cowplot plot_grid
#' ##' @author Guangchuang Yu
#' 
#' msea_plot.metPathMSEA <-
#'   function (x,
#'             metabolite_set_idx = 1,
#'             by = "all",
#'             title = "",
#'             color = 'black',
#'             color.line = "green",
#'             color.vline = "#FA5860",
#'             ...) {
#'     if (is.null(x)) {
#'       return(NULL)
#'     }
#'     by <- match.arg(by, c("running_score", "preranked", "all"))
#'     gs_data <- get_gs_info(x, metabolite_set_idx)
#'     
#'     p <- ggplot(gs_data, aes_(x = ~ x)) +
#'       theme_dose() +
#'       xlab("Position in the ranked list of features")
#'     
#'     if (by == "running_score" || by == "all") {
#'       p.res <-
#'         p + geom_linerange(aes_(ymin =  ~ ymin, ymax =  ~ ymax), color = color)
#'       p.res <-
#'         p.res + geom_line(aes_(y = ~ running_score),
#'                           color = color.line,
#'                           size =
#'                             1)
#'       enrichmentScore <-
#'         x@result[metabolite_set_idx, "enrichmentScore"]
#'       
#'       es.df <-
#'         data.frame(es = which.min(abs(
#'           p$data$running_score - enrichmentScore
#'         )))
#'       
#'       p.res <-
#'         p.res + geom_vline(
#'           data = es.df,
#'           aes_(xintercept = ~ es),
#'           colour = color.vline,
#'           linetype = "dashed"
#'         )
#'       p.res <- p.res + ylab("Running enrichment score")
#'       p.res <- p.res + geom_hline(yintercept = 0)
#'     }
#'     
#'     if (by == "preranked" || by == "all") {
#'       df2 <- data.frame(x = which(p$data$position == 1))
#'       df2$y <- p$data$feature_list[df2$x]
#'       p.pos <-
#'         p + geom_segment(data = df2,
#'                          aes_(
#'                            x =  ~ x,
#'                            xend =  ~ x,
#'                            y =  ~ y,
#'                            yend = 0
#'                          ),
#'                          color = color)
#'       p.pos <-
#'         p.pos + ylab("Ranked list metric") + xlim(0, length(p$data$feature_list))
#'     }
#'     if (by == "running_score")
#'       return(p.res + ggtitle(title))
#'     if (by == "preranked")
#'       return(p.pos + ggtitle(title))
#'     
#'     p.pos <-
#'       p.pos + xlab(NULL) + theme(axis.text.x = element_blank(),
#'                                  axis.ticks.x = element_blank())
#'     p.pos <- p.pos + ggtitle(title) +
#'       theme(plot.title = element_text(hjust = 0.5, size = rel(2)))
#'     cowplot::plot_grid(p.pos, p.res, ncol = 1, align = "v")
#'   }
#' 
#' 
#' 
#' get_gs_info <- function(object,
#'                         metabolite_set_idx = 1) {
#'   feature_list <- object@feature_list
#'   # if (is.numeric(metabolite_set_idx))
#'   #   metabolite_set_idx <- object@result[metabolite_set_idx, "ID"]
#'   metabolite_set <- object@metabolite_set[[metabolite_set_idx]]
#'   exponent <- object@params[["exponent"]]
#'   df <-
#'     get_msea_score(feature_list, metabolite_set, exponent, fortify = TRUE)
#'   df$ymin = 0
#'   df$ymax = 0
#'   pos <- df$position == 1
#'   h <- diff(range(df$running_score)) / 20
#'   df$ymin[pos] <- -h
#'   df$ymax[pos] <- h
#'   df$feature_list <- feature_list
#'   return(df)
#' }
#' 
#' 
#' ##' ggplot theme of DOSE
#' ##'
#' ##' @title theme_dose
#' ##' @param font.size font size
#' ##' @return ggplot theme
#' ##' @importFrom ggplot2 theme_bw
#' ##' @importFrom ggplot2 theme
#' ##' @importFrom ggplot2 element_text
#' ##' @importFrom ggplot2 margin
#' ##' @examples
#' ##' library(ggplot2)
#' ##' qplot(1:10) + theme_dose()
#' ##' @export
#' theme_dose <- function(font.size = 13) {
#'   theme_bw() +
#'     theme(
#'       axis.text.x = element_text(
#'         colour = "black",
#'         size = 12,
#'         vjust = 1
#'       ),
#'       axis.text.y = element_text(
#'         colour = "black",
#'         size = 12,
#'         hjust = 1
#'       ),
#'       axis.title = element_text(
#'         margin = margin(10, 5, 0, 0),
#'         color = "black",
#'         size = 13
#'       ),
#'       axis.title.y = element_text(angle = 90)
#'     )
#' }
