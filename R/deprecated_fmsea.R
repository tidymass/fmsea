#' ##------------------------------------------------------------------------------
#' #' @title feature MSEA
#' #' @description Install all packages in tidymass.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param from From github or gitee, if you are in China, try to set this as "gitee".
#' #' @param force Force installation, even if the remote state has not changed since the previous install.
#' #' @param upgrade One of "default", "ask", "always", or "never".
#' #' "default" respects the value of the R_REMOTES_UPGRADE environment variable if set,
#' #' and falls back to "ask" if unset. "ask" prompts the user for which out of date
#' #' packages to upgrade. For non-interactive sessions "ask" is equivalent to "always".
#' #' TRUE and FALSE are also accepted and correspond to "always" and "never" respectively.
#' #' @param dependencies Which dependencies do you want to check? Can be a character vector
#' #' (selecting from "Depends", "Imports", "LinkingTo", "Suggests", or "Enhances"),
#' #' or a logical vector.TRUE is shorthand for "Depends", "Imports", "LinkingTo"
#' #' and "Suggests". NA is shorthand for "Depends", "Imports" and "LinkingTo"
#' #' and is the default. FALSE is shorthand for no dependencies
#' #' (i.e. just check this package, not its dependencies).
#' #' The value "soft" means the same as TRUE, "hard" means the same as NA.
#' #' @param demo_data Install demo_data package or not.
#' #' @param which_package What packages you want to install? Default is all. You can set it as a character vector.
#' #' @param ... Other parameters from devtools::install_github() or devtools::install_git()
#' #' @importFrom devtools install_github
#' #' @importFrom devtools install_git
#' #' @export
#' 
#' 
#' ###
#' no_source()
#' sxtTools::setwd_project()
#' rm(list = ls())
#' source("R/annotate_isotope.R")
#' source("R/utils.R")
#' source("R/annotate_feature_table.R")
#' 
#' setwd("demo_data/denmark_project/metabolome/peaks/")
#' 
#' library(tidyverse)
#' 
#' ####load the correlation with GA
#' load("p_cor")
#' load("variable_info")
#' 
#' ##variable_info is the annotation information for
#' ##each peak according to MS2
#' feature_table <-
#'   data.frame(variable_info[, c(1:3)], p_cor[, -1], stringsAsFactors = FALSE) %>%
#'   dplyr::mutate(polarity = case_when(stringr::str_detect(name, "POS") ~ "positive",
#'                                      TRUE ~ "negative")) %>%
#'   dplyr::select(-p) %>%
#'   dplyr::rename(condition = cor)
#' 
#' head(feature_table)
#' 
#' ###get the annotation table
#' load("keggMS1database")
#' 
#' # annotation_table =
#' #   annotate_feature_table(
#' #     feature_table = feature_table,
#' #     feature_ms2_file = NULl,
#' #     column = "rp",
#' #     metabolite_database = keggMS1database,
#' #     ms1_match_ppm = 15,
#' #     peak_cluster_rt_tol = 5,
#' #     isotope_number = 3,
#' #     path = "."
#' #   )
#' # 
#' # save(annotation_table, file = "annotation_table")
#' load("annotation_table")
#' 
#' ####here, the pathway_database we must use the kegg_database
#' library(metPath)
#' load("kegg_hsa_pathway.rda")
#' pathway_database = kegg_hsa_pathway
#' min_size = 10
#' max_size = 100
#' p_value_cutoff = 0.05
#' p_adjust_method = "BH"
#' path = "."
#' 
#' fmse <-
#'   function(annotation_table,
#'            pathway_database,
#'            # organism = "hsa",
#'            exponent = 1,
#'            perm_num = 1000,
#'            min_size = 5,
#'            max_size = 1000,
#'            pvalue_cutoff = 0.05,
#'            p_adjust_method = c("holm",
#'                                "hochberg",
#'                                "hommel",
#'                                "bonferroni",
#'                                "BH",
#'                                "BY",
#'                                "fdr",
#'                                "none"),
#'            path = ".",
#'            seed = FALSE,
#'            verbose = TRUE,
#'            threads = 3,
#'            ...) {
#'     ###all the result should output to "Result" folder
#'     dir.create(file.path(path, "Result"))
#'     dir.create(file.path(path, "Result/intermediate_data"))
#'     
#'     ####---------------------------------------------------------------------------
#'     ####Metabolite set enrichment
#'     annotation_table <-
#'       annotation_table %>%
#'       dplyr::arrange(dplyr::desc(condition))
#'     
#'     purrr::map(1:length(pathway_database@compound_list), function(i) {
#'       ###Glycolysis / Gluconeogenesis pathway as an example
#'       metabolite_set = pathway_database@compound_list[[i]]
#'       result =
#'         fmse_single_set(
#'           annotation_table,
#'           metabolite_set,
#'           exponent = 1,
#'           perm_num = 1000,
#'           min_size = 5,
#'           max_size = 1000,
#'           pvalue_cutoff = 0.2,
#'           p_adjust_method = "fdr",
#'           seed = FALSE,
#'           verbose = TRUE,
#'         )
#'       
#'     })
#'     
#'     feature_table = annotation_table
#'     head(feature_table)
#'     metabolite_set = set1
#'     
#'     ###this is the weight for runing enrichment score
#'     exponent = 1
#'     perm_num = 1000
#'     min_size = 5
#'     max_size = 1000
#'     pvalue_cutoff = 0.2
#'     p_adjust_method = "fdr"
#'     seed = FALSE
#'     verbose = TRUE
#'     
#'     result <- do_msea2(
#'       feature_table = feature_table,
#'       metabolite_set = metabolite_set,
#'       exponent = 1,
#'       perm_num = 1000,
#'       min_size = 5,
#'       max_size = 1000,
#'       pvalue_cutoff = 0.25,
#'       p_adjust_method = "fdr",
#'       seed = TRUE,
#'       verbose = TRUE
#'     )
#'     
#'     msea_plot(x = result, title = result@result$Description[1])
#'     
#'     msea_plot(x = head_result[[idx]],
#'               title = paste(head_result[[idx]]@result$p_adjust))
#'     
#'     
#'     msea_plot(x = tail_result[[idx]],
#'               title = paste(tail_result[[idx]]@result$p_adjust))
#'     
#'     
#'     
#'     
#'   }
