##------------------------------------------------------------------------------
#' @title annotate_feature_table
#' @description annotate_feature_table
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param from From github or gitee, if you are in China, try to set this as "gitee".
#' @param force Force installation, even if the remote state has not changed since the previous install.
#' @param upgrade One of "default", "ask", "always", or "never".
#' "default" respects the value of the R_REMOTES_UPGRADE environment variable if set,
#' and falls back to "ask" if unset. "ask" prompts the user for which out of date
#' packages to upgrade. For non-interactive sessions "ask" is equivalent to "always".
#' TRUE and FALSE are also accepted and correspond to "always" and "never" respectively.
#' @param dependencies Which dependencies do you want to check? Can be a character vector
#' (selecting from "Depends", "Imports", "LinkingTo", "Suggests", or "Enhances"),
#' or a logical vector.TRUE is shorthand for "Depends", "Imports", "LinkingTo"
#' and "Suggests". NA is shorthand for "Depends", "Imports" and "LinkingTo"
#' and is the default. FALSE is shorthand for no dependencies
#' (i.e. just check this package, not its dependencies).
#' The value "soft" means the same as TRUE, "hard" means the same as NA.
#' @param demo_data Install demo_data package or not.
#' @param which_package What packages you want to install? Default is all. You can set it as a character vector.
#' @param ... Other parameters from devtools::install_github() or devtools::install_git()
#' @importFrom devtools install_github
#' @importFrom devtools install_git
#' @export


# ###
# no_source()
# sxtTools::setwd_project()
# rm(list = ls())
# source("R/annotate_isotope.R")
# source("R/utils.R")
#
# setwd("demo_data/denmark_project/metabolome/peaks/")
#
# library(tidyverse)
#
# ####load the correlation with GA
# load("p_cor")
# load("variable_info")
#
# ##variable_info is the annotation information for
# ##each peak according to MS2
# feature_table <-
#   data.frame(variable_info[, c(1:3)], p_cor[, -1], stringsAsFactors = FALSE) %>%
#   dplyr::mutate(polarity = case_when(stringr::str_detect(name, "POS") ~ "positive",
#                                      TRUE ~ "negative")) %>%
#   dplyr::select(-p) %>%
#   dplyr::rename(condition = cor)
#
# head(feature_table)
#
# load("keggMS1database")
# metabolite_database = keggMS1database

annotate_feature_table <-
  function(feature_table,
           feature_ms2_file,
           column = c("rp", "hilic"),
           metabolite_database,
           ms1_match_ppm = 15,
           peak_cluster_rt_tol = 5,
           isotope_number = 3,
           path = ".",
           ...) {
    column = match.arg(column)
    ###all the result should output to "Result" folder
    dir.create(file.path(path, "Result"))
    dir.create(file.path(path, "Result/intermediate_data"))
    
    ###should check feature_table here
    
    ms1_data_pos <-
      feature_table %>%
      dplyr::filter(polarity == "positive")
    
    ms1_data_neg <-
      feature_table %>%
      dplyr::filter(polarity == "negative")
    
    save(
      metabolite_database,
      file = file.path(path, "Result/intermediate_data/metabolite_database")
    )
    
    ############metabolite database matching for feature_table
    cat("Database matching for positive mode...\n")
    
    if (nrow(ms1_data_pos) == 0) {
      annotation_table_pos <- NULL
    } else{
      write.csv(
        ms1_data_pos,
        file = file.path(path, "Result/intermediate_data/ms1_data_pos.csv"),
        row.names = FALSE
      )
      
      annotation_table_pos <-
        metID::identify_metabolites(
          ms1.data = "ms1_data_pos.csv",
          ms1.match.ppm = ms1_match_ppm,
          column = column,
          polarity = "positive",
          database = "metabolite_database",
          candidate.num = 1000,
          path = file.path(path, "Result/intermediate_data")
        )
      annotation_table_pos <-
        metID::get_identification_table(annotation_table_pos,
                                        candidate.num = 10000,
                                        type = "new")
      annotation_table_pos$polarity <- "positive"
    }
    
    cat("Database matching for negative mode...\n")
    
    if (nrow(ms1_data_neg) == 0) {
      annotation_table_neg <- NULL
    } else{
      write.csv(
        ms1_data_neg,
        file = file.path(path, "Result/intermediate_data/ms1_data_neg.csv"),
        row.names = FALSE
      )
      
      annotation_table_neg <-
        metID::identify_metabolites(
          ms1.data = "ms1_data_neg.csv",
          ms1.match.ppm = ms1_match_ppm,
          column = column,
          polarity = "negative",
          database = "metabolite_database",
          candidate.num = 1000,
          path = file.path(path, "Result/intermediate_data")
        )
      annotation_table_neg <-
        metID::get_identification_table(annotation_table_neg,
                                        candidate.num = 10000,
                                        type = "new")
      
      annotation_table_neg$polarity <- "negative"
      
    }
    
    annotation_table =
      rbind(annotation_table_pos, annotation_table_neg) %>%
      as.data.frame()
    
    ###this should be remove when packged
    library(plyr)
    
    ####make sure that there is not NA in the annotation table
    annotation_table =
      annotation_table %>%
      plyr::dlply(.variables = .(name)) %>%
      purrr::map(function(x) {
        x$mz = x$mz[1]
        x$rt = x$rt[1]
        x$condition = stringr::str_trim(x$condition[1], side = "both") %>% as.numeric()
        x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    ###add Compound formula to annotation table
    annotation_table <-
      annotation_table %>%
      dplyr::arrange(desc(condition)) %>%
      dplyr::left_join(metabolite_database@spectra.info[, c("KEGG.ID", "Formula")], by = "KEGG.ID") %>%
      dplyr::filter(!is.na(KEGG.ID)) %>%
      dplyr::select(-c(MS2.spectra.name, Candidate.number, CE, SS, Total.score)) %>%
      dplyr::mutate(isotope = "[M]")
    
    ##------------------------------------------------------------------------------
    ####isotope annotation
    annotation_table_pos =
      annotation_table %>%
      dplyr::filter(polarity == "positive")
    
    annotation_table_neg =
      annotation_table %>%
      dplyr::filter(polarity == "negative")
    
    ###positive
    library(future)
    library(furrr)
    future::plan(multisession, workers = 8)
    
    cat("Isotope annotation for feature table positive mode...\n")
    
    system.time(
      isotope_pos <-
        furrr::future_map(
          .x = as.data.frame(t(annotation_table_pos)),
          .f = function(x) {
            adduct <-
              stringr::str_extract(x[11], "\\(.+\\)") %>%
              stringr::str_replace("\\(", "") %>%
              stringr::str_replace("\\)", "")
            
            mz = as.numeric(stringr::str_trim(x[2], side = "both"))
            rt = as.numeric(stringr::str_trim(x[3], side = "both"))
            
            temp_iso <- try(annotate_isotope(
              formula = stringr::str_trim(x[17], side = "both"),
              adduct = adduct,
              mz = mz,
              rt = rt,
              peak.mz = ms1_data_pos$mz,
              peak.rt = ms1_data_pos$rt,
              rt.tol = peak_cluster_rt_tol,
              mz.tol = ms1.match.ppm,
              max.isotope = 3
            ),
            silent = TRUE)
            
            if (class(temp_iso) == "try-error") {
              return(NULL)
            }
            
            if (is.null(temp_iso)) {
              return(NULL)
            }
            
            temp_iso <-
              cbind(ms1_data_pos[temp_iso$peakIndex, ], temp_iso) %>%
              dplyr::select(-c(peakIndex))
            
            colnames(temp_iso) <- c("name",
                                    "mz",
                                    "rt",
                                    "condition",
                                    "polarity",
                                    "mz.error",
                                    "isotope",
                                    "RT.error")
            
            x <- matrix(x, nrow = 1) %>% as.data.frame()
            colnames(x) <- colnames(annotation_table_pos)
            
            temp_iso$Compound.name <- x$Compound.name
            temp_iso$Lab.ID <- x$Lab.ID
            temp_iso$Adduct <- x$Adduct
            temp_iso$Formula <- x$Formula
            temp_iso$rt <- x$rt
            temp_iso$CAS.ID <- x$CAS.ID
            temp_iso$HMDB.ID <- x$HMDB.ID
            temp_iso$KEGG.ID <- x$KEGG.ID
            temp_iso$Database <- x$Database
            
            temp_iso$mz.match.score <-
              (ms1.match.ppm - temp_iso$mz.error) / ms1.match.ppm
            temp_iso$RT.match.score <-
              (peak_cluster_rt_tol - temp_iso$RT.error) / peak_cluster_rt_tol
            temp_iso %>%
              dplyr::select(colnames(x))
          },
          .progress = TRUE
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
    )
    
    annotation_table_pos =
      rbind(annotation_table_pos, isotope_pos) %>%
      as.data.frame() %>%
      dplyr::arrange(Compound.name, isotope, rt)
    
    ###negative
    library(future)
    library(furrr)
    future::plan(multisession, workers = 8)
    
    cat("Isotope annotation for feature table negative mode...\n")
    system.time(
      isotope_neg <-
        furrr::future_map(
          .x = as.data.frame(t(annotation_table_neg)),
          .f = function(x) {
            adduct <-
              stringr::str_extract(x[11], "\\(.+\\)") %>%
              stringr::str_replace("\\(", "") %>%
              stringr::str_replace("\\)", "")
            
            mz = as.numeric(stringr::str_trim(x[2], side = "both"))
            rt = as.numeric(stringr::str_trim(x[3], side = "both"))
            
            temp_iso <- try(annotate_isotope(
              formula = stringr::str_trim(x[17], side = "both"),
              adduct = adduct,
              mz = mz,
              rt = rt,
              peak.mz = ms1_data_neg$mz,
              peak.rt = ms1_data_neg$rt,
              rt.tol = peak_cluster_rt_tol,
              mz.tol = ms1.match.ppm,
              max.isotope = 3
            ),
            silent = TRUE)
            if (class(temp_iso) == "try-error") {
              return(NULL)
            }
            
            if (is.null(temp_iso)) {
              return(NULL)
            }
            
            temp_iso <-
              cbind(ms1_data_neg[temp_iso$peakIndex, ], temp_iso) %>%
              dplyr::select(-c(peakIndex))
            
            colnames(temp_iso) <- c("name",
                                    "mz",
                                    "rt",
                                    "condition",
                                    "polarity",
                                    "mz.error",
                                    "isotope",
                                    "RT.error")
            
            x <- matrix(x, nrow = 1) %>% as.data.frame()
            colnames(x) <- colnames(annotation_table_neg)
            
            temp_iso$Compound.name <- x$Compound.name
            temp_iso$Lab.ID <- x$Lab.ID
            temp_iso$Adduct <- x$Adduct
            temp_iso$Formula <- x$Formula
            temp_iso$rt <- x$rt
            temp_iso$CAS.ID <- x$CAS.ID
            temp_iso$HMDB.ID <- x$HMDB.ID
            temp_iso$KEGG.ID <- x$KEGG.ID
            temp_iso$Database <- x$Database
            
            temp_iso$mz.match.score <-
              (ms1.match.ppm - temp_iso$mz.error) / ms1.match.ppm
            temp_iso$RT.match.score <-
              (peak_cluster_rt_tol - temp_iso$RT.error) / peak_cluster_rt_tol
            temp_iso %>%
              dplyr::select(colnames(x))
          },
          .progress = TRUE
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
    )
    
    annotation_table_neg =
      rbind(annotation_table_neg, isotope_neg) %>%
      as.data.frame() %>%
      dplyr::arrange(Compound.name, isotope, rt)
    
    rm(list = c("isotope_pos", "isotope_neg"))
    
    annotation_table =
      rbind(annotation_table_pos, annotation_table_neg) %>%
      dplyr::arrange(Lab.ID, rt)
    
    ####combine peaks as group according to compound and rt
    library(plyr)
    
    annotation_table <-
      unique(annotation_table$Lab.ID) %>%
      furrr::future_map(
        .f = function(temp_id) {
          x = annotation_table %>%
            dplyr::filter(Lab.ID == temp_id)
          
          x <- x %>%
            dplyr::mutate(rt = as.numeric(rt)) %>%
            dplyr::arrange(rt)
          
          rt_class <-
            group_peaks_rt(rt = x$rt, rt.tol = peak_cluster_rt_tol) %>%
            dplyr::arrange(rt)
          
          rt_class <- paste(x$Lab.ID[1], rt_class$class, sep = "_")
          
          x =
            data.frame(x,
                       compound_class = rt_class,
                       stringsAsFactors = FALSE)
          
          x =
            unique(x$compound_class) %>%
            purrr::map(function(y) {
              z =
                x[x$compound_class == y, , drop = FALSE]
              score <- score_mfcp(z)
              z = data.frame(z, score, stringsAsFactors = FALSE)
              z
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame()
          x
        },
        .progress = TRUE
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    # head(annotation_table)
    
    # #####---------------------------------------------------------------------------
    # ###score distributation
    # temp_data <-
    #   annotation_table %>%
    #   dplyr::select(compound_class, score) %>%
    #   dplyr::distinct(compound_class, score) %>%
    #   dplyr::pull(score) %>%
    #   table() %>%
    #   as.data.frame()
    #
    # colnames(temp_data) <- c("Score", "Freq")
    #
    # temp_data =
    #   temp_data %>%
    #   dplyr::mutate(Score = as.numeric(as.character(Score)))
    #
    # plot <-
    #   temp_data %>%
    #   ggplot(aes(Score, Freq)) +
    #   geom_bar(stat = "identity",
    #            aes(x = Score,
    #                y = Freq,
    #                fill = Score),
    #            show.legend = FALSE) +
    #   # scale_fill_manual(
    #   #   values = c(
    #   #     "20" = "#3F4041FF",
    #   #     "30" = "#3F4041FF",
    #   #     "40" = "#84D7E1FF",
    #   #     "50" = "#84D7E1FF",
    #   #     "60" = "#84D7E1FF",
    #   #     "70" = "#84D7E1FF",
    #   #     "80" = "#84D7E1FF",
    #   #     "90" = "#84D7E1FF",
    #   #     "100" = "#008EA0FF",
    # #     "110" = "#008EA0FF",
    # #     "120" = "#008EA0FF",
    # #     "130" = "#008EA0FF",
    # #     "140" = "#008EA0FF",
    # #     "150" = "#008EA0FF",
    # #     "160" = "#008EA0FF",
    # #     "170" = "#008EA0FF",
    # #     "180" = "#008EA0FF",
    # #     "190" = "#008EA0FF",
    # #     "200" = "#C71000FF"
    # #   )
    # # ) +
    # theme_classic() +
    #   scale_y_continuous(expand = c(0, 10)) +
    #   scale_x_continuous(breaks = seq(0, 200, by = 20),
    #                      name = seq(0, 200, by = 20)) +
    #   theme(
    #     axis.title = element_text(size = 13),
    #     axis.text = element_text(size = 12),
    #     plot.background  = element_rect(fill = "transparent",
    #                                     color = NA),
    #     panel.background = element_rect(fill = "transparent",
    #                                     color = NA),
    #     legend.background = element_rect(fill = "transparent",
    #                                      color = NA),
    #     axis.text.x = element_text(
    #       size = 12,
    #       angle = 45,
    #       hjust = 1,
    #       vjust = 1
    #     )
    #   ) +
    #   labs(x = "Score") +
    #   ggforce::facet_zoom(
    #     xlim = c(110, 200),
    #     ylim = c(0, 400),
    #     horizontal = FALSE,
    #     zoom.size = 1
    #   )
    #
    # plot
    #
    # ggsave(
    #   plot,
    #   filename = "score_distributation.pdf",
    #   width = 8,
    #   height = 7,
    #   bg = "transparent"
    # )
    
    
    ######--------------------------------------------------------------------------
    ######gold standard from Level 1
    library(plyr)
    
    # golden_standard =
    #   variable_info %>%
    #   dplyr::filter(!is.na(Level)) %>%
    #   dplyr::filter(Level == 1 & SS > 0.7) %>%
    #   dplyr::select(
    #     name,
    #     Compound.name2 = Compound.name,
    #     CAS.ID2 = CAS.ID,
    #     HMDB.ID2 = HMDB.ID,
    #     KEGG.ID2 = KEGG.ID,
    #     Adduct2 = Adduct,
    #     mz.error2 = mz.error,
    #     RT.error2 = RT.error,
    #     SS2 = SS,
    #     Total.score2 = Total.score,
    #     Database2 = Database
    #   )
    #
    # temp_data <-
    #   golden_standard %>%
    #   dplyr::left_join(annotation_table, by = "name") %>%
    #   dplyr::select(
    #     c(
    #       compound_class,
    #       name,
    #       Compound.name,
    #       Compound.name2,
    #       CAS.ID,
    #       CAS.ID2,
    #       HMDB.ID,
    #       HMDB.ID2,
    #       KEGG.ID,
    #       KEGG.ID2,
    #       Adduct,
    #       Adduct2,
    #       Database,
    #       Database2,
    #       Formula,
    #       isotope,
    #       score
    #     )
    #   ) %>%
    #   plyr::dlply(.variables = .(name)) %>%
    #   purrr::map(function(x) {
    #     x %>%
    #       dplyr::arrange(desc(score), compound_class)
    #   })
    #
    # length(temp_data)
    #
    # temp_data[[2]]
    #
    # temp_data %>% lapply(function(x)
    #   max(x$score))  %>% unlist() %>% `>`(180) %>% which()
    #
    # temp_data[[196]]
    
    
    # #score distribution
    # score_rule <-
    #   data.frame(
    #     rule =
    #       c(
    #         "Adduct M+H",
    #         "Isotope (Adduct M+H)",
    #         "Other positive adduct",
    #         "Isotope (Other positive adduct)",
    #         "Adduct M-H",
    #         "Isotope (Adduct M-H)",
    #         "Other negative adduct",
    #         "Isotope (Other negative adduct)"
    #       ),
    #     score = c(50, 20, 20, 10, 50, 20, 20, 10),
    #     stringsAsFactors = FALSE
    #   )
    #
    # 3 * 3 * 3 * 3
    #
    # test <- list(c(1, 1), c(1, 0), c(0, 0))
    #
    # # comb <- vector(mode = "list", length = 3^4)
    # comb <- NULL
    # for (i in 1:3) {
    #   for (j in 1:3) {
    #     for (k in 1:3) {
    #       for (z in 1:3) {
    #         comb <- c(comb, list(c(test[[i]], test[[j]], test[[k]], test[[z]])))
    #       }
    #     }
    #   }
    # }
    #
    #
    # comb <-
    #   comb %>%
    #   do.call(cbind, .)
    #
    # remove_idx <-
    #   apply(comb, 2, function(x) {
    #     all(x == 0)
    #   }) %>%
    #   which()
    #
    # comb <- comb[, -remove_idx]
    #
    # score <- score_rule$score * comb
    # score <- score %>% colSums()
    #
    # idx <- order(score, decreasing = TRUE)
    #
    # score <- score[idx]
    #
    # comb <- comb[, idx]
    #
    # score_rule <-
    #   data.frame(score_rule, comb, stringsAsFactors = FALSE)
    #
    #
    # plot1 <-
    #   data.frame(index = factor(paste('X', 1:length(score), sep = ""),
    #                             levels = paste('X', 1:length(score), sep = "")),
    #              score,
    #              stringsAsFactors = TRUE) %>%
    #   ggplot(aes(index, score)) +
    #   geom_point(
    #     stat = "identity",
    #     aes(color = as.character(score)),
    #     shape = 16,
    #     show.legend = FALSE
    #   ) +
    #   geom_segment(aes(
    #     x = index,
    #     y = 0,
    #     xend = index,
    #     yend = score,
    #     color = as.character(score)
    #   ),
    #   show.legend = FALSE) +
    #   theme_classic() +
    #   scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    #   labs(x = "", y = "Score") +
    #   theme(
    #     axis.title.x = element_blank(),
    #     axis.text.x = element_blank(),
    #     axis.ticks.x = element_blank(),
    #     axis.title.y = element_text(size = 10),
    #     axis.text.y = element_text(size = 10),
    #     plot.background = element_rect(fill = "transparent", color = NA),
    #     panel.background = element_rect(fill = "transparent", color = NA),
    #     legend.background = element_rect(fill = "transparent", color = NA),
    #     plot.margin = unit(c(0, 0, 0, 0), "pt")
    #   )
    #
    # plot1
    #
    # plot2 <-
    #   score_rule %>%
    #   dplyr::select(-score) %>%
    #   tidyr::pivot_longer(cols = -rule,
    #                       names_to = "index",
    #                       values_to = "value") %>%
    #   dplyr::mutate(rule = factor(
    #     rule,
    #     levels =
    #       c(
    #         "Adduct M+H",
    #         "Isotope (Adduct M+H)",
    #         "Other positive adduct",
    #         "Isotope (Other positive adduct)",
    #         "Adduct M-H",
    #         "Isotope (Adduct M-H)",
    #         "Other negative adduct",
    #         "Isotope (Other negative adduct)"
    #       ) %>% rev()
    #   )) %>%
    #   ggplot(aes(index, rule)) +
    #   geom_tile(aes(fill = as.character(value)),
    #             color = "#FF6F00FF",
    #             show.legend = FALSE) +
    #   theme_bw() +
    #   labs(x = "", y = "") +
    #   scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    #   scale_fill_manual(values = c("0" = "white", "1" = "#008EA0FF")) +
    #   theme(
    #     panel.grid = element_blank(),
    #     axis.title.x = element_blank(),
    #     axis.text.x = element_blank(),
    #     axis.ticks.x = element_blank(),
    #     axis.text.y = element_text(size = 10),
    #     plot.margin = unit(c(0, 0, 0, 0), "pt")
    #   )
    #
    # plot2
    #
    # library(patchwork)
    #
    # plot <-
    #   plot1 + plot2 + patchwork::plot_layout(ncol = 1, heights = c(1, 3))
    #
    #
    # plot
    #
    #
    # ggsave(
    #   plot,
    #   filename = "score_rule.pdf",
    #   width = 12,
    #   height = 7,
    #   bg = "transparent"
    # )
    
    ###remove redundant annotation according to compound
    library(plyr)
    ###redundancy is
    calculate_redundance(annotation_table = annotation_table)
    
    cat("Remove redundancy...\n")
    annotation_table =
      remove_redundancy(annotation_table = annotation_table)
    
    calculate_redundance(annotation_table = annotation_table)
    
    
    ####---------------------------------------------------------------------------
    ####Metabolite set enrichment
    annotation_table <-
      annotation_table %>%
      dplyr::arrange(dplyr::desc(condition))
    cat("All done!\n")
    unlink(file.path(path, "Result"), recursive = TRUE)
    return(annotation_table)
  }
