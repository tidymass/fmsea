#' Annotate Metabolite Features Using MS1 and MS2 Data
#'
#' Performs metabolite annotation for LC-MS feature tables using MS1 data, with optional integration of MS2-based annotations.
#' The function supports polarity-specific annotation, adduct handling, isotope identification, and redundant annotation reduction.
#'
#' @param feature_table A data frame with metabolite features. Must include columns:
#'   `variable_id`, `mz`, `rt`, `condition`, `polarity`, and `mean_intensity`.
#' @param annotation_table_ms2 (Optional) A data frame containing MS2-based annotations.
#'   Must include `variable_id`, `mz`, `rt`, `KEGG.ID`, and other identification information.
#' @param column Character. Chromatographic column type: either `"rp"` (reverse-phase) or `"hilic"`. Default is `"rp"`.
#' @param metabolite_database A metabolite database object (from `metid`) used for annotation.
#' @param ms1_match_ppm Numeric. m/z tolerance (in ppm) for MS1 matching. Default is 15.
#' @param mfc_rt_tol Numeric. Retention time tolerance (in seconds) for metabolite feature cluster scoring. Default is 5.
#' @param isotope_number Integer. Maximum number of isotopes to annotate. Default is 3.
#'
#' @return A data frame containing the final annotated metabolite feature table. It includes matched annotations,
#'   isotope annotations, compound-level metadata, and scoring metrics.
#'
#' @details
#' The function includes the following steps:
#' \itemize{
#'   \item MS1-based metabolite annotation (positive and negative modes separately)
#'   \item Optional inclusion of MS2-based annotations
#'   \item Isotope matching and scoring
#'   \item Feature merging into metabolite feature clusters (MFCs)
#'   \item Annotation quality scoring and redundancy calculation
#' }
#'
#' Parallel computing is enabled for isotope annotation using `future` and `furrr` packages.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' data("example_feature_table")
#' data("example_metabolite_database")
#'
#' annotated <- annotate_feature_table(
#'   feature_table = example_feature_table,
#'   metabolite_database = example_metabolite_database,
#'   column = "rp",
#'   ms1_match_ppm = 15,
#'   mfc_rt_tol = 5
#' )
#' }
#'
#' @importFrom dplyr filter select mutate left_join arrange desc pull
#' @importFrom stringr str_extract str_replace
#' @importFrom future plan multisession
#' @importFrom furrr future_map
#' @importFrom data.table as.data.table rbindlist
#'
#' @export

annotate_feature_table <-
  function(feature_table,
           annotation_table_ms2,
           column = c("rp", "hilic"),
           metabolite_database,
           ms1_match_ppm = 15,
           mfc_rt_tol = 5,
           isotope_number = 3) {
    column = match.arg(column)
    ####feature_table is required
    ###check feature_table
    if (missing(feature_table)) {
      stop("feature_table is required.")
    } else{
      check_feature_table(feature_table)
    }
    
    ######annotation_table_ms2 is optional
    ##check annotation_table_ms2
    if (!missing(annotation_table_ms2)) {
      check_annotation_table_ms2(annotation_table_ms2)
      annotation_table_ms2 <-
        annotation_table_ms2 %>%
        dplyr::filter(!is.na(KEGG.ID)) %>%
        dplyr::mutate(isotope = "[M]")
      if (nrow(annotation_table_ms2) == 0) {
        annotation_table_ms2 <- NULL
      }
    } else{
      annotation_table_ms2 <- NULL
    }
    
    ##metabolite_database is required
    if (missing(metabolite_database)) {
      stop("metabolite_database is required.")
    }
    
    feature_table_pos <-
      feature_table %>%
      dplyr::filter(polarity == "positive")
    
    feature_table_neg <-
      feature_table %>%
      dplyr::filter(polarity == "negative")
    
    
    ############metabolite database matching for feature_table
    message("Annotating features positive mode...\n")
    
    if (nrow(feature_table_pos) == 0) {
      annotation_table_ms1_pos <- NULL
    } else{
      expression_data <-
        data.frame(sample_1 = rep(1, nrow(feature_table_pos)))
      sample_info <-
        data.frame(sample_id = "sample_1", class = "Subject")
      rownames(expression_data) <- feature_table_pos$variable_id
      object_pos <-
        massdataset::create_mass_dataset(
          expression_data = expression_data,
          sample_info = sample_info,
          variable_info = feature_table_pos
        )
      
      object_pos <-
        metid::annotate_metabolites(
          object = object_pos,
          based_on = "ms1",
          ms1.match.ppm = ms1_match_ppm,
          column = column,
          polarity = "positive",
          database = metabolite_database,
          candidate.num = 1000
        )
      
      annotation_table_ms1_pos <-
        object_pos@annotation_table
      annotation_table_ms1_pos$polarity <- "positive"
    }
    
    message("Annotating features negative mode...\n")
    
    if (nrow(feature_table_neg) == 0) {
      annotation_table_ms1_neg <- NULL
    } else{
      expression_data <-
        data.frame(sample_1 = rep(1, nrow(feature_table_neg)))
      sample_info = data.frame(sample_id = "sample_1", class = "Subject")
      rownames(expression_data) <- feature_table_neg$variable_id
      
      object_neg <-
        massdataset::create_mass_dataset(
          expression_data = expression_data,
          sample_info = sample_info,
          variable_info = feature_table_neg
        )
      
      object_neg <-
        metid::annotate_metabolites(
          object = object_neg,
          based_on = "ms1",
          ms1.match.ppm = ms1_match_ppm,
          column = column,
          polarity = "negative",
          database = metabolite_database,
          candidate.num = 1000
        )
      
      annotation_table_ms1_neg <-
        object_neg@annotation_table
      annotation_table_ms1_neg$polarity <- "negative"
    }
    
    annotation_table_ms1 <-
      rbind(annotation_table_ms1_pos, annotation_table_ms1_neg) %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(KEGG.ID)) %>%
      dplyr::mutate(isotope = "[M]")
    
    ###add Compound formula to annotation table
    # annotation_table_ms1 <-
    #   annotation_table_ms1 %>%
    # dplyr::arrange(desc(condition)) %>%
    # dplyr::left_join(metabolite_database@spectra.info[, c("KEGG.ID", "Formula")], by = "KEGG.ID") %>%
    # dplyr::filter(!is.na(KEGG.ID)) %>%
    # dplyr::select(-c(ms2_spectrum_id, CE)) %>%
    # dplyr::mutate(isotope = "[M]")
    
    ######add annotation information from MS2 results
    if (!is.null(annotation_table_ms2)) {
      annotation_table_ms1 <-
        annotation_table_ms1 %>%
        dplyr::filter(!variable_id %in% annotation_table_ms2$variable_id)
      
      ###rbind annotation_table_ms1 and annotation_table_ms2
      intersect_column_names <-
        intersect(colnames(annotation_table_ms1),
                  colnames(annotation_table_ms2))
      
      annotation_table_final <-
        rbind(annotation_table_ms1[, intersect_column_names], annotation_table_ms2[, intersect_column_names])
      
    } else{
      annotation_table_final <-
        annotation_table_ms1 %>%
        dplyr::select(
          "variable_id",
          "Compound.name",
          "CAS.ID",
          "HMDB.ID",
          "KEGG.ID",
          "Lab.ID",
          "Adduct",
          "mz.error",
          "mz.match.score",
          "RT.error",
          "RT.match.score",
          "CE",
          "SS",
          "Total.score",
          "Database",
          "Level",
          "polarity",
          "isotope"
        )
    }
    
    #####add formula to annotation table
    annotation_table_final <-
      annotation_table_final %>%
      dplyr::left_join(metabolite_database@spectra.info[, c("KEGG.ID", "Formula")], by = "KEGG.ID") %>%
      dplyr::filter(!is.na(KEGG.ID)) %>%
      dplyr::filter(KEGG.ID != "")
    
    ####add condition and other information
    annotation_table_final <-
      annotation_table_final %>%
      dplyr::left_join(feature_table %>% dplyr::select(-polarity), by = "variable_id")
    
    ##------------------------------------------------------------------------------
    ####isotope annotation
    annotation_table_final_pos <-
      annotation_table_final %>%
      dplyr::filter(polarity == "positive")
    
    annotation_table_final_neg =
      annotation_table_final %>%
      dplyr::filter(polarity == "negative")
    
    ###positive
    library(future)
    library(furrr)
    future::plan(multisession, workers = 8)
    
    message("Isotope annotation for feature table positive mode...\n")
    
    if (nrow(annotation_table_final_pos) == 0) {
      isotope_pos <- NULL
    } else{
      system.time(
        isotope_pos <-
          furrr::future_map(
            .x = 1:nrow(annotation_table_final_pos),
            .f = function(x) {
              # cat(x, " ")
              adduct <-
                stringr::str_extract(annotation_table_final_pos$Adduct[x], "\\(.+\\)") %>%
                stringr::str_replace("\\(", "") %>%
                stringr::str_replace("\\)", "")
              
              temp_iso <- try(annotate_isotope(
                formula = annotation_table_final_pos$Formula[x],
                adduct = adduct,
                mz = annotation_table_final_pos$mz[x],
                rt = annotation_table_final_pos$rt[x],
                mean_intensity = annotation_table_final_pos$mean_intensity[x],
                peak.mz = feature_table_pos$mz,
                peak.rt = feature_table_pos$rt,
                peak.int = feature_table_pos$mean_intensity,
                rt.tol = mfc_rt_tol,
                mz.tol = ms1.match.ppm,
                int.tol = 0.3,
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
                cbind(feature_table_pos[temp_iso$peakIndex, ], temp_iso) %>%
                dplyr::select(-c(peakIndex))
              
              colnames(temp_iso) <- c(
                "variable_id",
                "mz",
                "rt",
                "condition",
                "polarity",
                "mean_intensity",
                "mz.error",
                "isotope",
                "RT.error"
              )
              temp_iso$Compound.name <-
                annotation_table_final_pos$Compound.name[x]
              temp_iso$Lab.ID <- annotation_table_final_pos$Lab.ID[x]
              temp_iso$Adduct <- annotation_table_final_pos$Adduct[x]
              temp_iso$Formula <- annotation_table_final_pos$Formula[x]
              temp_iso$rt <- annotation_table_final_pos$rt[x]
              temp_iso$CAS.ID <- annotation_table_final_pos$CAS.ID[x]
              temp_iso$HMDB.ID <- annotation_table_final_pos$HMDB.ID[x]
              temp_iso$KEGG.ID <- annotation_table_final_pos$KEGG.ID[x]
              temp_iso$Database <- annotation_table_final_pos$Database[x]
              temp_iso$SS <- annotation_table_final_pos$SS[x]
              temp_iso$Total.score <- annotation_table_final_pos$Total.score[x]
              temp_iso$Level <- annotation_table_final_pos$Level[x]
              temp_iso$CE <- annotation_table_final_pos$CE[x]
              
              temp_iso$mz.match.score <-
                (ms1.match.ppm - temp_iso$mz.error) / ms1.match.ppm
              temp_iso$RT.match.score <-
                (mfc_rt_tol - temp_iso$RT.error) / mfc_rt_tol
              temp_iso %>%
                dplyr::select(colnames(annotation_table_final_pos))
            },
            .progress = TRUE
          ) %>%
          do.call(rbind, .) %>%
          as.data.frame()
      )
    }
    
    annotation_table_final_pos <-
      rbind(annotation_table_final_pos, isotope_pos) %>%
      as.data.frame() %>%
      dplyr::arrange(Compound.name, isotope, rt)
    
    message("Isotope annotation for feature table negative mode...\n")
    
    system.time(
      isotope_neg <-
        furrr::future_map(
          .x = 1:nrow(annotation_table_final_neg),
          .f = function(x) {
            adduct <-
              stringr::str_extract(annotation_table_final_neg$Adduct[x], "\\(.+\\)") %>%
              stringr::str_replace("\\(", "") %>%
              stringr::str_replace("\\)", "")
            
            temp_iso <- try(annotate_isotope(
              formula = annotation_table_final_neg$Formula[x],
              adduct = adduct,
              mz = annotation_table_final_neg$mz[x],
              rt = annotation_table_final_neg$rt[x],
              mean_intensity = annotation_table_final_neg$mean_intensity[x],
              peak.mz = feature_table_neg$mz,
              peak.rt = feature_table_neg$rt,
              peak.int = feature_table_neg$mean_intensity,
              rt.tol = mfc_rt_tol,
              mz.tol = ms1.match.ppm,
              int.tol = 0.3,
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
              cbind(feature_table_neg[temp_iso$peakIndex, ], temp_iso) %>%
              dplyr::select(-c(peakIndex))
            
            colnames(temp_iso) <- c(
              "variable_id",
              "mz",
              "rt",
              "condition",
              "polarity",
              "mean_intensity",
              "mz.error",
              "isotope",
              "RT.error"
            )
            temp_iso$Compound.name <-
              annotation_table_final_neg$Compound.name[x]
            temp_iso$Lab.ID <- annotation_table_final_neg$Lab.ID[x]
            temp_iso$Adduct <- annotation_table_final_neg$Adduct[x]
            temp_iso$Formula <- annotation_table_final_neg$Formula[x]
            temp_iso$rt <- annotation_table_final_neg$rt[x]
            temp_iso$CAS.ID <- annotation_table_final_neg$CAS.ID[x]
            temp_iso$HMDB.ID <- annotation_table_final_neg$HMDB.ID[x]
            temp_iso$KEGG.ID <- annotation_table_final_neg$KEGG.ID[x]
            temp_iso$Database <- annotation_table_final_neg$Database[x]
            temp_iso$SS <- annotation_table_final_neg$SS[x]
            temp_iso$Total.score <- annotation_table_final_neg$Total.score[x]
            temp_iso$Level <- annotation_table_final_neg$Level[x]
            temp_iso$CE <- annotation_table_final_neg$CE[x]
            
            temp_iso$mz.match.score <-
              (ms1.match.ppm - temp_iso$mz.error) / ms1.match.ppm
            temp_iso$RT.match.score <-
              (mfc_rt_tol - temp_iso$RT.error) / mfc_rt_tol
            temp_iso %>%
              dplyr::select(colnames(annotation_table_final_neg))
          },
          .progress = TRUE
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
    )
    
    annotation_table_final_neg <-
      rbind(annotation_table_final_neg, isotope_neg) %>%
      as.data.frame() %>%
      dplyr::arrange(Compound.name, isotope, rt)
    
    annotation_table_final <-
      rbind(annotation_table_final_pos, annotation_table_final_neg) %>%
      dplyr::arrange(Lab.ID, rt)
    
    ####combine peaks to metabolite feature cluster (MFC)
    message("Scoring metabolite feature clusters...\n")
    annotation_table_final <-
      score_annotation_table(annotation_table = annotation_table_final, mfc_rt_tol = mfc_rt_tol)
    
    
    ###remove redundant annotation according to metabolite
    # library(plyr)
    ###redundancy is
    # browser()
    calculate_redundancy(annotation_table = annotation_table_final)
    
    # message("Removing redundancy...\n")
    #
    # annotation_table_final <-
    #   remove_redundancy(annotation_table = annotation_table_final)
    #
    # calculate_redundancy(annotation_table = annotation_table_final)
    
    
    ####---------------------------------------------------------------------------
    ####Metabolite set enrichment
    # annotation_table_final <-
    #   annotation_table_final %>%
    #   dplyr::arrange(dplyr::desc(condition))
    message("All done!\n")
    # unlink(file.path(path, "Result"), recursive = TRUE)
    return(annotation_table_final)
  }




#' Score Metabolite Feature Clusters (MFCs) Based on RT Grouping
#'
#' Assigns metabolite feature clusters (MFCs) within each compound (`Lab.ID`) by grouping features based on retention time (RT)
#' and calculating a cluster-level score using `score_mfc()`. Each feature is tagged with a `metabolite_feature_cluster` identifier and a computed `score`.
#'
#' @param annotation_table A data frame containing annotated metabolite features. Must include at least the columns
#'   `Lab.ID` and `rt` (retention time).
#' @param mfc_rt_tol Numeric. Retention time tolerance (in seconds) used for grouping features into the same cluster. Default is 10.
#'
#' @return A data frame with an additional `metabolite_feature_cluster` column (unique cluster ID) and `score` column for each feature.
#'
#' @details
#' \itemize{
#'   \item Features are grouped by `Lab.ID`, then clustered into MFCs based on RT proximity.
#'   \item Each MFC is scored via the `score_mfc()` function (assumed to be user-defined).
#'   \item Parallel computation is used via `furrr::future_map()` for faster processing.
#' }
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' data("example_annotation_table")
#' scored_table <- score_annotation_table(annotation_table = example_annotation_table, mfc_rt_tol = 5)
#' head(scored_table)
#' }
#'
#' @importFrom dplyr filter arrange
#' @importFrom furrr future_map
#' @importFrom purrr map
#' @importFrom future plan multisession
#'
#' @export

score_annotation_table <-
  function(annotation_table, mfc_rt_tol = 10) {
    ##calculate the score for each metabolite feature cluster
    annotation_table <-
      unique(annotation_table$Lab.ID) %>%
      furrr::future_map(
        .f = function(temp_id) {
          x <- annotation_table %>%
            dplyr::filter(Lab.ID == temp_id)
          
          x <- x %>%
            dplyr::arrange(rt)
          
          rt_class <-
            group_peaks_rt(rt = x$rt, rt.tol = mfc_rt_tol) %>%
            dplyr::arrange(rt)
          
          rt_class <- paste(x$Lab.ID[1], rt_class$class, sep = "@")
          
          x$metabolite_feature_cluster <-
            rt_class
          
          x <-
            unique(x$metabolite_feature_cluster) %>%
            purrr::map(function(y) {
              z =
                x[x$metabolite_feature_cluster == y, , drop = FALSE]
              score <- score_mfc(z)
              z$score <- score
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
    annotation_table <-
      annotation_table %>%
      dplyr::arrange(metabolite_feature_cluster, rt)
    return(annotation_table)
  }
