#' @title calculate_redundance
#' @description calculate_redundance
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param annotation_table annotation table
#' @importFrom data.table rbindlist .N
#' @return redundancy
#' @export

calculate_redundance <-
  function(annotation_table) {
    if (is(annotation_table, "list")) {
      annotation_table <-
        annotation_table %>%
        data.table::rbindlist()
    }
    
    ##redundancy1 means one compound contains how many compound class
    annotation_table <-
      data.table::as.data.table(annotation_table)
    
    redundancy1 <-
      annotation_table[, .(number = length(unique(metabolite_feature_cluster))), by = Lab.ID]
    
    redundancy1 <-
      mean(redundancy1$number)
    
    ##redundancy2 means one peak matches how many compound
    redundancy2 <-
      annotation_table[, .N, by = variable_id] %>%
      pull(N) %>%
      mean()
    
    c("mfcs/metabolite" = redundancy1, "metabolites/feature" = redundancy2)
  }



remove_redundancy <-
  function(annotation_table) {
    ##for one compound, if one compound class with score > 100, then remove the
    ##compound class with score <= 20
    
    redundancy_diff = c(-1, -1)
    
    while (any(redundancy_diff < 0)) {
      # cat("i", " ")
      before_redundancy =
        calculate_redundance(annotation_table = annotation_table)
      
      annotation_table <-
        annotation_table %>%
        dplyr::group_by(Lab.ID) %>%
        dplyr::filter(if (any(score > 100)) {
          score > 20
        } else{
          score > 0
        }) %>%
        dplyr::ungroup()
      
      ##for one peak, if it has a annotation with score > 100, then remove other
      ##annotations with score <= 20
      annotation_table <-
        annotation_table %>%
        dplyr::group_by(name) %>%
        dplyr::filter(if (any(score > 100)) {
          score > 20
        } else{
          score > 0
        }) %>%
        dplyr::ungroup()
      
      
      ###re-calculated confidence score for each compound class.
      annotation_table =
        annotation_table %>%
        plyr::dlply(.variables = .(metabolite_feature_cluster)) %>%
        purrr::map(function(x) {
          score <- score_mfc(x)
          x$score = score
          x
        }) %>%
        dplyr::bind_rows()
      
      after_redundancy =
        calculate_redundance(annotation_table = annotation_table)
      
      redundancy_diff = after_redundancy - before_redundancy
    }
    
    annotation_table
  }
