#' @title calculate_redundancy
#' @description calculate_redundancy
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param annotation_table annotation table
#' @importFrom data.table rbindlist .N
#' @return redundancy
#' @export

calculate_redundancy <-
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
      annotation_table[, .(number = length(unique(compound_class))), by = Lab.ID]
    
    redundancy1 <-
      mean(redundancy1$number)
    
    ##redundancy2 means one peak matches how many compound
    redundancy2 <-
      annotation_table[, .N, by = variable_id] %>%
      pull(N) %>%
      mean()
    
    c(redundancy1 = redundancy1,
      redundancy2 = redundancy2)
  }