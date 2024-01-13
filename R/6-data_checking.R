#' @title Check data format for fMSEA
#' @description Check data format for fMSEA
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object_list object list, the items should be mass_dataset class from massdataset
#' @return Checking results
#' @export

check4fmsea <- function(object_list) {
  if (!is.list(object_list)) {
    stop("object_list should be a list.")
  }
  
  check_is_mass_dataset <-
    object_list %>%
    lapply(function(x) {
      is(x, class2 = "mass_dataset")
    }) %>%
    unlist()
  
  if (any(!check_is_mass_dataset)) {
    stop(paste(which(!check_is_mass_dataset), collapse = ","),
         " in object_list is (are) not mass_dataset.")
  }
  
  column_polarity_check <-
    object_list %>%
    purrr::map(
      .f = function(x) {
        variable_info <- slot(x, "variable_info")
        if (all(colnames(variable_info) != "column")) {
          error1 <- c("no column")
          column <- NA
        } else{
          error1 <- c("no")
          column <- unique(variable_info$column)
        }
        
        if (all(colnames(variable_info) != "polarity")) {
          error2 <- c("no polarity")
          polarity <- NA
        } else{
          error2 <- c("no")
          polarity <- unique(variable_info$polarity)
        }
        
        list(
          error = c(error1, error2),
          column = column,
          polarity = polarity
        )
      }
    )
  
  error <-
    column_polarity_check %>%
    purrr::map(function(x) {
      x$error
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  seq_len(nrow(error)) %>%
    lapply(function(i) {
      temp_error <-
        as.character(error[i, ])
      
      if (any(temp_error != "no")) {
        temp_error[temp_error != "no"] %>%
          purrr::walk(function(x) {
            stop("object ", i , " has ", x)
          })
      }
    })
  
  column <-
    column_polarity_check %>%
    purrr::map(function(x) {
      x$column
    })
  
  seq_along(column) %>%
    purrr::walk(function(i) {
      if (length(column[[i]]) > 1) {
        stop("object ",
             i,
             " column are not unique: ",
             paste(column[[i]], collapse = ";"))
      }
      if (!column[[i]] %in% c("RPLC", "HILIC")) {
        stop("object ", i, " column should be RPLC or HILIC: ", column[[i]])
      }
    })
  
  
  polarity <-
    column_polarity_check %>%
    purrr::map(function(x) {
      x$polarity
    })
  
  seq_along(polarity) %>%
    purrr::walk(function(i) {
      if (length(polarity[[i]]) > 1) {
        stop("object ",
             i,
             " polarity are not unique: ",
             paste(polarity[[i]], collapse = ";"))
      }
      if (!polarity[[i]] %in% c("positive", "negative")) {
        stop("object ",
             i,
             " polarity should be positive or negative: ",
             polarity[[i]])
      }
    })
  
}