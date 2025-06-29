#' @title group peaks according to RT
#' @description group peaks according to RT
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param rt rt
#' @param rt.tol rt.tol
#' @export
## group peaks according to RT
group_peaks_rt <- function(rt, rt.tol = 10) {
  rt_class <-
    lapply(rt, function(x) {
      which(rt >= x & rt < x + rt.tol)
    })
  
  rt_class <-
    lapply(seq_along(rt_class)[-1], function(i) {
      setdiff(rt_class[[i]], unlist(rt_class[1:(i - 1)]) %>% unique())
    }) %>%
    `c`(rt_class[1], .)
  
  rt_class <-
    rt_class[which(lapply(rt_class, length) != 0)]
  
  names(rt_class) <- seq_along(rt_class)
  
  rt_class <-
    purrr::map2(
      .x = rt_class,
      .y = names(rt_class),
      .f = function(x, y) {
        data.frame(rt = rt[x],
                   class = y,
                   stringsAsFactors = FALSE)
      }
    ) %>%
    do.call(rbind, .)
  rownames(rt_class) <- NULL
  return(rt_class)
}