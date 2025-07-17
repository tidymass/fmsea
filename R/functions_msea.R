





##' generic function for gene set enrichment analysis
##'
##'
##' @title MSEA_internal
##' @param annotation_table order ranked annotation_table
##' @param exponent weight of each step
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param eps This parameter sets the boundary for calculating the p value.
##' @param pvalueCutoff p value Cutoff
##' @param pAdjustMethod p value adjustment method
##' @param verbose print message or not
##' @param seed set seed inside the function to make result reproducible. FALSE by default.
##' @param USER_DATA annotation data
##' @param by one of 'fmsea' or 'DOSE'
##' @param ... other parameters
##' @return mseaResult object
##' @author Yu Guangchuang
MSEA_internal <- function(annotation_table,
                          exponent,
                          minGSSize,
                          maxGSSize,
                          eps,
                          pvalueCutoff,
                          pAdjustMethod,
                          verbose,
                          seed = FALSE,
                          USER_DATA,
                          by = "fmsea",
                          ...) {
  by <- match.arg(by, c("fmsea", "DOSE"))
  if (!is.sorted(annotation_table))
    stop("annotation_table should be a decreasing sorted vector...")
  if (by == 'fmsea') {
    .MSEA <- MSEA_fmsea
  } else {
    .MSEA <- MSEA_DOSE
  }
  
  res <- .MSEA(
    annotation_table          = annotation_table,
    exponent          = exponent,
    minGSSize         = minGSSize,
    maxGSSize         = maxGSSize,
    eps               = eps,
    pvalueCutoff      = pvalueCutoff,
    pAdjustMethod     = pAdjustMethod,
    verbose           = verbose,
    seed              = seed,
    USER_DATA         = USER_DATA,
    ...
  )
  
  res@organism <- "UNKNOWN"
  res@setType <- "UNKNOWN"
  res@keytype <- "UNKNOWN"
  return(res)
}


leading_edge <- function(observed_info) {
  core_enrichment <- lapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    running_enrichment_score <-
      running_enrichment_score[running_enrichment_score$position == 1, ]
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      leading_gene <- running_enrichment_score$gene[1:i]
    } else {
      i <- which.min(running_enrichment_score$running_score)
      leading_gene <- running_enrichment_score$gene[-c(1:(i - 1))]
    }
    return(leading_gene)
  })
  
  rank <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      rr <- which.max(running_enrichment_score$running_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      rr <- nrow(running_enrichment_score) - i + 1
    }
    return(rr)
  })
  
  tags <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    running_enrichment_score <-
      running_enrichment_score[running_enrichment_score$position == 1, ]
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      res <- i / nrow(running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      res <-
        (nrow(running_enrichment_score) - i + 1) / nrow(running_enrichment_score)
    }
    return(res)
  })
  
  ll <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      res <- i / nrow(running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      res <-
        (nrow(running_enrichment_score) - i + 1) / nrow(running_enrichment_score)
    }
    return(res)
  })
  
  N <- nrow(observed_info[[1]]$running_enrichment_score)
  setSize <-
    sapply(observed_info, function(x)
      sum(x$running_enrichment_score$position))
  signal <- tags * (1 - ll) * (N / (N - setSize))
  
  tags <- paste0(round(tags * 100), "%")
  ll <- paste0(round(ll * 100), "%")
  signal <- paste0(round(signal * 100), "%")
  leading_edge <-
    paste0('tags=', tags, ", list=", ll, ", signal=", signal)
  
  res <- list(
    rank = rank,
    tags = tags,
    list = ll,
    signal = signal,
    leading_edge = leading_edge,
    core_enrichment = core_enrichment
  )
  return(res)
}

permute_annotation_table <-
  function(annotation_table) {
    permute_idx <- sample.int(nrow(annotation_table))
    permute_annotation_table <- annotation_table
    permute_annotation_table$condition <-
      annotation_table$condition[permute_idx]
    
    permute_annotation_table =
      permute_annotation_table %>%
      dplyr::arrange(desc(condition))
    return(permute_annotation_table)
  }

# annotation_table = annotation_table
# metabolite_set = metabolite_set
# exponent = exponent

permute_msea_enrichment_score <-
  function(annotation_table,
           metabolite_set,
           exponent = 1) {
    random_annotation_table <-
      permute_annotation_table(annotation_table)
    
    res =
      get_msea_score(
        metabolite_set = metabolite_set,
        annotation_table = random_annotation_table,
        exponent = exponent
      )
    
    return(res$enrichment_score)
  }


geneSet_filter <-
  function(metabolite_set,
           annotation_table,
           minGSSize,
           maxGSSize) {
    metabolite_set <-
      sapply(metabolite_set, intersect, names(annotation_table))
    
    gs.idx <-
      get_geneSet_index(metabolite_set, minGSSize, maxGSSize)
    nGeneSet <- sum(gs.idx)
    
    if (nGeneSet == 0) {
      msg <-
        paste0("No gene set have size between [",
               minGSSize,
               ", ",
               maxGSSize,
               "]...")
      message(msg)
      message("--> return NULL...")
      return(NULL)
    }
    metabolite_set[gs.idx]
  }






# feature_list <- feature_list
# feature_set <- kegg_hsa_compound_pathway
# exponent = 1
# perm_num = 1000
# min_size = 5
# max_size = 1000
# pvalue_cutoff = 0.2
# p_adjust_method = "fdr"
# seed = FALSE
# verbose = TRUE
# library(tidyverse)
# result <-do_msea(feature_list = feature_list,
#                  feature_set = feature_set,
#                  exponent = 1, perm_num = 1000,
#                  min_size = 5, max_size = 1000,
#                  pvalue_cutoff = 0.2,
#                  p_adjust_method = "fdr",
#                  seed = TRUE,
#                  verbose = TRUE)


##' @name metPathMSEA-class
##' @aliases metPathMSEA-class
##'   show,metPathMSEA-method summary,metPathMSEA-method
##'
##' @docType class
##' @slot result MSEA anaysis result.
##' @slot metabolite_set metabolite_set
##' @slot feature_list order rank feature_list
##' @slot permute_scores permutation scores
##' @slot params parameters
##' @exportClass metPathMSEA
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @keywords classes
setClass(
  "metPathMSEA",
  representation   = representation(
    result          = "data.frame",
    feature_set        = "list",
    feature_list        = "numeric",
    perm_scores      = "matrix",
    params          = "list"
  )
)













###----------------------------------------------------------------------------
# feature_set <- kegg_hsa_compound_pathway
# feature_list <- feature_list

filter_feature_set <- function(feature_set,
                               feature_list,
                               min_size = 5,
                               max_size = 1000) {
  len <- lapply(feature_set, function(x) {
    length(intersect(x, names(feature_list)))
  }) %>%
    unlist() %>%
    unname()
  
  # len <- lapply(feature_set, length) %>%
  #   unlist()
  
  remain_idx <-
    which(len > min_size & len < max_size)
  
  feature_set <- feature_set[remain_idx]
  return(feature_set)
}





######credit to Guangchuang Yu
#annotation_table is the annotation table
#metabolite_set is the
get_msea_score <- function(annotation_table,
                           metabolite_set,
                           exponent = 1,
                           fortify = FALSE) {
  ###################################################################
  ##    annotation_table                                           ##
  ##                                                               ##
  ## 1. Rank order the N genes in D to form L = { g_1, ... , g_N}  ##
  ##    according to the correlation, r(g_j)=r_j,                  ##
  ##    of their expression profiles with C.                       ##
  ##                                                               ##
  ###################################################################
  
  ###################################################################
  ##    exponent                                                   ##
  ##                                                               ##
  ## An exponent p to control the weight of the step.              ##
  ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
  ##   the standard Kolmogorov-Smirnov statistic.                  ##
  ##   When p = 1, we are weighting the genes in S                 ##
  ##   by their correlation with C normalized                      ##
  ##   by the sum of the correlations over all of the genes in S.  ##
  ##                                                               ##
  ###################################################################
  
  ## metabolites defined in metabolite_set should appear in annotation_table.
  metabolite_set <-
    intersect(metabolite_set, annotation_table$Lab.ID)
  annotation_table =
    annotation_table %>%
    dplyr::arrange(desc(condition))
  annotation_table$order =
    1:nrow(annotation_table)
  ##N the length of annotation_table
  number_annotation_table <- nrow(annotation_table)
  ##number_metabolite_set the length of metabolite set
  number_metabolite_set <- length(metabolite_set)
  
  position_hit <-
    position_miss <- numeric(length = number_annotation_table)
  hits <- annotation_table$Lab.ID %in% metabolite_set ## logical
  
  position_hit[hits] <- abs(annotation_table$order[hits])^exponent
  NR <- sum(position_hit)
  position_hit <- cumsum(position_hit / NR)
  
  position_miss[!hits] <-
    1 / (number_annotation_table - number_metabolite_set)
  position_miss <- cumsum(position_miss)
  
  running_enrichment_score <- position_hit - position_miss
  
  # data.frame(index = 1:length(position_hit),
  #            position_hit, position_miss, running_enrichment_score,
  #            stringsAsFactors = FALSE) %>%
  #   tidyr::pivot_longer(cols = -index, names_to = "class", values_to = "value") %>%
  #   ggplot(aes(x = index, y = value, color = class)) +
  #   geom_point()
  
  ## ES is the maximum deviation from zero of position_hit-position_miss
  max_enrichment_score <- max(running_enrichment_score)
  min_enrichment_score <- min(running_enrichment_score)
  
  if (abs(max_enrichment_score) > abs(min_enrichment_score)) {
    enrichment_score <- max_enrichment_score
  } else {
    enrichment_score <- min_enrichment_score
  }
  
  df <- data.frame(
    x = seq_along(running_enrichment_score),
    running_score = running_enrichment_score,
    position = as.integer(hits)
  )
  
  # df %>%
  #   ggplot(aes(x, running_score)) +
  #   geom_point(aes(color = as.character(position)))
  
  if (fortify == TRUE) {
    return(df)
  }
  
  df$Lab.ID = annotation_table$Lab.ID
  df$name = annotation_table$name
  res <-
    list(enrichment_score = enrichment_score,
         running_enrichment_score = df)
  return(res)
}

perm_feature_list <- function(feature_list) {
  permute_idx <- sample.int(length(feature_list))
  perm_feature_list <- feature_list
  names(perm_feature_list) <- names(feature_list)[permute_idx]
  return(perm_feature_list)
}


perm_msea_es <- function(feature_list, feature_set, exponent = 1) {
  feature_list <- perm_feature_list(feature_list)
  res <- sapply(1:length(feature_set), function(i)
    get_msea_score(
      feature_list = feature_list,
      feature_set = feature_set[[i]],
      exponent = exponent
    )$enrichment_score)
  return(res)
}




calculate_qvalue <- function(pvals) {
  if (length(pvals) == 0)
    return(numeric(0))
  
  qobj <- tryCatch(
    qvalue::qvalue(pvals, lambda = 0.05, pi0.method = "bootstrap"),
    error = function(e)
      NULL
  )
  
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  return(qvalues)
}


###----------------------------------------------------------------------------
get_leading_edge <- function(observed_info) {
  core_enrichment <- lapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    running_enrichment_score <-
      running_enrichment_score[running_enrichment_score$position == 1, ]
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      leading_gene <- running_enrichment_score$gene[1:i]
    } else {
      i <- which.min(running_enrichment_score$running_score)
      leading_gene <- running_enrichment_score$gene[-c(1:(i - 1))]
    }
    return(leading_gene)
  })
  
  rank <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      rr <- which.max(running_enrichment_score$running_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      rr <- nrow(running_enrichment_score) - i + 1
    }
    return(rr)
  })
  
  tags <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    running_enrichment_score <-
      running_enrichment_score[running_enrichment_score$position == 1, ]
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      res <- i / nrow(running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      res <-
        (nrow(running_enrichment_score) - i + 1) / nrow(running_enrichment_score)
    }
    return(res)
  })
  
  ll <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      res <- i / nrow(running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      res <-
        (nrow(running_enrichment_score) - i + 1) / nrow(running_enrichment_score)
    }
    return(res)
  })
  
  N <- nrow(observed_info[[1]]$running_enrichment_score)
  setSize <-
    sapply(observed_info, function(x)
      sum(x$running_enrichment_score$position))
  signal <- tags * (1 - ll) * (N / (N - setSize))
  
  tags <- paste0(round(tags * 100), "%")
  ll <- paste0(round(ll * 100), "%")
  signal <- paste0(round(signal * 100), "%")
  leading_edge <-
    paste0('tags=', tags, ", list=", ll, ", signal=", signal)
  
  res <- list(
    rank = rank,
    tags = tags,
    list = ll,
    signal = signal,
    leading_edge = leading_edge,
    core_enrichment = core_enrichment
  )
  return(res)
}



##' show method for \code{metPathMSEA} instance
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @return message
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("show", signature(object = "metPathMSEA"), function (object) {
  params <- object@params
  cat("#\n# Gene Set Enrichment Analysis\n#\n")
  cat("#...@feature_list", "\t")
  str(object@feature_list)
  cat("#...nPerm", "\t", params$perm_num, "\n")
  cat(
    "#...pvalues adjusted by",
    paste0("'", params$p_adjust_method, "'"),
    paste0("with cutoff <", params$pvalueCutoff),
    "\n"
  )
  cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
  str(object@result)
})



##' visualize analyzing result of MSEA
##'
##' plotting function for mseaResult
##' @title msea_plot
##' @rdname mseaplot
##' @param x object of msea result
##' @param geneSetID geneSet ID
##' @param by one of "running_score" or "position"
##' @param title plot title
##' @param ... additional parameters
##' @return ggplot2 object
##' @export
##' @examples
##' library(DOSE)
##' data(annotation_table)
##' x <- gseDO(annotation_table)
##' mseaplot(x, geneSetID=1)
setGeneric(name = "msea_plot", function(x,
                                        feature_set_idx = 1,
                                        by = "all",
                                        title = "",
                                        color = 'black',
                                        color.line = "#8DD3C7",
                                        color.vline = "#FB8072",
                                        ...) {
  standardGeneric("msea_plot")
})


##' @rdname msea_plot
##' @exportMethod msea_plot
setMethod(f = "msea_plot", signature(x = "metPathMSEA"), function (x,
                                                                   feature_set_idx = 1,
                                                                   by = "all",
                                                                   title = "",
                                                                   color = 'black',
                                                                   color.line = "#8DD3C7",
                                                                   color.vline = "#FB8072",
                                                                   ...) {
  msea_plot.metPathMSEA(
    x,
    feature_set_idx = feature_set_idx,
    by = by,
    title = title,
    color = color,
    color.line = color.line,
    color.vline = color.vline,
    ...
  )
})

##' @rdname msea_plot
##' @param color color of line segments
##' @param color.line color of running enrichment score line
##' @param color.vline color of vertical line which indicating the maximum/minimal running enrichment score
##' @return ggplot2 object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_linerange
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 ggplotGrob
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 ggplot_gtable
##' @importFrom ggplot2 ggplot_build
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 rel
##' @importFrom cowplot plot_grid
##' @author Guangchuang Yu

msea_plot.metPathMSEA <-
  function (x,
            feature_set_idx = 1,
            by = "all",
            title = "",
            color = 'black',
            color.line = "green",
            color.vline = "#FA5860",
            ...) {
    if (is.null(x)) {
      return(NULL)
    }
    by <- match.arg(by, c("running_score", "preranked", "all"))
    gs_data <- get_gs_info(x, feature_set_idx)
    
    p <- ggplot(gs_data, aes_(x = ~ x)) +
      theme_dose() +
      xlab("Position in the ranked list of features")
    
    if (by == "running_score" || by == "all") {
      p.res <-
        p + geom_linerange(aes_(ymin =  ~ ymin, ymax =  ~ ymax), color = color)
      p.res <-
        p.res + geom_line(aes_(y = ~ running_score),
                          color = color.line,
                          size =
                            1)
      enrichmentScore <-
        x@result[feature_set_idx, "enrichmentScore"]
      
      es.df <-
        data.frame(es = which.min(abs(
          p$data$running_score - enrichmentScore
        )))
      
      p.res <-
        p.res + geom_vline(
          data = es.df,
          aes_(xintercept = ~ es),
          colour = color.vline,
          linetype = "dashed"
        )
      p.res <- p.res + ylab("Running enrichment score")
      p.res <- p.res + geom_hline(yintercept = 0)
    }
    
    if (by == "preranked" || by == "all") {
      df2 <- data.frame(x = which(p$data$position == 1))
      df2$y <- p$data$feature_list[df2$x]
      p.pos <-
        p + geom_segment(data = df2,
                         aes_(
                           x =  ~ x,
                           xend =  ~ x,
                           y =  ~ y,
                           yend = 0
                         ),
                         color = color)
      p.pos <-
        p.pos + ylab("Ranked list metric") + xlim(0, length(p$data$feature_list))
    }
    if (by == "running_score")
      return(p.res + ggtitle(title))
    if (by == "preranked")
      return(p.pos + ggtitle(title))
    
    p.pos <-
      p.pos + xlab(NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p.pos <- p.pos + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = rel(2)))
    cowplot::plot_grid(p.pos, p.res, ncol = 1, align = "v")
  }



get_gs_info <- function(object, feature_set_idx = 1) {
  feature_list <- object@feature_list
  # if (is.numeric(feature_set_idx))
  #   feature_set_idx <- object@result[feature_set_idx, "ID"]
  feature_set <- object@feature_set[[feature_set_idx]]
  exponent <- object@params[["exponent"]]
  df <-
    get_msea_score(feature_list, feature_set, exponent, fortify = TRUE)
  df$ymin = 0
  df$ymax = 0
  pos <- df$position == 1
  h <- diff(range(df$running_score)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$feature_list <- feature_list
  return(df)
}


##' ggplot theme of DOSE
##'
##' @title theme_dose
##' @param font.size font size
##' @return ggplot theme
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @examples
##' library(ggplot2)
##' qplot(1:10) + theme_dose()
##' @export
theme_dose <- function(font.size = 13) {
  theme_bw() +
    theme(
      axis.text.x = element_text(
        colour = "black",
        size = 12,
        vjust = 1
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 12,
        hjust = 1
      ),
      axis.title = element_text(
        margin = margin(10, 5, 0, 0),
        color = "black",
        size = 13
      ),
      axis.title.y = element_text(angle = 90)
    )
}









##' @name metPathMSEA-class
##' @aliases metPathMSEA-class
##'   show,metPathMSEA-method summary,metPathMSEA-method
##'
##' @docType class
##' @slot result MSEA anaysis result.
##' @slot metabolite_set metabolite_set
##' @slot feature_list order rank feature_list
##' @slot permute_scores permutation scores
##' @slot params parameters
##' @exportClass metPathMSEA
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @keywords classes
setClass(
  "metPathMSEA",
  representation   = representation(
    result          = "data.frame",
    feature_set        = "list",
    feature_list        = "numeric",
    perm_scores      = "matrix",
    params          = "list"
  )
)













###----------------------------------------------------------------------------
# feature_set <- kegg_hsa_compound_pathway
# feature_list <- feature_list

filter_feature_set2 <- function(feature_set,
                                feature_table,
                                min_size = 5,
                                max_size = 1000) {
  len <- lapply(feature_set, function(x) {
    length(intersect(x, feature_table$Lab.ID))
  }) %>%
    unlist() %>%
    unname()
  
  remain_idx <-
    which(len > min_size & len < max_size)
  
  feature_set <- feature_set[remain_idx]
  return(feature_set)
}



perm_feature_table <- function(feature_table) {
  permute_idx <- sample.int(nrow(feature_table))
  perm_feature_table <- feature_table
  perm_feature_table$Lab.ID <-
    perm_feature_table$Lab.ID[permute_idx]
  return(perm_feature_table)
}











calculate_qvalue <- function(pvals) {
  if (length(pvals) == 0)
    return(numeric(0))
  
  qobj <- tryCatch(
    qvalue::qvalue(pvals, lambda = 0.05, pi0.method = "bootstrap"),
    error = function(e)
      NULL
  )
  
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  return(qvalues)
}


###----------------------------------------------------------------------------
get_leading_edge <- function(observed_info) {
  core_enrichment <- lapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    running_enrichment_score <-
      running_enrichment_score[running_enrichment_score$position == 1, ]
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      leading_gene <- running_enrichment_score$gene[1:i]
    } else {
      i <- which.min(running_enrichment_score$running_score)
      leading_gene <- running_enrichment_score$gene[-c(1:(i - 1))]
    }
    return(leading_gene)
  })
  
  rank <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      rr <- which.max(running_enrichment_score$running_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      rr <- nrow(running_enrichment_score) - i + 1
    }
    return(rr)
  })
  
  tags <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    running_enrichment_score <-
      running_enrichment_score[running_enrichment_score$position == 1, ]
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      res <- i / nrow(running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      res <-
        (nrow(running_enrichment_score) - i + 1) / nrow(running_enrichment_score)
    }
    return(res)
  })
  
  ll <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_score)
      res <- i / nrow(running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_score)
      res <-
        (nrow(running_enrichment_score) - i + 1) / nrow(running_enrichment_score)
    }
    return(res)
  })
  
  N <- nrow(observed_info[[1]]$running_enrichment_score)
  setSize <-
    sapply(observed_info, function(x)
      sum(x$running_enrichment_score$position))
  signal <- tags * (1 - ll) * (N / (N - setSize))
  
  tags <- paste0(round(tags * 100), "%")
  ll <- paste0(round(ll * 100), "%")
  signal <- paste0(round(signal * 100), "%")
  leading_edge <-
    paste0('tags=', tags, ", list=", ll, ", signal=", signal)
  
  res <- list(
    rank = rank,
    tags = tags,
    list = ll,
    signal = signal,
    leading_edge = leading_edge,
    core_enrichment = core_enrichment
  )
  return(res)
}



##' show method for \code{metPathMSEA} instance
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @return message
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("show", signature(object = "metPathMSEA"), function (object) {
  params <- object@params
  cat("#\n# Gene Set Enrichment Analysis\n#\n")
  cat("#...@feature_list", "\t")
  str(object@feature_list)
  cat("#...nPerm", "\t", params$perm_num, "\n")
  cat(
    "#...pvalues adjusted by",
    paste0("'", params$p_adjust_method, "'"),
    paste0("with cutoff <", params$pvalueCutoff),
    "\n"
  )
  cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
  str(object@result)
})



##' visualize analyzing result of MSEA
##'
##' plotting function for mseaResult
##' @title msea_plot
##' @rdname mseaplot
##' @param x object of msea result
##' @param geneSetID geneSet ID
##' @param by one of "running_score" or "position"
##' @param title plot title
##' @param ... additional parameters
##' @return ggplot2 object
##' @export
##' @examples
##' library(DOSE)
##' data(annotation_table)
##' x <- gseDO(annotation_table)
##' mseaplot(x, geneSetID=1)
setGeneric(name = "msea_plot", function(x,
                                        feature_set_idx = 1,
                                        by = "all",
                                        title = "",
                                        color = 'black',
                                        color.line = "#8DD3C7",
                                        color.vline = "#FB8072",
                                        ...) {
  standardGeneric("msea_plot")
})


##' @rdname msea_plot
##' @exportMethod msea_plot
setMethod(f = "msea_plot", signature(x = "metPathMSEA"), function (x,
                                                                   feature_set_idx = 1,
                                                                   by = "all",
                                                                   title = "",
                                                                   color = 'black',
                                                                   color.line = "#8DD3C7",
                                                                   color.vline = "#FB8072",
                                                                   ...) {
  msea_plot.metPathMSEA(
    x,
    feature_set_idx = feature_set_idx,
    by = by,
    title = title,
    color = color,
    color.line = color.line,
    color.vline = color.vline,
    ...
  )
})

##' @rdname msea_plot
##' @param color color of line segments
##' @param color.line color of running enrichment score line
##' @param color.vline color of vertical line which indicating the maximum/minimal running enrichment score
##' @return ggplot2 object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_linerange
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 ggplotGrob
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 ggplot_gtable
##' @importFrom ggplot2 ggplot_build
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 rel
##' @importFrom cowplot plot_grid
##' @author Guangchuang Yu

msea_plot.metPathMSEA <-
  function (x,
            feature_set_idx = 1,
            by = "all",
            title = "",
            color = 'black',
            color.line = "green",
            color.vline = "#FA5860",
            ...) {
    if (is.null(x)) {
      return(NULL)
    }
    by <- match.arg(by, c("running_score", "preranked", "all"))
    gs_data <- get_gs_info(x, feature_set_idx)
    
    p <- ggplot(gs_data, aes_(x = ~ x)) +
      theme_dose() +
      xlab("Position in the ranked list of features")
    
    if (by == "running_score" || by == "all") {
      p.res <-
        p + geom_linerange(aes_(ymin =  ~ ymin, ymax =  ~ ymax), color = color)
      p.res <-
        p.res + geom_line(aes_(y = ~ running_score),
                          color = color.line,
                          size =
                            1)
      enrichmentScore <-
        x@result[feature_set_idx, "enrichmentScore"]
      
      es.df <-
        data.frame(es = which.min(abs(
          p$data$running_score - enrichmentScore
        )))
      
      p.res <-
        p.res + geom_vline(
          data = es.df,
          aes_(xintercept = ~ es),
          colour = color.vline,
          linetype = "dashed"
        )
      p.res <- p.res + ylab("Running enrichment score")
      p.res <- p.res + geom_hline(yintercept = 0)
    }
    
    if (by == "preranked" || by == "all") {
      df2 <- data.frame(x = which(p$data$position == 1))
      df2$y <- p$data$feature_list[df2$x]
      p.pos <-
        p + geom_segment(data = df2,
                         aes_(
                           x =  ~ x,
                           xend =  ~ x,
                           y =  ~ y,
                           yend = 0
                         ),
                         color = color)
      p.pos <-
        p.pos + ylab("Ranked list metric") + xlim(0, length(p$data$feature_list))
    }
    if (by == "running_score")
      return(p.res + ggtitle(title))
    if (by == "preranked")
      return(p.pos + ggtitle(title))
    
    p.pos <-
      p.pos + xlab(NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p.pos <- p.pos + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = rel(2)))
    cowplot::plot_grid(p.pos, p.res, ncol = 1, align = "v")
  }



get_gs_info <- function(object, feature_set_idx = 1) {
  feature_list <- object@feature_list
  # if (is.numeric(feature_set_idx))
  #   feature_set_idx <- object@result[feature_set_idx, "ID"]
  feature_set <- object@feature_set[[feature_set_idx]]
  exponent <- object@params[["exponent"]]
  df <-
    get_msea_score(feature_list, feature_set, exponent, fortify = TRUE)
  df$ymin = 0
  df$ymax = 0
  pos <- df$position == 1
  h <- diff(range(df$running_score)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$feature_list <- feature_list
  return(df)
}


##' ggplot theme of DOSE
##'
##' @title theme_dose
##' @param font.size font size
##' @return ggplot theme
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @examples
##' library(ggplot2)
##' qplot(1:10) + theme_dose()
##' @export
theme_dose <- function(font.size = 13) {
  theme_bw() +
    theme(
      axis.text.x = element_text(
        colour = "black",
        size = 12,
        vjust = 1
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 12,
        hjust = 1
      ),
      axis.title = element_text(
        margin = margin(10, 5, 0, 0),
        color = "black",
        size = 13
      ),
      axis.title.y = element_text(angle = 90)
    )
}






###----------------------------------------------------------------------------
# metabolite_set <- kegg_hsa_compound_pathway
# feature_list <- feature_list

filter_metabolite_set <- function(metabolite_set,
                                  annotation_table,
                                  min_size = 5,
                                  max_size = 1000) {
  len = length(metabolite_set)
  # len <- lapply(metabolite_set, function(x) {
  #   length(intersect(x, annotation_table$Lab.ID))
  # }) %>%
  #   unlist() %>%
  #   unname()
  
  if (len > min_size & len < max_size) {
    return(metabolite_set)
  } else{
    return(NULL)
  }
  
  # remain_idx <-
  #   which(len > min_size & len < max_size)
  #
  # metabolite_set <- metabolite_set[remain_idx]
  # return(metabolite_set)
}


perm_annotation_table <- function(annotation_table) {
  permute_idx <- sample.int(nrow(annotation_table))
  perm_annotation_table <- annotation_table
  perm_annotation_table$Lab.ID <-
    perm_annotation_table$Lab.ID[permute_idx]
  return(perm_annotation_table)
}
