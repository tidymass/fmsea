#' Annotate Isotopic Peaks for a Given Compound Formula
#'
#' Predicts theoretical isotopic peaks for a given compound formula and adduct, then matches them to observed peaks based on m/z,
#' retention time (RT), and intensity within user-defined tolerances. Supports up to `[M+n]` isotope matching.
#'
#' @param formula Character. The molecular formula (e.g., `"C9H14N2O12P2"`).
#' @param adduct Character. The adduct form (e.g., `"M-H"`, `"M+H"`).
#' @param charge Integer. The ion charge. Default is `1`.
#' @param mz Numeric. The observed m/z of the monoisotopic peak `[M]`.
#' @param rt Numeric. The retention time (in seconds) of the monoisotopic peak.
#' @param mean_intensity Numeric. The mean intensity of the monoisotopic peak.
#' @param peak.mz Numeric vector. Observed m/z values of all detected peaks.
#' @param peak.rt Numeric vector. Observed retention times of all detected peaks.
#' @param peak.int Numeric vector. Observed intensities of all detected peaks.
#' @param rt.tol Numeric. RT tolerance (in seconds) for matching isotopic peaks. Default is `5`.
#' @param mz.tol Numeric. m/z tolerance (in ppm) for matching. Default is `15`.
#' @param int.tol Numeric. Intensity ratio tolerance. Default is `0.3`.
#' @param max.isotope Integer. Maximum number of isotopes to predict (e.g., 3 includes `[M+1]`, `[M+2]`, `[M+3]`). Default is `3`.
#'
#' @return A data frame of matched isotopic peaks with the following columns:
#' \describe{
#'   \item{`peakIndex`}{Index of matched peak in the original peak list}
#'   \item{`mzError.ppm`}{m/z error in ppm}
#'   \item{`isotopes`}{Isotope label (e.g., `[M+1]`, `[M+2]`)}
#'   \item{`rtError.s`}{RT error in seconds}
#' }
#' Returns `NULL` if no `[M+1]` isotopic peak is matched.
#'
#' @details
#' The function uses the `Rdisop` package to predict isotope distribution, then filters observed peaks
#' by RT, m/z, and intensity match to detect isotopic features. It prioritizes matches with minimal intensity error.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' annotate_isotope(
#'   formula = "C9H14N2O12P2",
#'   adduct = "M-H",
#'   mz = 402.9998,
#'   rt = 823.1462,
#'   mean_intensity = 1000,
#'   peak.mz = c(403.0, 404.003, 405.008),
#'   peak.rt = c(823.2, 823.3, 823.5),
#'   peak.int = c(1000, 500, 200),
#'   rt.tol = 5,
#'   mz.tol = 15,
#'   int.tol = 0.3,
#'   max.isotope = 2
#' )
#' }
#'
#' @importFrom Rdisop getMolecule getIsotope
#' @importFrom masstools sum_formula
#' @importFrom purrr map
#' @importFrom dplyr mutate
#'
#' @export

annotate_isotope2 <-
  function(formula = "C9H14N2O12P2",
           adduct = "M-H",
           charge = 1,
           mz = 402.9998,
           rt = 823.1462,
           mean_intensity = 1000,
           ## peak information
           peak.mz = feature_table_pos$mz,
           peak.rt = feature_table_pos$rt,
           peak.int = feature_table_pos$mean_intensity,
           ## other parameters
           rt.tol = 5,
           mz.tol = 15,
           int.tol = 0.3,
           max.isotope = 3) {
    
    formula1 <-
      masstools::sum_formula(formula = formula, adduct = adduct)
    ###should be fix latter
    if (is.na(formula1)) {
      formula1 <- formula
    }
    
    molecule <- Rdisop::getMolecule(
      formula = formula1,
      # z = charge,
      maxisotopes = max.isotope + 1
    )
    
    isotopes <- t(Rdisop::getIsotope(molecule = molecule)[[1]])
    
    rownames(isotopes) <-
      c("[M]", paste("[M", "+", c(1:(nrow(
        isotopes
      ) - 1)), "]", sep = ""))
    
    isotopes <- data.frame(isotopes, rownames(isotopes), stringsAsFactors = FALSE)
    colnames(isotopes) <- c("mz", "intensity", "isotope")
    accurate.mz <- mz
    
    isotopes <- isotopes[-1, , drop = FALSE]
    isotopes$intensity <-
      isotopes$intensity * mean_intensity
    isotopes$rt <- rt
    
    ###rt filtering
    rt.error <- abs(rt - peak.rt)
    index1 <- which(rt.error <= rt.tol)
    
    if (length(index1) == 0) {
      return(NULL)
    }
    
    # iso.info <-
    #   purrr::map(
    #     isotopes,
    #     .f = function(x) {
    #       ###mz error
    #       mz.error <-
    #         abs(as.numeric(x[1]) - peak.mz) * 10^6 / ifelse(as.numeric(x[1]) >= 400, as.numeric(x[1]), 400)
    #       
    #       intensity.error <-
    #         abs(as.numeric(x[2]) - peak.int) / peak.int
    #       
    #       idx <- which(mz.error <= mz.tol &
    #                      rt.error < rt.tol & intensity.error < int.tol)
    #       
    #       if (length(idx) == 0) {
    #         return(NULL)
    #       }
    
    iso.info <- apply(isotopes, 1, function(x) {
      mz.error <- abs(as.numeric(x[1]) - peak.mz) * 1e6 / ifelse(as.numeric(x[1]) >= 400, as.numeric(x[1]), 400)
      intensity.error <- abs(as.numeric(x[2]) - peak.int) / peak.int
      
      idx <- which(mz.error <= mz.tol & rt.error < rt.tol & intensity.error < int.tol)
      
      if (length(idx) == 0) {
        return(NULL)
      }
      
      idx <- idx[which.min(intensity.error[idx])]
      
      iso_table <- data.frame(
        peakIndex = idx,
        mzError.ppm = mz.error[idx],
        isotopes = x[3],
        stringsAsFactors = FALSE
      )
      
      return(iso_table)
      
    }) 
    
    # remove NULL entries
    iso.info <- iso.info[!vapply(iso.info, is.null, logical(1))]
    
    # if no isotopic peaks are matched, return NULL
    if (length(iso.info) == 0) {
      return(NULL)
    }
    
    # merge the list of data frames into a single data frame
    iso.info <- do.call(rbind, iso.info) %>%
      as.data.frame()
    
    # if no [M+1] isotopic peak is matched, return NULL
    if (!"[M+1]" %in% iso.info$isotopes) {
      return(NULL)
    }
    

    iso.info$`rtError.s` <-
      rt.error[as.numeric(iso.info$peakIndex)]
    
    iso.info <-
      iso.info %>%
      dplyr::mutate(
        peakIndex = as.numeric(peakIndex),
        mzError.ppm = as.numeric(mzError.ppm),
        rtError.s = as.numeric(rtError.s)
      )
    
    # if the is no M+2, M+3 will be removed
    if (!"[M+2]" %in% iso.info$isotopes) {
      iso.info <- iso.info[1, ]
    }
    
    # rm(list = c("peak.mz", "peak.rt", "peak.int", "cor"))
    # gc()
    return(iso.info)
  }
