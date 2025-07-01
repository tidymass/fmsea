annotate_isotope <-
  function(formula = "C9H14N2O12P2",
           adduct = "M-H",
           charge = 1,
           mz = 402.9998,
           rt = 823.1462,
           mean_intensity = 1000,
           ## peak information
           peak.mz,
           peak.rt,
           peak.int,
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
    
    isotopes <- t(Rdisop::getIsotope(molecule = molecule))
    
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
    
    iso.info <-
      purrr::map(
        as.data.frame(t(isotopes)),
        .f = function(x) {
          ###mz error
          mz.error <-
            abs(as.numeric(x[1]) - peak.mz) * 10^6 / ifelse(as.numeric(x[1]) >= 400, as.numeric(x[1]), 400)
          
          intensity.error <-
            abs(as.numeric(x[2]) - peak.int) / peak.int
          
          idx <- which(mz.error <= mz.tol &
                         rt.error < rt.tol & intensity.error < int.tol)
          
          if (length(idx) == 0) {
            return(NULL)
          }
          
          ## more than one matched, see the intensity ratio error
          if (length(idx) > 0) {
            idx <- idx[which.min(intensity.error[idx])]
            return(c(idx, mz.error[idx], x[3]))
          }
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    
    if (nrow(iso.info) == 0) {
      return(NULL)
    }
    
    if (!"[M+1]" %in% iso.info$V3) {
      return(NULL)
    }
    
    colnames(iso.info) <- c("peakIndex", "mzError.ppm", "isotopes")
    
    if (!"[M+2]" %in% iso.info$isotopes) {
      iso.info <- iso.info[1, , drop = FALSE]
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
    
    if (!"[M+2]" %in% iso.info$isotopes) {
      iso.info <- iso.info[1, ]
    }
    
    # rm(list = c("peak.mz", "peak.rt", "peak.int", "cor"))
    # gc()
    return(iso.info)
  }
