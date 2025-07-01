###MS2 matching score 150
###positive
### [M+H] 50
### [M+H] isotope 20
### other adduct 20
### other adduct isotope 10

###negative
### [M-H] 50
### [M-H] isotope 20
### other adduct 20
### other adduct isotope 10

####The max score should be 500
#####The min score is 20

#' @title score_mfc
#' @description score_mfc
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param mfc mfc
#' @export
score_mfc <- function(mfc) {
  score <- 0
  
  ###
  if (any(!is.na(mfc$SS))) {
    score <-
      score + 150
  }
  
  ##This should be optimized in the future
  ###if positive adduct are +H, score add 50
  if (any(mfc$Adduct == "(M+H)+")) {
    score <- score + 50
  }
  
  if (any(mfc$Adduct == "(M+H)+" &
          mfc$isotope != "[M]")) {
    score <- score + 20
  }
  
  if (any(mfc$polarity == "positive" &
          mfc$Adduct != "(M+H)+")) {
    score <- score + 20
  }
  
  if (any(mfc$polarity == "positive" & mfc$Adduct != "(M+H)+" &
          mfc$isotope != "[M]")) {
    score <- score + 10
  }
  
  if (any(mfc$Adduct == "(M-H)-")) {
    score <- score + 50
  }
  
  if (any(mfc$Adduct == "(M-H)-" &
          mfc$isotope != "[M]")) {
    score <- score + 20
  }
  
  if (any(mfc$polarity == "negative" &
          mfc$Adduct != "(M-H)-")) {
    score <- score + 20
  }
  
  if (any(mfc$polarity == "negative" & mfc$Adduct != "(M-H)-" &
          mfc$isotope != "[M]")) {
    score <- score + 10
  }
  
  return(score)
}