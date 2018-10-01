#' Title
#'
#' @param ThetaEst Estimated adjacency matrix
#' @param Theta  True adjacency matrix
#'
#' @return
#' @export
#'
#' @examples
Diagnostic = function(ThetaEst,Theta){

  # Kuismin, M., & Sillanpää, M. J. (2016). Use of Wishart prior and simple extensions for
  # sparse precision matrix estimation. PloS one, 11(2), e0148171.

  TN = ifelse(Theta[upper.tri(Theta)] == 0 & ThetaEst[upper.tri(ThetaEst)] == 0, 1, 0); TN = sum(TN) # True Negative

  FP = ifelse(Theta[upper.tri(Theta)] == 0 & ThetaEst[upper.tri(ThetaEst)] != 0, 1, 0); FP = sum(FP) # False Positive

  TP = ifelse(Theta[upper.tri(Theta)] != 0 & ThetaEst[upper.tri(ThetaEst)] != 0, 1, 0); TP = sum(TP) # True Positive

  FN = ifelse(Theta[upper.tri(Theta)] != 0 & ThetaEst[upper.tri(ThetaEst)] == 0, 1, 0); FN = sum(FN) # False Negatives

  Specificity = TN/(TN + FP) # aka True Negative rate

  Sensitivity = TP/(TP + FN) # aka True Positive rate

  Fallout = FP/(TN + FP) # aka False Positive rate

  Precision = TP/(TP + FP) # aka Positive predictive value

  MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  lista = list("Sp" = Specificity, "Sen" = Sensitivity, "Fall" = Fallout, "Pre" = Precision ,"MCC" = MCC)

  return(lista)

}
