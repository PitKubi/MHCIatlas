#' @title SigCorrsMouse
#' @description Retrieve significantly correlating proteins and their gene names
#' @param Mouse_ProtCorrelations Output list from function MakeCorrMouse
#' @param pValue_Threshold maximum pValue at or below which a correlation is considered significant, Default: 0.01
#' @return Data frame containing significantly correlating genes and the mean correlation pValue
#' @details Used to get an overview of singificantly proteins that can be used for later gene ontology analysis etc.
#' @examples
#' \dontrun{
#' sig_prots_m<-SigCorrsHuman(corrs_m)
#' }
#' @rdname SigCorrsMouse
#' @export
SigCorrsMouse<- function(Mouse_ProtCorrelations, pValue_Threshold = 0.01){
  df.proteome<- Mouse_ProtCorrelations[[1]]
  prots_sigm<- subset(df.proteome,Corr.pValue<= pValue_Threshold)[c(19,20)]
  prots_sigm$Gene_names<- apply(prots_sigm[c(1:2)],1,function(x) ifelse(is.na(x[1])==T,x[2],x[1]))
  row.names(prots_sigm)<- seq(1,nrow(prots_sigm),1)
  return(prots_sigm)
}
