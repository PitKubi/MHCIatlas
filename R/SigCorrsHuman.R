#' @title SigCorrsHuman
#' @description Retrieve significantly correlating proteins and their gene names
#' @param Human_ProtCorrelations Output list from function MakeCorrHuman
#' @param RankSig_Threshold Min number of occurrences where correlation was found significant, Default: 2
#' @return Data frame containing significantly correlating genes and the mean correlation pValue
#' @details Used to get an overview of singificantly proteins that can be used for later gene ontology analysis etc.
#' @examples
#' \dontrun{
#' sig_prots_h<-SigCorrsHuman(corrs_h)
#' }
#' @rdname SigCorrsHuman
#' @export
SigCorrsHuman<- function(Human_ProtCorrelations,RankSig_Threshold=2){
  return( subset(Human_ProtCorrelations[[1]],RankSig>=RankSig_Threshold)[c('Gene_names','MeanCorr_pValue')]   )
}
