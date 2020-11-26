#' @title AnalCorrHuman
#' @description Function to extract proteins which significanlty correlate with MHCI peptide counts across tissues and donors.
#' @param DfList_Human_ProteinCorrelations Output from 'MakeCorrHuman' with parameter 'runAnalCorrHuman=FALSE'
#' @param pValue_Threshold Protein correlations with a pValue smaller than this threshold are considered significant, Default: 0.05
#' @return Data frame with all human proteins subjected to correlations and annotated significance across donors
#' @details This function should always be used as part of the function 'MakeCorrHuman', only to be run seperately to save on computation time/space
#' @examples
#' \dontrun{
#'  df_corrhuman <- AnalCorrHuman(MakeCorrHumanRaw,pValue_Threshold = 0.05)
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{scale_colour_brewer}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{facet_grid}},\code{\link[ggplot2]{lims}}
#' @rdname AnalCorrHuman
#' @export
#' @importFrom ggplot2 ggplot aes geom_line scale_color_brewer geom_vline theme element_blank element_rect labs facet_grid element_text xlab ylab ylim
AnalCorrHuman<- function(DfList_Human_ProteinCorrelations,pValue_Threshold = 0.05){
  #Internal Functions:
  getIndices <- function(dataframe,string){ which(grepl(string,tolower(names(dataframe))))}
  #Code
  df_prots<- DfList_Human_ProteinCorrelations[[1]]
  Variables_Main<- unique(as.character(DfList_Human_ProteinCorrelations[[2]][,1]))
  df_corr<-cbind(df_prots[1],df_prots[getIndices(df_prots, 'corr_' )[1]:ncol(df_prots)])
  df_corr[,order(colnames(df_corr))]->df_corr
  indices_corr<- getIndices(df_corr,'corr')
  indices_n<- getIndices(df_corr,'n_')
  indices_pvalues<- getIndices(df_corr,'pvalue')
  indices_slopes<- getIndices(df_corr,'slope')
  df_corr<- cbind(df_corr['Gene_names'],df_corr[c(indices_corr)],df_corr[c(indices_n)],df_corr[c(indices_slopes)],df_corr[c(indices_pvalues)])
  indices_corr<- getIndices(df_corr,'corr_')
  indices_n<- getIndices(df_corr,'n_')
  indices_pvalues<- getIndices(df_corr,'pvalue')
  indices_slopes<- getIndices(df_corr,'slope')
  df_corr$MaxCorr<- suppressWarnings(apply(df_corr[c(indices_corr)],1, function(x) max(x,na.rm=T)  ))
  df_corr$MaxCorr[df_corr$MaxCorr=='-Inf']<- NA

  df_corr$MeanCorr_SigVars<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'MeanCorr_SigVars'] <- tryCatch(mean(as.numeric(df_corr[i,which(df_corr[i,c(indices_pvalues)] <= pValue_Threshold)+(indices_corr[1]-1)]),na.rm=T),error=function(e) {NA} )
  }
  df_corr$MeanCorr_pValue<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'MeanCorr_pValue'] <- suppressWarnings(ifelse( is.na(mean(as.numeric(df_corr[i,c(indices_pvalues)]),na.rm=T)) == T,NA,  tryCatch(mean(as.numeric(df_corr[i,which(df_corr[i,c(indices_pvalues)] <= pValue_Threshold)+(indices_pvalues[1]-1)]),na.rm=T),error=function(e) {NA} )) )
  }
  df_corr$MeanCorr_n<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'MeanCorr_n'] <- suppressWarnings(ifelse( is.na(mean(as.numeric(df_corr[i,c(indices_pvalues)]),na.rm=T)) == T,NA, tryCatch(mean(as.numeric(df_corr[i,which(df_corr[i,c(indices_pvalues)] <= pValue_Threshold)+(indices_n[1]-1)]),na.rm=T),error=function(e) {NA} )) )
  }
  df_corr$MeanCorr_slope<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'MeanCorr_slope'] <- suppressWarnings(ifelse( is.na(mean(as.numeric(df_corr[i,c(indices_pvalues)]),na.rm=T)) == T,NA,  tryCatch(mean(as.numeric(df_corr[i,which(df_corr[i,c(indices_pvalues)] <= pValue_Threshold)+(indices_slopes[1]-1)]),na.rm=T),error=function(e) {NA} )) )
  }
  df_corr$SigCorr_Vars<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'SigCorr_Vars'] <-  tryCatch(as.character(paste(do.call(rbind,(strsplit(names(df_corr)[which(df_corr[i,c(indices_pvalues)] <= pValue_Threshold)+1],"_")))[,2],sep='_',collapse='_') ),error=function(e) {NA} )
  }

  df_corr$RankSig<- apply(df_corr[c(indices_pvalues)],1,function(x) length(which(x<=pValue_Threshold)) )

  df_corr$Corr_BestVar<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'Corr_BestVar'] <- suppressWarnings( ifelse( is.na(mean(as.numeric(df_corr[i,c(indices_pvalues)]),na.rm=T)) == T,NA,tryCatch(as.numeric(df_corr[i, as.numeric(which.min(df_corr[i,c(indices_pvalues)])) +(indices_corr[1]-1) ]),error=function(e) {NA} )  ) )
  }
  df_corr$pValue_BestVar<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'pValue_BestVar'] <- suppressWarnings( ifelse( is.na(mean(as.numeric(df_corr[i,c(indices_pvalues)]),na.rm=T)) == T,NA,tryCatch(as.numeric(df_corr[i, as.numeric(which.min(df_corr[i,c(indices_pvalues)])) +(indices_pvalues[1]-1) ]),error=function(e) {NA} )  ))
  }
  df_corr$n_BestVar<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'n_BestVar'] <-suppressWarnings( ifelse( is.na(mean(as.numeric(df_corr[i,c(indices_pvalues)]),na.rm=T)) == T,NA,tryCatch(as.numeric(df_corr[i, as.numeric(which.min(df_corr[i,c(indices_pvalues)])) +(indices_n[1]-1) ]),error=function(e) {NA} )  ))
  }
  df_corr$Slope_BestVar<- NA
  for (i in 1:nrow(df_corr)){
    df_corr[i,'Slope_BestVar'] <-suppressWarnings( ifelse( is.na(mean(as.numeric(df_corr[i,c(indices_pvalues)]),na.rm=T)) == T,NA,  tryCatch(as.numeric(df_corr[i, as.numeric(which.min(df_corr[i,c(indices_pvalues)]))+(indices_slopes[1]-1)]) ,error=function(e) {NA} )))
  }

  df_corr$Title<- "Protein intensities vs. # HLAI peptides/tissue (Human)"
  plot1<- ggplot2::ggplot(subset(df_corr,!is.na(n_BestVar) ),ggplot2::aes(x=pValue_BestVar,y=Corr_BestVar,color=factor(n_BestVar)))+ggplot2::geom_line(size=0.7)+ggplot2::scale_color_brewer(palette='Set1')+ggplot2::geom_vline(xintercept = pValue_Threshold)+
    ggplot2::theme(legend.position = c(0.8,0.7),panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill=NA))+ ggplot2::labs(color = "# of datapoints")+ ggplot2::facet_grid(. ~ Title) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill='darkgrey'),strip.text = ggplot2::element_text(size=15, colour="white"))+ggplot2::xlab('pValue')+ggplot2::ylab('R-squared')+ggplot2::ylim(0,0.85)
  tryCatch( print(plot1),error=function(e) {'No significant correlations found'} )
  ind_maxcorr<- getIndices(df_corr,'maxcorr')
  df_proteome_correllations<- cbind(df_prots,df_corr[ind_maxcorr:ncol(df_corr)])
  return(df_proteome_correllations)
}
