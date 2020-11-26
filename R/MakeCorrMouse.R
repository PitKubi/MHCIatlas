#' @title MakeCorrMouse
#' @description Compute correlations between tissue dependent intensities of proteins and the tissue dependent number of MHCI peptides
#' @param df_mouse Input MHCI tissue draft from GetMouseMHCIdata
#' @param pValue_Threshold Maximum pValue at which a correlation will be considered significant, Default: 0.01
#' @param useSILAC Logical, use SILAC fold change values instead of normalized raw intensities of protein expression data, Default: F
#' @return List of dataframes with results
#' @details Mouse protein expression data are retireved from this publication: Geiger et al. Molecular & Cellular Proteomics June 1, 2013, 12 (6) 1709-1722; PMID: 30777892 https://www.mcponline.org/content/12/6/1709
#' @examples
#' \dontrun{
#' corrs_m<-MakeCorrMouse(df_mouse)
#' }
#' @seealso
#'  \code{\link[dplyr]{tally}}
#'  \code{\link[stats]{lm}}
#'  \code{\link[broom]{reexports}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{scale_colour_brewer}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{facet_grid}},\code{\link[ggplot2]{lims}}
#' @rdname MakeCorrMouse
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr count
#' @importFrom stats lm
#' @importFrom broom glance
#' @importFrom ggplot2 ggplot aes geom_line scale_color_brewer geom_vline theme element_blank element_rect labs facet_grid element_text xlab ylab ylim
MakeCorrMouse<- function(df_mouse,pValue_Threshold=0.01,useSILAC=F){
  #Internal functions:
  `%>%` <- magrittr::`%>%`
  capstr<- function(x){paste(toupper(substring(x, 1,1)), substring(x, 2), sep="", collapse=" ")}
  #Code main function:
  cat('Using mouse protein data: \nDownloaded from * Geiger et al. Molecular & Cellular Proteomics June 1, 2013, 12 (6) 1709-1722; PMID: 30777892 * \nhttps://www.mcponline.org/content/12/6/1709 -> Figures & Data -> Suppl Table 1 - Protein table  -> Tab labelled "Suppl Table S1\n"')
  file<- system.file("extdata", "Suppl_Table_S1_all_proteins_Geiger_et_al_2013.txt", package = "MHCIatlas")
  table_s1_geiger <- read.csv(file,sep='\t')
  df.proteome<- table_s1_geiger
  if(useSILAC==T){
    df.proteome<- cbind(df.proteome[2:4],df.proteome[9],df.proteome[12],df.proteome[17],df.proteome[22],df.proteome[24:25],df.proteome[27:29],df.proteome[32:33],df.proteome[35:38])
    df.proteome[5:18][df.proteome[5:18]==0]<- NA
    df.proteome[5:18]<- log2(df.proteome[5:18])
    }else{
    df.proteome<- cbind(df.proteome[2:4],df.proteome[9],df.proteome[40],df.proteome[45],df.proteome[50],df.proteome[52:53],df.proteome[55:57],df.proteome[60:61],df.proteome[63:66])
    df.proteome[5:18][df.proteome[5:18]==0]<- NA
    df.proteome[5:18]<- log10(df.proteome[5:18])
      }

  colnames(df.proteome)<- c(names(df.proteome)[1:3],'MW',"adrenal",  "colon",    "heart"  ,'smallintestine' , "kidney" ,  "liver"  ,  "lung"   ,  "brain", "ovary"  ,  "pancreas" ,"spleen"  ,"stomach",  "thymus"  , "uterus"  )
  df.proteome$Gene_names<-apply(df.proteome[3],1,function(x) do.call(rbind,strsplit(x, split=';', fixed=TRUE))[1] )
  df.proteome$UniprotAcc<- apply(df.proteome[1],1,function(x) do.call(rbind,strsplit(x, split=';', fixed=TRUE))[1] )
  df.proteome$RankN<-rowSums(!is.na(df.proteome[,5:18]))
  df.proteome<- subset(df.proteome,RankN>9)

  if(useSILAC!=T){
    #Normalize raw intensities (If  no SILAC data are used)
    prot.mean<-mean(as.matrix(df.proteome[,5:18]),na.rm=T)
    for(i in 5:18){
      v<-NULL
      v2<-NULL
      v<-df.proteome[,i]
      v2 =  prot.mean- mean(as.matrix(df.proteome[,i]),na.rm=T)
      df.proteome[,i]<- v2+v
    }
  }


  mhc.counts <- data.frame(df_mouse['Tissue'] %>% dplyr::count(Tissue))
  df.proteome$Corr.rsq<- NA
  df.proteome$Corr.slope<- NA
  df.proteome$Corr.n<- NA
  df.proteome$Corr.pValue<- NA
  df.protl<-reshape(df.proteome,idvar="UniprotAcc",direction="long",varying=c(list(5:18)),times=c(names(df.proteome)[5:18]),timevar="Tissue",v.names=c("Int"))
  df.protl$Tissue<- strsplit(capstr(df.protl$Tissue),' ')[[1]]

  for (i in 1:length(df.proteome$UniprotAcc)) {
    prot<- subset(df.protl,UniprotAcc==df.proteome[i,20])[12:13]
    na.omit(as.matrix(merge(mhc.counts,prot,by='Tissue')[2:3]))->pro
    pro[,2]->y
    pro[,1]->x
    if(nrow(pro)>9){
      df.proteome[i,paste0('Corr',".rsq")]<- summary(stats::lm(y~x))[[8]]
      df.proteome[i,paste0('Corr',".slope")]<-summary(stats::lm(y~x))[[4]][2,1]
      df.proteome[i,paste0('Corr',".n")]<- nrow(pro)
      df.proteome[i,paste0('Corr',".pValue")]<-signif(data.frame(broom::glance(stats::lm(y~x))[5])[1,1],digits = 3)}else{
        df.proteome[i,paste0('Corr',".rsq")]<- NA
        df.proteome[i,paste0('Corr',".slope")]<-NA
        df.proteome[i,paste0('Corr',".n")]<- NA
        df.proteome[i,paste0('Corr',".pValue")]<-NA
      }
    ProgressBar(i,len_i=length(df.proteome$UniprotAcc) )
  }
  df.proteome$SigCorr<- apply(df.proteome['Corr.pValue'],1,function(x) {ifelse(x<= pValue_Threshold,1,NA)} )
  if(useSILAC==T){
    df.proteome$Title<- "Protein Super-SILAC H/L vs. # MHCI peptides/tissue (Mouse)"
  }else{
    df.proteome$Title<- "Protein intensities vs. # MHCI peptides/tissue (Mouse)"
  }

  df.proteome.m.annot<- df.proteome
  plot1<-ggplot2::ggplot(subset(df.proteome.m.annot,!is.na(Corr.n) ),ggplot2::aes(x=Corr.pValue,y=Corr.rsq,color=factor(Corr.n)))+ggplot2::geom_line(size=0.7)+ggplot2::scale_color_brewer(palette='Set1')+ggplot2::geom_vline(xintercept = pValue_Threshold)+
    ggplot2::theme(legend.position = c(0.8,0.8),panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill=NA))+ ggplot2::labs(color = "# of datapoints")+ ggplot2::facet_grid(. ~ Title) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill='darkgrey'),
          strip.text = ggplot2::element_text(size=15, colour="white"))+ggplot2::xlab('pValue')+ggplot2::ylab('R-squared')+ggplot2::ylim(0,0.85)
  suppressWarnings(print(plot1))
  l.return <- list(df.proteome,mhc.counts,table_s1_geiger,pValue_Threshold,useSILAC)
  names(l.return)<- c('MouseProteomeCorr','PeptideTissueCounts','Table_S1_Geiger2013','pValue_Threshold','useSILAC')
  return( l.return )
}
