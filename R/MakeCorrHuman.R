#' @title MakeCorrHuman
#' @description Compute correlations between tissue dependent intensities of proteins and the tissue dependent number of MHCI peptides
#' @param df_human Input MHCI tissue draft from GetHumanMHCIdata
#' @param Donors Vector with names of donors for which correlations are to be computed, Default: c("all")
#' @param deconvolute_byHLAGene Logical, can be used to deconvolute by allele instead of by donor (not recommended), Default: FALSE
#' @param pValue_Threshold Maximum pValue at which a correlation will be considered significant, Default: 0.05
#' @param runAnalCorrHuman Logical, directly runs the function AnalCorrHuman to retrieve significantly correlating proteins, Default: TRUE
#' @return List of dataframes with results
#' @details Running this function for all donors might take a while (around 15 minutes). Default pValue threshold is higher than for mouse since by default a pValue smaller the set value must be found in at least two donors/alleles
#' Human protein expression data were retrieved from: Wang et al.; Mol Syst Biol. 2019 Feb 18;15(2); PMID: 30777892 https://www.embopress.org/doi/full/10.15252/msb.20188503
#' @examples
#' \dontrun{
#' corrs_h<-MakeCorrHuman(df_human,Donors=c('AUT01-DN11'))
#' }
#' @seealso
#'  \code{\link[dplyr]{group_by}},\code{\link[dplyr]{summarise}},\code{\link[dplyr]{setops}}
#'  \code{\link[stats]{lm}}
#'  \code{\link[broom]{reexports}}
#' @rdname MakeCorrHuman
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr group_by summarise intersect
#' @importFrom stats lm
#' @importFrom broom glance
MakeCorrHuman<- function(df_human,Donors=c('all'),deconvolute_byHLAGene=FALSE,pValue_Threshold=0.05,runAnalCorrHuman=TRUE){
  #Internal Functions:
  `%>%` <- magrittr::`%>%`
  getIndices <- function(dataframe,string){ which(grepl(string,tolower(names(dataframe))))}
  #Load human protein data (Downloaded from * Wang et al.; Mol Syst Biol. 2019 Feb 18;15(2); PMID: 30777892 *, Suppl. materials -> Table_EV1.xlsx -> Tab 'Genes')
  cat('Using human protein data: \nDownloaded from * Wang et al.; Mol Syst Biol. 2019 Feb 18;15(2); PMID: 30777892 * \nhttps://www.embopress.org/doi/full/10.15252/msb.20188503 -> Supporting Information -> Table_EV1.xlsx -> Tab labelled "Genes"')
  file<- system.file("extdata", "Table_EV1_Genes.csv", package = "MHCIatlas")
  df_proteome<-read.csv(file)
  df.proteome<- cbind(df_proteome[2:31])
  colnames(df.proteome)<- c('Gene_names','AdrenalGland','Appendix','Brain','Colon','Duodenum','Endometrium','Esophargus','FallopianTube','Fat','Gallbladder','Heart','Kidney','Liver','Lung','LymphNode','Ovary','Pancreas','Placenta','Prostate','Rectum','SalivaryGland','SmallIntestine','Muscle','Spleen','Stomach','Testis','Thyroid','Tonsil','Bladder')
  df.proteome[df.proteome==0]<- NA
  df.proteome$RankN<-rowSums(!is.na(df.proteome[,2:30]))
  df.proteome<- subset(df.proteome,RankN>9)
  df.proteome<- df.proteome[order(df.proteome$Gene_names,df.proteome$RankN,decreasing = T),]
  df.proteome <- df.proteome[!duplicated(df.proteome$Gene_names),]
  df.proteomel<-reshape(df.proteome,idvar="Gene_names",direction="long",varying=c(list(2:30)),times=c(names(df.proteome)[2:30]),timevar="Tissue",v.names=c("Int"))
  df_prots<- df.proteome[1:31]

  cat('\nUsing the Human MHCI tissue atlas in the form of user input \'',deparse(substitute(df_human)),'\' \n' )
  if(deconvolute_byHLAGene==T){
    cat('Deconvoluting by HLA Gene (HLA-A, HLA-B, HLA-C)\n')
    df_human$Tissue_Donor<- paste0(df_human$Tissue,'_',df_human$Donor)
    dfh<- cbind(df_human[c('sequence','Tissue','Tissue_Donor')],substring(df_human[,'Best_Patient.Allele'],1,1) )
    colnames(dfh)<- c('sequence','Tissue','Tissue_Donor','HLA_Gene')
    dfhg<- dfh %>%
      dplyr::group_by(HLA_Gene, Tissue_Donor) %>%
      dplyr::summarise(n = n())
    dfhg<- data.frame(dfhg)
    dfhg$Tissue<- do.call(rbind,strsplit(dfhg$Tissue_Donor,'_'))[,1]
    dfhg<- dfhg[c(1,3,4)]
    dfhg<- dfhg %>%
      dplyr::group_by(HLA_Gene,Tissue) %>%
      dplyr::summarise(Mean = mean(n, na.rm=TRUE),sd = sd(n, na.rm=TRUE))
    dfhg<- data.frame(dfhg)
    dfhg<- dfhg[order(dfhg$HLA_Gene,dfhg$Mean,decreasing = F),]
    dfhg$Tissue<- factor(dfhg$Tissue, levels = unique(dfhg$Tissue) )
    colnames(dfhg)<- c('var','Tissue','mean_n','sd')
    Donors <- c('all')
  }else{
    cat('Deconvoluting by Donors\n')
    dfh<- df_human[c('sequence','Tissue','Donor')]
    colnames(dfh)<- c('sequence','Tissue','Donor')
    dfhg<- dfh %>%
      dplyr::group_by(Donor, Tissue) %>%
      dplyr::summarise(n = n())
    dfhg<- data.frame(dfhg)
    colnames(dfhg)<- c('var','Tissue','n')
  }

  Variables_Main<- as.character(unique(dfhg[,1]))

  if(Donors[1] =='all'){
    Variables_Main<- as.character(unique(dfhg[,1]))
  }else{
    Variables_Main<- as.character(Donors)
  }

  for (j in 1:length(Variables_Main)){
    var_main<- as.character(Variables_Main[j])
    peptide.counts<- subset(dfhg,var==var_main)[2:3]
    df_prots[,paste0("Corr_",var_main)]<- NA
    df_prots[,paste0("Slope_",var_main)]<- NA
    df_prots[,paste0("pValue_",var_main)]<- NA
    df_prots[,paste0("n_",var_main)]<- NA

    as.character(subset(dfhg,var==var_main)$Tissue)

    if(   length(dplyr::intersect(unique(df.proteomel$Tissue),as.character(subset(dfhg,var==var_main)$Tissue)))   <10){
      next
    }
    for (i in 1:length(df.proteome$Gene_names)) {
      prot<- subset(df.proteomel,Gene_names==df_prots[i,1])[3:4]
      na.omit(as.matrix(merge(peptide.counts,prot,by='Tissue')[2:3]))->pro
      log10(pro[,2])->y
      pro[,1]->x
      if(nrow(pro)>9){
        df_prots[i,paste0("Corr_",var_main)]<- summary(stats::lm(y~x))[[8]]
        df_prots[i,paste0("Slope_",var_main)]<-summary(stats::lm(y~x))[[4]][2,1]
        df_prots[i,paste0("pValue_",var_main)]<-signif(data.frame(broom::glance(stats::lm(y~x))[5])[1,1],digits = 3)
        df_prots[i,paste0("n_",var_main)]<- nrow(pro)
        }
      else{
        df_prots[i,paste0("Corr_",var_main)]<- NA
        df_prots[i,paste0("Slope_",var_main)]<- NA
        df_prots[i,paste0("pValue_",var_main)]<- NA
        df_prots[i,paste0("n_",var_main)]<- NA
        }
      ProgressBar(i,len_i = length(df.proteome$Gene_names),j, len_j = length(Variables_Main) )
    }

  }
  ProgressBar(i=100)
  cat('\n')

  if(runAnalCorrHuman==T){
    df_protcorr<- AnalCorrHuman(list(df_prots,dfhg),pValue_Threshold = pValue_Threshold)
    l.toreturn<- list(df_protcorr,dfhg,df.proteome)
    names(l.toreturn)<- c('HumanProtCorr','PeptideTissueCounts','Table_EV1_Wang2019')
  return( l.toreturn )}else{
    l.toreturn<- list(df_prots,dfhg,df.proteome)
    names(l.toreturn)<- c('HumanProtCorrRaw','PeptideTissueCounts','Table_EV1_Wang2019')
    return( l.toreturn )
  }
}
