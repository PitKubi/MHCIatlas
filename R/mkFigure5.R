#' @title mkFigure5
#' @description Function to plot Figure 4 of the manuscript
#' @param df_human Human MHCI tissue draft data frame from GetHumanMHCI data
#' @param HumanHousekeepers List output from function HousekeepersHuman
#' @param df_mouse Mouse MHCI tissue draft data frame from GetMouseMHCI data
#' @param useDefaultCons Uses conservatio scores calculated with defulat values as in manuscript, Default: TRUE
#' @param ConsMouse Output list from ConservationMouse (Needed if useDefaultCons is set to FLASE), Default: NA
#' @param ConsHuman Output list from ConservationHuman (Needed if useDefaultCons is set to FLASE), Default: NA
#' @return List with ggplot2 objects and data to generate the figure
#' @details DETAILS
#' @examples
#' \dontrun{
#' F5<- mkFigure5(df_human,cons_h,df_mouse)
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{scale_colour_brewer}},\code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{facet_grid}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{geom_bar}},\code{\link[ggplot2]{lims}},\code{\link[ggplot2]{geom_boxplot}},\code{\link[ggplot2]{ggtheme}}
#'  \code{\link[cowplot]{ggdraw}},\code{\link[cowplot]{draw_plot}},\code{\link[cowplot]{draw_plot_label}}
#'  \code{\link[ggplotify]{as.ggplot}}
#'  \code{\link[org.Mm.eg.db]{org.Mm.eg.db}}
#' @rdname mkFigure5
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_brewer geom_line facet_grid theme element_rect element_text element_blank labs geom_bar xlim geom_boxplot theme_void
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
#' @importFrom ggplotify as.ggplot
#' @importFrom AnnotationDbi select
#' @importFrom org.Mm.eg.db org.Mm.eg.db
mkFigure5<- function(df_human,HumanHousekeepers,df_mouse,useDefaultCons=TRUE,ConsMouse=NA,ConsHuman=NA){
  #Internal Functions:
  getIndices <- function(dataframe,string){ which(grepl(string,tolower(names(dataframe))))}
  #Main:
  #Generate panels C and D of Figure 5:
  if(useDefaultCons==T){
    filem<- system.file("extdata", "mouse_conservation_default.csv", package = "MHCIatlas")
    df_cons <- read.csv(filem,header=T)
    ptmp<-ggplot2::ggplot(df_cons,ggplot2::aes(x=as.numeric(phastCons),y=relative,color=factor(GeneSet)))+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::geom_line()+ggplot2::facet_grid(~Type)+ggplot2::theme(panel.border = ggplot2::element_rect(fill=NA),legend.position = c(0.8, 0.3),axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_rect(fill = 'white'))+
      ggplot2::labs(color = "Mouse Genes:",x='phastCons',y='Cumulative Frequency')
    pval_exon<-paste('Wilcoxon rank sum test, p-value: ',signif(wilcox.test(subset(df_cons,Type=='exons'&GeneSet=='Tissue Specific')$relative,subset(df_cons,Type=='exons'&GeneSet=='Universal')$relative)[[3]],digits=4))
    pval_prom<-paste('Wilcoxon rank sum test, p-value: ',signif(wilcox.test(subset(df_cons,Type=='promoters'&GeneSet=='Tissue Specific')$relative,subset(df_cons,Type=='promoters'&GeneSet=='Universal')$relative)[[3]],digits=4))
    p3<- cowplot::ggdraw() +cowplot::draw_plot(ggplotify::as.ggplot(ptmp), x = 0, y = 0, width = 1, height = 1)+cowplot::draw_plot_label(label = c(pval_exon, pval_prom), size = 9,x = c(0.0,0.45), y = c(0.9, 0.9))
    # p3 in conservation analysis mouse
    fileh<- system.file("extdata", "human_conservation_default.csv", package = "MHCIatlas")
    df_cons <- read.csv(fileh,header=T)
    ptmp<-ggplot2::ggplot(df_cons,ggplot2::aes(x=as.numeric(phastCons),y=relative,color=factor(GeneSet)))+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::geom_line()+ggplot2::facet_grid(~Type)+ggplot2::theme(panel.border = ggplot2::element_rect(fill=NA),legend.position = c(0.8, 0.3),axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_rect(fill = 'white'))+
      ggplot2::labs(color = "Human Genes:",x='phastCons',y='Cumulative Frequency')
    pval_exon<-paste('Wilcoxon rank sum test, p-value: ',signif(wilcox.test(subset(df_cons,Type=='exons'&GeneSet=='Tissue Specific')$relative,subset(df_cons,Type=='exons'&GeneSet=='Universal')$relative)[[3]],digits=4))
    pval_prom<-paste('Wilcoxon rank sum test, p-value: ',signif(wilcox.test(subset(df_cons,Type=='promoters'&GeneSet=='Tissue Specific')$relative,subset(df_cons,Type=='promoters'&GeneSet=='Universal')$relative)[[3]],digits=4))
    p4<- cowplot::ggdraw() +cowplot::draw_plot(ggplotify::as.ggplot(ptmp), x = 0, y = 0, width = 1, height = 1)+cowplot::draw_plot_label(label = c(pval_exon, pval_prom), size = 9,x = c(0.0,0.45), y = c(0.9, 0.9))
    # p3 in conservation analysis human
  }else{
    if(is.list(ConsMouse)==T){df_cons<- data.frame(ConsMouse[[1]]) }else{df_cons<- ConsMouse}
    ptmp<-ggplot2::ggplot(df_cons,ggplot2::aes(x=as.numeric(phastCons),y=relative,color=factor(GeneSet)))+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::geom_line()+ggplot2::facet_grid(~Type)+ggplot2::theme(panel.border = ggplot2::element_rect(fill=NA),legend.position = c(0.8, 0.3),axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_rect(fill = 'white'))+
      ggplot2::labs(color = "Mouse Genes:",x='phastCons',y='Cumulative Frequency')
    p3<- cowplot::ggdraw() +cowplot::draw_plot(ggplotify::as.ggplot(ptmp), x = 0, y = 0, width = 1, height = 1)+cowplot::draw_plot_label(label = c(pval_exon, pval_prom), size = 9,x = c(0.0,0.45), y = c(0.9, 0.9))
    # p3 in conservation analysis mouse
    if(is.list(ConsHuman)==T){df_cons<- data.frame(ConsHuman[[1]]) }else{df_cons<- ConsHuman}
    ptmp<-ggplot2::ggplot(df_cons,ggplot2::aes(x=as.numeric(phastCons),y=relative,color=factor(GeneSet)))+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::geom_line()+ggplot2::facet_grid(~Type)+ggplot2::theme(panel.border = ggplot2::element_rect(fill=NA),legend.position = c(0.8, 0.3),axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_rect(fill = 'white'))+
      ggplot2::labs(color = "Human Genes:",x='phastCons',y='Cumulative Frequency')
    p4<- cowplot::ggdraw() +cowplot::draw_plot(ggplotify::as.ggplot(ptmp), x = 0, y = 0, width = 1, height = 1)+cowplot::draw_plot_label(label = c(pval_exon, pval_prom), size = 9,x = c(0.0,0.45), y = c(0.9, 0.9))
    # p3 in conservation analysis human
  }

  ProgressBar(i=1,len_i=10)
  #Generate panels A and B of Figure 5:
  df_allgenes<- HumanHousekeepers[[2]]
  df_allgenes$Housekeeping<-1
  file <- system.file("extdata", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", package = "MHCIatlas")
  df_mRNA<-read.csv(file,header=T,sep='\t',skip=2)
  names_rna<- gsub('\\.\\.','_',names(df_mRNA))
  names_rna<- gsub('\\_\\.','_',names_rna)
  names_rna<- gsub('\\.','_',names_rna)
  colnames(df_mRNA)<- names_rna
  tissues_mRNA_h<-c('Adrenal_Gland','Artery_Aorta','Bladder','Brain_Cerebellum','Brain_Substantia_nigra','Colon_Sigmoid','Colon_Transverse','Esophagus_Gastroesophageal_Junction','Esophagus_Mucosa','Esophagus_Muscularis','Heart_Atrial_Appendage','Heart_Left_Ventricle','Kidney_Cortex','Kidney_Medulla','Liver','Lung','Muscle_Skeletal','Ovary','Pancreas','Prostate','Skin_Not_Sun_Exposed_Suprapubic_','Skin_Sun_Exposed_Lower_leg_','Spleen','Stomach','Testis','Thyroid','Uterus')
  df_rnaseq<- subset(df_mRNA,select=c('Name','Description',tissues_mRNA_h))
  df_rnaseq$AvmRNAExpr<- rowMeans(df_rnaseq[3:29])
  df_rnaseq$Ensembl<- gsub('\\..','',df_rnaseq$Name)
  #Merge mRNA and housekeeping peptides/genes based on ENSEMBL gene ID:
  df_housekeeping<- merge(df_allgenes,cbind(df_rnaseq[30:31]),by='Ensembl',all=T)
  df_housekeeping<- df_housekeeping[order(df_housekeeping$Ensembl,decreasing = T),]
  df_housekeeping<- df_housekeeping[!duplicated(df_housekeeping$Ensembl),]
  df_housekeeping<- df_housekeeping[order(df_housekeeping$AvmRNAExpr,decreasing = T),]
  df_housekeeping$GeneRank<- seq(1,nrow(df_housekeeping),1)
  df_housekeeping$Housekeeping[is.na(df_housekeeping$Housekeeping)]<- 0
  df_housekeeping<- df_housekeeping[order(df_housekeeping$RankN,df_housekeeping$Gene,decreasing = T),]
  ProgressBar(i=3,len_i=10)
  p2.1<-ggplot2::ggplot(df_housekeeping,ggplot2::aes(x=GeneRank,y= log2(AvmRNAExpr),fill=factor(Housekeeping)))+ggplot2::geom_bar(stat='identity')+ ggplot2::theme(legend.position = 'none',panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_blank() )+ggplot2::xlim(0,44000)
  p2.2<-ggplot2::ggplot(df_housekeeping,ggplot2::aes(x=factor(Housekeeping),y= log2(AvmRNAExpr),fill=factor(Housekeeping) ))+ggplot2::geom_boxplot(outlier.shape = NA)+ ggplot2::theme_void() + ggplot2::theme(legend.position="none")
  p2 <- cowplot::ggdraw() +cowplot::draw_plot(ggplotify::as.ggplot(p2.1), x = 0, y = 0, width = 0.9, height = 1) +cowplot::draw_plot(p2.2, x = 0.85, y = 0, width = 0.1, height = 1)
  ProgressBar(i=5,len_i=10)

  #Mouse :
  dfm<- df_mouse
  dfmr<- reshape(dfm,direction="wide",timevar="Tissue",idvar=c("Peptide"))
  dfmr<- dfmr[,order(names(dfmr))]
  dfmr$Av.Int<-rowMeans(subset(dfmr, select = getIndices(dfmr,'area') ), na.rm = TRUE)
  dfmr$Rank<- apply(dfmr[,min(getIndices(dfmr,'rank')):max(getIndices(dfmr,'rank'))],1,max,na.rm=T)
  dfmr2<- cbind(dfmr['Peptide'],dfmr[getIndices(dfmr,'area')],dfmr['Av.Int'],dfmr['Rank'])
  dfmr2$RankN<-rowSums(!is.na(dfmr2[,getIndices(dfmr2,'area')]))
  df_m<- df_mouse[c(1,10)]
  df_m<- df_m[order(df_m$Peptide),]
  df_m<- df_m[!duplicated(df_m$Peptide),]
  dfmr3<- merge(dfmr2,df_m,by='Peptide',all=F)
  subset(dfmr3,dfmr3$RankN>=17)-> df.housekeepers
  ENSEMBL_house<-suppressMessages( na.omit(data.frame(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=as.character(df.housekeepers$Uniprot), columns=c('ENSEMBL'), keytype='UNIPROT'))))
  ENSEMBL_house$Housekeeping<- 1
  #Mouse rna data:
  file<-system.file("extdata", "mouse_tpm.txt", package = "MHCIatlas")
  read.csv(file,header=T,sep='\t')->df.RNAm
  df.RNAm$ensembl<- rownames(df.RNAm)
  df.RNAm$RNAseq<- 1
  tissues_mouse<- c("Colon",    "Heart"  ,'Jejunum' , "Kidney" ,  "Liver"  ,  "Brain" ,  "Pancreas"  ,"Stomach",  "Thymus")
  dfrna<- cbind(df.RNAm[40:41])
  for (tissue in tissues_mouse ) {
    dfrna[,paste0(tissue)]<- rowMeans(df.RNAm[grep(tissue,names(df.RNAm))],na.rm=T)
  }
  colnames(dfrna)<- gsub('Jejunum','Smallintestine',names(dfrna))
  dfrna$AvmRNAExpr<- rowMeans(dfrna[3:11])
  subset(merge(ENSEMBL_house,dfrna,by.x='ENSEMBL',by.y='ensembl',suffixes = c('.Peptide','.RNAseq'),all=T),RNAseq==1)-> dfrna.housem
  dfrna.housem<-dfrna.housem[order(dfrna.housem$AvmRNAExpr,decreasing = T),]
  dfrna.housem$GeneRank<- seq(1,nrow(dfrna.housem),1)
  dfrna.housem$Housekeeping[is.na(dfrna.housem$Housekeeping)]<- 0
  ProgressBar(i=7,len_i=10)
  p1.1<- ggplot2::ggplot(dfrna.housem,ggplot2::aes(x=GeneRank,y= log2(AvmRNAExpr),fill=factor(Housekeeping)))+ggplot2::geom_bar(stat='identity')+ ggplot2::theme(legend.position = 'none',panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_blank() )+ggplot2::xlim(0,32000)
  p1.2<-ggplot2::ggplot(dfrna.housem,ggplot2::aes(x=factor(Housekeeping),y= log2(AvmRNAExpr),fill=factor(Housekeeping) ))+ggplot2::geom_boxplot(outlier.shape = NA)+ ggplot2::theme_void() + ggplot2::theme(legend.position="none")
  p1 <- cowplot::ggdraw() +cowplot::draw_plot(ggplotify::as.ggplot(p1.1), x = 0, y = 0, width = 0.9, height = 1) +cowplot::draw_plot(ggplotify::as.ggplot(p1.2), x = 0.85, y = 0, width = 0.1, height = 1)
  ProgressBar(i=9,len_i=10)
  ##Plot Figure
  plot<-suppressWarnings( cowplot::ggdraw() +
    cowplot::draw_plot(ggplotify::as.ggplot(p1), x = 0, y = 0.5, width = 0.5, height = 0.5)+
    cowplot::draw_plot(ggplotify::as.ggplot(p2), x = 0.5, y = 0.5, width = 0.5, height = 0.5)+
    cowplot::draw_plot(ggplotify::as.ggplot(p3), x = 0, y = 0, width = 0.5, height = 0.5)+
    cowplot::draw_plot(ggplotify::as.ggplot(p4), x = 0.5, y = 0, width = 0.5, height = 0.5)+
    cowplot::draw_plot_label(label = c('A','B','C','D'), size = 15,x = c(0,0.5,0,0.5), y = c(1, 1,0.52,0.52)) )
  ProgressBar(i=10,len_i=10)
  cat('\n')
  suppressWarnings( print(plot) )
  suppressWarnings( return(plot) )
}
