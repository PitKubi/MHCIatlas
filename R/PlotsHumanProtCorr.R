#' @title PlotsHumanProtCorr
#' @description Function to visualize correlations of protein intensities and MHCI peptide counts across tissues and donors in human
#' @param Human_ProtCorrelations List which is output object from the function MakeCorrHuman
#' @param SigGene_names List of genes of interest, Default: NULL
#' @param allSigprots Logical to generate plots for all genes (TRUE) or for a specified list of genes (FALSE), Default: TRUE
#' @param RankSigThreshold Number of donors where a gene has to be found with significant correlation, Default: 2
#' @param path_filename path and filename where a pdf file containing all plots is to be written to (optional), Default: NA
#' @param return_list_of_plots Logical, can return a list with plots for each gene, Default: TRUE
#' @return list of plots for each gene as ggplot2 object
#' @examples
#' \dontrun{
#' plots_h<-PlotsHumanProtCorr(corrs_h,SigGene_names=c('HLA-B'),
#' allSigprots=F,path_filename='~/HLAB.pdf')
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_smooth}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{scale_colour_brewer}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}}
#'  \code{\link[ggpmisc]{stat_poly_eq}}
#'  \code{\link[ggrepel]{geom_label_repel}}
#' @rdname PlotsHumanProtCorr
#' @export
#' @importFrom ggplot2 ggplot aes geom_smooth geom_point xlab ylab ggtitle scale_color_brewer theme element_rect element_text element_blank labs
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom ggrepel geom_text_repel
PlotsHumanProtCorr<- function(Human_ProtCorrelations,SigGene_names=NULL,allSigprots=TRUE,RankSigThreshold=2,path_filename=NA,return_list_of_plots=TRUE){
  #Internal Functions:
  getIndices <- function(dataframe,string){ which(grepl(string,tolower(names(dataframe))))}
  GG_save_pdf = function(list, path_filename) {
    pdf(path_filename)
    for (p in list) {
      print(p)
    }
    dev.off()
    invisible(NULL)
  }

  #Code
  if(allSigprots==T){
    df_corrs<-subset(Human_ProtCorrelations[[1]],RankSig>=RankSigThreshold)
  }else{
    Gene_names<- SigGene_names
    df_corrs<- merge(Human_ProtCorrelations[[1]],data.frame(Gene_names),by='Gene_names',all=F)
  }

  df_plotprot<- cbind(df_corrs[1],df_corrs[c(getIndices(df_corrs,'maxcorr'):ncol(df_corrs))])
  dfhg<- Human_ProtCorrelations[[2]]
  df.proteome<- Human_ProtCorrelations[[3]]
  df.proteomelh<-reshape(df.proteome,idvar="Gene_names",direction="long",varying=c(list(2:30)),times=c(names(df.proteome)[2:30]),timevar="Tissue",v.names=c("Int"))

  ploth<- function(gene_in){
    #gene_in<-deparse(substitute(gene_in))
    vars_toplot<- unlist(strsplit( subset(df_plotprot,Gene_names==gene_in)[,'SigCorr_Vars'],'_'))
    peptide.counts<-dfhg[which(dfhg[,1] %in% vars_toplot),][c(1:3)]
    prot<- subset(df.proteomelh,Gene_names==gene_in)[3:4]
    na.omit(merge(peptide.counts,prot,by='Tissue'))->pro
    my.formula <- y ~ x
    colnames(pro)<-c('Tissue','var','n','Int')
    ggplot2::ggplot(data = pro,ggplot2::aes(x=n,y=log10(Int),label =Tissue,color=factor(var) )) +
      ggplot2::geom_smooth(method = "lm", se=F, formula = my.formula) +
      ggpmisc::stat_poly_eq(formula = my.formula,
                   eq.with.lhs = "italic(hat(y))~`=`~",
                   ggplot2::aes(label = paste(..rr.label.., sep = "~")),label.x='right',label.y = 'bottom',
                   parse = TRUE)+
      ggplot2::geom_point(ggplot2::aes(color=factor(var)) )+ggrepel::geom_text_repel()+ggplot2::xlab("# of HLAI Peptides")+ggplot2::ylab("Protein Intensity, log10")+ggplot2::ggtitle(paste0(gene_in,' -- ','Human'))+ggplot2::scale_color_brewer(palette = 'Set1')+
      ggplot2::theme(panel.border = ggplot2::element_rect(fill=NA),axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_rect(fill = 'white'))+ggplot2::labs(color='Variable')
  }

  l.plots<- list()
  for (i in 1:nrow(df_plotprot)) {
    l.plots[[i]]<-ploth(df_plotprot[i,1])
  }
  names(l.plots)<- as.character(df_plotprot[,1])

  if(is.na(path_filename)==F ){
    cat('\nWriting .pdf file to',path_filename,'\n')
    GG_save_pdf(l.plots,path_filename)
  }

  if(return_list_of_plots ==T){
    return(l.plots)
  }

}
