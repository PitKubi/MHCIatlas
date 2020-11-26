#' @title PlotsMouseProtCorr
#' @description Function to visualize correlations of protein intensities and MHCI peptide counts across tissues in mouse
#' @param Mouse_ProtCorrelations List which is output object from the function MakeCorrMouse
#' @param pValue_Threshold Max pValue below which correlations are considered signifcant, Default: 0.01
#' @param SigGene_names List of genes of interest, Default: NULL, Default: NULL
#' @param path_filename path and filename where a pdf file containing all plots is to be written to (optional), Default: NA
#' @param return_list_of_plots Logical, can return a list with plots for each gene, Default: TRUE
#' @return list of plots for each gene as ggplot2 object
#' @details ggplot2 plots can directly be accessed from the list of plots
#' @examples
#' \dontrun{
#' plots_m<-PlotsMouseProtCorr(corrs_m,SigGene_names=c('Erap1','Psmb1'),
#' allSigprots=F,path_filename='~/Erap1andPsmb1.pdf')
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_smooth}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{scale_colour_brewer}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}}
#'  \code{\link[ggpmisc]{stat_poly_eq}}
#'  \code{\link[ggrepel]{geom_label_repel}}
#' @rdname PlotsMouseProtCorr
#' @export
#' @importFrom ggplot2 ggplot aes geom_smooth xlab ylab ggtitle scale_color_brewer theme element_rect element_text element_blank
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom ggrepel geom_text_repel
PlotsMouseProtCorr <- function(Mouse_ProtCorrelations,pValue_Threshold =0.01,SigGene_names=NULL,path_filename=NA,return_list_of_plots=TRUE){
  #Internal Functions:
  capstr<- function(x){paste(toupper(substring(x, 1,1)), substring(x, 2), sep="", collapse=" ")}
  GG_save_pdf = function(list, path_filename) {
    pdf(path_filename)
    for (p in list) {
      print(p)
    }
    dev.off()
    invisible(NULL)
  }
  #Main Function
  df.proteome<- Mouse_ProtCorrelations[[1]]
  df.protl<-reshape(df.proteome,idvar="UniprotAcc",direction="long",varying=c(list(5:18)),times=c(names(df.proteome)[5:18]),timevar="Tissue",v.names=c("Int"))
  df.protl$Tissue<- strsplit(capstr(df.protl$Tissue),' ')[[1]]
  mhc.counts <-Mouse_ProtCorrelations[[2]]

  plotm<- function(genem_in){
    genem_in<- as.character(genem_in)
    title<-subset(df.protl,UniprotAcc== genem_in)[1,5]
    if(is.na(title)==T){ title<- genem_in}
    prot<- subset(df.protl,UniprotAcc== genem_in)[14:15]
    na.omit(merge(mhc.counts,prot,by='Tissue'))->pro
    pro$Tissue <-apply(pro[1],1, capstr)
    my.formula <- y ~ x
    ggplot2::ggplot(data = pro,ggplot2::aes(x=n,y=Int,label =Tissue)) +
      ggplot2::geom_smooth(method = "lm", se=T, formula = my.formula,color='black') +
      ggpmisc::stat_poly_eq(formula = my.formula,
                   eq.with.lhs = "italic(hat(y))~`=`~",
                   ggplot2::aes(label = paste(..rr.label.., sep = "~")),label.x='right',label.y = 'bottom',
                   parse = TRUE) +geom_point(aes(color='black' ))+ggrepel::geom_text_repel()+ggplot2::xlab("# of MHCI Peptides")+ggplot2::ylab(ifelse(Mouse_ProtCorrelations[[5]]==T,"Protein Super-SILAC Ratio H/L, log2","Protein Intensity, log10")[[1]])+ggplot2::ggtitle(paste0(title,' -- ','Mouse'))+ggplot2::scale_color_brewer(palette = 'Set1')+
      ggplot2::theme(panel.border = ggplot2::element_rect(fill=NA),legend.position = 'none',axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_rect(fill = 'white'))

  }

  prots_sigm<- subset(df.proteome,Corr.pValue<= pValue_Threshold)[c(19,20)]
  prots_sigm$Gene_names<- apply(prots_sigm[c(1:2)],1,function(x) ifelse(is.na(x[1])==T,x[2],x[1]))
  l.plotsm<- list()
  for (i in 1:nrow(prots_sigm)) {
    l.plotsm[[i]]<-plotm(prots_sigm[i,2])
  }
  names(l.plotsm)<- as.character(prots_sigm[,1])

  if(is.na(path_filename)==F ){
    cat('\nWriting .pdf file to',path_filename,'\n')
    GG_save_pdf(l.plotsm,path_filename)
  }

  if(return_list_of_plots ==T){
    return(l.plotsm)
  }
}
