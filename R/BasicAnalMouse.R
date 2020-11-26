#' @title BasicAnalMouse
#' @description Function to provide basic analysis and visualization of mouse MHCI tissue draft (MHCI peptides per tissue and principal component analysis)
#' @param df_mouse Mouse MHCI data frame, usually an output from the function 'GetMouseMHCIdata'
#' @param returnPlotsOnly Only plots will be returned if TRUE and plots and data will be returned as a list if set to FALSE, Default: FALSE
#' @return List of plots and data
#' @details Returned data frame is a matrix of peptide intensities across tissues
#' @examples
#' \dontrun{
#' BasicAnalMouse(df_mouse)
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_bar}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{scale_colour_brewer}}
#'  \code{\link[FactoMineR]{PCA}}
#'  \code{\link[factoextra]{get_pca}},\code{\link[factoextra]{fviz_pca}}
#'  \code{\link[cowplot]{ggdraw}},\code{\link[cowplot]{draw_plot}},\code{\link[cowplot]{draw_plot_label}}
#' @rdname BasicAnalMouse
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar theme labs scale_fill_brewer
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_pca_var fviz_pca_var
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
BasicAnalMouse<- function(df_mouse,returnPlotsOnly=FALSE){
  #Internal functions:
  getIndices <- function(dataframe,string){ which(grepl(string,tolower(names(dataframe))))}
  dfm<- df_mouse
  tablem <- within(dfm, Tissue <- factor(Tissue,levels=names(sort(table(Tissue),decreasing=TRUE))))
  p1<- ggplot2::ggplot(tablem,ggplot2::aes(x=Tissue,fill=factor(MHC)))+ggplot2::geom_bar(position="dodge")+ggplot2::theme(legend.position = c(0.8, 0.5),axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'))+ ggplot2::labs(fill = "MHC")+ggplot2::scale_fill_brewer(palette="Set1")
  #PCR mouse:
  dfmr<- reshape(dfm,direction="wide",timevar="Tissue",idvar=c("Peptide"))
  dfmr<- dfmr[,order(names(dfmr))]
  dfmr$Av.Int<-rowMeans(subset(dfmr, select = getIndices(dfmr,'area') ), na.rm = TRUE)
  dfmr$Rank<- apply(dfmr[,min(getIndices(dfmr,'rank')):max(getIndices(dfmr,'rank'))],1,max,na.rm=T)
  dfmr2<- cbind(dfmr['Peptide'],dfmr[getIndices(dfmr,'area')],dfmr['Av.Int'],dfmr['Rank'])
  dfmr2$RankN<-rowSums(!is.na(dfmr2[,getIndices(dfmr2,'area')]))
  df.pcrm<- subset(dfmr2,RankN>0)
  TissueInd <- getIndices(df.pcrm,'area.')
  pcrm.names<- do.call(rbind,strsplit(names(df.pcrm), split='.', fixed=TRUE))[,2]
  colnames(df.pcrm)<-pcrm.names
  df.pcrm<- log2(df.pcrm[TissueInd])
  m.pcr<- as.matrix(df.pcrm)
  rownames(m.pcr)<- subset(dfmr2,RankN>0)[,'Peptide']
  m.pcr[is.na(m.pcr)]<- 0
  m.pcr[m.pcr== -Inf]<- 0
  m.pcr->m.pca
  m.pca.scale_m <- FactoMineR::PCA(m.pca, scale.unit = F, ncp = 23, graph = F)
  varm <- factoextra::get_pca_var(m.pca.scale_m)
  set.seed(123)
  res.km <- kmeans(varm$coord, centers = 2, nstart = 19)
  grpm <- as.factor(res.km$cluster)
  p2<- factoextra::fviz_pca_var(m.pca.scale_m,repel=T, col.var = grpm,palette = "Set1",legend.title = "Cluster", geom = c("point", "text"), addEllipses=TRUE,pointsize=4,labelsize=4)+ggplot2::theme(legend.position = c(0.1,0.1))
  print(cowplot::ggdraw() +
    cowplot::draw_plot(p1+ ggplot2::theme(legend.position = c(0.7, 0.6)), x = 0, y = 0, width = .5, height = 1) +
    cowplot::draw_plot(p2, x = 0.5, y = 0, width = .5, height = 1) +
    cowplot::draw_plot_label(label = c("A", "B"), size = 15,
                    x = c(0, 0.5), y = c(1, 1)))
  if(returnPlotsOnly==T){
    return(list(p1,p2))
  }else{
  return(dfmr2)
  }
}
