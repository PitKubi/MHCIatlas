#' @title BasicAnalHuman
#' @description Function to provide basic analysis and visualization of human MHCI tissue draft (MHCI peptides per tissue and principal component analysis)
#' @param df_human Human MHCI data frame, usually an output from the function 'GetHumanMHCIdata'
#' @param df_mouse Optional parameter; Mouse MHCI data frame, usually an output from the function 'GetMouseMHCIdata', Default: NULL
#' @return List of plots of the analysis
#' @details If df_mouse is specified, mouse and human MHCI tissue drafts are directly compared
#' @examples
#' \dontrun{
#' BasicAnalHuman(df_human)
#' }
#' @seealso
#'  \code{\link[dplyr]{tally}},\code{\link[dplyr]{group_by}},\code{\link[dplyr]{summarise_all}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_boxplot}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{geom_smooth}},\code{\link[ggplot2]{geom_point}}
#'  \code{\link[forcats]{fct_reorder}}
#'  \code{\link[FactoMineR]{PCA}}
#'  \code{\link[factoextra]{get_pca}},\code{\link[factoextra]{fviz_pca}}
#'  \code{\link[ggpmisc]{stat_poly_eq}}
#'  \code{\link[ggrepel]{geom_label_repel}}
#'  \code{\link[cowplot]{ggdraw}},\code{\link[cowplot]{draw_plot}},\code{\link[cowplot]{draw_plot_label}}
#' @rdname BasicAnalHuman
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr count group_by summarise_all
#' @importFrom ggplot2 ggplot aes geom_boxplot theme element_text element_blank ylab xlab geom_smooth geom_point
#' @importFrom forcats fct_reorder
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_pca_var fviz_pca_var
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
BasicAnalHuman<-function(df_human, df_mouse=NULL){

############### ToDo: make output good for future functions (Figures etc)
  `%>%` <- magrittr::`%>%`
  #Tissue peptide counts:
  dfh<- df_human
  tissues<-data.frame(dfh["Donor_Tissue" ] %>% dplyr::count(Donor_Tissue))
  tissues2 <- data.frame(cbind(unlist(do.call(rbind,strsplit(tissues$Donor_Tissue,split='_'))[,2]),tissues$n))
  colnames(tissues2)<- c("Tissue",'freq')
  tissues2$freq<- as.numeric(as.character(tissues2$freq))
  p1<-ggplot2::ggplot(tissues2,ggplot2::aes(x=forcats::fct_reorder(Tissue,freq, .desc=TRUE),y=freq))+ggplot2::geom_boxplot()+ ggplot2::theme(legend.position = 'none',axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'))+ ggplot2::theme(legend.title = ggplot2::element_blank())+ggplot2::ylab('HLA type I peptide counts')+ggplot2::xlab('')

  #PCA:
  dfh$Donor_sequence<- paste0(dfh$Donor,'_',dfh$sequence)
  dfh_r<-reshape(cbind(dfh[21],dfh[8],dfh[6]),direction="wide",timevar="Tissue",idvar=c("Donor_sequence"))
  dfh_r$sequence<-unlist(do.call(rbind,strsplit(dfh_r$Donor_sequence,split='_'))[,2])
  dfh_r <- cbind(dfh_r[,ncol(dfh_r)],dfh_r[,2:(ncol(dfh_r)-1)])
  colnames(dfh_r)<- c('sequence',names(dfh_r)[2:ncol(dfh_r)])
  dfh_rs<- data.frame(dfh_r %>% dplyr::group_by(sequence) %>% dplyr::summarise_all( list(~ mean(., na.rm = TRUE)) ) )
  dfh_rs$RankN<-rowSums(!is.na(dfh_rs[,2:ncol(dfh_rs)]))
  dfh_rs$Av.Int<- rowMeans(dfh_rs[,2:(ncol(dfh_rs)-1)],na.rm = T )
  df.pcr<- subset(dfh_rs,RankN>0)
  pcr.names<- gsub( pattern = 'median_int.','',names(df.pcr))
  colnames(df.pcr)<-pcr.names
  df.pcr<- df.pcr[2:(ncol(df.pcr)-2)]
  m.pcr<- as.matrix(df.pcr)
  rownames(m.pcr)<- subset(dfh_rs,RankN>0)[,1]
  m.pcr[is.na(m.pcr)]<- 0
  m.pcr[m.pcr== -Inf]<- 0
  m.pcr->m.pca
  FactoMineR::PCA(m.pca, scale.unit = F, ncp = 3, graph = F)->m.pca.scale
  var <- factoextra::get_pca_var(m.pca.scale)
  set.seed(123)
  res.km <- kmeans(var$coord, centers = 3, nstart = 4)
  grp <- as.factor(res.km$cluster)
  p2 <- factoextra::fviz_pca_var(m.pca.scale,repel=T, col.var = grp,palette = 'Set1',legend.title = "Cluster", geom = c("point", "text"), addEllipses=TRUE,pointsize=4,labelsize=3)+ggplot2::theme(legend.position = c(0.1,0.1))+ggplot2::theme(legend.position = c(0.1,0.1))

  #### Human versus mouse peptides per tissue:
  if( is.null(df_mouse)==F){
  pt_h<- data.frame(tissues2 %>% dplyr::group_by(Tissue) %>% dplyr::summarise_all( list(~ max(., na.rm = TRUE)) ) )
  pt_m<-data.frame(df_mouse['Tissue'] %>% dplyr::count(Tissue))
  pt_m<- data.frame(pt_m)
  colnames(pt_m)<- c('Tissue','freq')
  pt<-merge(pt_h,pt_m,by='Tissue',suffixes = c('.human','.mouse'))
  my.formula <- y ~ x
  p3<-ggplot2::ggplot(data = pt,ggplot2::aes(x=freq.mouse,y=freq.human,label =Tissue )) +
    ggplot2::geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    ggpmisc::stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 ggplot2::aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE) +ggplot2::geom_point(color='red')+ggrepel::geom_text_repel()+ggplot2::xlab("# of HLA Peptides/Tissue in Mouse")+ggplot2::ylab("# of HLA Peptides/Tissue in Human")
  suppressMessages( print(cowplot::ggdraw() +
          cowplot::draw_plot(p1+ ggplot2::theme(legend.position = c(0.7, 0.6)), x = 0, y = 0.5, width = .5, height = 0.5) +
          cowplot::draw_plot(p2, x = 0.5, y = 0.5, width = .5, height = 0.5) +
          cowplot::draw_plot(p3, x = 0.0, y = 0, width = .5, height = 0.5) +
          cowplot::draw_plot_label(label = c("A", "B"), size = 15,
                                   x = c(0, 0.5), y = c(1, 1))) )
  return(list(p1,p2,p3))
  }else{
    suppressMessages(print(cowplot::ggdraw() +
          cowplot::draw_plot(p1+ ggplot2::theme(legend.position = c(0.7, 0.6)), x = 0, y = 0, width = .5, height = 1) +
          cowplot::draw_plot(p2, x = 0.5, y = 0, width = .5, height = 1) +
          cowplot::draw_plot_label(label = c("A", "B"), size = 15,
                                   x = c(0, 0.5), y = c(1, 1))))
  }
}
