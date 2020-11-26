#' @title mkFigure2
#' @description Function to plot Figure 2 of the manuscript
#' @param df_human Human MHCI tissue draft data frame from GetHumanMHCI data
#' @return List with ggplot2 objects and data to generate the figure
#' @details DETAILS
#' @examples
#' \dontrun{
#' F2 <- mkFigure2(df_human)
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_bar}},\code{\link[ggplot2]{facet_wrap}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{coord_flip}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{scale_colour_brewer}}
#'  \code{\link[cowplot]{ggdraw}},\code{\link[cowplot]{draw_plot}},\code{\link[cowplot]{draw_plot_label}}
#'  \code{\link[ggplotify]{as.ggplot}}
#' @rdname mkFigure2
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap theme element_text coord_flip labs scale_fill_brewer xlab ylab
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
#' @importFrom ggplotify as.ggplot
mkFigure2<- function(df_human){

p2<- ggplot2::ggplot(subset(df_human,Donor=='AUT01-DN11'),ggplot2::aes(x=factor(Tissue),fill=factor(Best_Donor.Allele)))+ggplot2::geom_bar(position='fill')+ggplot2::facet_wrap(~Donor,scales='free')+ ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),axis.text.y=ggplot2::element_text(size=7,colour = 'black'))+ggplot2::coord_flip()+ ggplot2::labs(fill = "Alleles")+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::xlab('Tissues (Donor AUT01-DN11)')+ggplot2::ylab('Relative Proportion')
p3<-ggplot2::ggplot(subset(df_human,Donor=='AUT01-DN13'),ggplot2::aes(x=factor(Tissue),fill=factor(Best_Donor.Allele)))+ggplot2::geom_bar(position='fill')+ggplot2::facet_wrap(~Donor,scales='free')+ ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),axis.text.y=ggplot2::element_text(size=7,colour = 'black'))+ggplot2::coord_flip()+ ggplot2::labs(fill = "Alleles")+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::xlab('Tissues (Donor AUT01-DN13)')+ggplot2::ylab('Relative Proportion')
p4<-ggplot2::ggplot(subset(df_human,Donor=='AUT01-DN12'),ggplot2::aes(x=factor(Tissue),fill=factor(Best_Donor.Allele)))+ggplot2::geom_bar(position='fill')+ggplot2::facet_wrap(~Donor,scales='free')+ ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),axis.text.y=ggplot2::element_text(size=7,colour = 'black'))+ggplot2::coord_flip()+ ggplot2::labs(fill = "Alleles")+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::xlab('Tissues (Donor AUT01-DN12)')+ggplot2::ylab('Relative Proportion')
p1<-ggplot2::ggplot(df_human,ggplot2::aes(x=factor(Tissue),fill=factor(gsub("\\:","",gsub("[0-9]","",Best_Donor.Allele))  )))+ggplot2::geom_bar(position='fill')+ ggplot2::theme(legend.position = 'right',axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,size=8,colour = 'black'),axis.text.y=ggplot2::element_text(size=7,colour = 'black'))+ ggplot2::labs(fill = "Alleles")+ggplot2::xlab('Tissues (All Subjects)')+ggplot2::ylab('Relative Proportion')
HumanConnectMap<- mkHumanConnectivityMap(df_human)
p5<- HumanConnectMap[[1]]
Figure2<- cowplot::ggdraw() +
  cowplot::draw_plot(p1, x = 0, y = .7, width = .5, height = .3) +
  cowplot::draw_plot(p2, x = .5, y = .7, width = .5, height = .3) +
  cowplot::draw_plot(p3, x = 0, y =0.4, width = .5, height =.3) +
  cowplot::draw_plot(p4, x = .5, y = 0.4, width = .5, height = .3) +
  cowplot::draw_plot(ggplotify::as.ggplot(p5), x = 0.05, y = 0, width = 0.9, height = .4) +
  cowplot::draw_plot_label(label = c("A", "B", "C","D","E"), size = 15,
                  x = c(0, 0.5, 0,0.5,0), y = c(1, 1, 0.71,0.71,0.4))

plot.new()
print(Figure2)
return.list<- list(HumanConnectMap,Figure2)
return(return.list)
cat('\n')
}
