#' @title mkFigure4
#' @description Function to plot Figure 4 of the manuscript
#' @param df_mouse Mouse MHCI tissue draft data frame from GetMouseMHCI data
#' @return List with ggplot2 objects and data to generate the figure
#' @details Mouse mRNA expression data were retrieved from: Söllner, J., Leparc, G., Hildebrandt, T. et al. An RNA-Seq atlas of gene expression in mouse and rat normal tissues. Sci Data 4, 170185 (2017). https://doi.org/10.1038/sdata.2017.185
#' @examples
#' \dontrun{
#' F4 <- mkFigure4(df_human,df_mouse)
#' }
#' @seealso
#'  \code{\link[org.Mm.eg.db]{org.Mm.eg.db}}
#'  \code{\link[dplyr]{group_by}},\code{\link[dplyr]{summarise_all}},\code{\link[dplyr]{tally}}
#'  \code{\link[pheatmap]{pheatmap}}
#'  \code{\link[RColorBrewer]{RColorBrewer}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_violin}},\code{\link[ggplot2]{facet_wrap}},\code{\link[ggplot2]{geom_boxplot}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{scale_colour_brewer}}
#'  \code{\link[webr]{PieDonut}}
#'  \code{\link[cowplot]{ggdraw}},\code{\link[cowplot]{draw_plot}},\code{\link[cowplot]{draw_plot_label}}
#'  \code{\link[ggplotify]{as.ggplot}}
#' @rdname mkFigure4
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom AnnotationDbi select
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom dplyr group_by summarise_all count
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes geom_violin facet_wrap geom_boxplot theme element_text scale_color_brewer scale_fill_brewer
#' @importFrom webr PieDonut
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
#' @importFrom ggplotify as.ggplot
mkFigure4<- function(df_mouse){

`%>%` <- magrittr::`%>%`

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
##################### prepare peptides:
df.mhc<- df_mouse
df.mhc.r<-reshape(df.mhc,direction="wide",timevar="Tissue",idvar=c("Peptide"))
df.mhc.r[,order(names(df.mhc.r))]->df.mhc.r
df.mhc.r$AvInt<-rowMeans(subset(df.mhc.r, select =grep('Area',names(df.mhc.r) ) ), na.rm = TRUE)
df.mhc.r$Rank<- apply(df.mhc.r[,grep('Rank',names(df.mhc.r) )],1,max,na.rm=T)
df.mhc.r$fasta_header<- apply(df.mhc.r[,grep('Accession',names(df.mhc.r) )],1,function(x)( x[!is.na(x)][1]))
df.mhc.r$MHC<- apply(df.mhc.r[,grep('MHC',names(df.mhc.r) )],1,function(x)( x[!is.na(x)][1]))
df.mhc.2r<- cbind(df.mhc.r['Peptide'],df.mhc.r[grep('Area',names(df.mhc.r))],df.mhc.r[154:157])
df.mhc.2r$RankN<-rowSums(!is.na(df.mhc.2r[,2:20]))
colnames(df.mhc.2r)<- do.call(rbind,strsplit(names(df.mhc.2r), split='.', fixed=TRUE))[,2]
df.mhc3 <- cbind(df.mhc.2r[1],df.mhc.2r[21:24],df.mhc.2r[ , names(dfrna)[3:11] ])
as.character(apply(df.mhc3["fasta_header"],1, function(x) do.call(rbind,strsplit(x, split='|', fixed=TRUE))[length(do.call(rbind,strsplit(x, split='|', fixed=TRUE)))-1])) ->df.mhc3$Uniprot
#########


#keytypes(org.Mm.eg.db::org.Mm.eg.db)
Uniprot_m <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=df.mhc3$Uniprot, columns='ENSEMBL', keytype='UNIPROT')
colnames(Uniprot_m)<- c('Uniprot','ensembl')
Uniprots<- unique(na.omit(merge(dfrna["ensembl"],Uniprot_m,by='ensembl',all=T,suffixes=c('.rna','.pep'))))
df.mhc3m<- merge(df.mhc3,Uniprots,by="Uniprot",all.x=T,all.y=F)
df.mhc3<-df.mhc3m[!duplicated(df.mhc3m$Peptide),]
###########
df.mhc3$Uniprot_Peptide<-paste(df.mhc3$Uniprot, df.mhc3$Peptide, sep="_")
df.mhc3$RankN<- rowSums(!is.na(df.mhc3[,7:15]))
df.mhc4<- subset(df.mhc3,RankN==1)[7:17]
t.names<- names(df.mhc4[1:9])
for (i in 1:9){
  df.mhc4[i]<- apply(  df.mhc4[i] ,1,function(x) ifelse( is.na(x),NA,t.names[i]))
}
df.mhc4$U.Tissue<- apply(df.mhc4[,1:(ncol(df.mhc4)-2)],1,function(x)  x[!is.na(x)==T] )
df.mhc5<- df.mhc4[10:12]
subset(merge(df.mhc5,dfrna,by='ensembl',suffixes = c('.Peptide','.RNAseq')),RNAseq==1)-> dfrna.m
ma_rna<- as.matrix(dfrna.m[5:13])
ma_rna.norm<- round(t(apply(ma_rna, 1, function(x)(x)/sum(x) )),digits=3)
dfrna.m[5:13]<- ma_rna.norm
rnaseq <- cbind(dfrna.m[3],dfrna.m[5:13])
unique(rnaseq$U.Tissue)
X<-data.frame(rnaseq %>% dplyr::group_by(U.Tissue) %>% dplyr::summarise_all( list(~ mean(., na.rm = TRUE)) ) )
ma_x<- as.matrix(X[,2:ncol(X)])
rownames(ma_x)<- data.frame(X)[,1]
ma_x[order(rownames(ma_x),decreasing = F),]->ma_y
ma_y[,order(colnames(ma_y),decreasing = F)]->ma_y
p9<-pheatmap::pheatmap(ma_y,scale='row',cluster_rows = F,cluster_cols = F,color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name ="RdYlBu")))(50))

##################### Stopped here: Make violin plots of abundance of tissue spec genes and mhc generating genes and all genes.
######### HLA peptides in mRNA data
df.mhc$Uniprot_Tissue<- paste( tolower(df.mhc$UniprotAcc),df.mhc$Tissue,sep='_')
dfrna$id<- seq(1,nrow(dfrna),1)

dfrna.l<- reshape(dfrna,idvar="ensembl",direction="long",varying=c(list(3:11)),times=c(names(dfrna)[3:11]),timevar="Tissue",v.names=c("IntmRNA"))
dfpepm<- subset(cbind(df.mhc3[18],df.mhc3[16]),RankN>0)
dfpepm$T.spec<- apply(dfpepm[1],1,function(x) ifelse(x==1,1,0))
df_rnahla<- merge(dfrna.l,dfpepm,by='ensembl',suffixes = c('.mRNA','.MHCI'),all=T)
df_rnahla<- subset(df_rnahla,is.na(Tissue)==F)
df_rnahla$T.spec[is.na(df_rnahla$T.spec)]<- 2
p10<-ggplot2::ggplot(df_rnahla,ggplot2::aes(x=factor(T.spec),y=log10(IntmRNA),fill=factor(T.spec)))+ggplot2::geom_violin()+ggplot2::facet_wrap(~Tissue,scales='free',nrow=2)+ggplot2::geom_boxplot(width=0.2,outlier.shape = NA)+ggplot2::theme(legend.position = c(0.9, 0.2),axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'))+ggplot2::scale_color_brewer(palette = 'Set1')


pie<-data.frame(rnaseq["U.Tissue"] %>% dplyr::count(U.Tissue ))
colnames(pie)<- c('Tissue','freq')
pie<-pie[order(pie$freq,decreasing = T),]
pie<-within(pie, Tissue <- factor(Tissue,levels=Tissue))


p11<-suppressMessages(webr::PieDonut(pie,ggplot2::aes(Tissue,count=freq),r0=0.4,start=0,labelposition=1)+ggplot2::scale_fill_brewer(palette="Set1"))
## Plot Figure 3
Figure4<- cowplot::ggdraw() +
  cowplot::draw_plot(ggplotify::as.ggplot(p10), x = 0, y = .5, width = 1, height = .5) +
  cowplot::draw_plot(ggplotify::as.ggplot(p9), x = 0.5, y = 0, width = 0.5, height = .48) +
  cowplot::draw_plot(ggplotify::as.ggplot(p11), x = 0, y = 0.05, width = 0.5, height = .5) +
  cowplot::draw_plot_label(label = c("A", "B","C"), size = 15,
                  x = c(0, 0,0.5), y = c(1, 0.51,0.51))

suppressWarnings( plot.new() )
suppressWarnings( print(Figure4) )
suppressWarnings(return(list(p10,p9,p11,Figure4) ) ) 
}
