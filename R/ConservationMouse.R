#' @title ConservationMouse
#' @description Function to calculate conservation rates of genes representing tissue specific MHCI peptides versus genes representing universal (found across all or most tissues) MHCI peptides
#' @param df_mouse Mouse MHCI dataset, usually output from the function 'GetMouseMHCIdata'
#' @param pathBW_mouse path to the BigWig file mm10.60way.phastCons.bw downloadable from https://genome.ucsc.edu/, Default: '~/Downloads/mm10.60way.phastCons.bw'
#' @param samplesize Number of tissue specific genes that will be randomly selected to compute conservation, Default: 2000
#' @param returnplots Logical value to return or not return visualization object in output list, Default: TRUE
#' @return List containing ggplot2 objects (plots) and output data to create plots
#' @details Output data using the default values are provided with the package and can be accessed throug running the function 'mkFigure5'
#' @examples
#' \dontrun{
#' cons_m<-ConservationMouse(df_mouse,samplesize = 200)
#' }
#' @seealso
#'  \code{\link[TxDb.Mmusculus.UCSC.mm10.knownGene]{TxDb.Mmusculus.UCSC.mm10.knownGene}}
#'  \code{\link[org.Mm.eg.db]{org.Mm.eg.db}}
#'  \code{\link[GenomicFeatures]{transcripts}},\code{\link[GenomicFeatures]{transcriptsBy}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{scale_colour_brewer}},\code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{facet_grid}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{labs}}
#'  \code{\link[cowplot]{ggdraw}}
#'  \code{\link[ggplotify]{as.ggplot}}
#' @rdname ConservationMouse
#' @export
#' @importFrom TxDb.Mmusculus.UCSC.mm10.knownGene TxDb.Mmusculus.UCSC.mm10.knownGene
#' @importFrom rtracklayer BigWigFile import
#' @importFrom AnnotationDbi select
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom GenomicFeatures promoters genes exonsBy
#' @importFrom zoo rollapply zoo
#' @importFrom ggplot2 ggplot aes scale_fill_brewer geom_line facet_grid theme element_rect element_text element_blank labs
#' @importFrom cowplot ggdraw
#' @importFrom ggplotify as.ggplot
ConservationMouse<- function(df_mouse,pathBW_mouse = "~/Downloads/mm10.60way.phastCons.bw",samplesize = 2000, returnplots=TRUE){
###Internal functions:
getIndices <- function(dataframe,string){ which(grepl(string,tolower(names(dataframe))))}
####Main:
#Preparation of data:
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
if(file.exists(pathBW_mouse)==F){stop('BigWig file "mm10.60way.phastCons.bw" not found, correct path specified?; BigWig file can be downloaded at:
https://hgdownload.soe.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw \n') }
bw.cons = rtracklayer::BigWigFile(pathBW_mouse)
promoters_txdb <- suppressWarnings( GenomicFeatures::promoters( txdb,columns=c('gene_id') ))
df_promoters<- data.frame(promoters_txdb)
GRanges_tmp<- GenomicFeatures::exonsBy(txdb, "gene")

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

############ Genes represented by Tissue specific:
subset(dfmr3,dfmr3$RankN==1)-> df.tissuespec
df.tissuespec<- subset(df.tissuespec,Av.Int>0) 
ENTREZID<-suppressMessages( data.frame(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=as.character(df.tissuespec$Uniprot), columns=c('ENTREZID'), keytype='UNIPROT')))
df_in<- na.omit(ENTREZID)
sample<- seq(1,nrow(df_in),1)
if(nrow(df_in)<samplesize) {message('Samplesize is larger than number of genes, all genes were used to compute conservation rates of tissue-specific genes')
  sample<-sample
}else{
  sample<- sample(sample,samplesize)
}

l.exon_scores<- list()
l.promoter_scores<- list()
suppressWarnings(
  for (i in 1:length(sample) ) {
    tryCatch({l.exon_scores[[i]]<-data.frame(zoo::rollapply(zoo::zoo( unlist(rtracklayer::import(bw.cons, as="NumericList", selection=GRanges_tmp[[ as.character(df_in[sample[i],2]) ]]  )) ),width=12,by=12,FUN=max,align='left'))},error = function(e) e)
    tryCatch({l.promoter_scores[[i]] <-data.frame(zoo::rollapply(zoo::zoo( unlist(rtracklayer::import(bw.cons, as="NumericList", selection=promoters_txdb[ as.numeric(   row.names( df_promoters[which(df_promoters$gene_id %in% df_in[sample[i],2] ),])    ) ] ) ) ),width=12,by=12,FUN=max,align='left'))},error = function(e) e)
    ProgressBar(i,len_i=length(sample) ,j=1,len_j = 2)
  } )

exons_ts<- data.frame(do.call(rbind,l.exon_scores))
colnames(exons_ts)<- c('phastCons')
exons_ts<- data.frame(exons_ts[order(exons_ts$phastCons),])
colnames(exons_ts)<- c('phastCons')
df_exonstsm<- data.frame(cbind( Freq=table(exons_ts$phastCons), Cumul=cumsum(table(exons_ts$phastCons)), relative=cumsum(prop.table(table(exons_ts$phastCons)))))
df_exonstsm$phastCons<- rownames(df_exonstsm)

promoters_ts<- data.frame(do.call(rbind,l.promoter_scores))
colnames(promoters_ts)<- c('phastCons')
promoters_ts<- data.frame(promoters_ts[order(promoters_ts$phastCons),])
colnames(promoters_ts)<- c('phastCons')
df_promoterstsm<- data.frame(cbind( Freq=table(promoters_ts$phastCons), Cumul=cumsum(table(promoters_ts$phastCons)), relative=cumsum(prop.table(table(promoters_ts$phastCons)))))
df_promoterstsm$phastCons<- rownames(df_promoterstsm)

df_exonstsm$Type<- 'exons'
df_promoterstsm$Type<- 'promoters'
df_ts.m<- rbind(df_exonstsm,df_promoterstsm)
df_ts.m$GeneSet<- 'Tissue Specific'

############ Genes represented by Universal peptides:
subset(dfmr3,dfmr3$RankN>=19)-> df.housekeepers
ENTREZID<-suppressMessages( data.frame(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=as.character(df.housekeepers$Uniprot), columns=c('ENTREZID'), keytype='UNIPROT')))
df_in<- na.omit(ENTREZID)
sample<- seq(1,nrow(df_in),1)
l.exon_scores<- list()
l.promoter_scores<- list()
suppressWarnings(
  for (i in 1:length(sample) ) {
    tryCatch({l.exon_scores[[i]]<-data.frame(zoo::rollapply(zoo::zoo( unlist(rtracklayer::import(bw.cons, as="NumericList", selection=GRanges_tmp[[ as.character(df_in[sample[i],2]) ]]  )) ),width=12,by=12,FUN=max,align='left'))},error = function(e) e)
    tryCatch({l.promoter_scores[[i]] <-data.frame(zoo::rollapply(zoo::zoo( unlist(rtracklayer::import(bw.cons, as="NumericList", selection=promoters_txdb[ as.numeric(   row.names( df_promoters[which(df_promoters$gene_id %in% df_in[sample[i],2] ),])    ) ] ) ) ),width=12,by=12,FUN=max,align='left'))},error = function(e) e)
    ProgressBar(i,len_i=length(sample) ,j=2,len_j = 2)
  } )

exons_house<- data.frame(do.call(rbind,l.exon_scores))
colnames(exons_house)<- c('phastCons')
exons_house<- data.frame(exons_house[order(exons_house$phastCons),])
colnames(exons_house)<- c('phastCons')
df_exonsm<- data.frame(cbind( Freq=table(exons_house$phastCons), Cumul=cumsum(table(exons_house$phastCons)), relative=cumsum(prop.table(table(exons_house$phastCons)))))
df_exonsm$phastCons<- rownames(df_exonsm)

promoters_house<- data.frame(do.call(rbind,l.promoter_scores))
colnames(promoters_house)<- c('phastCons')
promoters_house<- data.frame(promoters_house[order(promoters_house$phastCons),])
colnames(promoters_house)<- c('phastCons')
df_promotersm<- data.frame(cbind( Freq=table(promoters_house$phastCons), Cumul=cumsum(table(promoters_house$phastCons)), relative=cumsum(prop.table(table(promoters_house$phastCons)))))
df_promotersm$phastCons<- rownames(df_promotersm)

df_exonsm$Type<- 'exons'
df_promotersm$Type<- 'promoters'
df_hk.m<- rbind(df_exonsm,df_promotersm)
df_hk.m$GeneSet<- 'Universal'

df_cons.mouse<- rbind(df_hk.m,df_ts.m)
pval_exon<-paste('Wilcoxon rank sum test, p-value: ',signif(wilcox.test(df_exonstsm$relative,df_exonsm$relative)[[3]],digits=4))
pval_prom<-paste('Wilcoxon rank sum test, p-value: ',signif(wilcox.test(df_promoterstsm$relative,df_promotersm$relative)[[3]],digits=4))

if(returnplots==T){
p1<-ggplot2::ggplot(df_cons.mouse,ggplot2::aes(x=as.numeric(phastCons),y=relative,color=factor(GeneSet)))+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::geom_line()+ggplot2::facet_grid(~Type)+ggplot2::theme(panel.border = ggplot2::element_rect(fill=NA),legend.position = c(0.8, 0.3),axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_rect(fill = 'white'))+
  ggplot2::labs(color = "Mouse Genes:",x='phastCons',y='Cumulative Frequency')
p2<- cowplot::ggdraw() +cowplot::draw_plot(ggplotify::as.ggplot(p1), x = 0, y = 0, width = 1, height = 1)+cowplot::draw_plot_label(label = c(pval_exon, pval_prom), size = 9,x = c(0.0,0.45), y = c(0.9, 0.9))
print(p2)

return(list(df_cons.mouse,p2))
}else{
  return(df_cons.mouse)
}

}
