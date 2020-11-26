#' @title ConservationHuman
#' @description Function to calculate conservation rates of genes representing tissue specific MHCI peptides versus genes representing universal (found across all or most tissues/donors) MHCI peptides
#' @param df_human Human MHCI dataset, usually output from the function GetHumanMHCIdata
#' @param HumanHousekeepers List object from the function HousekeepersHuman
#' @param pathBW_human Path to bigwig file containing human PhastCons conservation scores, Default: '~/Downloads/hg38.phastCons100way.bw'
#' @param samplesize Number of tissue specific genes that will be randomly selected to compute conservation, Default: 2000
#' @param quantile Intensity quantile above which tissue specific peptides will be chosen (4 means that only the 25 percent most intense peptides are considered, 1 means that all peptides are considered), Default: 4
#' @param ts_DonorSpecific Logical value to set if tissue specific (ts) peptides are filtered donor by donor (TRUE) or across all donors (FALSE), Default: FALSE
#' @param MinTissuesPerDonor Peptide must be present in samples from a donor for which the specified minimum number of tissues was measured, Default: 15
#' @param returnplots Logical value to return or not return visualization object in output list, Default: TRUE
#' @return List containing ggplot2 objects (plots) and output data to create plots
#' @details Output data using the default values are provided with the package and can be accessed throug running the function 'mkFigure5'
#' @examples 
#' \dontrun{
#' cons_h<-ConservationHuman(df_human,HousekeepersHuman(df_human),samplesize=200)
#' }
#' @seealso 
#'  \code{\link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}}
#'  \code{\link[dplyr]{count}}
#'  \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#'  \code{\link[GenomicFeatures]{transcripts}},\code{\link[GenomicFeatures]{transcriptsBy}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{scale_colour_brewer}},\code{\link[ggplot2]{geom_path}},\code{\link[ggplot2]{facet_grid}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{labs}}
#'  \code{\link[cowplot]{ggdraw}},\code{\link[cowplot]{draw_plot}},\code{\link[cowplot]{draw_plot_label}}
#'  \code{\link[ggplotify]{as.ggplot}}
#' @rdname ConservationHuman
#' @export 
#' @importFrom magrittr `%>%`
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom rtracklayer BigWigFile import
#' @importFrom dplyr count
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom GenomicFeatures promoters exonsBy
#' @importFrom zoo rollapply zoo
#' @importFrom ggplot2 ggplot aes scale_fill_brewer geom_line facet_grid theme element_rect element_text element_blank labs
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
#' @importFrom ggplotify as.ggplot
ConservationHuman<- function(df_human,HumanHousekeepers,pathBW_human = "~/Downloads/hg38.phastCons100way.bw",samplesize=2000,quantile=4, ts_DonorSpecific=FALSE,MinTissuesPerDonor=15, returnplots=TRUE){
  ###Internal functions:
  getIndices <- function(dataframe,string){ which(grepl(string,tolower(names(dataframe))))}
  `%>%` <- magrittr::`%>%`
  ####Main:
  #Preparation of data:
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  if(file.exists(pathBW_human)==F){stop('BigWig file "hg38.phastCons100way.bw" not found, correct path specified?; BigWig file can be downloaded at:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw \n') }
  bw.cons = rtracklayer::BigWigFile(pathBW_human)
  dfh<- df_human[c(1,2,6,7,8,15,20)]
  tissues<-data.frame(dfh["Donor_Tissue" ] %>% dplyr::count(Donor_Tissue))
  df_Donor_Tissues <- data.frame(table(unlist(do.call(rbind,strsplit(tissues$Donor_Tissue,split='_'))[,1])))
  ##### Make vector of patients with more than 15 tissues count for further selection of data.
  v_Dn<- as.character(subset(df_Donor_Tissues, Freq>= MinTissuesPerDonor)$Var1)
  int_cutoff<- as.numeric(quantile(dfh$median_int,na.rm=T)[quantile])

  if(ts_DonorSpecific==F){
  dfh_s<- dfh[which(dfh$Donor %in% v_Dn),]
  dfh_s<- subset(dfh_s,median_int> int_cutoff)
  df_di<- dfh_s
  df_acc<- unique(df_di[1:2])
  df_h<- cbind(df_di[1],df_di[7],df_di[5])
  df_hw<- reshape(df_h,direction='wide',timevar='Donor_Tissue',idvar='sequence')
  df_hw$RankN<- rowSums(!is.na(df_hw[2:ncol(df_hw)]))
  df_hws<- subset(df_hw,RankN==1)
  df_hwsm<- merge(df_hws,df_acc,by='sequence',all=F)
  df_hwsm<-df_hwsm[order(df_hwsm$sequence,df_hwsm$accession,decreasing=T),]
  df_hwsm<- df_hwsm[!duplicated(df_hwsm$sequence),]
  df_hwsm$Donor<- 'all'
  df.tissuespec<- cbind(df_hwsm[1],df_hwsm[(ncol(df_hwsm)-2):ncol(df_hwsm)])
  a <- apply(df.tissuespec[3], 1, function(x) as.character(unlist(do.call(rbind,strsplit(as.character(unlist(do.call(rbind,strsplit(x,split='_'))[,1])),split='\\|'))[,2]))  )
  df.tissuespec$Uniprot<- a
  df_hshort<- df.tissuespec[4:5]#unique(cbind(df.tissuespec[4:5]))
  df_hg<-data.frame(df_hshort["Uniprot" ] %>% dplyr::count(Uniprot))
  df_hg<-df_hg[order(df_hg$n,decreasing = T),]
    if(nrow(df_hg)<samplesize) {message('Samplesize is larger than number of genes, all genes were used to compute conservation rates of tissue-specific genes')
      df.tissuespec<- df_hg
    }else{
      df.tissuespec<- df_hg[sample(seq(1,nrow(df_hg)),samplesize),]
      }
    #p0<-ggplot2::ggplot(df_hg,ggplot2::aes(x=factor(n),fill=factor(ifelse(n>=10,'Used for analysis','Discarded')) ))+ggplot2::geom_bar()+ggplot2::scale_y_log10()+ggplot2::xlab("Number of tissue-specific MHCI peptides per protein")+ggplot2::labs(fill='')

  }else{
  tissue_specific<- list()
  for (i in 1:length(v_Dn) ){
    df_di<- subset(dfh,Donor== v_Dn[i]&median_int>int_cutoff )
    df_acc<- unique(df_di[1:2])
    df_h<- cbind(df_di[1],df_di[3],df_di[5])
    df_hw<- reshape(df_h,direction='wide',timevar='Tissue',idvar='sequence')
    df_hw$RankN<- rowSums(!is.na(df_hw[2:ncol(df_hw)]))
    df_hws<- subset(df_hw,RankN==1)
    df_hwsm<- merge(df_hws,df_acc,by='sequence',all=F)
    df_hwsm<-df_hwsm[order(df_hwsm$sequence,df_hwsm$accession,decreasing=T),]
    df_hwsm<- df_hwsm[!duplicated(df_hwsm$sequence),]
    if(length(df_hwsm[,1])>0){
      df_hwsm$Donor<- v_Dn[i]
      tissue_specific[[length(tissue_specific)+1]]<- cbind(df_hwsm[1],df_hwsm[(ncol(df_hwsm)-2):ncol(df_hwsm)])
    }
  }
  df.tissuespec<-do.call(rbind,tissue_specific)
  a <- apply(df.tissuespec[3], 1, function(x) as.character(unlist(do.call(rbind,strsplit(as.character(unlist(do.call(rbind,strsplit(x,split='_'))[,1])),split='\\|'))[,2]))  )
  df.tissuespec$Uniprot<- a
  df_hshort<- unique(cbind(df.tissuespec[4:5]))
  df_hg<-data.frame(df_hshort["Uniprot" ] %>% dplyr::count(Uniprot))
  df_hg<-df_hg[order(df_hg$n,decreasing = T),]
  if(nrow(df_hg)<samplesize) {warning('samplesize is larger than number of genes, all genes were used to compute conservation rates of tissue-specific genes')
    df.tissuespec<- df_hg
  }else{
    df.tissuespec<- df_hg[sample(seq(1,nrow(df_hg)),samplesize),]
  }
  #p0<-ggplot2::ggplot(df_hg,ggplot2::aes(x=factor(n),fill=factor(ifelse(n>=10,'Used for analysis','Discarded')) ))+ggplot2::geom_bar()+ggplot2::scale_y_log10()+ggplot2::xlab("Number of tissue-specific MHCI peptides per protein")+ggplot2::labs(fill='')
}
  ############ Genes represented by Tissue specific:
  ENTREZID<-suppressMessages( data.frame(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=as.character(df.tissuespec$Uniprot), columns=c('ENTREZID'), keytype='UNIPROT')))
  df_in<- na.omit(ENTREZID)

  promoters_txdb <- suppressWarnings( GenomicFeatures::promoters( txdb,columns=c('gene_id') ))
  df_promoters<- data.frame(promoters_txdb)
  l.exon_scores<- list()
  l.promoter_scores<- list()
  GRanges_tmp<- GenomicFeatures::exonsBy(txdb, "gene")
  suppressWarnings(
  for (i in 1:nrow(df_in) ) {
    tryCatch({l.exon_scores[[i]]<-data.frame(zoo::rollapply(zoo::zoo( unlist(rtracklayer::import(bw.cons, as="NumericList", selection=GRanges_tmp[[ as.character(df_in[i,2]) ]]  )) ),width=12,by=12,FUN=max,align='left'))},error = function(e) e)
    tryCatch({l.promoter_scores[[i]] <-data.frame(zoo::rollapply(zoo::zoo( unlist(rtracklayer::import(bw.cons, as="NumericList", selection=promoters_txdb[ as.numeric(   row.names( df_promoters[which(df_promoters$gene_id %in% df_in[i,2] ),])    ) ] ) ) ),width=12,by=12,FUN=max,align='left'))},error = function(e) e)
    ProgressBar(i,len_i=nrow(df_in) ,j=1,len_j = 2)
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
  HumanHousekeepers[[2]] -> df.housekeepers
  ENTREZID<-suppressMessages( data.frame(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=as.character(df.housekeepers$Uniprot), columns=c('ENTREZID'), keytype='UNIPROT')))
  df_in<- na.omit(ENTREZID)
  l.exon_scores<- list()
  l.promoter_scores<- list()
  suppressWarnings(
    for (i in 1:nrow(df_in) ) {
      tryCatch({l.exon_scores[[i]]<-data.frame(zoo::rollapply(zoo::zoo( unlist(rtracklayer::import(bw.cons, as="NumericList", selection=GRanges_tmp[[ as.character(df_in[i,2]) ]]  )) ),width=12,by=12,FUN=max,align='left'))},error = function(e) e)
      tryCatch({l.promoter_scores[[i]] <-data.frame(zoo::rollapply(zoo::zoo( unlist(rtracklayer::import(bw.cons, as="NumericList", selection=promoters_txdb[ as.numeric(   row.names( df_promoters[which(df_promoters$gene_id %in% df_in[i,2] ),])    ) ] ) ) ),width=12,by=12,FUN=max,align='left'))},error = function(e) e)
      ProgressBar(i,len_i=nrow(df_in) ,j=2,len_j = 2)
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
  df_ts.m<- rbind(df_exonstsm,df_promoterstsm)
  df_ts.m$GeneSet<- 'Tissue Specific'
  df_cons<- rbind(df_hk.m,df_ts.m)
  pval_exon<-paste('Wilcoxon rank sum test, p-value: ',signif(wilcox.test(df_exonstsm$relative,df_exonsm$relative)[[3]],digits=4))
  pval_prom<-paste('Wilcoxon rank sum test, p-value: ',signif(wilcox.test(df_promoterstsm$relative,df_promotersm$relative)[[3]],digits=4))

  if(returnplots==T){
    p1<-ggplot2::ggplot(df_cons,ggplot2::aes(x=as.numeric(phastCons),y=relative,color=factor(GeneSet)))+ggplot2::scale_fill_brewer(palette="Set1")+ggplot2::geom_line()+ggplot2::facet_grid(~Type)+ggplot2::theme(panel.border = ggplot2::element_rect(fill=NA),legend.position = c(0.8, 0.3),axis.text.x=ggplot2::element_text(colour = 'black'),axis.text.y=ggplot2::element_text(colour = 'black'),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_rect(fill = 'white'))+
      ggplot2::labs(color = "Human Genes:",x='phastCons',y='Cumulative Frequency')
    p2<- cowplot::ggdraw() +cowplot::draw_plot(ggplotify::as.ggplot(p1), x = 0, y = 0, width = 1, height = 1)+cowplot::draw_plot_label(label = c(pval_exon, pval_prom), size = 9,x = c(0.0,0.45), y = c(0.9, 0.9))
    print(p2)

    return(list(df_cons,p2))#p0 not included
  }else{
    return(df_cons)
  }

}
