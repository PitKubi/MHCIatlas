#' @title mkHumanConnectivityMap
#' @description Function to generate human connectivity map and allele enrichment map
#' @param df_human Data frame containing MHCI tissue draft data from GetHumanMHCIdata
#' @param returnPlotsOnly Loical, Default: FALSE
#' @param AlleleEnrichThr Min fold change value above which an allele is considered enriched within a given tissue, Default: 1.5
#' @return list with ggplot2 objects and results data
#' @examples
#' \dontrun{
#' mkHumanConnectivityMap(df_human, returnPlotsOnly=F,AlleleEnrichThr=1.5)
#' }
#' @seealso
#'  \code{\link[dplyr]{group_by}},\code{\link[dplyr]{summarise}},\code{\link[dplyr]{n}},\code{\link[dplyr]{mutate}},\code{\link[dplyr]{tally}},\code{\link[dplyr]{mutate_all}}
#'  \code{\link[tidyr]{spread}}
#'  \code{\link[reshape2]{melt}},\code{\link[reshape2]{cast}}
#'  \code{\link[RColorBrewer]{RColorBrewer}}
#'  \code{\link[pheatmap]{pheatmap}}
#'  \code{\link[Matrix]{forceSymmetric}}
#' @rdname mkHumanConnectivityMap
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr group_by summarise n mutate count mutate_if
#' @importFrom tidyr spread
#' @importFrom reshape2 melt dcast
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom Matrix forceSymmetric
mkHumanConnectivityMap<- function(df_human, returnPlotsOnly=FALSE,AlleleEnrichThr=1.5){
  #Internal:
  `%>%` <- magrittr::`%>%`
  #Main:
  df4<- data.frame(unique(df_human$Donor_Tissue))
  df4<- data.frame(do.call(rbind,strsplit(as.character(df4$unique.df_human.Donor_Tissue.), split='_', fixed=TRUE)))
  colnames(df4)<- c('Donor','Tissue')
  Donors<- sort(as.character(unique(df4$Donor)))
  Donors<- Donors[grep('AUT01',Donors)]
  Donors<- Donors[Donors!='AUT01-DN16']
  dft<- subset(df_human,Donor==Donors[1])
  dft_s<- cbind(dft[6],dft[11])
  data1<- dft_s %>%
    dplyr::group_by(Tissue, Best_Donor.Allele) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(freq = n / sum(n)*100)
  data1<- data.frame(data1)
  data1_wide <- tidyr::spread(cbind(data1[1:2],data1[4]), Tissue, freq)
  data1_wide$Av<- rowMeans(data1_wide[2:ncol(data1_wide)])
  data1_wide$Max<- apply(data1_wide[2:ncol(data1_wide)],1,function(x) max(x))
  data1_wide$Enrich<- round(data1_wide$Max/data1_wide$Av,digits=2)
  data1_wide$Donor<-Donors[1]
  Enriched_Allele<- cbind(subset(data1_wide,Enrich>AlleleEnrichThr)[1],subset(data1_wide,Enrich>AlleleEnrichThr)[(ncol(data1_wide)-1):ncol(data1_wide)])
  df_en<- as.matrix(data1_wide[2:(ncol(data1_wide)-4)])
  rownames(df_en)<- data1_wide[,1]
  ma_en<-sweep(df_en, 1,data1_wide$Av , `/`)
  which(ma_en >=0, arr.ind = TRUE)-> df_index
  alleles<-rownames(ma_en)[df_index[1,1]]
  tissues<- colnames(ma_en)[df_index[1,2]]
  enrichment<-  ma_en[df_index[1,1],df_index[1,2]]
  alleles_heat<- data.frame(tissues,alleles,enrichment)
  for (l in 1:nrow(df_index)) {
    alleles<-rownames(ma_en)[df_index[l,1]]
    tissues<- colnames(ma_en)[df_index[l,2]]
    enrichment<-  ma_en[df_index[l,1],df_index[l,2]]
    alleles_heat_tmp<- data.frame(tissues,alleles,enrichment)
    alleles_heat<- rbind(alleles_heat,alleles_heat_tmp)
  }
  dft_r<-reshape(cbind(dft[1],dft[6],dft[11]),direction="wide",timevar="Tissue",idvar=c("sequence"))
  t.names<- do.call(rbind,strsplit(names(dft_r), split='.', fixed=TRUE))[,3]
  for (j in 2:ncol(dft_r)){
    dft_r[j]<- apply(  dft_r[j] ,1,function(x) ifelse( is.na(x),NA,t.names[j]))
  }
  dfss<- dft_r[2:ncol(dft_r)]
  colnames(dfss) <- do.call(rbind,strsplit(names(dfss), split='.', fixed=TRUE))[,3]
  rownames(dfss)<- seq(1,nrow(dfss),1)
  dfss$RankN<-rowSums(!is.na(dfss[,1:ncol(dfss)]))
  dfc<- cbind(dfss[ncol(dfss)],dfss[1:(ncol(dfss)-1)])
  rownames(dfc)<- seq(1,nrow(dfc),1)
  dfcc<- dfc[2:ncol(dfc)]
  names.dfcc<- names(dfcc)
  df.forself<- subset(dfc,RankN==1)[2:ncol(dfc)]
  tissues.list=list()
  imax<-length(names(dfcc))
  for (i in 1:imax) {
    if (i<imax){
      df.tissues<- na.omit(reshape2::melt(dfcc[i:imax],id=names.dfcc[i]))
      df.tissues<- cbind(df.tissues[1],df.tissues[3])
      colnames(df.tissues)<- c('T1','T2')
      df.self<-df.forself[i]
      df.self$T1<- names.dfcc[i]
      df.self<- cbind(df.self[2],df.self[1])
      colnames(df.self)<- c('T1',"T2")
      df.tissues<-na.omit(rbind(df.tissues,df.self))
      tissues.list[[i]]<- df.tissues
    }
    else{
      df.self<-df.forself[i]
      df.self$T1<- names.dfcc[i]
      df.self<- cbind(df.self[2],df.self[1])
      colnames(df.self)<- c('T1',"T2")
      df.tissues<-na.omit(df.self)
      tissues.list[[i]]<- df.tissues
    }
  }
  df.tissues_combos_h = na.omit(do.call(rbind, tissues.list))
  df.tissues_combos_h$Donor<- Donors[1]
  ##################################################-- Loop over Donors
  for (k in 2:length(Donors)) {
    dft<- subset(df_human,Donor==Donors[k])
    dft_s<- cbind(dft[6],dft[11])
    data1<- dft_s %>%
      dplyr::group_by(Tissue, Best_Donor.Allele) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::mutate(freq = n / sum(n)*100)
    data1<- data.frame(data1)
    data1_wide <- tidyr::spread(cbind(data1[1:2],data1[4]), Tissue, freq)
    data1_wide$Av<- rowMeans(data1_wide[2:ncol(data1_wide)])
    data1_wide$Max<- apply(data1_wide[2:ncol(data1_wide)],1,function(x) max(x))
    data1_wide$Enrich<- round(data1_wide$Max/data1_wide$Av,digits=1)
    data1_wide$Donor<-Donors[k]
    Enriched_Allele_tmp<- cbind(subset(data1_wide,Enrich>AlleleEnrichThr)[1],subset(data1_wide,Enrich>AlleleEnrichThr)[(ncol(data1_wide)-1):ncol(data1_wide)])
    Enriched_Allele<- rbind(Enriched_Allele,Enriched_Allele_tmp)
    df_en<- as.matrix(data1_wide[2:(ncol(data1_wide)-4)])
    rownames(df_en)<- data1_wide[,1]
    ma_en<-sweep(df_en, 1,data1_wide$Av , `/`)
    which(ma_en >=0, arr.ind = TRUE)-> df_index
    for (l in 1:nrow(df_index)) {
      alleles<-rownames(ma_en)[df_index[l,1]]
      tissues<- colnames(ma_en)[df_index[l,2]]
      enrichment<-  ma_en[df_index[l,1],df_index[l,2]]
      alleles_heat_tmp<- data.frame(tissues,alleles,enrichment)
      alleles_heat<- rbind(alleles_heat,alleles_heat_tmp)
    }
    dft_r<-reshape(cbind(dft[1],dft[6],dft[11]),direction="wide",timevar="Tissue",idvar=c("sequence"))
    t.names<- do.call(rbind,strsplit(names(dft_r), split='.', fixed=TRUE))[,3]
    for (j in 2:ncol(dft_r)){
      dft_r[j]<- apply(  dft_r[j] ,1,function(x) ifelse( is.na(x),NA,t.names[j]))
    }
    dfss<- dft_r[2:ncol(dft_r)]
    colnames(dfss) <- do.call(rbind,strsplit(names(dfss), split='.', fixed=TRUE))[,3]
    rownames(dfss)<- seq(1,nrow(dfss),1)
    dfss$RankN<-rowSums(!is.na(dfss[,1:ncol(dfss)]))
    dfc<- cbind(dfss[ncol(dfss)],dfss[1:(ncol(dfss)-1)])
    rownames(dfc)<- seq(1,nrow(dfc),1)
    dfcc<- dfc[2:ncol(dfc)]
    names.dfcc<- names(dfcc)
    df.forself<- subset(dfc,RankN==1)[2:ncol(dfc)]
    tissues.list=list()
    imax<-length(names(dfcc))
    for (i in 1:imax) {
      if (i<imax){
        df.tissues<- na.omit(reshape2::melt(dfcc[i:imax],id=names.dfcc[i]))
        df.tissues<- cbind(df.tissues[1],df.tissues[3])
        colnames(df.tissues)<- c('T1','T2')
        df.self<-df.forself[i]
        df.self$T1<- names.dfcc[i]
        df.self<- cbind(df.self[2],df.self[1])
        colnames(df.self)<- c('T1',"T2")
        df.tissues<-na.omit(rbind(df.tissues,df.self))
        tissues.list[[i]]<- df.tissues
        ProgressBar(i,len_i=imax, j=k,len_j = length(Donors))
      }
      else{
        df.self<-df.forself[i]
        df.self$T1<- names.dfcc[i]
        df.self<- cbind(df.self[2],df.self[1])
        colnames(df.self)<- c('T1',"T2")
        df.tissues<-na.omit(df.self)
        tissues.list[[i]]<- df.tissues
        ProgressBar(i,len_i=imax, j=k,len_j = length(Donors))
      }
    }
    df.tissues_combos_tmp = na.omit(do.call(rbind, tissues.list))
    df.tissues_combos_tmp$Donor<- Donors[k]
    df.tissues_combos_h <- rbind(df.tissues_combos_h,df.tissues_combos_tmp)
  }

  #Allele enrichment map:
  alleles_heat_all<- alleles_heat %>%
    dplyr::group_by(tissues,alleles) %>%
    dplyr::summarise(mm = median(enrichment, na.rm=TRUE))

  alleles_heat_all<- subset(data.frame(alleles_heat_all),mm>=AlleleEnrichThr)
  reshape2::dcast(alleles_heat_all,tissues~alleles,value.var = 'mm')-> ma_heat
  ma_heat[is.na(ma_heat)]<- 0
  mma_heat<- as.matrix(ma_heat[2:ncol(ma_heat)])
  rownames(mma_heat)<- ma_heat[,1]
  heatmap(mma_heat,symm = F,col =RColorBrewer::brewer.pal(8,"Reds"))
  mma_heat->test
  paletteLength <- 50
  myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
  myBreaks <- c(seq(min(test), 0.1, length.out=ceiling(paletteLength/2) + 1),
                seq(max(test)/paletteLength, max(test)*0.5, length.out=floor(paletteLength/2)))
  p1<-pheatmap::pheatmap(test[,sort(colnames(test))],cluster_cols = F, treeheight_row = 0, treeheight_col = 0, color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(50), breaks=myBreaks,fontsize_row=9,fontsize_col=9)


  #Connectivity map:

  df.tissues_combos_sort<- t(apply(df.tissues_combos_h[1:2], 1, function(x) sort(x) ))
  data.frame(df.tissues_combos_sort)-> df.tissues_combos_sort
  df.tissues_combos_sort$Donor <- df.tissues_combos_h$Donor
  colnames(df.tissues_combos_sort)<- c('T1','T2','Donor')
  df.tissues_combos_h<- df.tissues_combos_sort
  dfcount<- data.frame(df.tissues_combos_h %>% dplyr::count(Donor,T1,T2))
  dfcount_s<- cbind(dfcount[2:3],dfcount[4])
  df_count_all<- dfcount_s %>%
  dplyr::group_by(T1,T2) %>%
  dplyr::summarise(mm = mean(n, na.rm=TRUE))
  dat.sort = t(apply(df_count_all, 1, sort))
  dfcount1<- data.frame(dat.sort[!duplicated(dat.sort),])
  dfcount<-cbind(dfcount1[2],dfcount1[3],dfcount1[1])
  colnames(dfcount)<- c('a1','a2','freq')
  dfcount$freq<- as.numeric(as.character(dfcount$freq))
  data.frame(setdiff(unique(df_human$Tissue),unique(dfcount$a1)))-> df.nonshared
  colnames(df.nonshared)<- c('a1')
  df.nonshared$a2<- df.nonshared$a1
  df.nonshared$freq<- 0
  rbind(dfcount,df.nonshared)->dfcount
  dfcount %>% dplyr::mutate_if(is.factor, as.character) -> dfcount
  dfcount<-transform(dfcount, freq = as.numeric(freq))
  freqpairs<- dfcount
  reshape2::dcast(freqpairs,a1~a2,value.var = 'freq')-> df.adjacent
  suppressWarnings(data.matrix(df.adjacent) )-> m.adjacent
  row.names(m.adjacent)<- df.adjacent[,1]
  m.adjacent[is.na(m.adjacent)]<- 0
  m.adjacent<- m.adjacent[,2:ncol(m.adjacent)]
  as.matrix(Matrix::forceSymmetric(m.adjacent))->m.symadjacent_h

  p2<-pheatmap::pheatmap(m.symadjacent_h,color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),breaks =exp(seq(log(20), log(2100), length.out = 101)) , treeheight_row = 0, treeheight_col = 0,fontsize_row=9,fontsize_col=9)
  results<- list(p1,p2,test,m.symadjacent_h)
  names(results)<- c('AllelePlot','HumanConnectivityPlot','AlleleEnrichMatrix','HumanConnectivityMatrix')
  return(results)
}
