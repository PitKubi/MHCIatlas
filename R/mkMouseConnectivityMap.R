#' @title mkMouseConnectivityMap
#' @description unction to generate mouse connectivity map
#' @param df_mouse Data frame containing MHCI tissue draft data from GetMouseMHCIdata
#' @return list with ggplot2 objects and results data
#' @examples
#' \dontrun{
#' mkMouseConnectivityMap(df_mouse)
#' }
#' @seealso
#'  \code{\link[reshape2]{melt}},\code{\link[reshape2]{cast}}
#'  \code{\link[dplyr]{group_by}},\code{\link[dplyr]{summarise}},\code{\link[dplyr]{n}},\code{\link[dplyr]{mutate_all}}
#'  \code{\link[Matrix]{forceSymmetric}}
#'  \code{\link[pheatmap]{pheatmap}}
#'  \code{\link[RColorBrewer]{RColorBrewer}}
#' @rdname mkMouseConnectivityMap
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom reshape2 melt dcast
#' @importFrom dplyr group_by summarise n mutate_if
#' @importFrom Matrix forceSymmetric
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
mkMouseConnectivityMap<- function(df_mouse){
  #Internal functions:
  `%>%` <- magrittr::`%>%`
  #Main
  df2r<- BasicAnalMouse(df_mouse)
  dfss<- df2r[2:20]
  t.names<- do.call(rbind,strsplit(names(dfss), split='.', fixed=TRUE))[,2]
  for (i in 1:19){
    dfss[i]<- apply(  dfss[i] ,1,function(x) ifelse( is.na(x),NA,t.names[i]))
  }

  rownames(dfss)<- seq(1,nrow(dfss),1)
  dfc<- cbind(df2r[1],dfss[1:19])
  rownames(dfc)<- seq(1,nrow(dfc),1)
  dfc$RankN<-rowSums(!is.na(dfc[,2:20]))
  t.names<- do.call(rbind,strsplit(names(dfc), split='.', fixed=TRUE))[,2]
  colnames(dfc)<- t.names
  dfcc<- dfc[2:20]
  names.dfcc<- names(dfcc)
  df.forself<- subset(dfc,RankN==1)[2:20]
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
      ProgressBar(i,len_i=imax)
    }
    else{
      df.self<-df.forself[i]
      df.self$T1<- names.dfcc[i]
      df.self<- cbind(df.self[2],df.self[1])
      colnames(df.self)<- c('T1',"T2")
      df.tissues<-na.omit(df.self)
      tissues.list[[i]]<- df.tissues
      ProgressBar(imax,len_i=imax)
      cat('\n')
    }
  }
  df.tissues_combos = na.omit(do.call(rbind, tissues.list))


  rownames(df.tissues_combos)<- seq(1,nrow(df.tissues_combos),1)

  dfcount<- df.tissues_combos %>%
    dplyr::group_by(T1, T2) %>%
    dplyr::summarise(n = dplyr::n())
  dfcount<- data.frame(dfcount)
  dat.sort = t(apply(dfcount, 1, sort))
  dfcount1<- data.frame(dat.sort[!duplicated(dat.sort),])
  dfcount<-cbind(dfcount1[2],dfcount1[3],dfcount1[1])
  colnames(dfcount)<- c('a1','a2','freq')

  dfcount %>% dplyr::mutate_if(is.factor, as.character) -> dfcount
  dfcount<-transform(dfcount, freq = as.numeric(freq))
  freqpairs<- dfcount
  reshape2::dcast(freqpairs,a1~a2,value.var = 'freq')-> df.adjacent
  suppressWarnings( data.matrix(df.adjacent) )-> m.adjacent
  row.names(m.adjacent)<- df.adjacent[,1]
  m.adjacent[is.na(m.adjacent)]<- 0
  m.adjacent<- m.adjacent[,2:20]
  as.matrix(Matrix::forceSymmetric(m.adjacent))->m.symadjacent_m
  p<-pheatmap::pheatmap(m.symadjacent_m,color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),breaks =exp(seq(log(20), log(2100), length.out = 101)) , treeheight_row = 0, treeheight_col = 0,fontsize_row=9,fontsize_col=9)
  return(p)
}
