#' @title HousekeepersHuman
#' @description Functon to retrieve genes of MHCI peptides represented across most tissues/donors/alleles
#' @param df_human Input data frame usually from GetHumanMHCI data
#' @param NumDonorsperGene Minimum number of donors where a gene was found to be represented by universal (housekeeping) MHCI peptides, Default: 2
#' @param NumAllelesperGene Minimum number of alleles across all donors where a gene was found to be represented by universal (housekeeping) MHCI peptides, Default: 1
#' @param topSeq Selection of the top x peptides measured the most frequently as univeral (housekeeping) from which genes will be determined as universal, Default: 100
#' @param minNumTissues Min number of tissues that must have been sampled for a donor to be considered, Default: 14
#' @return List containing ggplot2 figures and output data frame
#' @details The term universal peptide is used synonymously to housekeeping peptide
#' @examples
#' \dontrun{
#' HousekeepersHuman(df_human)
#' }
#' @seealso
#'  \code{\link[dplyr]{count}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_bar}},\code{\link[ggplot2]{coord_flip}},\code{\link[ggplot2]{scale_continuous}}
#'  \code{\link[limma]{venn}}
#'  \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#' @rdname HousekeepersHuman
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr count
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip scale_y_log10
#' @importFrom limma vennCounts vennDiagram
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
HousekeepersHuman<- function(df_human,NumDonorsperGene=2,NumAllelesperGene=1,topSeq=100,minNumTissues = 14){
  `%>%` <- magrittr::`%>%`

  dfh<- df_human[c(1,2,6,7,8,15,20)]
  tissues<-data.frame(dfh["Donor_Tissue" ] %>% dplyr::count(Donor_Tissue))
  tissues2 <- data.frame(cbind(unlist(do.call(rbind,strsplit(tissues$Donor_Tissue,split='_'))[,2]),tissues$n))
  colnames(tissues2)<- c("Tissue",'freq')
  tissues2$freq<- as.numeric(as.character(tissues2$freq))

  df_Donor_Tissues <- data.frame(table(unlist(do.call(rbind,strsplit(tissues$Donor_Tissue,split='_'))[,1])))
  ##### Make vector of patients with more than 15 tissues count for further selection of data.
  v_Dn<- as.character(subset(df_Donor_Tissues, Freq>minNumTissues)$Var1)
  housekeepers<- list()
  for (i in 1:length(v_Dn) ){
    df_di<- subset(dfh,Donor== v_Dn[i])
    df_acc<- unique(df_di[1:2])
    df_h<- cbind(df_di[1],df_di[3],df_di[5])
    df_hw<- reshape(df_h,direction='wide',timevar='Tissue',idvar='sequence')
    df_hw$RankN<- rowSums(is.na(df_hw[2:ncol(df_hw)]))
    df_hws<- subset(df_hw,RankN<1)
    df_hwsm<- merge(df_hws,df_acc,by='sequence',all=F)
    df_hwsm<-df_hwsm[order(df_hwsm$sequence,df_hwsm$accession,decreasing=T),]
    df_hwsm<- df_hwsm[!duplicated(df_hwsm$sequence),]
    assign(paste0('dfhouse_',v_Dn[i]),df_hwsm)
    if(length(df_hwsm[,1])>0){
      df_hwsm$Donor<- v_Dn[i]
      housekeepers[[length(housekeepers)+1]]<- cbind(df_hwsm[1],df_hwsm[(ncol(df_hwsm)-2):ncol(df_hwsm)])
    }
    ProgressBar(i,len_i = length(v_Dn),j=1,len_j = 3)
  }
  df_housekeeper<-do.call(rbind,housekeepers)
  p1<-ggplot2::ggplot(df_housekeeper,ggplot2::aes(x=factor(Donor)))+ggplot2::geom_bar()+ggplot2::coord_flip()

  a <- apply(df_housekeeper[3], 1, function(x) as.character(unlist(do.call(rbind,strsplit(as.character(unlist(do.call(rbind,strsplit(x,split='_'))[,1])),split='\\|'))[,2]))  )
  df_housekeeper$Gene<- a
  housepep_patient<-df_housekeeper
  df_hshort<- unique(cbind(df_housekeeper[4:5]))
  df_hg<-data.frame(df_hshort["Gene" ] %>% dplyr::count(Gene))
  df_hg<-df_hg[order(df_hg$n,decreasing = T),]
  df.hgdonor<- subset(df_hg,n>=NumDonorsperGene)
  ###########################
  ##########
  ############################
  ############# Another way to select housekeeping peptides would be to look at the problem in an allele specific way:
  dfh<- df_human[c(1,2,6,7,8,11,15,20)]
  dfh$Allele_Tissue<- paste0(dfh$Best_Donor.Allele,'_',dfh$Tissue)
  dfh<- data.frame(unique(dfh$Allele_Tissue))
  dfh<- data.frame(do.call(rbind,strsplit(as.character(dfh$unique.dfh.Allele_Tissue.), split='_', fixed=TRUE)))
  colnames(dfh)<- c('Allele','Tissue')
  df_Allele_Tissues <-  data.frame(dfh["Allele" ] %>% dplyr::count(Allele))

  ##### Make vector of patients with more than 15 tissues count for further selection of data.
  dfh<- df_human[c(1,2,6,7,8,11,15,20)]
  v_Dn<- as.character(subset(df_Allele_Tissues, n>minNumTissues)$Allele)
  housekeepers<- list()
  for (i in 1:length(v_Dn) ){
    df_pi<- subset(dfh,Best_Donor.Allele== v_Dn[i])
    df_acc<- unique(df_pi[1:2])
    df_h<- cbind(df_pi[1],df_pi[5],df_pi[8])
    df_hw<- reshape(df_h,direction='wide',timevar='Donor_Tissue',idvar='sequence')
    df_hw$RankN<- rowSums(is.na(df_hw[2:ncol(df_hw)]))
    df_hws<- subset(df_hw,RankN==0)
    df_hwsm<- merge(df_hws,df_acc,by='sequence',all=F)
    df_hwsm<-df_hwsm[order(df_hwsm$sequence,df_hwsm$accession,decreasing=T),]
    df_hwsm <- df_hwsm[!duplicated(df_hwsm$sequence),]
    assign(paste0('dfhouse_',v_Dn[i]),df_hwsm)
    if(length(df_hwsm[,1])>0){
      df_hwsm$Allele<- v_Dn[i]
      housekeepers[[length(housekeepers)+1]]<- cbind(df_hwsm[1],df_hwsm[(ncol(df_hwsm)-2):ncol(df_hwsm)])
    }
    ProgressBar(i,len_i = length(v_Dn),j=2,len_j = 3)
  }
  df_housekeeper<-do.call(rbind,housekeepers)
  #How many allele specific housekeeping pepetides are there?
  p2<- ggplot2::ggplot(df_housekeeper,ggplot2::aes(x=factor(Allele)))+ggplot2::geom_bar()+ggplot2::coord_flip()
  a <- apply(df_housekeeper[3], 1, function(x) as.character(unlist(do.call(rbind,strsplit(as.character(unlist(do.call(rbind,strsplit(x,split='_'))[,1])),split='\\|'))[,2]))  )
  df_housekeeper$Gene<- a
  housepep_allele<-df_housekeeper
  df_hshort<- unique(cbind(df_housekeeper[4:5]))
  df_hg<-data.frame(df_hshort["Gene" ] %>% dplyr::count(Gene))
  df_hg<-df_hg[order(df_hg$n,decreasing = T),]
  df.hgallele<- subset(df_hg,n>= NumAllelesperGene)


  seq.per.gene<- data.frame(unique(cbind(df_housekeeper[1],df_housekeeper[5])))
  spg_m<- merge(seq.per.gene,df.hgallele,by="Gene",all=F)[1:2]
  count.gene<- data.frame(spg_m["Gene" ] %>% dplyr::count(Gene))
  count.gene<- count.gene[order(count.gene$n,decreasing = T),]
  p3<-ggplot2::ggplot(count.gene,ggplot2::aes(x=factor(n) ))+ggplot2::geom_bar()
  ######################
  #################
  ################# Okay and now the approach to just take the top 100 pepetides that are measured the most often


  dfh_top<- dfh[c(1,5,8)]
  count.seq<- data.frame(dfh_top["sequence" ] %>% dplyr::count(sequence))
  count.seq<- count.seq[order(count.seq$n,decreasing = T),]
  count.seq_top<- count.seq[1:topSeq,]
  p4<- ggplot2::ggplot(count.seq,ggplot2::aes(x=n ))+ggplot2::geom_bar()+ggplot2::scale_y_log10()

  df_hgs<- merge(count.seq_top,unique(dfh[1:2]),all=F)
  df_hgs<-df_hgs[order(df_hgs$n,decreasing=T),]
  a <- apply(df_hgs[3], 1, function(x) as.character(unlist(do.call(rbind,strsplit(as.character(unlist(do.call(rbind,strsplit(x,split='_'))[,1])),split='\\|'))[,2]))  )
  df_hgs$Gene<- a
  housepep_top<- df_hgs
  df_tophouse<-data.frame(df_hgs["Gene" ] %>% dplyr::count(Gene))
  df_tophouse<-df_tophouse[order(df_tophouse$n,decreasing=T),]
  colnames(df_tophouse)<- c("Gene","n")
  df_tophouse$Top100<- 1
  df.hgallele$Allele<- 1
  df.hgdonor$Donor<- 1
  df_allgenes<- merge(df_tophouse,df.hgallele,by='Gene',all=T)
  df_allgenes<- merge(df_allgenes,df.hgdonor,by='Gene',all=T)
  df_allgenes<- cbind(df_allgenes[1],df_allgenes[3],df_allgenes[5],df_allgenes[7])
  df_allgenes$RankN<- rowSums(!is.na(df_allgenes[2:ncol(df_allgenes)]))
  df_allgenes[is.na(df_allgenes)]<- 0
  df_allgenes<- df_allgenes[order(df_allgenes$RankN,decreasing = T),]
  VC<- limma::vennCounts(df_allgenes[2:4])
  limma::vennDiagram(VC)

 df_allgenes$Gene<- as.character(df_allgenes$Gene)
 df_allgenes$Gene[df_allgenes$Gene=='P04222'|df_allgenes$Gene=="P30504"]<- c('P10321')
 Uniprot_all <- apply(dfh[2], 1, function(x) as.character(unlist(do.call(rbind,strsplit(as.character(unlist(do.call(rbind,strsplit(x,split='_'))[,1])),split='\\|'))[,2]))  )
 Gene_names_all<-apply(dfh[2], 1, function(x) as.character(unlist(do.call(rbind,strsplit(as.character(unlist(do.call(rbind,strsplit(x,split='_'))[,1])),split='\\|'))[,3]))  )
 proteins_all<- unique(data.frame(Gene_names_all,Uniprot_all))
 colnames(proteins_all)<- c('Gene','Uniprot')
 colnames(df_allgenes)<- c('Uniprot',names(df_allgenes[2:5]))
 df_allgenes_a<-merge(df_allgenes,proteins_all,by='Uniprot')

 Protein_annotations<-suppressMessages( data.frame(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=as.character(df_allgenes_a$Uniprot), columns=c('ENSEMBL','ENTREZID'), keytype='UNIPROT')) )
 colnames(Protein_annotations)<- c('Uniprot','Ensembl',"ENTREZID")
 housekeeping_genes<- merge(df_allgenes_a,Protein_annotations,by='Uniprot')
 housekeeping_genes<- housekeeping_genes[order(housekeeping_genes$Uniprot),]
 housekeeping_genes<- housekeeping_genes[!duplicated(housekeeping_genes$Uniprot),]
 ProgressBar(1,len_i = 1,j=3,len_j = 3)
 cat('\n')
 plots<- list(p1,p2,p3)
 list_return<- list(plots,housekeeping_genes)
 return(list_return)

}

