#' @title GetMouseMHCIdata
#' @description Function to retrieve the mouse MHCI tissue draft dataset
#' @param NetMHC_Rank_Threshold NetMHCpan4.0 rank threshold value (peptides with binding scores below or equal to this value are selected), Default: 2
#' @param return_all_rawData Logical, TRUE returns the raw dataset while FALSE returns a cleaned up version that is to be used with most functions of the MHCIatlas Rpackage, Default: FALSE
#' @return Data frame of the human MHCI peptide tissue draft
#' @details More about the NetMHCpan-4.0 scoring system and MHCI binding prediction can be found here:
#'     NetMHCpan-4.0: Improved Peptide–MHC Class I Interaction Predictions Integrating Eluted Ligand and Peptide
#'     Binding Affinity Data
#'     Vanessa Jurtz, Sinu Paul, Massimo Andreatta, Paolo Marcatili, Bjoern Peters and Morten Nielsen
#'     The Journal of Immunology (2017) ji1700893; DOI: 10.4049/jimmunol.1700893 
#' @examples
#' \dontrun{
#' df_mouse<- GetMouseMHCIdata()
#' }
#' @rdname GetMouseMHCIdata
#' @export
GetMouseMHCIdata<- function(NetMHC_Rank_Threshold = 2,return_all_rawData = FALSE){
# Internal Functions:
  capstr<- function(x){paste(toupper(substring(x, 1,1)), substring(x, 2), sep="", collapse=" ")}
#Code
  file<-system.file("extdata", "2020-04-08_files.combined_5percFDR_Mouse.txt", package = "MHCIatlas")
  read.csv(file,header=T)->dfm
  dfm$Tissue <- gsub('spleen', 'Spleen', dfm$Tissue)
  dfm$Type[dfm$Type=='LINEAR']<- NA
  dfm$Hybrid<- apply(dfm[13],1,function(x) ifelse(!is.na(x),1,0))
  dfm$Tissue_Peptide<-paste(dfm$Peptide, dfm$Tissue, sep="_")
  subset(dfm,grepl('qa',dfm$Cluster)==F)->dfm
  dfm<-dfm[order(dfm$Tissue_Peptide, dfm$Area, decreasing=TRUE),]
  dfm<-dfm[!duplicated(dfm$Tissue_Peptide),]
  apply(dfm[26],1, function(x) do.call(rbind,strsplit(do.call(rbind,strsplit(x, split='|', fixed=TRUE))[length(do.call(rbind,strsplit(x, split='|', fixed=TRUE)))],split='_',fixed=T))[1])->dfm$Gene_name
  apply(dfm[26],1,function(x) getAccessions(x,1))-> dfm$UniprotAcc
  if (return_all_rawData == TRUE){
    dfm<- cbind(dfm[c(1,8,15:23,27,29,32:33)])
    cat("Returning Mouse MHCI rawData")
    return(dfm)
  }
  else{
    dfm<- subset(dfm,Rank <= NetMHC_Rank_Threshold)
    dfm<- subset(dfm,Gene_name!='denovo')
    rows_to_delete = which(dfm$Tissue =='melanoma'|dfm$Tissue =='Lymphoma'|dfm$Tissue =='malignantglioma'|dfm$Tissue =='Lewislungcarcinoma')
    dfm<- dfm[-(rows_to_delete),]
    dfm$Tissue<- tolower(dfm$Tissue)
    dfm$Tissue<- apply(dfm['Tissue'],1, capstr)
    dfm<- dfm[c(1,8,16,20,26,32,33,37:39)]
    cat('Returning cleaned up Mouse MHCI rawData\n')
    return(dfm)
  }
}
