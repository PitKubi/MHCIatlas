#' @title GetHumanMHCIdata
#' @description Function to retrieve the human MHCI tissue draft dataset
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
#' df_human<- GetHumanMHCIdata()
#' }
#' @rdname GetHumanMHCIdata
#' @export
GetHumanMHCIdata<- function(NetMHC_Rank_Threshold = 2,return_all_rawData = FALSE){
  files<- system.file("extdata", c("HumanHLAatlas_reprocessed1.csv","HumanHLAatlas_reprocessed2.csv","HumanHLAatlas_reprocessed3.csv","HumanHLAatlas_reprocessed4.csv","HumanHLAatlas_reprocessed5.csv"), package = "MHCIatlas")
  dfh<- list()
  for (i in 1:length(files)) {
    dfh[[i]] <- read.csv(files[i],header=T)
    ProgressBar(i,len_i=length(files))
  }
  cat('\n')
  dfh <- do.call(rbind,dfh)
  colnames(dfh)<- gsub('Patient','Donor',names(dfh))
  if (return_all_rawData == TRUE){
    cat("Returning Human-MHCI-tissue-atlas rawData\n")
    return(dfh)
  }
  else{
    dfh<- subset(dfh,Best_Donor.Rank <= NetMHC_Rank_Threshold)
    dfh<- dfh[c(1,3,5:9,13:15,72:80)]
    dfh$Donor_Tissue <- apply(cbind(dfh['Donor'],dfh['Tissue']),1,function(x){ paste0(as.character(x[1]),'_',x[2] )} )
    cat('Returning cleaned up Human-MHCI-tissue-atlas rawData\n')
    return(dfh)
    }
  }
