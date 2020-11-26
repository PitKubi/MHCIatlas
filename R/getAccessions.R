#' @title getAccessions
#' @description Function to extract entries from a UniprotKB/Swissprot fasta header
#' @param x fasta header
#' @param y position of string separated by '|' that is to be eaxtracted (From the right starting at position 0)
#' @return Substring
#' @details Used withing the package MHCIatlas to aid data analysis
#' @examples 
#' \dontrun{
#' getAccessions('sp|P55012|S12A2_MOUSE',0)
#' > "S12A2"
#' }
#' @rdname getAccessions
#' @export 
getAccessions<-function(x,y){
  do.call(rbind,strsplit(do.call(rbind,strsplit(x, split='|', fixed=TRUE))[length(do.call(rbind,strsplit(x, split='|', fixed=TRUE))) - y],split='_',fixed=T))[,1]
}
