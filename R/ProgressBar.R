#' @title ProgressBar
#' @description Draws progress bar, useful to show progress in loops and nested loops
#' @param i integer i, Default: 1
#' @param len_i maximum value of i, Default: 100
#' @param j integer j, Default: 1
#' @param len_j maximum value of j, Default: 1
#' @return A progress bar
#' @examples
#' \dontrun{
#' ProgressBar(1,3)
#' }
#' @rdname ProgressBar
#' @export
ProgressBar<- function(i=1,len_i=100,j=1,len_j=1){
  len <- len_i*len_j
  k <- i+(j-1)*len_i
    width <- (round(k/len,digits = 1) *10)
    widthc <- ifelse(width==0,1,width )
    cat('\rProgress:', paste0((round(k/len,digits = 2) *100),'% (',paste0(rep('=',widthc)[1:width],collapse=''),
                              paste0( ifelse( is.na(rep(' ',10-widthc)[1:(10-widthc)])==T,'',rep(' ',10-widthc)[1:(10-widthc)] ) ,collapse=''),')'))
}
