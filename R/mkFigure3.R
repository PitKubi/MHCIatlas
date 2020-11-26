#' @title mkFigure3
#' @description Function to plot Figure 3 of the manuscript
#' @param df_human Human MHCI tissue draft data frame from GetHumanMHCI data
#' @param df_mouse Mouse MHCI tissue draft data frame from GetMouseMHCI data
#' @return List with ggplot2 objects and data to generate the figure
#' @details DETAILS
#' @examples
#' \dontrun{
#' F3 <- mkFigure3(df_human,df_mouse)
#' }
#' @seealso
#'  \code{\link[cowplot]{ggdraw}},\code{\link[cowplot]{draw_plot}},\code{\link[cowplot]{draw_plot_label}}
#'  \code{\link[ggplot2]{theme}}
#'  \code{\link[ggplotify]{as.ggplot}}
#' @rdname mkFigure3
#' @export
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
#' @importFrom ggplot2 theme
#' @importFrom ggplotify as.ggplot
mkFigure3<- function(df_human,df_mouse){

plots_m<- BasicAnalMouse(df_mouse,returnPlotsOnly = T)
plots_h<- BasicAnalHuman(df_human,df_mouse)
connectivitymap_mouse<- mkMouseConnectivityMap(df_mouse)

suppressMessages(
print(
  cowplot::ggdraw() +
    cowplot::draw_plot(plots_m[[1]]+ ggplot2::theme(legend.position = c(0.8, 0.5)), x = 0, y = .66, width = .5, height = .33) +
    cowplot::draw_plot(plots_h[[1]], x = .5, y = .66, width = .5, height = .33) +
    cowplot::draw_plot(plots_h[[3]], x = 0, y = 0.33, width = .5, height = 0.33) +
    cowplot::draw_plot(plots_m[[2]], x = .5, y = 0.33, width = .5, height = 0.33) +
    cowplot::draw_plot(plots_h[[2]], x = 0.0, y = 0, width = .45, height = 0.33) +
    cowplot::draw_plot(ggplotify::as.ggplot(connectivitymap_mouse), x = .45, y = 0, width = .55, height = 0.33) +
    cowplot::draw_plot_label(label = c("A", "B", "C","D","E","F"), size = 15,
                  x = c(0, 0.5, 0,0.5,0,0.45), y = c(1, 1, 0.66,0.66,0.35,0.35))   ) )

}
