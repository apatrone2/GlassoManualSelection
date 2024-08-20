


#' glasso_manual
#'
#' Creates an visualization of glasso solution.
#' @param data Dataframe
#' @param real_network Adjadency matrix of real network, if available.
#' @param file Path to the file, where solution is to be saved.
#' @examples
#' data <- huge.generator(n = 75, d = 100, graph = "scale-free", verbose = FALSE)
#' app <- glasso_manual(data = data$data,  real_network = data$theta)
#' shiny::runApp(app)
#' @export
glasso_manual <- function(data, file = "./manual_selection.RData", real_network = NULL){
  if(is.null(real_network)){
    ManualSelection(data = data, output_file = file)
  } else{
    ManualSelection_known_network(data = data, output_file = file, real_network = real_network)
  }
}
