



#' ManualSelection
#'
#' Creates visualization tool for glasso. Allows user to choose GLASSO-solution and visualize it immediately.
#' @param data dataframe of observations with column-names or without.
#' @param output_file Name of the file, where user selected solution in stored
#'
#' @examples
#' data <- huge.generator(n = 75, d = 100, graph = "scale-free", verbose = FALSE)
#' app <- ManualSelection(data = data$data)
#' shiny::runApp(app)
ManualSelection <- function(data, output_file = "./rho_value.RData") {

  # number of columns and rows in data
  n <- nrow(data)
  p <- ncol(data)

  # set up for min. and max. of the available regularization to choose form
  min <- 0.1 # min of regularization
  max <- 1 # maxmimum of regularization. Often, little as 0.6 will create empty network.
  step <- 0.01 # fineness of adjustment

  # Create solution path for glasso in HUGE
  solution_path <- huge(x = data, method = "glasso", lambda = seq(min, max, by = step))

  # Function to convert adjacency-matrix to a tibble format for shiny
  adjacency_to_tibble <- function(adj_matrix) {
    if (sum(adj_matrix) == 0) {
      edge_list <- tibble(from = integer(0), to = integer(0)) # Tyhjä edge-lista
    } else {
      edges <- which(adj_matrix == 1, arr.ind = TRUE)
      edge_list <- tibble(from = edges[, 1] - 1, to = edges[, 2] - 1) # Muunna 0-indeksointiin
    }
    return(edge_list)
  }

  # Function to check if node names given in data
  NodeID_identity <- function(data){
    if(is.null(colnames(data))){
      # Use the variable numbers as nodeID
      return(data.frame(NodeID = 0:(p-1)))
    }else{
      return(data.frame(NodeID = colnames(data)))
    }
  }

  nodes <- NodeID_identity(data)


  # Luo Shiny-sovellus
  ui <- fluidPage(
    titlePanel("Interactive Glasso Network Visualization"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("rho", "Choose regularization parameter:", min = min, max = max, value = 0.5, step = step),
        actionButton("saveButton", "Save choise")
      ),
      mainPanel(
        forceNetworkOutput("network")
      )
    )
  )

  server <- function(input, output, session) {
    edges <- reactive({
      # choose the nearest solution of choise.
      rho <- input$rho
      index <- which.min(abs(solution_path$lambda - rho))
      adj_matrix <- adjacency_matrix(solution_path$path[[index]])
      adjacency_to_tibble(adj_matrix)
    })

    output$network <- renderForceNetwork({
      links <- as.data.frame(edges())
      if (nrow(links) == 0) {
        links <- data.frame(matrix(0, nrow = p, ncol = 2)) # Create an empty edge-list.
        colnames(links) <- c("from", "to")
      }
      forceNetwork(
        Links = links,
        Nodes = nodes,
        Group = 1,
        Source = "from",
        Target = "to",
        NodeID = "NodeID",
        opacity = 0.8,
        zoom = TRUE
      )
    })

    observeEvent(input$saveButton, {
      rho <- input$rho
      index <- which.min(abs(solution_path$lambda - rho))
      selected_solution <- solution_path$path[[index]]
      save(rho, selected_solution, file = output_file)
      stopApp()  # Stop Shiny-application after model is saved
    })
  }

  return(list(ui = ui, server = server))
}





