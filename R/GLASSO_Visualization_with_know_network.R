



#' ManualSelection_known_network
#'
#' Creates an visualization for GLASSO-solution path in a simulated setup.
#' @param data Sampled data in a data.frame form.
#' @param output_file The name of the file if user decides to save the solution chosen.
#' @param real_network Real network in either matrix or sparse-matrix form
#' @examples
#' data <- huge.generator(n = 75, d = 100, graph = "scale-free", verbose = FALSE)
#' app <- ManualSelection_known_network(data = data$data,  real_network = data$theta)
#' shiny::runApp(app)
ManualSelection_known_network <- function(data, output_file = "./rho_value.RData", real_network) {

  # Parameters
  n <- nrow(data)
  p <- ncol(data)

  min <- 0.1
  max <- 1
  step <- 0.01
  path <- seq(min, max, by = step)

  # Create solution path using the huge package
  solution_path <- huge(data, method = "glasso", lambda = path)

  # Create a matrix to store scores for each solution in the path
  scores <- matrix(nrow = length(path), ncol = 14)
  for (i in seq_along(path)) {
    cm <- conf_matrix(adjacency_matrix(solution_path$path[[i]]), real_network)
    scores[i, ] <- unpack_score_list(calculate_scores(cm))
  }

  # Compute clustering coefficients for each path
  ccs <- matrix(ncol = p, nrow = length(path))
  for (i in seq_along(path)) {
    cc <- clustcoef_auto(adjacency_matrix(solution_path$path[[i]]), thresholdWS = 0, thresholdON = 0)
    ccs[i, ] <- as.matrix(round(cc["clustWS"], 3))
  }

  # Function to convert adjacency matrix to tibble format
  adjacency_to_tibble <- function(adj_matrix) {
    if (sum(adj_matrix) == 0) {
      edge_list <- tibble(from = integer(0), to = integer(0)) # Empty edge list
    } else {
      edges <- which(adj_matrix == 1, arr.ind = TRUE)
      edge_list <- tibble(from = edges[, 1] - 1, to = edges[, 2] - 1) # Convert to 0-indexing
    }
    return(edge_list)
  }

  # Shiny application
  ui <- fluidPage(
    titlePanel("Interactive Glasso Network Visualization"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("rho", "Choose regularization parameter:", min = min, max = max, value = 0.5, step = step),
        actionButton("saveButton", "Save choise"),
        checkboxInput("checkbox", "Show wrong connections", value = TRUE),
        radioButtons("nodeOption", "Choose source of node ID:",
                     choices = c("Show clustering coefficent" = "default",
                                 "Show Walktrap-glusters" = "cluster_tags"),
                     selected = "default")  # Add radio-buttons
      ),
      sidebarPanel(
        uiOutput("uu")
      )
    ),
    mainPanel(
      forceNetworkOutput("network")
    )
  )

  server <- function(input, output, session) {

    edges <- reactive({
      # Select the closest solution from the solution path.
      rho <- input$rho
      index <- which.min(abs(solution_path$lambda - rho))
      adj_matrix <- adjacency_matrix(solution_path$path[[index]])
      adjacency_to_tibble(adj_matrix)
    })

    edge_colors <- reactive({
      rho <- input$rho
      index <- which.min(abs(solution_path$lambda - rho))
      adj_matrix <- adjacency_matrix(solution_path$path[[index]])

      # Access real network correctly
      real_network <- as.matrix(real_network)

      # Create a matrix to store colors
      colors <- matrix("gray", nrow = nrow(real_network), ncol = ncol(real_network))

      if (input$checkbox) {
        # Mark correct edges as green and incorrect edges as red
        correct_edges <- (adj_matrix == real_network)
        colors[correct_edges] <- "green"
        colors[!correct_edges] <- "red"
      }

      # Convert matrix to edge list format
      edge_list <- adjacency_to_tibble(adj_matrix)
      edge_color_list <- apply(edge_list, 1, function(edge) {
        colors[edge[1] + 1, edge[2] + 1] # Adjust for 0-indexing in edge list
      })

      return(edge_color_list)
    })

    nodes <- reactive({
      rho <- input$rho
      index <- which.min(abs(solution_path$lambda - rho))

      # Define source of node ID's source
      # Either walktrap clustering or
      if (input$nodeOption == "default") {
        NodeID <- ccs[index, ]  # Use clustering coefficent
      } else if (input$nodeOption == "cluster_tags") {
        NodeID <- cluster_tags(solution_path$path[[index]])  # Use cluster_tags-function
      }

      tibble(NodeID = NodeID)
    })

    output$uu <- renderUI({
      rho <- input$rho
      index <- which.min(abs(solution_path$lambda - rho))
      paste("MCC:", round(scores[index, 3], 2), "\nACC:", round(scores[index, 1], 2))  # Display the third column score
    })

    output$network <- renderForceNetwork({
      links <- as.data.frame(edges())
      node_data <- as.data.frame(nodes())
      link_colors <- edge_colors()

      if (nrow(links) == 0) {
        links <- data.frame(matrix(0, nrow = p, ncol = 2)) # Empty edge list
        colnames(links) <- c("from", "to")
        link_colors <- rep("green", nrow(links))
      }

      forceNetwork(
        Links = links,
        Nodes = node_data,
        Group = 1,
        Source = "from",
        Target = "to",
        NodeID = "NodeID",
        opacity = 1,
        zoom = TRUE,
        linkColour = link_colors, # Assign edge colors here
        bounded = FALSE,
        opacityNoHover = 1
      )
    })

    observeEvent(input$saveButton, {
      rho <- input$rho
      index <- which.min(abs(solution_path$lambda - rho))
      selected_solution <- solution_path$path[[index]]
      save(rho, selected_solution, file = output_file)
      stopApp()  # Stop the Shiny app
    })
  }

  return(list(ui = ui, server = server))
}
