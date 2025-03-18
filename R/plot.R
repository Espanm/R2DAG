#library(igraph)

dag_plot <- function(table,
                     name="DAG",
                     vertex_label_size=1,  # Vertex label size
                     arrow_size=1,  # Arrow size for directed edges
                     edge_label_size = 2,  # Edge label size
                     label_round = 3  # Decimals to round for label
) {

  nodes <- colnames(table)

  # Create edge list from matrix
  edge_list <- which(table != 0, arr.ind = TRUE)

  # Build a data frame for edges (source, target, weight)
  edges <- data.frame(
    from = nodes[edge_list[, 2]],  # Column is the source (j)
    to = nodes[edge_list[, 1]],    # Row is the target (i)
    weight = table[edge_list]
  )

  # Create igraph object
  g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

  # Set vertex attributes
  V(g)$size <- 15
  V(g)$label.cex <- vertex_label_size  # Customizable label size
  V(g)$label.color <- "black"
  V(g)$frame.color <- NA

  # Set edge attributes
  E(g)$width <- (E(g)$weight / max(E(g)$weight)) * 5  # Scale edge width
  E(g)$color <- "gray50"
  E(g)$arrow.size <- arrow_size  # Customizable arrow size

  # Display edge weights as labels
  E(g)$label <- round(E(g)$weight, label_round)  # Round to 2 decimal places
  E(g)$label.cex <- edge_label_size  # Larger edge label size
  E(g)$label.color <- "blue"  # Edge label color

  # Apply curvature parameter
  E(g)$curved <- 0

  # Use a force-directed layout
  set.seed(123)
  layout_pos <- layout_with_fr(g)

  # Plot the graph
  plot(
    g,
    layout = layout_pos,
    vertex.label = V(g)$name,
    vertex.label.family = "Helvetica",
    edge.curved = FALSE,  # Disable curvature globally (individual edges still use E(g)$curved)
    edge.label = E(g)$label,  # Show edge weights
    edge.label.cex = E(g)$label.cex,  # Larger edge label size
    edge.label.color = E(g)$label.color,  # Edge label color
    main = name
  )
}
