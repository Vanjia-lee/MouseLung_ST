#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Set colour
#' @param The number of color to plot
#' @return Related colors 
#' 
CustomCol <- function(n){
  my_palette <- c(
      "#1C79B3", "#A5CFE3", "#FF8000", "#FFBD6F", "#31A229", 
      "#B4DF8C", "#E31A1C", "#F89B99", "#683E9A", "#DBB9EC", 
      "#A95526", "#A5771B", "#C75DAA", "#DEA9CC", "#006027", 
      "#7CC57D", "#BBB000", "#EEF05C", "#009CA0", "#7EC6C8", 
      "#7D0112")
  return(my_palette[n])
}

#' Identify PCA dims
#' 
#'
#'
FindElbow <- function(plot_data) {
  # Prepare all required parameters from plot_data.
  dimensions <- plot_data$dims
  total_dims <- length(dimensions)
  stdev <- plot_data$stdev
  slopes_so_far <- c() # Will keep track of all the slopes (between two consecutive points) calculated.
  last_dim <- length(dimensions)

  # Calculates slopes between every pair of consecutive points.
  for (dim in 2:last_dim) {
    slopes_so_far <- c(slopes_so_far, (stdev[dim] - stdev[dim - 1]))
  }

  # Default dimensionality to return is 1.
  dimensionality <- 1
  nth_point <- 2
  for (slope in slopes_so_far) {
    # After trial and error, these conditions seem to work well to determine the elbow.
    # Note: slopes are negative in this plot, so using max(slopes_so_far).
    # The first condition looks to see if we are approaching a lower limit, but this wouldn't be enough to
    # confidently say that it is, so we need the second condition.
    # The second condition makes sure that the point used to determine the slope is very close to the
    # smallest value, which means that it increases the confidence that this slope is indeed approaching
    # a lower limit.
    if (slope > 10 * max(slopes_so_far) && stdev[nth_point] < min(stdev) + ((max(stdev) - min(stdev)) * 0.05)) {
      dimensionality <- nth_point
      return(dimensionality - 1) # Subtract 1 because we want the point that leads into the lower limit or the
      # flat part of the graph.
    }
    nth_point <- nth_point + 1
  }

  # The function could not determine dimensionality or elbow from the plot.
  return(dimensionality)
}

#' stack barplot for cluster components
#'
PlotAbundances <- function(object, groupby, propby){
    mt <- object@meta.data[,c(groupby,propby)] %>%
        group_by_at(c(groupby,propby)) %>%
        dplyr::summarise(spot_number = n()) %>%
        mutate(Freq = (spot_number / sum(spot_number)) * 100)
    
    plot <- ggplot(mt, aes_string(groupby, "Freq", fill = propby)) + 
        geom_bar(stat = "identity",position = "fill") +
        theme_classic() +
        theme(
               text = element_text(size = 18),
               axis.text.x = element_text( angle = 270, hjust = 0, vjust = 0.5),
               axis.title.x = element_text( margin = margin(10,0,0,0))) +
        labs( y = "Frequency(%)", fill = "Slides")  
    
    return(plot)
}

