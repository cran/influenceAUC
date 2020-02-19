#' Internal function: Build an AUC ggplot2 object
#' @keywords internal
#' @param df A data frame consists of `x`, `y`, `color`, `shape`, `label`
#' @param threshold A numeric inherits `threshold` in `IAUC` and `LAUC`
#' @param title A character represents the plot's title
#' @param show.legend Logical which controls the illustration of the legends
#' @param ylimit A vector for lower and upper limits of y-axis
#' @return A ggplot2 object
build_AUC_plot <- function(df, threshold, title = "", show.legend = FALSE, ylimit = c(-1, 1), yintercept = NULL, is_two_sided = TRUE) {
  if(is.null(yintercept)) yintercept <- threshold

  # Set the default theme
  theme_set(
    theme_minimal() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5))
  )
  # df consists of `Index`, `unitslope`, `Outcome`, `slopepoten`, `label`
  result <- ggplot(df,
                   aes_string(x = "x",
                              y = "y",
                              color = "color",
                              shape = "color",
                              label = "label")) +
    geom_point() +
    scale_shape_manual(values = c(4, 16)) +
    scale_color_manual(values = c("#e63535", "#3d46eb")) +
    scale_y_continuous("Influence Value", limits = ylimit) +
    ggtitle(title)+
    geom_text_repel(show.legend = show.legend)

  if(is_two_sided)
    result <- result +
    geom_hline(yintercept = yintercept, linetype = "dashed") +
    geom_hline(yintercept = - yintercept, linetype = "dashed")
  else
    result <- result +
    geom_hline(yintercept = yintercept, linetype = "dashed")

  return(result)
}

#' Internal function: Create lift-chart ggplot2 object
#' @keywords internal
#' @param score A vector inherits `score` in `ICLC`
#' @param binary A vector inherits `binary` in `ICLC`
#' @param prop A numeric inherits `prop` in `ICLC`
#' @param xlab A character represents plot's x axis label
#' @param ylab A character represents plot's y axis label
#' @param  title A character represents plot's title
#' @return A ggplot2 object
create_lift_chart <- function(score, binary, prop, xlab, ylab, title) {
  pred <- prediction(score, binary)
  perf <- performance(pred, "lift", "rpp")
  temp <- unique(data.frame(binary, score) %>%
                   arrange(desc(score)))

  data <- data.frame(
    x = perf@x.values[[1]][-1],  #  rate of positive predictions
    y = perf@y.values[[1]][-1],  #  lift index
    temp)

  npoten <- data %>%
    filter(x <= prop, binary == 0) %>%
    mutate(y_zero = 0)

  # customize the ggplot setting
  theme_set(
    theme_minimal()+
      theme(plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size = rel(1.15)),
            axis.title.y = element_text(size = rel(1.15)))
  )

  # Create a cumulative lift chart
  clc <- ggplot(data, aes(x = .data$x, y = .data$y)) +
    geom_line() +
    labs(title = title,
         x = xlab,
         y = ylab) +
    geom_hline(yintercept = 1, linetype = "solid", col = "blue") +
    geom_segment(data = npoten,
                 mapping = aes_string(x = "x",
                                      y = "y_zero",
                                      xend = "x",
                                      yend= "y"),
                 size = 0.3,
                 color = "red",
                 linetype = 3)

  return(clc)
}

#' Internal function: Print output
#' @keywords internal
#' @param output A valid object (e.g., `IAUC`, `LAUC`)
#' @param attrs A vector of attribute names
print_output <- function(output, attrs) {
  if(sum(attrs %in% names(output)) == 0)
    cat("Wrong input object")
  else
    for(attr in attrs)
      if(attr %in% names(output) & !is.null(output[[attr]])) {
        cat(paste(attr, "is: \n", sep = " "))
        print(output[[attr]])
      }
}

#' Internal function: Plot sequentially
#' @keywords internal
#' @param output A valid object (e.g., `IAUC`, `LAUC`)
plot_sequentially <- function(objs){
  par(ask = TRUE)
  for(obj in objs) show(obj)
  par(ask = FALSE)
}
