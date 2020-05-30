#' @title Show IAUC results
#'
#' @description  Print IAUC output in detail
#'
#' @param x An IAUC class object for `print method
#' @param ... Not used directly
#' @seealso \link{IAUC}
#' @export
#' @method print IAUC
#' @examples
#' library(ROCR)
#' data("ROCR.simple")
#' Ioutput <- IAUC(ROCR.simple$predictions, ROCR.simple$labels)
#' print(Ioutput)

print.IAUC <- function(x, ...) {
  attrs <- c("output", "test_output")
  print_output(x, attrs)
}

#' @title Show LAUC results
#'
#' @description  Print LAUC output in detail
#'
#' @param x An LAUC class object for `print` method
#' @param ... Not used directly
#' @seealso \link{LAUC}
#' @export
#' @method print LAUC
#' @examples
#' library(ROCR)
#' data("ROCR.simple")
#' Loutput <- LAUC(ROCR.simple$predictions, ROCR.simple$labels)
#' print(Loutput)

print.LAUC <- function(x, ...) {
  print_output(x, "output")
}

#' @title Visualize IAUC result
#'
#' @description  Visualize IAUC output sequentially
#'
#' @param x An IAUC class object for `plot` method
#' @param ... Not used directly
#' @seealso \link{IAUC}
#' @export
#' @method plot IAUC
#' @examples
#' library(ROCR)
#' data("ROCR.simple")
#' Ioutput <- IAUC(ROCR.simple$predictions, ROCR.simple$labels)
#' plot(Ioutput)

plot.IAUC <- function(x, ...) {
  # Fetch the used dataframe
  df <- x[["rdata"]] %>%
    mutate(
      x = Index,
      y = SIF,
      color = Outcome,
      shape = Outcome,
      label = sifpoten
    )
  # Build AUC plots
  sif <- build_AUC_plot(
    df = df %>%
      select(x, y, color, shape, label),
    threshold = x[["threshold"]],
    title = "Sample Influence Function"
  )

  df$y <- df$DEIF
  df$label <- df$deifpoten
  deif <- build_AUC_plot(
    df = df %>%
      select(x, y, color, shape, label),
    threshold = x[["threshold"]],
    title = "Deleted Empirical Influence Function"
  )
  plot_objs <- list(sif = sif, deif = deif)
  # Produce a test plot if needed
  if (!is.null(x[["test_data"]])) {
    title <-
      paste("Sample Influence Function \n",
            "(Testing Statistic with d = ",
            x[["testdiff"]],
            ")",
            sep = "")
    test <- x[["test_data"]] %>%
      mutate(
        x = Index,
        y = pivot,
        color = Outcome,
        shape = Outcome,
        label = testpoten
      ) %>%
      select(x, y, color, shape, label)
    test_sif <- build_AUC_plot(
      df = test,
      threshold = x[["threshold"]],
      title = title,
      ylimit = NULL,
      yintercept = qnorm(1 - x[["alpha"]]),
      is_two_sided = FALSE
    )
    plot_objs[["test_sif"]] <- test_sif
  }

  plot_sequentially(plot_objs)
}

#' @title Visualize LAUC results
#'
#' @description  Visualize LAUC output sequentially
#' @param x An LAUC class object for `plot` method
#' @param ... Not used directly
#' @seealso \link{LAUC}
#' @export
#' @method plot LAUC
#' @examples
#' library(ROCR)
#' data("ROCR.simple")
#' Loutput <- LAUC(ROCR.simple$predictions, ROCR.simple$labels)
#' plot(Loutput)

plot.LAUC <- function(x, ...) {
  # Fetch the used dataframe
  df <- x[["rdata"]] %>%
    mutate(
      x = Index,
      y = unitslope,
      color = Outcome,
      shape = Outcome,
      label = slopepoten
    )
  # Build AUC plots
  slope <- build_AUC_plot(
    df = df %>%
      select(x, y, color, shape, label),
    threshold = x[["threshold"]],
    title = "Slope"
  )
  df$y <- df$unitcurvature
  df$label <- df$curvaturepoten
  curvature <- build_AUC_plot(
    df = df %>%
      select(x, y, color, shape, label),
    threshold = x[["threshold"]],
    title = "Curvature"
  )

  plot_sequentially(list(slope = slope, curvature = curvature))
}

#' @title Visualize ICLC results
#'
#' @description  Visualize ggplot2 objects in ICLC sequentially
#'
#' @param x An ICLC class object
#' @param ... Not used directly
#' @seealso \link{ICLC}
#' @export
#' @method plot ICLC
#' @examples
#' library(ROCR)
#' data("ROCR.simple")
#' Coutput <- ICLC(ROCR.simple$predictions, ROCR.simple$labels)
#' plot(Coutput)

plot.ICLC <- function(x, ...) {
  plot_sequentially(x)
}


#' @title Determine Identified Influential Cases
#'
#' @description Provide either mutually identified influential cases through IAUC and LAUC or compare with cumulative lift charts to determine which theoretical approach is more appropriate.
#'
#' @param inf_list An IAUC class object
#' @param local_list An LAUC class object
#' @seealso \link{IAUC} \link{LAUC}
#' @export
#' @references Ke, B. S., Chiang, A. J., & Chang, Y. C. I. (2018). Influence Analysis for the Area Under the Receiver Operating Characteristic Curve. Journal of biopharmaceutical statistics, 28(4), 722-734.
#' @examples
#'
#' library(ROCR)
#' data("ROCR.simple")
#' Ioutput <- IAUC(ROCR.simple$predictions, ROCR.simple$labels)
#' Loutput <- LAUC(ROCR.simple$predictions, ROCR.simple$labels)
#' pinpoint(Ioutput, Loutput)

pinpoint <- function(inf_list, local_list) {
  hampel <- fetch_output_indeces(inf_list)
  cook <- fetch_output_indeces(local_list)

  ratio <-
    length(intersect(hampel, cook)) / length(union(hampel, cook))

  if (ratio > 0.5) {
    cat("The possible influential cases are\n",
        sort(intersect(hampel, cook)),
        ".")
  } else {
    writeLines(
      "IAUC and LAUC reach an inconsistent conclusion.\n
      Select the one (IAUC or LAUC) that is more consistent with cumulative lift charts."
    )
  }
}
