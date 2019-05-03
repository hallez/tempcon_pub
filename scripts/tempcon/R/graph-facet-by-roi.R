#' Create graphs that are faceted into columns by ROIs
#'
#' This function takes a dataframe (input_data) and a title (as a character string).
#' To get graphs to print, must wrap function call in print(), e.g., `print(graph_facet_by_roi(cur_summarized_data,limits,cur_title,cur_pvals_vs_0))`
#' @param data: dataframe created by `calculate_summary.R` and filtered by current conditions of interest by `select_cur_data_and_limits.R`. Should include pvalues and limits
#' @param title: character string that will be the overall title of the plot

graph_facet_by_roi <- function(input_data, title){
  ggplot2::ggplot(data = input_data, ggplot2::aes(condition, mean_PS, fill=roi)) +
    ggplot2::geom_bar(position="dodge", stat="identity") +
    ggplot2::facet_grid(. ~ roi) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_errorbar(ggplot2::aes(ymax = upper, ymin = lower), width=0.25, color="black") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(size = 7),
                   axis.title.y=ggplot2::element_text(size=20),  legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_blank(),plot.title = ggplot2::element_text(size=20,vjust=2)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_text(data=input_data, ggplot2::aes(y = min(lower) - min(lower * 0.5),
                                                label = pval_vs_0_formatted,
                                                parse = TRUE,
                                                size = 0.5))
}

graph_facet_by_roi_xlabels <- function(input_data, title){
  ggplot2::ggplot(data = input_data, ggplot2::aes(condition, mean_PS, fill=roi)) +
    ggplot2::geom_bar(position="dodge", stat="identity") +
    ggplot2::facet_grid(. ~ roi) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_errorbar(ggplot2::aes(ymax = upper, ymin = lower), width=0.25, color="black") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, color = "black", angle = 45, vjust = 1, hjust = 1), axis.title.x = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(size = 7),
                   axis.title.y=ggplot2::element_text(size=20),  legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_blank(),plot.title = ggplot2::element_text(size=20,vjust=2)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_text(data=input_data, ggplot2::aes(y = min(lower) - min(lower * 0.5),
                                                     label = pval_vs_0_formatted,
                                                     parse = TRUE,
                                                     size = 0.5))
}

graph_facet_hemi_by_roi_xlabels <- function(input_data, title){
  ggplot2::ggplot(data = input_data, ggplot2::aes(condition, mean_PS, fill=roi)) +
    ggplot2::geom_bar(position="dodge", stat="identity") +
    ggplot2::facet_grid(roi ~ hemi) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_errorbar(ggplot2::aes(ymax = upper, ymin = lower), width=0.25, color="black") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, color = "black"), axis.title.x = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(size = 7),
                   axis.title.y=ggplot2::element_text(size=20),  legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_blank(),plot.title = ggplot2::element_text(size=20,vjust=2)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_text(data=input_data, ggplot2::aes(y = min(lower) - min(lower * 0.5),
                                                     label = pval_vs_0_formatted,
                                                     parse = TRUE,
                                                     size = 0.5))
}

graph_facet_hemi_by_roi_xlabels_no_pval <- function(input_data, title){
  ggplot2::ggplot(data = input_data, ggplot2::aes(condition, mean_PS, fill=roi)) +
    ggplot2::geom_bar(position="dodge", stat="identity") +
    ggplot2::facet_grid(roi ~ hemi) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_errorbar(ggplot2::aes(ymax = upper, ymin = lower), width=0.25, color="black") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, color = "black"), axis.title.x = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(size = 7),
                   axis.title.y=ggplot2::element_text(size=20),  legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_blank(),plot.title = ggplot2::element_text(size=20,vjust=2)) +
    ggplot2::theme(legend.position = "none")
}
