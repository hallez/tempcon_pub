#' Select current data of interest (and define limits for this subset of the data)
#'
#' This function filters a summarized dataframe (input_data) based on a ROI type (cur_filter; entire, anterior, or posterior), current condition of interest (cur_condition) and hemisphere of interest (cur_hemi).
#' The summarized dataframe should first be computed by using `calculate_summary_dplyr.R` from the `halle` package and will be called `summarized_data`.
#' It returns the data, filtered on these conditions, as the global variable `cur_summarized_data`.
#' It then filters `pvals_vs_0` based on the current condition of interest and returns these filtered p-values as the global variable `cur_pvals_vs_0`.
#' @param input_data: dataframe created by `calculate_summary.R` Should be called `summarized_data`
#' @param cur_condition: condition of interest. Should be the name of a column in `data`
#' @param cur_hemi: "left" or "right"
#' @param pval_lbl: the label of the colum in `pvals_vs_0` that corresponds to `cur_condition`
#' @param cur_filter: which ROIs are being pulled in. Should be one of "entire","head", "body", or "tail"

select_cur_data_and_limits <- function (input_data, cur_condition, cur_hemi, cur_filter) {
  # clean out prior variables
  cur_summarized_data <- NULL
  limits <- NULL

  # filter based on the current iteration
  # of for loop (entire, head, body, tail)
  if(cur_filter=="entire"){
    input_data %>%
      dplyr::filter(!grepl("head|body|tail",roi),
                    roi != "CA2",
                    condition == cur_condition,
                    hemi == cur_hemi) ->>
      cur_summarized_data_w_pvals


  } else if(cur_filter == "head" | cur_filter =="body" | cur_filter =="tail"){
    input_data %>%
      dplyr::filter(grepl(cur_filter,roi),
                    condition == cur_condition,
                    hemi == cur_hemi)->>
      cur_summarized_data_w_pvals
  }
}
