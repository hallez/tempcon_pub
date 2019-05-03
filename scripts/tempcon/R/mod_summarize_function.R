mod_summarize <- function(input_df) {
  out_df <<- NULL
  out_df <<- input_df %>%
    dplyr::select(starts_with("aov"), roi) %>%
    # NB: will get an error message about column names in anova that are not recognized
    # this is because `broom::tidy` is expecting an actual anova object and not a chisq'd log likelihood test
    dplyr::mutate_at(vars(starts_with("aov")), funs(purrr::map(., ~broom::tidy(.x)))) %>%
    tidyr::gather(key = "aov_type", value = "metrics", -roi) %>%
    tidyr::unnest()
}

mod_summarize_by_hemi_and_roi <- function(input_df) {
  out_df <<- NULL
  out_df <<- input_df %>%
    dplyr::select(starts_with("aov"), roi, hemi) %>%
    # NB: will get an error message about column names in anova that are not recognized
    # this is because `broom::tidy` is expecting an actual anova object and not a chisq'd log likelihood test
    dplyr::mutate_at(vars(starts_with("aov")), funs(purrr::map(., ~broom::tidy(.x)))) %>%
    tidyr::gather(key = "aov_type", value = "metrics", -roi, -hemi) %>%
    tidyr::unnest()
}
