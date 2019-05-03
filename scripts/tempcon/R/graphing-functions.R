dodged_boxplot_listXroi <- function(input_dat){
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = list_fact, y = r, fill = roi_fact_ordered)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
}

dodged_boxplot_listhalvesXquest <- function(input_dat){
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = list_halves_fact, y = r, fill = question_fact)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::facet_wrap(~roi_fact_ordered) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
}

dodged_boxplot_listSDxquest <- function(input_dat){
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = list_fact, y = r, fill = question_fact)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::facet_wrap(~roi_fact_ordered) +
    ggplot2::xlab("list")
}

tiled_by_enc_trial_id <- function(input_dat){
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = row_enc_trial_id, y = col_enc_trial_id, fill = mean_r)) +
    # formatting based on: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  midpoint = 0, space = "Lab",
                                  name="Pearson\nCorrelation") +
    ggplot2::facet_grid(roi_fact_ordered~hemi_fact) +
    ggplot2::xlab("encoding trial ID") +
    ggplot2::ylab("encoding trial ID")
}

subtraction_boxplot <- function(input_dat) {
  p <<- NULL
  p<<- input_dat %>%
    # remove subtraction notion for sake of plotting
    dplyr::mutate(condition = sub("_sameQ.diffQ", "", condition)) %>%
    ggplot2::ggplot(ggplot2::aes(x = condition, y = mean_r, fill = roi_fact_ordered)) +
    ggplot2::geom_boxplot() +
    ggplot2::ylab("mean PS: same question - different question") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
    ggplot2::facet_wrap(~hemi_fact)
}

sfn_2016_plot <- function(input_dat) {
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(data = ., ggplot2::aes(question_fact, mean_PS, fill=question_fact)) +
    ggplot2::geom_bar(position="dodge", stat="identity", color = "darkgreen", width = 0.7) +
    ggplot2::facet_grid(. ~ hemi_fact) +
    # ggplot2::geom_errorbar(ggplot2::aes(ymax = upper, ymin = lower), width=0.15, color="black") +
    ggplot2::ggtitle("Mean Pattern Similarity\n by Encoding Question") +
    ggplot2::ylab("mean pattern similarity") +
    ggplot2::scale_fill_manual(values = c("darkseagreen1", "springgreen3")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   strip.background = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 20),
                   strip.text.x = ggplot2::element_text(size = 15),
                   strip.text.y = ggplot2::element_text(size = 15),
                   # margin based on: http://stackoverflow.com/questions/14487188/increase-distance-between-text-and-title-on-the-y-axis
                   axis.title.y = ggplot2::element_text(size = 20, margin = ggplot2::margin(0,20,0,0)), axis.text.y = ggplot2::element_text(size = 20),
                   plot.title = ggplot2::element_text(size = 20, margin = ggplot2::margin(0, 0, 10, 0))) +
    ggplot2::ylim(-0.0015, 0.004)
}

thresholded_corrplot <- function(input_corr){

  cur_corr <- NULL
  cur_corr_signif <- NULL

  cur_corr <- cor(input_corr, use="pairwise.complete.obs")
  cur_corr_signif <- corrplot::cor.mtest(input_corr, conf.level = .95)

  corrplot::corrplot(cur_corr,
                     p.mat = cur_corr_signif$p, sig.level = .05,
                     type = "lower", method = "color", tl.col = "black")
}

tiled_by_enc_list <- function(input_dat){
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = row_list, y = col_list, fill = group_mean_r)) +
    # formatting based on: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  midpoint = 0, space = "Lab",
                                  name="Mean\nPearson\nCorrelation") +
    ggplot2::facet_grid(question~roi_formatted) +
    ggplot2::xlab("encoding list ID") +
    ggplot2::ylab("encoding list ID")
}

tiled_by_enc_list_single_roi_facet_hemi <- function(input_dat){
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = row_list, y = col_list, fill = group_mean_r)) +
    # formatting based on: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  midpoint = 0, space = "Lab",
                                  name="Mean\nPearson\nCorrelation") +
    ggplot2::facet_grid(question~hemi) +
    ggplot2::xlab("encoding list ID") +
    ggplot2::ylab("encoding list ID")
}

tiled_by_enc_list_no_quest <- function(input_dat){
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = row_list, y = col_list, fill = group_mean_r)) +
    # formatting based on: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "#fee5d9", high = "#a50f15", mid = "#fc9272",
                                  midpoint = 0, space = "Lab",
                                  name="mean pattern similarity (r)") +
    ggplot2::facet_grid(.~roi_formatted) +
    ggplot2::xlab("list") +
    ggplot2::ylab("list")
}

tiled_by_enc_list_no_quest_by_hemi <- function(input_dat){
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = row_list, y = col_list, fill = group_mean_r)) +
    # formatting based on: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "#fee5d9", high = "#a50f15", mid = "#fc9272",
                                  midpoint = 0, space = "Lab",
                                  name="mean pattern similarity (r)") +
    ggplot2::facet_grid(roi_formatted ~ hemi) +
    ggplot2::xlab("list") +
    ggplot2::ylab("list")
}

sfn_2017_theme <- ggplot2::theme(axis.title.x = ggplot2::element_text(size = 30), axis.text.x = ggplot2::element_text(size = 30),
                                 strip.text.x = ggplot2::element_text(size = 28),
                                 # margin based on: http://stackoverflow.com/questions/14487188/increase-distance-between-text-and-title-on-the-y-axis
                                 axis.title.y = ggplot2::element_text(size = 30, margin = ggplot2::margin(0,20,0,0)), axis.text.y = ggplot2::element_text(size = 25),
                                 strip.text.y = ggplot2::element_text(size = 28),
                                 plot.title = ggplot2::element_text(size = 30))

mixed_models_corr_theme <- ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16, angle = 45, hjust = 1),
                                          axis.title.x = ggplot2::element_blank(),
                                          strip.text.x = ggplot2::element_text(size = 18),
                                          axis.text.y = ggplot2::element_text(size = 16),
                                          axis.title.y = ggplot2::element_blank(),
                                          strip.text.y = ggplot2::element_text(size = 18),
                                          legend.text = ggplot2::element_text(size = 16),
                                          legend.title = ggplot2::element_text(size = 16),
                                          plot.title = ggplot2::element_text(size = 30))

mixed_models_jitterbox_theme <- ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12),
                                          axis.title.x = ggplot2::element_text(size = 18),
                                          strip.text.x = ggplot2::element_text(size = 8),
                                          axis.text.y = ggplot2::element_text(size = 16),
                                          axis.title.y = ggplot2::element_text(size = 18),
                                          strip.text.y = ggplot2::element_text(size = 12),
                                          legend.text = ggplot2::element_text(size = 16),
                                          legend.title = ggplot2::element_text(size = 16),
                                          plot.title = ggplot2::element_text(size = 30))

barplot_by_quest <- function(input_dat){
  p <<- NULL
  p <<- ggplot2::ggplot(data = input_dat, ggplot2::aes(question_fact, mean_PS, fill=roi_formatted, alpha = question_fact)) +
    ggplot2::geom_bar(position="dodge", stat="identity", width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymax = upper, ymin = lower), width = 0.1, color = "dimgray", alpha = 1) +
    # alpha scaling based on: https://stackoverflow.com/questions/24800626/r-ggplot-transparency-alpha-values-conditional-on-other-variable
    ggplot2::scale_alpha_discrete(range = c(0.5, 0.9)) +
    ggplot2::facet_grid(~roi_formatted) +
    ggplot2::xlab("encoding question") +
    ggplot2::ylab("mean PS") +
    sfn_2017_theme +
    ggplot2::theme(legend.position="none")
}

barplot_by_quest_trialnums <- function(input_dat){
  p <<- NULL
  p <<- ggplot2::ggplot(data = input_dat, ggplot2::aes(question_fact, N, fill=roi_formatted, alpha = question_fact)) +
    ggplot2::geom_bar(position="dodge", stat="identity", width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymax = upper, ymin = lower), width = 0.1, color = "dimgray", alpha = 1) +
    # alpha scaling based on: https://stackoverflow.com/questions/24800626/r-ggplot-transparency-alpha-values-conditional-on-other-variable
    ggplot2::scale_alpha_discrete(range = c(0.5, 0.9)) +
    ggplot2::facet_grid(roi_formatted~subj) +
    ggplot2::xlab("encoding question") +
    ggplot2::ylab("mean PS") +
    ggplot2::theme(legend.position="none")
}

plot_permutation <- function(input_dat, chisq_perm, p95_thresh, chisq_observed, wrap_var){
  # calling column names based on: https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot (NOT WORKING)
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(data = ., ggplot2::aes_string(x = chisq_perm)) +
    ggplot2::geom_histogram() +
    # plot the statistical cutoff
    ggplot2::geom_vline(ggplot2::aes_string(xintercept = p95_thresh), linetype="dashed",
                        color = "red") +
    # plot the observed p value
    ggplot2::geom_vline(ggplot2::aes_string(xintercept = chisq_observed),
                        color = "blue") +
    ggplot2::facet_wrap(facets = ggplot2::aes_string(wrap_var)) +
    ggplot2::ggtitle("Observed chisq (blue), p < 0.05 cutoff (red)")
}

plot_trial_nums_interaction <- function(input_df, var1_name, var2_name, subj_varname, hemi_varname, roi_varname){
  # var1_name = label for x axis
  # var2_name = row variable for faceting

  trialnums_cutoff <- 100
  trialnums_cutoff_global <<- 100
  dodge_value <- 0.8

  p_summary <<- NULL
  p_summary <<- input_df %>%
    dplyr::group_by(!!sym(subj_varname), !!sym(hemi_varname), !!sym(roi_varname), !!sym(var1_name), !!sym(var2_name)) %>%
    # this will create a column called `n`
    dplyr::tally() %>%
    summary()

  # report subject ids who have low trial numbers
  p_summary_by_subj <- input_df %>%
    dplyr::group_by(!!sym(subj_varname), !!sym(hemi_varname), !!sym(roi_varname), !!sym(var1_name)) %>%
    # this will create a column called `n`
    dplyr::tally() %>%
    dplyr::filter(n <= trialnums_cutoff)

  subj_low_trialnums <<- NULL
  subj_low_trialnums <<- unique(p_summary_by_subj$subj)
  message(sprintf("There are %d subjects with less than %d trials.", length(subj_low_trialnums), trialnums_cutoff_global))
  if(length(subj_low_trialnums > 0)){
    message(sprintf(paste0(subj_low_trialnums, collapse = ', ')))
  }

  p <<- NULL
  p <<- input_df %>%
    dplyr::group_by(!!sym(subj_varname), !!sym(hemi_varname), !!sym(roi_varname), !!sym(var1_name), !!sym(var2_name)) %>%
    # this will create a column called `n`
    dplyr::tally() %>%
    dplyr::mutate(outlier_trialnums = ifelse(n < trialnums_cutoff, "TRUE", "FALSE")) %>%
    ggplot2::ggplot(ggplot2::aes(x = factor(!!sym(var1_name)), y = n, fill = factor(!!sym(hemi_varname)))) +
    # dodge values less than 1 decrease the space between the grouped bars
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(dodge_value)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(), ggplot2::aes(color = outlier_trialnums)) +
    # just use grayscale for hemisphere differences (reserve colors for subjects)
    # re-scale so that not using black/darkest gray end of the scale
    ggplot2::scale_fill_grey(start=0.8, end=0.4) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    # feeding in columns as characters based on: https://stackoverflow.com/questions/21588096/pass-string-to-facet-grid-ggplot2
    ggplot2::facet_grid(reformulate(var2_name, roi_varname)) +
    ggplot2::ggtitle(sprintf("red < %d trials", trialnums_cutoff)) +
    mixed_models_jitterbox_theme
}

plot_trial_nums_single_var <- function(input_df, var1_name, subj_varname, hemi_varname, roi_varname){
  # var1_name = label for x axis

  trialnums_cutoff <- 100
  trialnums_cutoff_global <<- trialnums_cutoff
  dodge_value <- 0.8

  p_summary <<- NULL
  p_summary <<- input_df %>%
    dplyr::group_by(!!sym(subj_varname), !!sym(hemi_varname), !!sym(roi_varname), !!sym(var1_name)) %>%
    # this will create a column called `n`
    dplyr::tally() %>%
    summary()

  # report subject ids who have low trial numbers
  p_summary_by_subj <- input_df %>%
    dplyr::group_by(!!sym(subj_varname), !!sym(hemi_varname), !!sym(roi_varname), !!sym(var1_name)) %>%
    # this will create a column called `n`
    dplyr::tally() %>%
    dplyr::filter(n <= trialnums_cutoff)

  subj_low_trialnums <<- NULL
  subj_low_trialnums <<- unique(p_summary_by_subj$subj)
  message(sprintf("There are %d subjects with less than %d trials.", length(subj_low_trialnums), trialnums_cutoff_global))
  if(length(subj_low_trialnums > 0)){
    message(sprintf(paste0(subj_low_trialnums, collapse = ', ')))
  }

  p <<- NULL
  p <<- input_df %>%
    dplyr::group_by(!!sym(subj_varname), !!sym(hemi_varname), !!sym(roi_varname), !!sym(var1_name)) %>%
    # this will create a column called `n`
    dplyr::tally() %>%
    dplyr::mutate(outlier_trialnums = ifelse(n < trialnums_cutoff, "TRUE", "FALSE")) %>%
    ggplot2::ggplot(ggplot2::aes(x = factor(!!sym(var1_name)), y = n, fill = factor(!!sym(hemi_varname)))) +
    # dodge values less than 1 decrease the space between the grouped bars
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(dodge_value)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(), ggplot2::aes(color = outlier_trialnums)) +
    # just use grayscale for hemisphere differences (reserve colors for subjects)
    # re-scale so that not using black/darkest gray end of the scale
    ggplot2::scale_fill_grey(start=0.8, end=0.4) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    # feeding in columns as characters based on: https://stackoverflow.com/questions/21588096/pass-string-to-facet-grid-ggplot2
    ggplot2::facet_grid(reformulate(roi_varname)) +
    ggplot2::ggtitle(sprintf("red < %d trials", trialnums_cutoff)) +
    mixed_models_jitterbox_theme
}

boxplot_facet_single_var <- function(input_df, x_var, y_var, color_var, facet_var) {
  p <<- NULL
  p <<- ggplot2::ggplot(data = input_df, ggplot2::aes(x = !!sym(x_var), y = !!sym(y_var))) +
    ggplot2::geom_boxplot(ggplot2::aes(color = !!sym(color_var))) +
    ggplot2::facet_wrap(reformulate(facet_var)) +
    mixed_models_jitterbox_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title.x = ggplot2::element_blank())
}

boxplot_no_facet <- function(input_df, x_var, y_var, color_var) {
  p <<- NULL
  p <<- ggplot2::ggplot(data = input_df, ggplot2::aes(x = !!sym(x_var), y = !!sym(y_var))) +
    ggplot2::geom_boxplot(ggplot2::aes(color = !!sym(color_var))) +
    mixed_models_jitterbox_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title.x = ggplot2::element_blank())
}

boxplot_facet_hemi_and_somethingelse <- function(input_df, x_var, y_var, color_var, hemi_facet_var, facet_var) {
  p <<- NULL
  p <<- ggplot2::ggplot(data = input_df, ggplot2::aes(x = !!sym(x_var), y = !!sym(y_var))) +
    ggplot2::geom_boxplot(ggplot2::aes(color = !!sym(color_var))) +
    ggplot2::facet_grid(reformulate(hemi_facet_var, facet_var)) +
    mixed_models_jitterbox_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title.x = ggplot2::element_blank())
}

sfn_2017_theme <- ggplot2::theme(axis.title.x = ggplot2::element_text(size = 30), axis.text.x = ggplot2::element_text(size = 30),
                                 strip.text.x = ggplot2::element_text(size = 28),
                                 # margin based on: http://stackoverflow.com/questions/14487188/increase-distance-between-text-and-title-on-the-y-axis
                                 axis.title.y = ggplot2::element_text(size = 30, margin = ggplot2::margin(0,20,0,0)), axis.text.y = ggplot2::element_text(size = 25),
                                 strip.text.y = ggplot2::element_text(size = 28),
                                 plot.title = ggplot2::element_text(size = 30))
