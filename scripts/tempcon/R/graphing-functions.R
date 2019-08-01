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

plot_permutation_single_roi <- function(input_dat, chisq_perm, p95_thresh, chisq_observed){
  # calling column names based on: https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot (NOT WORKING)
  p <<- NULL
  p <<- input_dat %>%
    ggplot2::ggplot(data = ., ggplot2::aes_string(x = chisq_perm)) +
    ggplot2::geom_histogram() +
    # combining string and non-string aes based on:
    # https://stackoverflow.com/questions/28777626/how-do-i-combine-aes-and-aes-string-options
    # plot the statistical cutoff
    ggplot2::geom_vline(ggplot2::aes_q(xintercept = as.name(p95_thresh), color = "red"), linetype="dashed") +
    # plot the observed p value
    ggplot2::geom_vline(ggplot2::aes_q(xintercept = as.name(chisq_observed), color = "blue")) +
    ggplot2::ggtitle("Permutation for Interaction Effect") +
    ggplot2::xlab("chisq value") +
    # adding a manual legend based on:
    # https://stackoverflow.com/questions/32490043/adding-manual-legend-in-ggplot?rq=1
    ggplot2::scale_color_identity("temporary title", labels = c("observed chisq", "p < 0.05\ncutoff"), guide = "legend") +
    # adding a dashed line in the legend based on:
    # https://stackoverflow.com/questions/31519202/dashed-line-in-ggplot-legend
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    publication_theme +
    # override the theme to get back our x-axis title
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 20),
                   legend.position = "right", legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 12))
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

# based on: https://stackoverflow.com/questions/32936506/store-custom-ggplot-styles-in-object
four_colors_two_groups <- c("#e66101", "#fdb863", "#5e3c99", "#b2abd2")
four_colors_two_groups_v2 <- c("#a6611a", "#dfc27d", "#018571", "#80cdc1")
three_colors_no_groups <- c("darkgoldenrod1", "#0571b0", "#ca0020")
three_colors_no_groups_v2 <- c("gold", "#0571b0", "violetred3")
four_colors_fill_theme <- list(ggplot2::scale_fill_manual(values = four_colors_two_groups),
                               ggplot2::scale_color_manual(values = four_colors_two_groups))
four_colors_v2_fill_theme <- list(ggplot2::scale_fill_manual(values = four_colors_two_groups_v2),
                               ggplot2::scale_color_manual(values = four_colors_two_groups_v2))
three_colors_fill_theme <- list(ggplot2::scale_fill_manual(values = three_colors_no_groups),
                               ggplot2::scale_color_manual(values = three_colors_no_groups))
three_colors_v2_fill_theme <- list(ggplot2::scale_fill_manual(values = three_colors_no_groups_v2),
                                ggplot2::scale_color_manual(values = three_colors_no_groups_v2))
publication_theme <- list(ggplot2::theme_gray(),
                          ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 18),
                                    axis.title.y = ggplot2::element_text(size = 20), axis.text.y = ggplot2::element_text(size = 15),
                                    legend.position = "none",
                                    plot.title = ggplot2::element_text(size = 20)))

# approach of means and individual subject points from here:
# https://drsimonj.svbtle.com/plotting-individual-observations-and-group-means-with-ggplot2
# approach of nested x axis labels from here (see answer 2):
# https://stackoverflow.com/questions/18165863/multirow-axis-labels-with-nested-grouping-variables
meanbars_w_subjscatter_facetted <- function(subj_data_df, group_means_df, x_var, facet_var, y_var, sem_min, sem_max, color_var) {
  p <<- NULL
  p <<- ggplot2::ggplot(data = subj_data_df,
                        ggplot2::aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(color_var))) +
    ggplot2::geom_point(ggplot2::aes(color = !!sym(color_var))) +
    ggplot2::geom_col(data = group_means_df, alpha = 0.5, width = 0.7) +
    ggplot2::geom_errorbar(data = group_means_df,
                           ggplot2::aes(ymin = !!sym(sem_min), ymax = !!sym(sem_max)),
                           width = 0.1) +
    ggplot2::facet_wrap(reformulate(facet_var), strip.position = "bottom") +
    publication_theme +
    ggplot2::theme(panel.spacing = ggplot2::unit(1, "lines"),
                   strip.background = ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.text.x = ggplot2::element_text(size = 14),
                   axis.text.x = ggplot2::element_text(size = 9))
}

meanbars_w_subjscatter <- function(subj_data_df, group_means_df, x_var, y_var, sem_min, sem_max, color_var) {
  p <<- NULL
  p <<- ggplot2::ggplot(data = subj_data_df,
                        ggplot2::aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(color_var))) +
    ggplot2::geom_point(ggplot2::aes(color = !!sym(color_var))) +
    ggplot2::geom_col(data = group_means_df, alpha = 0.5, width = 0.5) +
    ggplot2::geom_errorbar(data = group_means_df,
                           ggplot2::aes(ymin = !!sym(sem_min), ymax = !!sym(sem_max)),
                           width = 0.1) +
    publication_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 13))
}

meanbars_no_scatter <- function(group_means_df, x_var, y_var, sem_min, sem_max, color_var) {
  p <<- NULL
  p <<- ggplot2::ggplot(data = group_means_df,
                        ggplot2::aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(color_var))) +
    ggplot2::geom_col(alpha = 0.5, width = 0.5) +
    ggplot2::geom_errorbar(data = group_means_df,
                           ggplot2::aes(ymin = !!sym(sem_min), ymax = !!sym(sem_max)),
                           width = 0.1) +
    publication_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 13))
}

compute_cor_and_plot <- function(input_df, col1, col2, subj_col){
  pp <- NULL
  stat_vals <- NULL

  # calling column names in a function based on: https://stackoverflow.com/questions/2641653/pass-a-data-frame-column-name-to-a-function
  print(cor.test(input_df[[col1]], input_df[[col2]]))

  stat_vals <<- cor.test(input_df[[col1]], input_df[[col2]])
  # if set values in the plot as global variables, shit breaks
  tmp <- cor.test(input_df[[col1]], input_df[[col2]])
  stat_for_plot <- sprintf("t(%.f) = %0.3f,\nr = %0.3f, p = %s",
                           tmp$parameter, tmp$statistic,
                           tmp$estimate,
                           ifelse(tmp$p.value < 0.05, sprintf("%0.2f*", tmp$p.value), sprintf("%0.2f", tmp$p.value)))

  pp <<- input_df %>%
    ggplot2::ggplot(., ggplot2::aes_string(x = col1, y = col2)) +
    ggplot2::geom_point(ggplot2::aes_string(color = subj_col)) +
    ggplot2::annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = stat_for_plot, size = 7)
}
