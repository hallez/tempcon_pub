#' ---
#' title: Graphs of mixed model data
#' author: Halle R. Dimsdale-Zucker
#' output:
#'  html_document:
#'    toc: true
#'    toc_depth: 4
#'    toc_float:
#'      collapsed: false
#'      smooth_scroll: false
#'    number_sections: true
#'    theme: spacelab
#'    code_folding: "hide"
#' ---

time_script_start <- Sys.time()
tictoc::tic(msg = "Graphs of data from mixed models")

#' # Setup
#+ devtools, warning = FALSE, message = FALSE
devtools::load_all()
devtools::install_github("hallez/halle", subdir="halle")

#+ label = "library"
# add dplyr with library()
# NB: this is non-standard
# correct way would be to make this script into a function,
# store in the R/ directory,
# and use @importFrom dplyr "%>%",
# but this is incompatible with getting in-line results with code in knitr output
library(dplyr)

#+ label = "config file"
#' ## Load in config file
# use packge::function syntax to avoid confusion
config <- yaml::yaml.load_file("../../config.yml")

#+ label = "setup paths"
#' ## Setup paths
project_dir <- ("../../")
analyzed_mri_dir <- paste0(halle::ensure_trailing_slash(config$directories$analyzed_mri))
raw_behavioral_dir <- paste0(halle::ensure_trailing_slash(config$directories$raw_behavioral))
analyzed_behavioral_dir <- paste0(halle::ensure_trailing_slash(config$directories$analyzed_behavioral))
dropbox_dir <- halle::ensure_trailing_slash(config$directories$dropbox_tempcon)
graph_out_dir <- paste0(halle::ensure_trailing_slash(analyzed_behavioral_dir),
                        halle::ensure_trailing_slash("plots"))
dir.create(graph_out_dir, recursive=TRUE)

#' ## Flags
SAVE_GRAPH_FLAG <- 1

#' ## Set script-wise variables
rois_of_interest <- c("CA1", "CA23DG", "SUB", "PRC", "aPHC", "pPHC")
# can add right hemisphere in the future, but trying to just plot what we need for now
hemi_of_interest <- "left"

#' # Load data
load(file.path(analyzed_mri_dir, "z_mm_trial_pairs_all_subj.RData"))
# this file gets created in `rsa-tidy-data-btwn-runs.R`
load(file.path(analyzed_mri_dir, "excluded_subjects_variables.RData"))

#' ## Load all permutation output files
# annoyingly, most of these variables have the same name when loaded so immediately re-name and remove
# also, append a unique ID incase column names overlap between dataframes
# based on: https://stackoverflow.com/questions/35697940/append-suffix-to-colnames
# https://stackoverflow.com/questions/29948876/adding-prefix-or-suffix-to-most-data-frame-variable-names-in-piped-r-workflow?rq=1
load(file.path(analyzed_mri_dir, "perms_list-samediff-halves_quest-samediff_chisq-perms.RData"))
lsdh_qsd <- same_diff_list_halves_truncated %>%
  dplyr::rename_at(vars(-roi, -nperm),function(x) paste0(x,".lsdh_qsd"))
rm(same_diff_list_halves_truncated)

load(file.path(analyzed_mri_dir, "perms_list-samediff-halves_quest-samediff-no-hemi-in-models_chisq-perms.RData"))
lsdh_qsd_no_hemi <- same_diff_list_halves_truncated %>%
  dplyr::rename_at(vars(-roi, -nperm, -hemi),function(x) paste0(x,".lsdh_qsd_no_hemi"))
rm(same_diff_list_halves_truncated)

load(file.path(analyzed_mri_dir, "perms_list-samediff_quest-samediff_chisq-perms.RData"))
lsd_qsd <- same_diff_list_sd_truncated %>%
  dplyr::rename_at(vars(-roi, -nperm),function(x) paste0(x,".lsd_qsd"))
rm(same_diff_list_sd_truncated)

#' ## Merge permutation files
all_perm_chisq <- lsdh_qsd %>%
  dplyr::left_join(lsdh_qsd_no_hemi, by = c("nperm", "roi")) %>%
  dplyr::left_join(lsd_qsd, by = c("nperm", "roi")) %>%
  dplyr::mutate(roi_fact = as.factor(roi),
  # NB: if any level is missing, it will get filled with NA
  roi_fact_ordered = factor(roi_fact, levels = c("CA1_body", "CA2_3_DG_body", "ERC", "subiculum_body",
                                               "alEC", "pmEC", "PRC", "aPHC", "pPHC",
                                               "whole_hippo", "whole_hippo_head", "whole_hippo_body", "whole_hippo_tail")),
  # in order to preserve the ordering, have to again use `levels` (i.e., it doesn't persist from `roi_fact_ordered`)
  roi_formatted = as.factor(car::recode(roi_fact_ordered, c("'CA1_body' = 'CA1';'CA2_3_DG_body' = 'CA23DG'; 'subiculum_body' = 'SUB'"),
                                        levels = c("CA1", "CA23DG", "SUB",
                                                   "alEC", "pmEC", "PRC", "aPHC", "pPHC",
                                                   "whole_hippo", "whole_hippo_head", "whole_hippo_body", "whole_hippo_tail", "ERC"))))

#' ## Filter z_mm based on minimum trial numbers requirement
z_mm_filt <- z_mm %>%
  dplyr::ungroup() %>%
  dplyr::filter(!subj_fact %in% exclude_subjects) %>%
  # this will filter out subjects who don't have enough trials, but won't fail if there are 0 subjects to exclude
  dplyr::filter(!subj_fact %in% low_trials_questXlist_subj) %>%
  dplyr::filter(item_mem == "rHit")

#' ## Quickly check the data that gets passed to analyses
message(sprintf("Now analyzing data from %d subjects after filtering on %d min trials.",
                length(unique(z_mm_filt$subj_fact)), min_trials))

#' ### Excluded subjects
message("information included in `z_mm_filt`")
message("excluded subjects: ")
exclude_subjects

#' ### Included subjects
message("included subjects:")
unique(z_mm_filt$subj)
length(unique(z_mm_filt$subj))

#' ### Column names, rois, and conditions
colnames(z_mm_filt)
unique(z_mm_filt$roi)
unique(z_mm_filt$condition)

#' # Plots for list halves (same/diff), question (same/diff)
for(iroi in 1:length(rois_of_interest)){
  cur_roi <- rois_of_interest[iroi]

  if(!is.element(cur_roi, unique(z_mm_filt$roi_formatted))) {
    message(sprintf("current roi %s not in dataframe, continuing.\n", cur_roi))
    next
  }

  # --- bar plots of means
  listXquest <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi_formatted == cur_roi) %>%
    # revalue labels for graphs
    dplyr::mutate(list_same_diff_half_recode = dplyr::recode(list_same_diff_half, "same_half" = "similar temporal\ncontext", "diff_half" = "different temporal\ncontext"),
                  question_recode = dplyr::recode_factor(question, "same" = "same cognitive context", "diff" = "different cognitive context",
                                                         .ordered = TRUE))

  listXquest_group_means <- listXquest %>%
    # first take the mean w/in subject
    dplyr::group_by(subj, list_same_diff_half_recode, question_recode, roi_formatted) %>%
    dplyr::summarise(subj_mean_r = mean(r, na.rm = TRUE)) %>%
    dplyr::group_by(list_same_diff_half_recode, question_recode, roi_formatted) %>%
    dplyr::summarise(mean_r = mean(subj_mean_r, na.rm = TRUE),
                     median_r = median(subj_mean_r, na.rm = TRUE),
                     num_subj = length(subj),
                     sd_r = sd(subj_mean_r, na.rm = TRUE),
                     sem_r = sd_r / sqrt(num_subj),
                     sem_ymin = mean_r - sem_r,
                     sem_ymax = mean_r + sem_r) %>%
    # for the sake of coloring the bars, create a new combination variable
    dplyr::mutate(list_quest = paste0(list_same_diff_half_recode, question_recode))

  listXquest_subj_means <- listXquest %>%
    dplyr::group_by(subj, list_same_diff_half_recode, question_recode, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE),
                     median_r = median(r, na.rm = TRUE)) %>%
    # for the sake of coloring the bars, create a new combination variable
    dplyr::mutate(list_quest = paste0(list_same_diff_half_recode, question_recode))

  meanbars_w_subjscatter_facetted(listXquest_subj_means, listXquest_group_means, "list_same_diff_half_recode", "question_recode", "mean_r", "sem_ymin", "sem_ymax", "list_quest")
  mean_w_scatt <- NULL
  mean_w_scatt <- p +
    ggplot2::ylab("mean pattern similarity") +
    ggplot2::ggtitle("Mean Pattern Similarity by Condition") +
    four_colors_fill_theme
  print(mean_w_scatt)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-w-subjjitt_listhalf-samediff_quest-samediff_%s.pdf", cur_roi)),
                    height = 4, width = 6)
  }

  # --- plot subtraction between conditions
  quest_diffs <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi_formatted == cur_roi) %>%
    # since each trial pair could only have be same or different questions,
    # first get an average for each condition
    dplyr::group_by(subj, question, roi_formatted) %>%
    dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
    # now, spread column of interest so can subtract
    tidyr::spread(key = "question", value = "mean_r") %>%
    # start subtracting conditions to get difference scores
    dplyr::mutate(samequest_minus_diffquest = same - diff)

  list_diffs <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi_formatted == cur_roi) %>%
    dplyr::group_by(subj, list_same_diff_half, roi_formatted) %>%
    dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
    tidyr::spread(key = "list_same_diff_half", value = "mean_r") %>%
    dplyr::mutate(samelisthalf_minus_difflisthalf = same_half - diff_half)

  listXquest_diffs <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi_formatted == cur_roi) %>%
    tidyr::unite(quest_list, question, list_same_diff_half) %>%
    dplyr::group_by(subj, quest_list, roi_formatted) %>%
    dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
    tidyr::spread(key = "quest_list", value = "mean_r") %>%
    dplyr::mutate(samequest.samehalf_minus_difftask.samehalf = same_same_half - diff_same_half,
                  samequest.diffhalf_minus_difftask.diffhalf = same_diff_half - diff_diff_half,
                  ssh.dsh_sdh.ddh = samequest.samehalf_minus_difftask.samehalf - samequest.diffhalf_minus_difftask.diffhalf)

  # now, lump all of these dataframes together so can be on same plot
  all_cond_diffs <- quest_diffs %>%
    dplyr::left_join(list_diffs, by = c("subj", "roi_formatted")) %>%
    dplyr::left_join(listXquest_diffs, by = c("subj", "roi_formatted")) %>%
    # rename so looks better on plots
    dplyr::rename(questME = samequest_minus_diffquest,
                  listhalfME = samelisthalf_minus_difflisthalf,
                  questXlisthalf = ssh.dsh_sdh.ddh) %>%
    # lump things we care about into a single column
    tidyr::gather(key = condition, value = diff_score, -subj, -roi_formatted) %>%
    # rename again
    dplyr::mutate(condition_recode = dplyr::recode_factor(condition, "listhalfME" = "temporal\ncontext\neffect",
                                                          "questME" = "cognitive\ncontext\neffect",
                                                          "questXlisthalf" = "temporal and cognitive\ncontext\ninteraction"),
                  .ordered = TRUE) %>%
    # only retain the conditions we want to plot
    dplyr::filter(condition_recode %in% c("cognitive\ncontext\neffect", "temporal\ncontext\neffect", "temporal and cognitive\ncontext\ninteraction")) %>%
    droplevels()

  # save out differences so can use in correlations
  save("all_cond_diffs", file = file.path(analyzed_mri_dir, sprintf("condition_diffs_list-samediff-halves_quest-samediff_%s.RData", cur_roi)))

  all_cond_diffs_group <- all_cond_diffs %>%
    # first take the mean w/in subject
    dplyr::group_by(subj, condition_recode, roi_formatted) %>%
    dplyr::summarise(subj_mean_diff = mean(diff_score, na.rm = TRUE)) %>%
    dplyr::group_by(condition_recode, roi_formatted) %>%
    dplyr::summarise(mean_diff = mean(subj_mean_diff, na.rm = TRUE),
                     median_diff = median(subj_mean_diff, na.rm = TRUE),
                     num_subj = length(subj),
                     sd_diff = sd(subj_mean_diff, na.rm = TRUE),
                     sem_diff = sd_diff / sqrt(num_subj),
                     sem_ymin = mean_diff - sem_diff,
                     sem_ymax = mean_diff + sem_diff)

  all_cond_diffs_subj <- all_cond_diffs %>%
    dplyr::group_by(subj, condition_recode, roi_formatted) %>%
    dplyr::summarize(mean_diff = mean(diff_score, na.rm = TRUE),
                     median_diff = median(diff_score, na.rm = TRUE))

  meanbars_w_subjscatter(all_cond_diffs_subj, all_cond_diffs_group, "condition_recode", "mean_diff", "sem_ymin", "sem_ymax", "condition_recode")
  diff_w_scatt <- NULL
  diff_w_scatt <- p +
    ggplot2::ylab("mean pattern similarity\ndifferences") +
    ggplot2::ggtitle("Differences Between Conditions") +
    three_colors_fill_theme
  print(diff_w_scatt)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-w-subjjitt_list-halves_condition-differences_%s.pdf", cur_roi)),
                    height = 4, width = 6)
  }

  meanbars_no_scatter(all_cond_diffs_group, "condition_recode", "mean_diff", "sem_ymin", "sem_ymax", "condition_recode")
  diff_no_scatt <- NULL
  diff_no_scatt <- p +
    ggplot2::ylab("mean pattern similarity\ndifferences") +
    ggplot2::ggtitle("Differences Between Conditions") +
    three_colors_fill_theme
  print(diff_no_scatt)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-no_subjjitt_list-halves_condition-differences_%s.pdf", cur_roi)),
                    height = 4, width = 6)
  }

  # --- plot permutation
  all_perm_chisq_roi <- all_perm_chisq %>%
    dplyr::filter(roi_formatted == cur_roi)

  plot_permutation_single_roi(all_perm_chisq_roi, "chisq_listhalvessamediff.listhalvesXsamediff_perm.lsdh_qsd", "chisq_listhalvessamediff.listhalvesXsamediff_perm_p95.lsdh_qsd", "chisq_listhalvesquest.listhalvesXquest.lsdh_qsd")
  perm_plot <- NULL
  perm_plot <- p
  print(perm_plot)

  # --- mush all the plots together into a figure
  # use `hjust` to ensure that labels don't overlap w/ figure titles
  fig_all_diff_no_scatt <- NULL
  fig_all_diff_no_scatt <- cowplot::plot_grid(mean_w_scatt, diff_no_scatt, perm_plot, align = "v", labels = "AUTO", hjust = -2)
  print(fig_all_diff_no_scatt)

  # save as a .png and also as a svg
  all_file_ext <- c("png", "svg")
  for(iext in 1:length(all_file_ext)){
    cur_file_ext <- all_file_ext[iext]
    fig_name <- sprintf("patt_sim_figure_diff-no-scatt_list-halves_%s.%s", cur_roi, cur_file_ext)
    cowplot::save_plot(file.path(graph_out_dir, fig_name),
                       fig_all_diff_no_scatt,
                       ncol = 2, # we're saving a grid plot of 2 columns
                       nrow = 2, # and 2 rows
                       base_aspect_ratio = 1.5)
  }

  # also save out a version w/ subject scatter
  fig_all_subj_scatt <- NULL
  fig_all_subj_scatt <- cowplot::plot_grid(mean_w_scatt, diff_w_scatt, perm_plot, align = "v", labels = "AUTO", hjust = -2)
  print(fig_all_subj_scatt)

  # save as a .png and also as a svg
  all_file_ext <- c("png", "svg")
  for(iext in 1:length(all_file_ext)){
    cur_file_ext <- all_file_ext[iext]
    fig_name <- sprintf("patt_sim_figure_diff-w-scatt_list-halves_%s.%s", cur_roi, cur_file_ext)
    cowplot::save_plot(file.path(graph_out_dir, fig_name),
                       fig_all_subj_scatt,
                       ncol = 2, # we're saving a grid plot of 2 columns
                       nrow = 2, # and 2 rows
                       base_aspect_ratio = 1.5)
  }

} # for(iroi

#' ## Split by hemisphere
for(iroi in 1:length(rois_of_interest)){
  for(ihemi in 1:length(hemi_of_interest)){
    cur_roi <- rois_of_interest[iroi]
    cur_hemi <- hemi_of_interest[ihemi]

    if(!is.element(cur_roi, unique(z_mm_filt$roi_formatted))) {
      message(sprintf("current roi %s not in dataframe, continuing.\n", cur_roi))
      next
    }

    # --- bar plots of means
    listXquest <- z_mm_filt %>%
      dplyr::filter(list_lag >= 0,
                    list_same_diff_half %in% c("same_half", "diff_half"),
                    question %in% c("same", "diff"),
                    hemi == cur_hemi) %>%
      dplyr::filter(roi_formatted == cur_roi) %>%
      # revalue labels for graphs
      dplyr::mutate(list_same_diff_half_recode = dplyr::recode(list_same_diff_half, "same_half" = "similar temporal\ncontext", "diff_half" = "different temporal\ncontext"),
                    question_recode = dplyr::recode_factor(question, "same" = "same cognitive context", "diff" = "different cognitive context",
                                                           .ordered = TRUE))

    listXquest_group_means <- listXquest %>%
      # first take the mean w/in subject
      dplyr::group_by(subj, list_same_diff_half_recode, question_recode, roi_formatted) %>%
      dplyr::summarise(subj_mean_r = mean(r, na.rm = TRUE)) %>%
      dplyr::group_by(list_same_diff_half_recode, question_recode, roi_formatted) %>%
      dplyr::summarise(mean_r = mean(subj_mean_r, na.rm = TRUE),
                       median_r = median(subj_mean_r, na.rm = TRUE),
                       num_subj = length(subj),
                       sd_r = sd(subj_mean_r, na.rm = TRUE),
                       sem_r = sd_r / sqrt(num_subj),
                       sem_ymin = mean_r - sem_r,
                       sem_ymax = mean_r + sem_r) %>%
      # for the sake of coloring the bars, create a new combination variable
      dplyr::mutate(list_quest = paste0(list_same_diff_half_recode, question_recode))

    listXquest_subj_means <- listXquest %>%
      dplyr::group_by(subj, list_same_diff_half_recode, question_recode, roi_formatted) %>%
      dplyr::summarize(mean_r = mean(r, na.rm = TRUE),
                       median_r = median(r, na.rm = TRUE)) %>%
      # for the sake of coloring the bars, create a new combination variable
      dplyr::mutate(list_quest = paste0(list_same_diff_half_recode, question_recode))

    meanbars_w_subjscatter_facetted(listXquest_subj_means, listXquest_group_means, "list_same_diff_half_recode", "question_recode", "mean_r", "sem_ymin", "sem_ymax", "list_quest")
    mean_w_scatt <- NULL
    mean_w_scatt <- p +
      ggplot2::ylab("mean pattern similarity") +
      ggplot2::ggtitle("Mean Pattern Similarity by Condition") +
      four_colors_fill_theme
    print(mean_w_scatt)

    if(SAVE_GRAPH_FLAG == 1){
      ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-w-subjjitt_listhalf-samediff_quest-samediff_%s_%s.pdf", cur_roi, cur_hemi)),
                      height = 4, width = 6)
    }

    # --- plot subtraction between conditions
    quest_diffs <- z_mm_filt %>%
      dplyr::filter(list_lag >= 0,
                    list_same_diff_half %in% c("same_half", "diff_half"),
                    question %in% c("same", "diff"),
                    hemi == cur_hemi) %>%
      dplyr::filter(roi_formatted == cur_roi) %>%
      # since each trial pair could only have be same or different questions,
      # first get an average for each condition
      dplyr::group_by(subj, question, roi_formatted) %>%
      dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
      # now, spread column of interest so can subtract
      tidyr::spread(key = "question", value = "mean_r") %>%
      # start subtracting conditions to get difference scores
      dplyr::mutate(samequest_minus_diffquest = same - diff)

    list_diffs <- z_mm_filt %>%
      dplyr::filter(list_lag >= 0,
                    list_same_diff_half %in% c("same_half", "diff_half"),
                    question %in% c("same", "diff"),
                    hemi == cur_hemi) %>%
      dplyr::filter(roi_formatted == cur_roi) %>%
      dplyr::group_by(subj, list_same_diff_half, roi_formatted) %>%
      dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
      tidyr::spread(key = "list_same_diff_half", value = "mean_r") %>%
      dplyr::mutate(samelisthalf_minus_difflisthalf = same_half - diff_half)

    listXquest_diffs <- z_mm_filt %>%
      dplyr::filter(list_lag >= 0,
                    list_same_diff_half %in% c("same_half", "diff_half"),
                    question %in% c("same", "diff"),
                    hemi == cur_hemi) %>%
      dplyr::filter(roi_formatted == cur_roi) %>%
      tidyr::unite(quest_list, question, list_same_diff_half) %>%
      dplyr::group_by(subj, quest_list, roi_formatted) %>%
      dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
      tidyr::spread(key = "quest_list", value = "mean_r") %>%
      dplyr::mutate(samequest.samehalf_minus_difftask.samehalf = same_same_half - diff_same_half,
                    samequest.diffhalf_minus_difftask.diffhalf = same_diff_half - diff_diff_half,
                    ssh.dsh_sdh.ddh = samequest.samehalf_minus_difftask.samehalf - samequest.diffhalf_minus_difftask.diffhalf)

    # now, lump all of these dataframes together so can be on same plot
    all_cond_diffs <- quest_diffs %>%
      dplyr::left_join(list_diffs, by = c("subj", "roi_formatted")) %>%
      dplyr::left_join(listXquest_diffs, by = c("subj", "roi_formatted")) %>%
      # rename so looks better on plots
      dplyr::rename(questME = samequest_minus_diffquest,
                    listhalfME = samelisthalf_minus_difflisthalf,
                    questXlisthalf = ssh.dsh_sdh.ddh) %>%
      # lump things we care about into a single column
      tidyr::gather(key = condition, value = diff_score, -subj, -roi_formatted) %>%
      # rename again
      dplyr::mutate(condition_recode = dplyr::recode_factor(condition, "listhalfME" = "temporal\ncontext\neffect",
                                                            "questME" = "cognitive\ncontext\neffect",
                                                            "questXlisthalf" = "temporal and cognitive\ncontext\ninteraction"),
                    .ordered = TRUE) %>%
      # only retain the conditions we want to plot
      dplyr::filter(condition_recode %in% c("cognitive\ncontext\neffect", "temporal\ncontext\neffect", "temporal and cognitive\ncontext\ninteraction")) %>%
      droplevels()

    # save out differences so can use in correlations
    save("all_cond_diffs", file = file.path(analyzed_mri_dir, sprintf("condition_diffs_list-samediff-halves_quest-samediff_%s_%s.RData", cur_roi, cur_hemi)))

    all_cond_diffs_group <- all_cond_diffs %>%
      # first take the mean w/in subject
      dplyr::group_by(subj, condition_recode, roi_formatted) %>%
      dplyr::summarise(subj_mean_diff = mean(diff_score, na.rm = TRUE)) %>%
      dplyr::group_by(condition_recode, roi_formatted) %>%
      dplyr::summarise(mean_diff = mean(subj_mean_diff, na.rm = TRUE),
                       median_diff = median(subj_mean_diff, na.rm = TRUE),
                       num_subj = length(subj),
                       sd_diff = sd(subj_mean_diff, na.rm = TRUE),
                       sem_diff = sd_diff / sqrt(num_subj),
                       sem_ymin = mean_diff - sem_diff,
                       sem_ymax = mean_diff + sem_diff)

    all_cond_diffs_subj <- all_cond_diffs %>%
      dplyr::group_by(subj, condition_recode, roi_formatted) %>%
      dplyr::summarize(mean_diff = mean(diff_score, na.rm = TRUE),
                       median_diff = median(diff_score, na.rm = TRUE))

    meanbars_w_subjscatter(all_cond_diffs_subj, all_cond_diffs_group, "condition_recode", "mean_diff", "sem_ymin", "sem_ymax", "condition_recode")
    diff_w_scatt <- NULL
    diff_w_scatt <- p +
      ggplot2::ylab("mean pattern similarity\ndifferences") +
      ggplot2::ggtitle("Differences Between Conditions") +
      three_colors_fill_theme
    print(diff_w_scatt)

    if(SAVE_GRAPH_FLAG == 1){
      ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-w-subjjitt_list-halves_condition-differences_%s_%s.pdf", cur_roi, cur_hemi)),
                      height = 4, width = 6)
    }

    meanbars_no_scatter(all_cond_diffs_group, "condition_recode", "mean_diff", "sem_ymin", "sem_ymax", "condition_recode")
    diff_no_scatt <- NULL
    diff_no_scatt <- p +
      ggplot2::ylab("mean pattern similarity\ndifferences") +
      ggplot2::ggtitle("Differences Between Conditions") +
      three_colors_fill_theme
    print(diff_no_scatt)

    if(SAVE_GRAPH_FLAG == 1){
      ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-no_subjjitt_list-halves_condition-differences_%s_%s.pdf", cur_roi, cur_hemi)),
                      height = 4, width = 6)
    }

    # --- plot permutation
    all_perm_chisq_roi <- all_perm_chisq %>%
      dplyr::filter(roi_formatted == cur_roi) %>%
      dplyr::filter(hemi == cur_hemi)

    plot_permutation_single_roi(all_perm_chisq_roi, "chisq_listhalvessamediff.listhalvesXsamediff_perm.lsdh_qsd", "chisq_listhalvessamediff.listhalvesXsamediff_perm_p95.lsdh_qsd", "chisq_listhalvesquest.listhalvesXquest.lsdh_qsd")
    perm_plot <- NULL
    perm_plot <- p
    print(perm_plot)

    # --- mush all the plots together into a figure
    # use `hjust` to ensure that labels don't overlap w/ figure titles
    fig_all_diff_no_scatt <- NULL
    fig_all_diff_no_scatt <- cowplot::plot_grid(mean_w_scatt, diff_no_scatt, perm_plot, align = "v", labels = "AUTO", hjust = -2)
    print(fig_all_diff_no_scatt)

    # save as a .png and also as a svg
    all_file_ext <- c("png", "svg")
    for(iext in 1:length(all_file_ext)){
      cur_file_ext <- all_file_ext[iext]
      fig_name <- sprintf("patt_sim_figure_diff-no-scatt_list-halves_%s_%s.%s", cur_roi, cur_hemi, cur_file_ext)
      cowplot::save_plot(file.path(graph_out_dir, fig_name),
                         fig_all_diff_no_scatt,
                         ncol = 2, # we're saving a grid plot of 2 columns
                         nrow = 2, # and 2 rows
                         base_aspect_ratio = 1.5)
    }

    # also save out a version w/ subject scatter
    fig_all_subj_scatt <- NULL
    fig_all_subj_scatt <- cowplot::plot_grid(mean_w_scatt, diff_w_scatt, perm_plot, align = "v", labels = "AUTO", hjust = -2)
    print(fig_all_subj_scatt)

    # save as a .png and also as a svg
    all_file_ext <- c("png", "svg")
    for(iext in 1:length(all_file_ext)){
      cur_file_ext <- all_file_ext[iext]
      fig_name <- sprintf("patt_sim_figure_diff-w-scatt_list-halves_%s_%s.%s", cur_roi, cur_hemi, cur_file_ext)
      cowplot::save_plot(file.path(graph_out_dir, fig_name),
                         fig_all_subj_scatt,
                         ncol = 2, # we're saving a grid plot of 2 columns
                         nrow = 2, # and 2 rows
                         base_aspect_ratio = 1.5)
    }
  }# for(ihemi
} # for(iroi

#' # Plots for list (same/diff), question (same/diff)
for(iroi in 1:length(rois_of_interest)){
  cur_roi <- rois_of_interest[iroi]

  if(!is.element(cur_roi, unique(z_mm_filt$roi_formatted))) {
    message(sprintf("current roi %s not in dataframe, continuing.\n", cur_roi))
    next
  }

  # --- bar plots of means
  listXquest <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff %in% c("same_list", "diff_list"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi_formatted == cur_roi) %>%
    # revalue labels for graphs
    dplyr::mutate(list_same_diff_recode = dplyr::recode(list_same_diff, "same_list" = "same temporal\ncontext", "diff_list" = "different temporal\ncontext"),
                  question_recode = dplyr::recode_factor(question, "same" = "same cognitive context", "diff" = "different cognitive context",
                                                         .ordered = TRUE))

  listXquest_group_means <- listXquest %>%
    # first take the mean w/in subject
    dplyr::group_by(subj, list_same_diff_recode, question_recode, roi_formatted) %>%
    dplyr::summarise(subj_mean_r = mean(r, na.rm = TRUE)) %>%
    dplyr::group_by(list_same_diff_recode, question_recode, roi_formatted) %>%
    dplyr::summarise(mean_r = mean(subj_mean_r, na.rm = TRUE),
                     median_r = median(subj_mean_r, na.rm = TRUE),
                     num_subj = length(subj),
                     sd_r = sd(subj_mean_r, na.rm = TRUE),
                     sem_r = sd_r / sqrt(num_subj),
                     sem_ymin = mean_r - sem_r,
                     sem_ymax = mean_r + sem_r) %>%
    # for the sake of coloring the bars, create a new combination variable
    dplyr::mutate(list_quest = paste0(list_same_diff_recode, question_recode))

  listXquest_subj_means <- listXquest %>%
    dplyr::group_by(subj, list_same_diff_recode, question_recode, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE),
                     median_r = median(r, na.rm = TRUE)) %>%
    # for the sake of coloring the bars, create a new combination variable
    dplyr::mutate(list_quest = paste0(list_same_diff_recode, question_recode))

  meanbars_w_subjscatter_facetted(listXquest_subj_means, listXquest_group_means, "list_same_diff_recode", "question_recode", "mean_r", "sem_ymin", "sem_ymax", "list_quest")
  mean_w_scatt <- NULL
  mean_w_scatt <- p +
    ggplot2::ylab("mean pattern similarity") +
    ggplot2::ggtitle("Mean Pattern Similarity by Condition") +
    four_colors_v2_fill_theme
  print(mean_w_scatt)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-w-subjjitt_list-samediff_quest-samediff_%s.pdf", cur_roi)),
                    height = 4, width = 6)
  }

  # --- plot subtraction between conditions
  quest_diffs <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff %in% c("same_list", "diff_list"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi_formatted == cur_roi) %>%
    # since each trial pair could only have be same or different questions,
    # first get an average for each condition
    dplyr::group_by(subj, question, roi_formatted) %>%
    dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
    # now, spread column of interest so can subtract
    tidyr::spread(key = "question", value = "mean_r") %>%
    # start subtracting conditions to get difference scores
    dplyr::mutate(samequest_minus_diffquest = same - diff)

  list_diffs <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff %in% c("same_list", "diff_list"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi_formatted == cur_roi) %>%
    dplyr::group_by(subj, list_same_diff, roi_formatted) %>%
    dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
    tidyr::spread(key = "list_same_diff", value = "mean_r") %>%
    dplyr::mutate(samelist_minus_difflist = same_list - diff_list)

  listXquest_diffs <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff %in% c("same_list", "diff_list"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi_formatted == cur_roi) %>%
    tidyr::unite(quest_list, question, list_same_diff) %>%
    dplyr::group_by(subj, quest_list, roi_formatted) %>%
    dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
    tidyr::spread(key = "quest_list", value = "mean_r") %>%
    dplyr::mutate(samequest.samelist_minus_difftask.samelist = same_same_list - diff_same_list,
                  samequest.difflist_minus_difftask.difflist = same_diff_list - diff_diff_list,
                  ssl.dsl_sdl.ddl = samequest.samelist_minus_difftask.samelist - samequest.difflist_minus_difftask.difflist)

  # now, lump all of these dataframes together so can be on same plot
  all_cond_diffs <- quest_diffs %>%
    dplyr::left_join(list_diffs, by = c("subj", "roi_formatted")) %>%
    dplyr::left_join(listXquest_diffs, by = c("subj", "roi_formatted")) %>%
    # rename so looks better on plots
    dplyr::rename(questME = samequest_minus_diffquest,
                  listME = samelist_minus_difflist,
                  questXlist = ssl.dsl_sdl.ddl) %>%
    # lump things we care about into a single column
    tidyr::gather(key = condition, value = diff_score, -subj, -roi_formatted) %>%
    # rename again
    dplyr::mutate(condition_recode = dplyr::recode_factor(condition, "listME" = "temporal\ncontext\neffect",
                                                          "questME" = "cognitive\ncontext\neffect",
                                                          "questXlist" = "temporal and cognitive\ncontext\ninteraction"),
                  .ordered = TRUE) %>%
    # only retain the conditions we want to plot
    dplyr::filter(condition_recode %in% c("cognitive\ncontext\neffect", "temporal\ncontext\neffect", "temporal and cognitive\ncontext\ninteraction")) %>%
    droplevels()

  # save out differences so can use in correlations
  save("all_cond_diffs", file = file.path(analyzed_mri_dir, sprintf("condition_diffs_list-samediff_quest-samediff_%s.RData", cur_roi)))

  all_cond_diffs_group <- all_cond_diffs %>%
    # first take the mean w/in subject
    dplyr::group_by(subj, condition_recode, roi_formatted) %>%
    dplyr::summarise(subj_mean_diff = mean(diff_score, na.rm = TRUE)) %>%
    dplyr::group_by(condition_recode, roi_formatted) %>%
    dplyr::summarise(mean_diff = mean(subj_mean_diff, na.rm = TRUE),
                     median_diff = median(subj_mean_diff, na.rm = TRUE),
                     num_subj = length(subj),
                     sd_diff = sd(subj_mean_diff, na.rm = TRUE),
                     sem_diff = sd_diff / sqrt(num_subj),
                     sem_ymin = mean_diff - sem_diff,
                     sem_ymax = mean_diff + sem_diff)

  all_cond_diffs_subj <- all_cond_diffs %>%
    dplyr::group_by(subj, condition_recode, roi_formatted) %>%
    dplyr::summarize(mean_diff = mean(diff_score, na.rm = TRUE),
                     median_diff = median(diff_score, na.rm = TRUE))

  meanbars_w_subjscatter(all_cond_diffs_subj, all_cond_diffs_group, "condition_recode", "mean_diff", "sem_ymin", "sem_ymax", "condition_recode")
  diff_w_scatt <- NULL
  diff_w_scatt <- p +
    ggplot2::ylab("mean pattern similarity\ndifferences") +
    ggplot2::ggtitle("Differences Between Conditions") +
    three_colors_v2_fill_theme
  print(diff_w_scatt)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-w-subjjitt_by-list_condition-differences_%s.pdf", cur_roi)),
                    height = 4, width = 6)
  }

  meanbars_no_scatter(all_cond_diffs_group, "condition_recode", "mean_diff", "sem_ymin", "sem_ymax", "condition_recode")
  diff_no_scatt <- NULL
  diff_no_scatt <- p +
    ggplot2::ylab("mean pattern similarity\ndifferences") +
    ggplot2::ggtitle("Differences Between Conditions") +
    three_colors_v2_fill_theme
  print(diff_no_scatt)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("mm-means-no_subjjitt_by-list_condition-differences_%s.pdf", cur_roi)),
                    height = 4, width = 6)
  }

  # --- plot permutation
  all_perm_chisq_roi <- all_perm_chisq %>%
    dplyr::filter(roi_formatted == cur_roi)

  plot_permutation_single_roi(all_perm_chisq_roi, "chisq_listsdsamediff.listsdXsamediff_perm.lsd_qsd", "chisq_listsdsamediff.listsdXsamediff_perm_p95.lsd_qsd", "chisq_listsdquest.listsdXquest.lsd_qsd")
  perm_plot <- NULL
  perm_plot <- p
  print(perm_plot)

  # --- mush all the plots together into a figure
  # use `hjust` to ensure that labels don't overlap w/ figure titles
  fig_all_diff_no_scatt <- NULL
  fig_all_diff_no_scatt <- cowplot::plot_grid(mean_w_scatt, diff_no_scatt, perm_plot, align = "v", labels = "AUTO", hjust = -2)
  print(fig_all_diff_no_scatt)

  # save as a .png and also as a svg
  all_file_ext <- c("png", "svg")
  for(iext in 1:length(all_file_ext)){
    cur_file_ext <- all_file_ext[iext]
    fig_name <- sprintf("patt_sim_figure_diff-no-scatt_by-list_%s.%s", cur_roi, cur_file_ext)
    cowplot::save_plot(file.path(graph_out_dir, fig_name),
                       fig_all_diff_no_scatt,
                       ncol = 2, # we're saving a grid plot of 2 columns
                       nrow = 2, # and 2 rows
                       base_aspect_ratio = 1.5)
  }

  # also save out a version w/ subject scatter
  fig_all_subj_scatt <- NULL
  fig_all_subj_scatt <- cowplot::plot_grid(mean_w_scatt, diff_w_scatt, perm_plot, align = "v", labels = "AUTO", hjust = -2)
  print(fig_all_subj_scatt)

  # save as a .png and also as a svg
  all_file_ext <- c("png", "svg")
  for(iext in 1:length(all_file_ext)){
    cur_file_ext <- all_file_ext[iext]
    fig_name <- sprintf("patt_sim_figure_diff-w-scatt_by-list_%s.%s", cur_roi, cur_file_ext)
    cowplot::save_plot(file.path(graph_out_dir, fig_name),
                       fig_all_subj_scatt,
                       ncol = 2, # we're saving a grid plot of 2 columns
                       nrow = 2, # and 2 rows
                       base_aspect_ratio = 1.5)
  }

} # for(iroi

#' # Timer
tictoc::toc()

#' # Print out R session info
sessionInfo()
sprintf("It took %.2f minutes to run all of this.", Sys.time() - time_script_start)
tictoc::toc()
