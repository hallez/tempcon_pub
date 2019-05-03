#' ---
#' title: Mixed models TempCon data &mdash; Between runs
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
numperm <- 2
alpha_value <- 0.05
hc_rois <- c("CA1_body", "CA2_3_DG_body", "subiculum_body")
mtlc_rois <- c("PRC", "aPHC", "pPHC")
rois_of_interest <- NULL
rois_of_interest[[1]] <- hc_rois
rois_of_interest[[2]] <- mtlc_rois

#' # Load data
load(file.path(analyzed_mri_dir, "z_mm_trial_pairs_all_subj.RData"))
# this file gets created in `rsa-tidy-data-btwn-runs.R`
load(file.path(analyzed_mri_dir, "excluded_subjects_variables.RData"))

# for now, also remove other subjects who are causing problems (s3, s15, s32); hopefully add them back in eventually
exclude_subjects <- c(exclude_subjects, "s3", "s15", "s32")

#' ## Filter z_mm based on minimum trial numbers requirement
z_mm_filt <- z_mm %>%
  dplyr::ungroup() %>%
  dplyr::filter(!subj_fact %in% exclude_subjects) %>%
  # this will filter out subjects who don't have enough trials, but won't fail if there are 0 subjects to exclude
  dplyr::filter(!subj_fact %in% low_trials_questXlist_subj) %>%
  dplyr::filter(item_mem == "rHit")

#' ## Separately filter so have rhits as well as fhits&misses
z_mm_RxFM <- z_mm %>%
  dplyr::ungroup() %>%
  dplyr::filter(!subj_fact %in% exclude_subjects) %>%
  dplyr::filter(!subj_fact %in% low_trials_questXlist_fHitANDmisses_subj) %>%
  dplyr::filter(item_mem %in% c("fHitANDmiss", "rHit"))

#' ## Quickly check the data that gets passed to analyses
message("information included in `z_mm_filt`")
message("excluded subjects: ")
exclude_subjects

message("included subjects:")
unique(z_mm_filt$subj)
length(unique(z_mm_filt$subj))

colnames(z_mm_filt)
unique(z_mm_filt$roi)
unique(z_mm_filt$condition)

message(sprintf("Now analyzing data from %d subjects after filtering on %d min trials.",
                length(unique(z_mm_filt$subj_fact)), min_trials))

#' # Notes about mixed models
# although some have recommended always fitting a "maximal" model
# (see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3881361/)
# however, more recent articles have argued that these models can be over-specified
# (see: https://arxiv.org/pdf/1506.04967.pdf)
# to evaluate model significance, actually performing a likelihood ratio test even though we are calling `anova`
# see: http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf (page 13)
# R tricks for speeding up model fitting with `lme4` and `purrr`
# trick to speed up based on: https://cran.r-project.org/web/packages/lme4/vignettes/lmerperf.html
# see `fitting the models` from: http://jofrhwld.github.io/blog/2016/05/01/many_models.html
# an approach for testing models from: https://stackoverflow.com/questions/50702152/compare-models-via-anova-with-purrr-or-dplyr
# notes on model convergence (hard to know if it's actually an issue or overly sensitive error reporting): https://rpubs.com/bbolker/lme4_convergence

#' # Notes about permutation tests
#' ## Toy example to check the code
test_permute <- dplyr::tibble(subj = rep(c("s1", "s2", "s3"), times = 5),
                              hemi = sample(c("left", "right"), size = length(subj), replace = TRUE),
                              source_mem = sample(c("rHit_questionLiberalMiss", "rHit_questionLiberalHit"),
                                                  size = length(subj), replace = TRUE),
                              z_r = rnorm(length(subj))) %>%
  dplyr::group_by(subj, hemi) %>%
  tidyr::nest(.key = "trials") %>%
  # permute the trials for each subject where what gets permuted is the column you are trying to test in the model
  dplyr::mutate(source_mem_perm = purrr::map(trials, ~sample(x = .$source_mem, size = length(.$source_mem), replace = FALSE)))
# grab 2 rows of data and check that they are in different orders
test_permute$trials[[1]]
test_permute$source_mem_perm[[1]]

#' ## Correction factor for permuted p-values
# numerator = (number of permutated statistics (e.g., chi squared value) that exceed the observed chi squared value) + 1
# denominator = (number of permutations) + 1
# this correction factor (+1) essentially assumes that if you ran ALL possible permutations
# that your observed data would have been one of these and so adds it to the total number
# of permutation iterations. NB: w/ small numbers of permutations, this correction factor
# mucks up the estimation of the p-values

#' ## Computed p-values
# based on idea that actual pvalues are superior (esp w/ lots of permutation tests)
# see: https://www.ncbi.nlm.nih.gov/pubmed/21044043, esp p. 6
# how many of the permutations exceed the observed value (w/ +1 added as a correction factor)
# [(# of permuted values > observed value) + 1] / [number of distinct permutations computed + 1]

#' ## Permutation p-values for interactions
# since what we care about is an interaction, need to permute the interaction column
# (or, this is conceptually the same as using the same permutation ordering for both columns of interest)
# seems like the best way to do this is to `tidyr::unite()` the two columns of interest,
# permute this united column, `tidyr::spread()`, and then run stats. this ensures that
# the relationship between the two columns remains the same, but their relationship to the
# dependent variable gets permuted. to check this, can plot a count of the observations in the original
# and permuted data

#' # Split by encoding list half (same/diff) *and* same/diff question
tictoc::tic(msg = "Split by encoding list half (same/diff) and same/diff question")

#' ## Graph
# repeat 2x so that can split plots by HC vs MTLc (since we'd expect pretty different scales)
for (iroi in 1:length(rois_of_interest)){
  cur_rois <- rois_of_interest[[iroi]]

  # figure out how to label the graph when it gets saved out
  if (any(is.element(hc_rois, cur_rois))) {
    roi_lbl <- "hc-rois"
  } else if (any(is.element(mtlc_rois, cur_rois))) {
    roi_lbl <- "mtlc-rois"
  }

  z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi %in% cur_rois) %>%
    dplyr::filter(!(is.na(roi_formatted))) %>%
    dplyr::group_by(subj, hemi, question, list_same_diff_half, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE)) %>%
    boxplot_facet_hemi_and_somethingelse(., "roi_formatted", "mean_r", "list_same_diff_half", "hemi", "question")
  print(p)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("list-same-diff-halves_quest-same-diff_corr_%s.pdf", roi_lbl)),
                    height = 6, width = 12)
  }

  z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi %in% cur_rois) %>%
    dplyr::filter(!(is.na(roi_formatted))) %>%
    tidyr::unite(quest_list, question, list_same_diff_half) %>%
    dplyr::group_by(subj, hemi, quest_list, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE)) %>%
    # revalue `quest_list` for plot formatting
    dplyr::mutate(quest_list_formatted = dplyr::recode(quest_list, "same_same_half" = "SQSH", "same_diff_half" = "SQDH",
                                                       "diff_same_half" = "DQSH", "diff_diff_half" = "DQDH"),
                  quest_list_formatted = factor(quest_list_formatted, levels = c("SQSH", "SQDH", "DQSH", "DQDH"))) %>%
    boxplot_facet_single_var(., "roi_formatted", "mean_r", "quest_list_formatted", "hemi")
  print(p)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("list-same-diff-halves_quest-same-diff_4bars_corr_%s.pdf", roi_lbl)),
                    height = 6, width = 12)
  }

  z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi %in% cur_rois) %>%
    dplyr::filter(!(is.na(roi_formatted))) %>%
    tidyr::unite(quest_list, question, list_same_diff_half) %>%
    dplyr::group_by(subj, quest_list, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE)) %>%
    # revalue `quest_list` for plot formatting
    dplyr::mutate(quest_list_formatted = dplyr::recode(quest_list, "same_same_half" = "SQSH", "same_diff_half" = "SQDH",
                                                       "diff_same_half" = "DQSH", "diff_diff_half" = "DQDH"),
                  quest_list_formatted = factor(quest_list_formatted, levels = c("SQSH", "SQDH", "DQSH", "DQDH"))) %>%
    boxplot_no_facet(., "roi_formatted", "mean_r", "quest_list_formatted")
  print(p)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("list-same-diff-halves_quest-same-diff_4bars_corr_%s_collapse-hemi.pdf", roi_lbl)),
                    height = 6, width = 12)
  }

  z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi %in% cur_rois) %>%
    dplyr::filter(!(is.na(roi_formatted))) %>%
    dplyr::group_by(subj, hemi, list_same_diff_half, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE)) %>%
    boxplot_facet_single_var(., "roi_formatted", "mean_r", "list_same_diff_half", "hemi")
  print(p)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("list-same-diff-halves_quest-same-diff_ME-list_corr_%s.pdf", roi_lbl)),
                    height = 6, width = 12)
  }

  z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi %in% cur_rois) %>%
    dplyr::filter(!(is.na(roi_formatted))) %>%
    dplyr::group_by(subj, list_same_diff_half, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE)) %>%
    boxplot_no_facet(., "roi_formatted", "mean_r", "list_same_diff_half")
  print(p)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("list-same-diff-halves_quest-same-diff_ME-list_corr_%s_collapse-hemi.pdf", roi_lbl)),
                    height = 6, width = 12)
  }

  z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi %in% cur_rois) %>%
    dplyr::filter(!(is.na(roi_formatted))) %>%
    dplyr::group_by(subj, hemi, question, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE)) %>%
    boxplot_facet_single_var(., "roi_formatted", "mean_r", "question", "hemi")
  print(p)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("list-same-diff-halves_quest-same-diff_ME-quest_corr_%s.pdf", roi_lbl)),
                    height = 6, width = 12)
  }

  z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  list_same_diff_half %in% c("same_half", "diff_half"),
                  question %in% c("same", "diff")) %>%
    dplyr::filter(roi %in% cur_rois) %>%
    dplyr::filter(!(is.na(roi_formatted))) %>%
    dplyr::group_by(subj, question, roi_formatted) %>%
    dplyr::summarize(mean_r = mean(r, na.rm = TRUE)) %>%
    boxplot_no_facet(., "roi_formatted", "mean_r", "question")
  print(p)

  if(SAVE_GRAPH_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("list-same-diff-halves_quest-same-diff_ME-quest_corr_%s_collapse-hemi.pdf", roi_lbl)),
                    height = 6, width = 12)
  }
}

#' ## Stats
list_sdhalves_quest_samediff_stats <- z_mm_filt %>%
  dplyr::filter(list_lag >= 0,
                list_same_diff_half %in% c("same_half", "diff_half"),
                question %in% c("same", "diff")) %>%
  dplyr::group_by(roi) %>%
  tidyr::nest(.key = "by_roi") %>%
  dplyr::mutate(mod_null = purrr::map(by_roi, ~lme4::lmer(z_r ~ 1 + (1 | subj), data = .,
                                                          REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_hemi = purrr::map(by_roi, ~lme4::lmer(z_r ~ hemi + (1 | subj), data = .,
                                                          REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_list_halves = purrr::map(by_roi, ~lme4::lmer(z_r ~ hemi + list_same_diff_half + (1 | subj), data = .,
                                                                 REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_quest = purrr::map(by_roi, ~lme4::lmer(z_r ~ hemi + question + (1 | subj), data = .,
                                                                 REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_list_halves_quest = purrr::map(by_roi, ~lme4::lmer(z_r ~ hemi + list_same_diff_half + question + (1 | subj), data = .,
                                                                       REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_list_halvesXquest = purrr::map(by_roi, ~lme4::lmer(z_r ~ hemi + (list_same_diff_half * question) + (1 | subj), data = .,
                                                                       REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_list_halvesXquestXhemi = purrr::map(by_roi, ~lme4::lmer(z_r ~ (list_same_diff_half * question * hemi) + (1 | subj), data = .,
                                                                       REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                aov_null.hemi = purrr::map2(mod_null, mod_hemi, anova),
                chisq_null.hemi = purrr::map_dbl(aov_null.hemi, ~.$Chisq[2]),
                aov_ME.listhalves = purrr::map2(mod_quest, mod_list_halves_quest, anova),
                chisq_ME.listhalves = purrr::map_dbl(aov_ME.listhalves, ~.$Chisq[2]),
                aov_ME.quest = purrr::map2(mod_list_halves, mod_list_halves_quest, anova),
                chisq_ME.quest = purrr::map_dbl(aov_ME.quest, ~.$Chisq[2]),
                aov_listhalvesquest.listhalvesXquest = purrr::map2(mod_list_halves_quest, mod_list_halvesXquest, anova),
                chisq_listhalvesquest.listhalvesXquest = purrr::map_dbl(aov_listhalvesquest.listhalvesXquest, ~.$Chisq[2]),
                aov_listhalvesXquest.listhalvesXquestXhemi = purrr::map2(mod_list_halvesXquest, mod_list_halvesXquestXhemi, anova),
                chisq_listhalvesXquest.listhalvesXquestXhemi = purrr::map_dbl(aov_listhalvesXquest.listhalvesXquestXhemi, ~.$Chisq[2]))

# save out for use in permutations
save(list_sdhalves_quest_samediff_stats, file = file.path(analyzed_mri_dir, "list_sdhalves_quest_samediff_stats.RData"))

#' ### Print model summaries
#+ warning = FALSE
mod_summarize(list_sdhalves_quest_samediff_stats)
out_df %>%
  dplyr::filter(term == ".y[[i]]") %>%
  dplyr::select(roi, aov_type, statistic, Chi.Df, p.value) %>%
  knitr::kable() %>%
  kableExtra::kable_styling(bootstrap_options = "striped")

#' ### Trial numbers
list_sdhalves_quest_samediff_stats %>%
  dplyr::select(roi, by_roi) %>%
  tidyr::unnest() %>%
  plot_trial_nums_interaction(., "list_same_diff_half", "question", "subj", "hemi", "roi")
p_summary
p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
print(p)

if(SAVE_GRAPH_FLAG == 1){
  ggplot2::ggsave(filename = file.path(graph_out_dir, "list-same-diff-halves_quest-same-diff_trialnums.pdf"),
                  height = 6, width = 12)
}

#' ## Stats: Split by hemi
list_sdhalves_quest_samediff_split_by_hemi_stats <- z_mm_filt %>%
  dplyr::filter(list_lag >= 0,
                list_same_diff_half %in% c("same_half", "diff_half"),
                question %in% c("same", "diff")) %>%
  dplyr::group_by(roi, hemi) %>%
  tidyr::nest(.key = "by_roi") %>%
  dplyr::mutate(mod_list_halves = purrr::map(by_roi, ~lme4::lmer(z_r ~ list_same_diff_half + (1 | subj), data = .,
                                                                 REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_quest = purrr::map(by_roi, ~lme4::lmer(z_r ~ question + (1 | subj), data = .,
                                                           REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_list_halves_quest = purrr::map(by_roi, ~lme4::lmer(z_r ~ list_same_diff_half + question + (1 | subj), data = .,
                                                                       REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                mod_list_halvesXquest = purrr::map(by_roi, ~lme4::lmer(z_r ~ (list_same_diff_half * question) + (1 | subj), data = .,
                                                                       REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                aov_ME.listhalves = purrr::map2(mod_quest, mod_list_halves_quest, anova),
                chisq_ME.listhalves = purrr::map_dbl(aov_ME.listhalves, ~.$Chisq[2]),
                aov_ME.quest = purrr::map2(mod_list_halves, mod_list_halves_quest, anova),
                chisq_ME.quest = purrr::map_dbl(aov_ME.quest, ~.$Chisq[2]),
                aov_listhalvesquest.listhalvesXquest = purrr::map2(mod_list_halves_quest, mod_list_halvesXquest, anova),
                chisq_listhalvesquest.listhalvesXquest = purrr::map_dbl(aov_listhalvesquest.listhalvesXquest, ~.$Chisq[2]))

# save out for use in permutations
save(list_sdhalves_quest_samediff_split_by_hemi_stats, file = file.path(analyzed_mri_dir, "list_sdhalves_quest_samediff_split_by_hemi_stats.RData"))

#' ### Print model summaries
#+ warning = FALSE
mod_summarize_by_hemi_and_roi(list_sdhalves_quest_samediff_split_by_hemi_stats)

#' #### Left
out_df %>%
  dplyr::filter(term == ".y[[i]]") %>%
  dplyr::select(roi, hemi, aov_type, statistic, Chi.Df, p.value) %>%
  dplyr::filter(hemi == "left") %>%
  knitr::kable() %>%
  kableExtra::kable_styling(bootstrap_options = "striped")

#' #### Right
out_df %>%
  dplyr::filter(term == ".y[[i]]") %>%
  dplyr::select(roi, hemi, aov_type, statistic, Chi.Df, p.value) %>%
  dplyr::filter(hemi == "right") %>%
  knitr::kable() %>%
  kableExtra::kable_styling(bootstrap_options = "striped")

#' ## Clean up
remove(list = c("list_sdhalves_quest_samediff_stats", "list_sdhalves_quest_samediff_split_by_hemi_stats"))

#' ## Timer
tictoc::toc()

#' # Print out R session info
sessionInfo()
sprintf("It took %.2f minutes to run all of this.", Sys.time() - time_script_start)
tictoc::toc()
