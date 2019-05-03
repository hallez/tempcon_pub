#' ---
#' title: Permutation test TempCon data
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
numperm <- 1000
alpha_value <- 0.05
# since permutations take ages to run, be selective about only including rois of interest
rois_of_interest <- c("CA1_body", "CA2_3_DG_body", "subiculum_body", "PRC", "aPHC", "pPHC")

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
  dplyr::filter(item_mem == "rHit") %>%
  dplyr::filter(roi %in% rois_of_interest)

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
message(sprintf("Analyzing ROIs of interest: "))
message(sprintf(paste0(rois_of_interest, collapse = ', ')))

#' # Split by encoding list half (same/diff) *and* same/diff question
tictoc::tic(msg = "Split by encoding list half (same/diff) and same/diff question")

#' ## Read in mixed models output
load(file = file.path(analyzed_mri_dir, "list_sdhalves_quest_samediff_split_by_hemi_stats.RData"))
load(file = file.path(analyzed_mri_dir, "list_sdhalves_quest_samediff_stats.RData"))

# filter so only have ROIs of interest
list_sdhalves_quest_samediff_split_by_hemi_stats <- list_sdhalves_quest_samediff_split_by_hemi_stats %>%
  dplyr::filter(roi %in% rois_of_interest)

head(list_sdhalves_quest_samediff_split_by_hemi_stats)

list_sdhalves_quest_samediff_stats <- list_sdhalves_quest_samediff_stats %>%
  dplyr::filter(roi %in% rois_of_interest)

head(list_sdhalves_quest_samediff_stats)

#' # Stats: Hemi permutation
qsd_halves_hemi_permute_anova <- data.frame(stringsAsFactors = FALSE)
for (iperm in 1:numperm){
  qsd_halves_hemi_permute <- z_mm_filt %>%
    dplyr::filter(list_lag >= 0,
                  question %in% c("same", "diff"),
                  list_same_diff_half %in% c("same_half", "diff_half")) %>%
    dplyr::group_by(subj, hemi, roi) %>%
    # since we're testing a 3-way interaction, permute all 3 relevant columns (but maintain relationships w/ one another)
    tidyr::unite(listhalves_qsamediff_hemi, list_same_diff_half, question, hemi, sep = ".") %>%
    tidyr::nest(.key = "trials") %>%
    dplyr::mutate(listhalves_qsamediff_hemi_perm = purrr::map(trials, ~sample(x = .$listhalves_qsamediff_hemi, size = length(.$listhalves_qsamediff_hemi),
                                                                         replace = FALSE))) %>%
    tidyr::unnest() %>%
    # use the `\\` to escape the period so it's read as a literal period (rather than regex)
    tidyr::separate(listhalves_qsamediff_hemi_perm, into = c("list_halves_perm", "question_perm", "hemi_perm"), sep = "\\.") %>%
    dplyr::group_by(roi) %>%
    tidyr::nest(.key = "by_roi") %>%
    dplyr::mutate(mod_list_halvesXquest_perm = purrr::map(by_roi, ~lme4::lmer(z_r ~ hemi_perm + (list_halves_perm * question_perm) + (1 | subj), data = .,
                                                                         REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                  mod_list_halvesXquestXhemi_perm = purrr::map(by_roi, ~lme4::lmer(z_r ~ (list_halves_perm * question_perm * hemi_perm) + (1 | subj), data = .,
                                                                              REML = FALSE, lme4::lmerControl(calc.derivs = FALSE))),
                  aov_listhalvesXquest.listhalvesXquestXhemi_perm = purrr::map2(mod_list_halvesXquest_perm, mod_list_halvesXquestXhemi_perm, anova),
                  chisq_listhalvesXquest.listhalvesXquestXhemi_perm = purrr::map_dbl(aov_listhalvesXquest.listhalvesXquestXhemi_perm, ~.$Chisq[2]),
                  nperm = as.numeric(iperm))

  print(message("permutation", as.numeric(iperm)))

  if(dim(qsd_halves_hemi_permute_anova)[1] == 0){
    qsd_halves_hemi_permute_anova <- qsd_halves_hemi_permute %>%
      dplyr::select(nperm, roi, starts_with("chisq"))
  } else {
    cur_dat <- qsd_halves_hemi_permute %>%
      dplyr::select(nperm, roi, starts_with("chisq"))

    qsd_halves_hemi_permute_anova <- dplyr::bind_rows(qsd_halves_hemi_permute_anova, cur_dat)
  }

  # write out progress to a file so that know script is running
  progress_fname <- file.path(analyzed_mri_dir,
                              "perms_list-samediff-halves_quest-samediff-hemi_interactions_progress.csv")
  write.csv(data.frame(perm_iteration = unique(qsd_halves_hemi_permute_anova$nperm)), file = progress_fname, append = TRUE, quote = FALSE, row.names = FALSE)

}

#' ## Compute permuted p-values
qsd_halves_hemi_permuted_pvals <- qsd_halves_hemi_permute_anova %>%
  dplyr::group_by(roi) %>%
  dplyr::mutate_at(vars(starts_with("chisq")),
                   funs(p95 = quantile(x = ., probs = 0.95),
                        p99 = quantile(x = ., probs = 0.99)))

#' ## Compute actual p values
qsd_halves_hemi_computed_pvals <- list_sdhalves_quest_samediff_stats %>%
  dplyr::left_join(., qsd_halves_hemi_permute_anova, by = c("roi")) %>%
  dplyr::mutate(listhalvesXsamediff.listhalvesXsamediffXhemi_perm_exceeds = case_when(chisq_listhalvesXquest.listhalvesXquestXhemi_perm > chisq_listhalvesXquest.listhalvesXquestXhemi ~ 1,
                                                                                      chisq_listhalvesXquest.listhalvesXquestXhemi_perm <= chisq_listhalvesXquest.listhalvesXquestXhemi ~ 0)) %>%
  dplyr::group_by(roi) %>%
  dplyr::summarise(chisq_listhalvesXquest.listhalvesXquestXhemi_num_signif_perms = sum(listhalvesXsamediff.listhalvesXsamediffXhemi_perm_exceeds)) %>%
  dplyr::mutate(chisq_listhalvesXquest.listhalvesXquestXhemi_computed_pval = (chisq_listhalvesXquest.listhalvesXquestXhemi_num_signif_perms + 1) / (numperm + 1))

#' ## Pretty print stat results
same_diff_list_halves_hemi_all_results <- list_sdhalves_quest_samediff_stats %>%
  dplyr::left_join(., qsd_halves_hemi_permute_anova, by = c("roi")) %>%
  dplyr::left_join(., qsd_halves_hemi_permuted_pvals, by = c("roi", "nperm", "chisq_listhalvesXquest.listhalvesXquestXhemi_perm")) %>%
  dplyr::left_join(., qsd_halves_hemi_computed_pvals, by = c("roi"))

same_diff_list_halves_hemi_all_results %>%
  dplyr::select(roi, contains("computed_pval"),
                "chisq_listhalvesXquest.listhalvesXquestXhemi") %>%
  dplyr::distinct() %>%
  dplyr::select(roi,
                chisq_listhalvesXquest.listhalvesXquestXhemi, chisq_listhalvesXquest.listhalvesXquestXhemi_computed_pval) %>%
  knitr::kable() %>%
  kableExtra::kable_styling(bootstrap_options = "striped")

#' ## Timer
tictoc::toc()
