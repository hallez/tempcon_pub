#' ---
#' title: Tidy between runs mixed models data
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
graph_out_dir <- paste0(halle::ensure_trailing_slash(dropbox_dir),
                        halle::ensure_trailing_slash("plots"))
dir.create(graph_out_dir, recursive=TRUE)

#+ label = "subjects"
#' ## Figure out subjects
# assume that subjects have folders in `raw_behavioral_dir` that start with `s` and are followed by one or two digits
subjects <- c(list.files(path=raw_behavioral_dir,pattern="^[s][(123456789)]$|^[s][(123456789)][(0123456789)]$"))
# chuck out any subjects who need to be excluded
exclude_subjects <- c("s1", "s2", "s3", "s7", "s11", "s13", "s14", "s24", "s26")
subjects <- subjects[!is.element(subjects, exclude_subjects)]
length(subjects)
subjects

#+ label = "setup vars"
#' ## Set other variables
roi_dirs <- c("left","right")
rois <- c("CA1","CA2","CA3","DG","ERC","subiculum",
          "35_36","CA2_3_DG","CA3_DG", "whole_hippo",
          "CA1_head", "CA1_body", "CA1_tail",
          "CA2_3_DG_head", "CA2_3_DG_body", "CA2_3_DG_tail",
          "whole_hippo_head", "whole_hippo_body", "whole_hippo_tail")

rois_of_most_interest <- c("subiculum_body",
                           "CA1_body", "CA2_3_DG_body",
                           "PRC_L", "PRC_R",
                           "aPHC_L", "aPHC_R", "pPHC_L", "pPHC_R")

#' ## Flags
SAVE_FLAG <- 1
SAVE_GRAPH_FLAG <- 1
EXCLUDE_SUBJ_FLAG <- 0
REDUCED_ROIS_FLAG <- 1
LOAD_EXISTING_GROUP_FILE_FLAG <- 0

if(REDUCED_ROIS_FLAG==1){
  rois <- rois_of_most_interest
}

#' # Load group trial pairs data
#+ label = "load data", message = TRUE
if(LOAD_EXISTING_GROUP_FILE_FLAG == 1){
  message(sprintf("\nLoad in existing group data file.\n"))
  # this file would have been created previously by `mixed-models-btwn-runs.R`
  load(file.path(analyzed_mri_dir, "z_mm_trial_pairs_all_subj.RData"))
} else {
  # these files are created by `rsa-load-data-btwn-runs.R`
  all_subj_mm_trial_pairs <- data.frame()
  for(isubj in subjects){
    message(sprintf("\nLoading data for %s\n", isubj))
    for(ihemi in 1:length(roi_dirs)) {
      for(iroi in 1:length(rois)) {
        hemi_label <- roi_dirs[ihemi]
        cur_roi <- rois[iroi]

        # the file is created in `rsa-load-data-btwn-runs.R`
        mm_filepath <- paste0(analyzed_mri_dir,
                              halle::ensure_trailing_slash(isubj),
                              sprintf("mm_trial_pairs_%s_%s.RData", hemi_label, cur_roi))
        if(!file.exists(mm_filepath)) {
          message(sprintf("\nMixed models file %s does not exist - skipping\n", mm_filepath))
          next
        } else {
          load(mm_filepath) # this loads a variable called `mm_trial_pairs`
        }

        if(dim(all_subj_mm_trial_pairs)[1]==0){
          all_subj_mm_trial_pairs <- mm_trial_pairs
        } else {
          all_subj_mm_trial_pairs <- dplyr::full_join(all_subj_mm_trial_pairs, mm_trial_pairs,
                                                      by=intersect(names(all_subj_mm_trial_pairs), names(mm_trial_pairs)))
          mm_trial_pairs <- NULL
        }
      } # iroi
    } # ihemi
    message(Sys.time() - time_script_start)
    message("save out group data file")
    save(all_subj_mm_trial_pairs, file = file.path(analyzed_mri_dir, "mm_trial_pairs_all_subj.RData"))
  } # for(isubj

  #+ label = "drop subjects"
  if(EXCLUDE_SUBJ_FLAG == 1){
    all_subj_mm_trial_pairs <-
      all_subj_mm_trial_pairs %>%
      dplyr::filter(!subj %in% c("s3", "s15", "s32"))
  }

  #+ label = "z-score data"
  #' # Fisher transform (for stats) and tidy up dataframe
  z_mm <-
    all_subj_mm_trial_pairs %>%
    dplyr::mutate(z_r = halle::fisher_r_to_z(r)) %>%
    dplyr::mutate(subj_fact = as.factor(subj),
                  roi = dplyr::recode(roi, "PRC_L" = "PRC", "PRC_R" = "PRC",
                                      "pPHC_L" = "pPHC", "pPHC_R" = "pPHC",
                                      "aPHC_L" = "aPHC", "aPHC_R" = "aPHC",
                                      "pmEC_L" = "pmEC", "pmEC_R" = "pmEC",
                                      "alEC_L" = "alEC", "alEC_R" = "alEC"),
                  roi_fact = as.factor(roi),
                  # NB: if any level is missing, it will get filled with NA
                  roi_fact_ordered = factor(roi_fact, levels = c("CA1_body", "CA2_3_DG_body", "ERC", "subiculum_body",
                                                                 "alEC", "pmEC", "PRC", "aPHC", "pPHC",
                                                                 "whole_hippo", "whole_hippo_head", "whole_hippo_body", "whole_hippo_tail")),
                  # in order to preserve the ordering, have to again use `levels` (i.e., it doesn't persist from `roi_fact_ordered`)
                  roi_formatted = as.factor(car::recode(roi_fact_ordered, c("'CA1_body' = 'CA1';'CA2_3_DG_body' = 'CA23DG'; 'subiculum_body' = 'SUB'"),
                                                        levels = c("CA1", "CA23DG", "SUB",
                                                                   "alEC", "pmEC", "PRC", "aPHC", "pPHC",
                                                                   "whole_hippo", "whole_hippo_head", "whole_hippo_body", "whole_hippo_tail", "ERC"))),
                  hemi_fact = as.factor(hemi),
                  item_mem_fact = as.factor(item_mem),
                  question_mem_fact = as.factor(question_mem),
                  list_mem_fact = as.factor(list_mem),
                  question_fact = factor(question, levels = c("same", "diff")),
                  question_formatted = as.factor(car::recode(question, c("'same' = 'same_task';'diff' = 'diff_task'"))),
                  question_formatted = factor(question_formatted, levels = c("same_task", "diff_task")),
                  list_fact = as.factor(list),
                  source_mem_fact = as.factor(source_mem),
                  condition_fact = as.factor(condition),
                  cond_short_fact = as.factor(car::recode(condition_fact, c("'sameList_sameQuestion' = 'SLSQ'; 'diffList_sameQuestion' = 'DLSQ';'sameList_diffQuestion' = 'SLDQ';'diffList_diffQuestion' = 'DLDQ'"))),
                  # NB: this is only partially correct - it only considers trial pairs from the same list
                  list_halves_fact = as.factor(ifelse(grepl("list[1-4]", list), "first_half",
                                                      ifelse(grepl("list[5-8]", list), "second_half",list))),
                  retrieval_lag = col_recog_trial_id - row_recog_trial_id) %>%
    # figure out encoding questions for trial pair
    # `question_id` only gives information when both trials were associated with the same question
    # mapping between question number and question type comes from `scripts/montchal/make_stimuli.m`
    dplyr::mutate(col_enctask_lbl = dplyr::recode(col_enctask, `1` = "fridge", `2` = "bathtub", `3` = "supermarket", `4` = "conveniencestore")) %>%
    dplyr::mutate(row_enctask_lbl = dplyr::recode(row_enctask, `1` = "fridge", `2` = "bathtub", `3` = "supermarket", `4` = "conveniencestore")) %>%
    tidyr::unite(c("col_enctask_lbl", "row_enctask_lbl"), col = "rowenc_colenc", sep = "_") %>%
    dplyr::mutate(encoding_quest = factor(x = rowenc_colenc, levels = c("bathtub_bathtub", "bathtub_fridge", "bathtub_conveniencestore", "bathtub_supermarket",
                                                                       "conveniencestore_conveniencestore", "conveniencestore_supermarket", "conveniencestore_bathtub", "conveniencestore_fridge",
                                                                       "fridge_fridge", "fridge_bathtub", "fridge_conveniencestore", "fridge_supermarket",
                                                                       "supermarket_supermarket", "supermarket_conveniencestore", "supermarket_bathtub", "supermarket_fridge"))) %>%
    dplyr::mutate(quest_fit_find = dplyr::recode(encoding_quest, "bathtub_bathtub" = "fit",
                                                 "conveniencestore_conveniencestore" = "find",
                                                 "fridge_fridge" = "fit",
                                                 "supermarket_supermarket" = "find",
                                                 "bathtub_fridge" = "fit",
                                                 "fridge_bathtub" = "fit",
                                                 "conveniencestore_supermarket" = "find",
                                                 "supermarket_conveniencestore" = "find")) %>%
    dplyr::group_by(list_fact) %>%
    dplyr::mutate(encoding_lag = col_enc_trial_id - row_enc_trial_id) %>%
    dplyr::ungroup() %>%
    # re-value same-list labels
    dplyr::mutate(list_fact_recode = dplyr::recode_factor(list_fact, list1 = "list1_list1", list2 = "list2_list2", list3 = "list3_list3", list4 = "list4_list4",
                                                          list5 = "list5_list5", list6 = "list6_list6", list7 = "list7_list7", list8 = "list8_list8")) %>%
    # compute list lag
    dplyr::mutate(list_fact_recode2 = list_fact_recode) %>%
    tidyr::separate(list_fact_recode2, into = c("row_list", "col_list")) %>%
    dplyr::mutate(row_list_num = as.numeric(sub("list", "", row_list)),
                  col_list_num = as.numeric(sub("list", "", col_list)),
                  list_lag = col_list_num - row_list_num) %>%
    # create a column that splits list by same/different half - NB this is different from early/late
    dplyr::mutate(list_same_diff_half = dplyr::recode_factor(list_fact_recode,
                                                             "list1_list1" = "same_half", "list1_list2" = "same_half", "list1_list3" = "same_half", "list1_list4" = "same_half",
                                                             "list2_list1" = "same_half", "list2_list2" = "same_half", "list2_list3" = "same_half", "list2_list4" = "same_half",
                                                             "list3_list1" = "same_half", "list3_list2" = "same_half", "list3_list3" = "same_half", "list3_list4" = "same_half",
                                                             "list4_list1" = "same_half", "list4_list2" = "same_half", "list4_list3" = "same_half", "list4_list4" = "same_half",
                                                             "list5_list5" = "same_half", "list5_list6" = "same_half", "list5_list7" = "same_half", "list5_list8" = "same_half",
                                                             "list6_list5" = "same_half", "list6_list6" = "same_half", "list6_list7" = "same_half", "list6_list8" = "same_half",
                                                             "list7_list5" = "same_half", "list7_list6" = "same_half", "list7_list7" = "same_half", "list7_list8" = "same_half",
                                                             "list8_list5" = "same_half", "list8_list6" = "same_half", "list8_list7" = "same_half", "list8_list8" = "same_half",
                                                             "list1_list5" = "diff_half", "list1_list6" = "diff_half", "list1_list7" = "diff_half", "list1_list8" = "diff_half",
                                                             "list2_list5" = "diff_half", "list2_list6" = "diff_half", "list2_list7" = "diff_half", "list2_list8" = "diff_half",
                                                             "list3_list5" = "diff_half", "list3_list6" = "diff_half", "list3_list7" = "diff_half", "list3_list8" = "diff_half",
                                                             "list4_list5" = "diff_half", "list4_list6" = "diff_half", "list4_list7" = "diff_half", "list4_list8" = "diff_half",
                                                             "list5_list1" = "diff_half", "list5_list2" = "diff_half", "list5_list3" = "diff_half", "list5_list4" = "diff_half",
                                                             "list6_list1" = "diff_half", "list6_list2" = "diff_half", "list6_list3" = "diff_half", "list6_list4" = "diff_half",
                                                             "list7_list1" = "diff_half", "list7_list2" = "diff_half", "list7_list3" = "diff_half", "list7_list4" = "diff_half",
                                                             "list8_list1" = "diff_half", "list8_list2" = "diff_half", "list8_list3" = "diff_half", "list8_list4" = "diff_half")) %>%
    # figure out if the list was in the first or second half, essentially split `same_half` from `list_same_diff_half`
    dplyr::mutate(list_first_second_half = dplyr::recode_factor(list_fact_recode,
                                                                "list1_list1" = "first_half", "list1_list2" = "first_half", "list1_list3" = "first_half", "list1_list4" = "first_half",
                                                                "list2_list1" = "first_half", "list2_list2" = "first_half", "list2_list3" = "first_half", "list2_list4" = "first_half",
                                                                "list3_list1" = "first_half", "list3_list2" = "first_half", "list3_list3" = "first_half", "list3_list4" = "first_half",
                                                                "list4_list1" = "first_half", "list4_list2" = "first_half", "list4_list3" = "first_half", "list4_list4" = "first_half",
                                                                "list5_list5" = "second_half", "list5_list6" = "second_half", "list5_list7" = "second_half", "list5_list8" = "second_half",
                                                                "list6_list5" = "second_half", "list6_list6" = "second_half", "list6_list7" = "second_half", "list6_list8" = "second_half",
                                                                "list7_list5" = "second_half", "list7_list6" = "second_half", "list7_list7" = "second_half", "list7_list8" = "second_half",
                                                                "list8_list5" = "second_half", "list8_list6" = "second_half", "list8_list7" = "second_half", "list8_list8" = "second_half",
                                                                "list1_list5" = "diff_half", "list1_list6" = "diff_half", "list1_list7" = "diff_half", "list1_list8" = "diff_half",
                                                                "list2_list5" = "diff_half", "list2_list6" = "diff_half", "list2_list7" = "diff_half", "list2_list8" = "diff_half",
                                                                "list3_list5" = "diff_half", "list3_list6" = "diff_half", "list3_list7" = "diff_half", "list3_list8" = "diff_half",
                                                                "list4_list5" = "diff_half", "list4_list6" = "diff_half", "list4_list7" = "diff_half", "list4_list8" = "diff_half",
                                                                "list5_list1" = "diff_half", "list5_list2" = "diff_half", "list5_list3" = "diff_half", "list5_list4" = "diff_half",
                                                                "list6_list1" = "diff_half", "list6_list2" = "diff_half", "list6_list3" = "diff_half", "list6_list4" = "diff_half",
                                                                "list7_list1" = "diff_half", "list7_list2" = "diff_half", "list7_list3" = "diff_half", "list7_list4" = "diff_half",
                                                                "list8_list1" = "diff_half", "list8_list2" = "diff_half", "list8_list3" = "diff_half", "list8_list4" = "diff_half")) %>%
    # split up `quest_fit_find` into same/different question stem - essentially, combine fit and find and label all other pairings as different
    dplyr::mutate(quest_same_diff_stem = dplyr::recode(quest_fit_find, "fit" = "same_quest_stem", "find" = "same_quest_stem",
                                                       "conveniencestore_bathtub" = "diff_quest_stem", "conveniencestore_fridge" = "diff_quest_stem",
                                                       "fridge_supermarket" = "diff_quest_stem", "fridge_conveniencestore" = "diff_quest_stem",
                                                       "supermarket_fridge" = "diff_quest_stem", "supermarket_bathtub" = "diff_quest_stem",
                                                       "bathtub_supermarket" = "diff_quest_stem", "bathtub_conveniencestore" = "diff_quest_stem"))
    print("column names in z_mm are: ")
    print(colnames(z_mm))

  #' # Write out group files for use elsewhere
  save(z_mm, file = file.path(analyzed_mri_dir, "z_mm_trial_pairs_all_subj.RData"))
} # if(LOAD_EXISTING_GROUP_FILE_FLAG

message(sprintf("\nFinished loading. Time so far: \n"))
message(Sys.time() - time_script_start)

#' ## List included subjects, conditions, ROIs
unique(z_mm$subj)
length(unique(z_mm$subj))
unique(z_mm$roi)
unique(z_mm$condition)

#+ label = "trial numbers"
#' # Trial numbers by subject
#' ## Set a minimum number of trials
min_trials <- 10
message(sprintf("Current minimum number of trials is %d", min_trials))

#' ## Graph
source_mem_filter <- "rHit"
z_mm %>%
  # to look at trial numbers, should just need a single ROI and hemi
  dplyr::filter(roi == "CA2_3_DG_body") %>%
  dplyr::filter(hemi == "right") %>%
  dplyr::filter(source_mem %in% source_mem_filter) %>%
  dplyr::filter(cond_short_fact %in% c("SLSQ", "DLSQ", "SLDQ", "DLDQ")) %>%
  dplyr::group_by(cond_short_fact, subj_fact) %>%
  dplyr::tally() %>%
  ggplot2::ggplot(ggplot2::aes(x = cond_short_fact, y = n)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_wrap(~subj_fact) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = min_trials, color = "red")) +
  ggplot2::ylab("count of trial pairs") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle=90, hjust=1), axis.title.x = ggplot2::element_blank(),
                 legend.position = "none")

if(SAVE_GRAPH_FLAG == 1){
  ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("trial-nums-%s-same-diff-list-by-question-%d-min-trials.pdf", source_mem_filter, min_trials)),
                  height = 10, width = 7)
}

z_mm %>%
  # to look at trial numbers, should just need a single ROI and hemi
  dplyr::filter(roi == "CA2_3_DG_body") %>%
  dplyr::filter(hemi == "right") %>%
  dplyr::filter(source_mem %in% source_mem_filter) %>%
  dplyr::filter(grepl("list[1-8]", list_fact)) %>%
  dplyr::group_by(list_fact, subj_fact) %>%
  dplyr::tally() %>%
  ggplot2::ggplot(ggplot2::aes(x = list_fact, y = n)) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(ggplot2::aes(color = subj_fact)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = min_trials)) +
  ggplot2::annotate("text", x = 10, y = min_trials, label = sprintf("Cutoff of %d trials", min_trials)) +
  ggplot2::ylab("count of trial pairs") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle=90, hjust=1), axis.title.x = ggplot2::element_blank())

if(SAVE_GRAPH_FLAG == 1){
  ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("trial-nums-%s-by-list-%d-min-trials.pdf", source_mem_filter, min_trials)),
                  height = 10, width = 7)
}

z_mm %>%
  # to look at trial numbers, should just need a single ROI and hemi
  dplyr::filter(roi == "CA2_3_DG_body") %>%
  dplyr::filter(hemi == "right") %>%
  dplyr::filter(source_mem == "fHitANDmiss") %>%
  dplyr::filter(grepl("list[1-8]", list_fact)) %>%
  dplyr::group_by(list_fact, subj_fact) %>%
  dplyr::tally() %>%
  ggplot2::ggplot(ggplot2::aes(x = list_fact, y = n)) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(ggplot2::aes(color = subj_fact)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = min_trials)) +
  ggplot2::annotate("text", x = 4, y = min_trials, label = sprintf("Cutoff of %d trials", min_trials)) +
  ggplot2::ylab("count of trial pairs") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle=90, hjust=1), axis.title.x = ggplot2::element_blank())

if(SAVE_GRAPH_FLAG == 1){
  ggplot2::ggsave(filename = file.path(graph_out_dir, sprintf("trial-nums-fHitANDmiss-by-list-%d-min-trials.pdf", min_trials)),
                  height = 10, width = 7)
}

#' ## figure out which subjects would be dropped w/ current trial number cutoff
low_trials_by_list_subj <- z_mm %>%
  # to look at trial numbers, should just need a single ROI and hemi
  dplyr::filter(roi == "CA2_3_DG_body") %>%
  dplyr::filter(hemi == "right") %>%
  dplyr::filter(source_mem %in% source_mem_filter) %>%
  dplyr::filter(grepl("list[1-8]", list_fact)) %>%
  dplyr::group_by(list_fact, subj_fact) %>%
  dplyr::tally() %>%
  dplyr::filter(n <= min_trials)

low_trials_by_list_subj_nums <- unique(low_trials_by_list_subj$subj_fact)
message(sprintf("There are %d subjects with less than %d trials when splitting by list1-8.", length(low_trials_by_list_subj_nums), min_trials))

low_trials_questXlist_subj <- z_mm %>%
  # to look at trial numbers, should just need a single ROI and hemi
  dplyr::filter(roi == "CA2_3_DG_body") %>%
  dplyr::filter(hemi == "right") %>%
  dplyr::filter(source_mem %in% source_mem_filter) %>%
  dplyr::filter(question_fact %in% c("same", "diff"),
                list_fact %in% c("same", "diff")) %>%
  dplyr::group_by(list_fact, question_fact, subj_fact) %>%
  dplyr::tally() %>%
  dplyr::filter(n <= min_trials)
low_trials_by_questXlist_subj_nums <- unique(low_trials_questXlist_subj$subj_fact)
message(sprintf("There are %d subjects with less than %d trials when splitting by same/diff question and list.", length(low_trials_by_questXlist_subj_nums), min_trials))

low_trials_questXlist_fHitANDmisses_subj <- z_mm %>%
  # to look at trial numbers, should just need a single ROI and hemi
  dplyr::filter(roi == "CA2_3_DG_body") %>%
  dplyr::filter(hemi == "right") %>%
  dplyr::filter(item_mem %in% c("rHit", "fHitANDmiss")) %>%
  dplyr::filter(question_fact %in% c("same", "diff"),
                list_fact %in% c("same", "diff")) %>%
  dplyr::group_by(item_mem, list_fact, question_fact, subj_fact) %>%
  dplyr::tally() %>%
  dplyr::filter(n <= min_trials)
low_trials_by_questXlist_fHitANDmisses_subj_nums <- unique(low_trials_questXlist_fHitANDmisses_subj$subj_fact)

message(sprintf("There are %d subjects with less than %d trials when splitting by same/diff question and list for fhit and miss trials.", length(low_trials_by_questXlist_fHitANDmisses_subj_nums), min_trials))

save("min_trials", "low_trials_questXlist_fHitANDmisses_subj", "low_trials_questXlist_subj", "low_trials_by_list_subj", "exclude_subjects",
     file = file.path(analyzed_mri_dir, "excluded_subjects_variables.RData"))

#' ## Filter z_mm based on minimum trial numbers requirement
z_mm_filt <- z_mm %>%
  # this will filter out subjects who don't have enough trials, but won't fail if there are 0 subjects to exclude
  dplyr::filter(!subj_fact %in% low_trials_questXlist_subj) %>%
  dplyr::filter(item_mem == "rHit") %>%
  dplyr::filter(source_mem == source_mem_filter)

message(sprintf("Now analyzing data from %d subjects after filtering on %d min trials and %s source memory.", length(unique(z_mm_filt$subj_fact)), min_trials, source_mem_filter))

#' ## Check that all subjects have all ROIs
z_mm_filt %>%
  dplyr::group_by(subj) %>%
  dplyr::summarise(num_rois = length(unique(roi_formatted))) %>%
  ggplot2::ggplot(ggplot2::aes(x = subj, y = num_rois)) +
  ggplot2::geom_point()

total_num_rois <- length(unique(z_mm_filt$roi_formatted))
subj_roi_counts <- z_mm_filt %>%
  dplyr::group_by(subj) %>%
  dplyr::summarise(num_rois = length(unique(roi_formatted)))
message(sprintf("The following subjects have less than the total possible number of ROIs (%d): ", total_num_rois))
subj_roi_counts %>%
  dplyr::filter(num_rois < total_num_rois)

#' ## Plot trial pair numbers
by_quest_by_subj <- z_mm_filt %>%
  dplyr::ungroup() %>%
  dplyr::filter(question %in% c("same", "diff")) %>%
  dplyr::group_by(subj, roi_formatted, question_fact) %>%
  dplyr::rename(value = r) %>%
  dplyr::summarise(mean_PS = mean(value, na.rm = TRUE),
                   sd_PS = sd(value, na.rm = TRUE),
                   N = length(value), # this will reflect the number of trial pairs, not subjects
                   se_PS = sd_PS / sqrt(N),
                   upper = (mean_PS + se_PS),
                   lower = (mean_PS - se_PS))

message("plot number of trial pairs (ideally, should have same number of trial pairs for all ROIs)")
by_quest_by_subj %>%
  dplyr::filter(roi_formatted %in% c("CA1", "CA23DG", "SUB")) %>%
  barplot_by_quest_trialnums(.)
print(p)
if(SAVE_GRAPH_FLAG == 1){
  ggplot2::ggsave(filename = file.path(graph_out_dir, "CA1-CA23DG-SUB_hemi-combo_same-diff-question_trialnums.pdf"),
                  width = 8, height = 6)
}

#' # Print out R session info
sessionInfo()
