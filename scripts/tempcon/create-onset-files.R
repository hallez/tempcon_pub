#' ---
#' title: TempCon onset files
#' author: Halle R. Dimsdale-Zucker
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    theme: spacelab
#' ---

#+ label="Setup stuff"
#+ initialize, warning = FALSE, message = FALSE
devtools::load_all()
SAVE_FLAG <- 1

# add dplyr with library()
# NB: this is non-standard
# correct way would be to make this script into a function,
# store in the R/ directory,
# and use @importFrom dplyr "%>%",
# but this is incompatible with getting in-line results with code in knitr output
library(dplyr)

#' ## Load in config file
# this won't knit
# use packge::function syntax to avoid confusion
config <- yaml::yaml.load_file("../../config.yml")

#' ## Setup paths
project_dir <- ("../../")
analyzed_mri_dir <- halle::ensure_trailing_slash(config$directories$analyzed_mri)
raw_behavioral_dir <- halle::ensure_trailing_slash(config$directories$raw_behavioral)
analyzed_behavioral_dir <- halle::ensure_trailing_slash(config$directories$analyzed_behavioral)

#' ## Figure out subjects
# assume that subjects have folders in `raw_behavioral_dir` that start with `s` and are followed by one or two digits
subjects <- c(list.files(path=raw_behavioral_dir,pattern="^[s][(123456789)]$|^[s][(123456789)][(0123456789)]$"))

# for subjects who were ran in the backroom, they have slightly different filenames
# subjects s11 and s12 were run in the backroom and have `backroom` appended onto both their folder AND their datafiles
# the rest of the backroom subjects just have `backroom` appended to their datafiles
backroom_subjects <- c(11, 12, 14:23, 26:28, 31)
backroom_subjects_formatted <- NULL
for(i in 1:length(backroom_subjects)){
  backroom_subjects_formatted[i] <- halle::prepend_s_to_subject_id(backroom_subjects[i])
}

#' ## Load excluded subjects
# this file gets created by `rsa-tidy-data-btwn-runs.R`
rsa_exclude_subjects_fpath <- file.path(analyzed_mri_dir, "excluded_subjects_variables.RData")

if(file.exists(rsa_exclude_subjects_fpath)) {
  load(file = rsa_exclude_subjects_fpath)
  exclude_subjects_formatted <- exclude_subjects
} else {
  exclude_subjects <- c(1, 2, 3, 11, 13)

  exclude_subjects_formatted <- NULL
  for(i in 1:length(exclude_subjects)){
    exclude_subjects_formatted[i] <- halle::prepend_s_to_subject_id(exclude_subjects[i])
  }
}

message(sprintf("Excluded subjects: %s", paste0(exclude_subjects, collapse = ', ')))

subjects <- subjects[!is.element(subjects, exclude_subjects_formatted)]

#' ## List included subjects
print(subjects)
length(subjects)

#' ## TRs by subject
expected_num_trs <- 268
tr_in_secs <- 2.01
tr_buffer <- 16 # could be more conservative, but 16 seconds after the last trial seems reasonable
subj_262_trs <- c("s1", "s2", "s3", "s4", "s5", "s10")
subj_267_trs <- c("s6", "s7", "s8", "s9")

#' # Load in data
# these are the files that are created by `load_data.R`
# this means that we only need to deal w/ the Presentation annoyances of file naming, etc. once
for(isubj in subjects){
  subj_analyzed_dir <- paste0(analyzed_behavioral_dir, halle::ensure_trailing_slash(isubj))

  # load doesn't like when you try to set it directly to a variable
  # so this will load in a variable called `cur_data_scored_for_onsets`
  # this file is created by `behav-tidy.R`
  cur_data_scored_for_onsets <- NULL
  load(paste0(subj_analyzed_dir,paste0(isubj,"_concat_for_onsets.Rdata")))

  # ---figure out onset times---
  # based on how MM does this in `onsets_rkn_hitmiss_0dur.m`
  run_onsets <- NULL
  run_onsets <-
    cur_data_scored_for_onsets %>%
    # just grab the trials that are the fixation cross at the start of each run
    dplyr::filter(listnum == 0) %>%
    # this is confusing b/c "oldnew" doesn't make sense as a column label,
    # but MM is pulling column 4 when computing run start time in this script
    dplyr::mutate(run_onset_time = oldnew / 10000) %>%
    # get rid of unnecessary columns
    dplyr::select(item_recog_run, subj, run_onset_time)

  # print this out to take a look
  run_onsets

  # calculate onset times for trials
  # again, this is based on how MM does it in `onsets_rkn_hitmiss_0dur.m`
  # a nice double check is that the `trial_onset_time_from_start` should be equal to `item_recog_jitter`
  # for the fixation cross (ie, first row w/in a run)
  for(irun in 1:6){
    cur_run <- dplyr::filter(cur_data_scored_for_onsets, item_recog_run == irun)
    cur_run_onset <- dplyr::filter(run_onsets, item_recog_run == irun)
    cur_run_w_onsets <-
      cur_run %>%
      dplyr::mutate(trial_onset_time_from_start = (trial_onset_time/10000) - cur_run_onset$run_onset_time)

    if(irun==1){
      cur_data_w_onsets <- NULL
      cur_data_w_onsets <- cur_run_w_onsets
    } else {
      cur_data_w_onsets <- dplyr::full_join(cur_data_w_onsets, cur_run_w_onsets, by = intersect(names(cur_data_w_onsets), names(cur_run_w_onsets)))
    }
  } #for(irun in 1:6)

  # check that the dimensions make sense once we've calculated trial onset times
  if(dim(cur_data_w_onsets)[1] != 6*62){
    print(sprintf("WARNING: %s does not have the correct number of rows in `cur_data_w_onsets`.", isubj))
  }

  # ---figure out subject-specific number of TRs---
  cur_num_trs <- NULL
  if(isubj %in% subj_262_trs){
    message(sprintf("%s has 262 TRs", isubj))
    cur_num_trs <- 262
  } else if(isubj %in% subj_267_trs){
    message(sprintf("%s has 267 TRs", isubj))
    cur_num_trs <- 267
  } else {
    message(sprintf("%s has 268 TRs (as expected", isubj))
    cur_num_trs <- expected_num_trs
  }

  # ---create regressors---
  # for now, just worry about item recognition hits
  # in the future, could also consider temporal source hits, question source hits, or some combination of all three
  # remember, `*_scored` columns are created in `behav-tidy.R`
  cur_regressors <- NULL
  cur_regressors <-
    cur_data_w_onsets %>%
    # add in information about number of TRs (so know if need to exlude trials at the end of runs)
    dplyr::mutate(num_trs = cur_num_trs,
                  trial_onset_plus_buffer = trial_onset_time_from_start + tr_buffer,
                  exclude_trial_exceeds_trs = ifelse(trial_onset_plus_buffer > (num_trs * tr_in_secs), "exclude", "keep")) %>%
    dplyr::mutate(itemHit_model = ifelse(exclude_trial_exceeds_trs == "exclude", 99,
                                         ifelse(item_recog_scored == "rhitNT", 1,
                                           ifelse(item_recog_scored == "rhitT", 1,
                                                  ifelse(item_recog_scored == "fhit", 2,
                                                         ifelse(item_recog_scored == "miss", 3,
                                                                ifelse(item_recog_scored == "fa", 4,
                                                                       ifelse(item_recog_scored == "cr", 5, 99))))))),
                  itemHit_key = stringr::str_replace_all(itemHit_model, "1", "rhit"),
                  itemHit_key = stringr::str_replace_all(itemHit_key, "2", "fhit"),
                  itemHit_key = stringr::str_replace_all(itemHit_key, "3", "miss"),
                  itemHit_key = stringr::str_replace_all(itemHit_key, "4", "fa"),
                  itemHit_key = stringr::str_replace_all(itemHit_key, "5", "cr"),
                  itemHit_key = stringr::str_replace_all(itemHit_key, "99", "excludeTrial"),
                  itemHitBYType_model = ifelse(exclude_trial_exceeds_trs == "exclude", 99,
                                               ifelse(item_recog_scored == "rhitNT", 1,
                                                 ifelse(item_recog_scored == "rhitT", 2,
                                                        ifelse(item_recog_scored == "fhit", 3,
                                                               ifelse(item_recog_scored == "miss", 4,
                                                                      ifelse(item_recog_scored == "fa", 5,
                                                                             ifelse(item_recog_scored == "cr", 6, 99))))))),
                  itemHitBYType_key = stringr::str_replace_all(itemHitBYType_model, "1", "rhitNT"),
                  itemHitBYType_key = stringr::str_replace_all(itemHitBYType_key, "2", "rhitT"),
                  itemHitBYType_key = stringr::str_replace_all(itemHitBYType_key, "3", "fhit"),
                  itemHitBYType_key = stringr::str_replace_all(itemHitBYType_key, "4", "miss"),
                  itemHitBYType_key = stringr::str_replace_all(itemHitBYType_key, "5", "fa"),
                  itemHitBYType_key = stringr::str_replace_all(itemHitBYType_key, "6", "cr"),
                  itemHitBYType_key = stringr::str_replace_all(itemHitBYType_key, "99", "excludeTrial"),
                  collapsedFhitMissFA_model = ifelse(exclude_trial_exceeds_trs == "exclude", 99,
                                                     ifelse(item_recog_scored == "rhitNT", 1,
                                                        ifelse(item_recog_scored == "rhitT", 1,
                                                               ifelse(item_recog_scored == "fhit", 2,
                                                                      ifelse(item_recog_scored == "miss", 2,
                                                                             ifelse(item_recog_scored == "fa", 2,
                                                                                    ifelse(item_recog_scored == "cr", 3, 99))))))),
                  collapsedFhitMissFA_key = stringr::str_replace_all(collapsedFhitMissFA_model, "1", "rhit"),
                  collapsedFhitMissFA_key = stringr::str_replace_all(collapsedFhitMissFA_key, "2", "FhitMissFA"),
                  collapsedFhitMissFA_key = stringr::str_replace_all(collapsedFhitMissFA_key, "2", "FhitMissFA"),
                  collapsedFhitMissFA_key = stringr::str_replace_all(collapsedFhitMissFA_key, "2", "FhitMissFA"),
                  collapsedFhitMissFA_key = stringr::str_replace_all(collapsedFhitMissFA_key, "3", "cr"),
                  collapsedFhitMissFA_key = stringr::str_replace_all(collapsedFhitMissFA_key, "99", "excludeTrial"))

  head(cur_regressors)

  # let's sanity check ourselves
  # sorting by `item_recog_scored` or by `itemHitBYType_model` should be the same
  cur_regressors %>%
    dplyr::group_by(item_recog_scored) %>%
    dplyr::summarise(count = n())

  cur_regressors %>%
    dplyr::group_by(itemHitBYType_model) %>%
    dplyr::summarise(count = n())

  # and if we sort by `itemHit_model`, we should collapse the `1` and `2` responses from `itemHitBYType_model`
  cur_regressors %>%
    dplyr::group_by(itemHit_model) %>%
    dplyr::summarise(count = n())

  # ---save out onsets file---
  cur_onsets <- NULL
  cur_onsets <-
    cur_regressors %>%
    dplyr::rename(onset = trial_onset_time_from_start) %>%
    dplyr::mutate(duration = 0) %>%
    # eliminate `NA` onset rows
    dplyr::filter(onset != "NA") %>%
    dplyr::select(subj, onset, duration, item_recog_run, objnum, item_recog_row_number,
                  temporal_scored, question_scored_exact, question_scored_liberal,
                  exclude_trial_exceeds_trs,
                  itemHit_model, itemHit_key,
                  itemHitBYType_model, itemHitBYType_key,
                  collapsedFhitMissFA_model, collapsedFhitMissFA_key)

  if(SAVE_FLAG == 1){
    onsets_fname <- paste0(analyzed_behavioral_dir, halle::ensure_trailing_slash(isubj), "onsets_file.csv")
    write.csv(cur_onsets, file = onsets_fname, row.names=FALSE, quote=FALSE)
  }

  #' ---create trial ID matrix---
  trial_info_mtx <- NULL
  trial_info_mtx <-
    cur_regressors %>%
    # because the first row is always the fixation cross, just subtract 1 from `item_recog_row_number`
    # this feels a bit risky, but I'm 99% sure this is always how the Presentation output was formatted
    #
    # for the other columns, 1 = true (based on column label), 0 = false
    # this should allow easy indexing when using this in matlab
    dplyr::mutate(recog_trial_id = item_recog_row_number - 1,
                  item_temporal_rhit = ifelse(item_recog_scored == "rhitT", 1, 0),
                  item_non_temporal_rhit = ifelse(item_recog_scored == "rhitNT", 1, 0),
                  item_fhit = ifelse(item_recog_scored == "fhit", 1, 0),
                  item_miss = ifelse(item_recog_scored == "miss", 1, 0),
                  temporal_exact_hit = ifelse(temporal_scored == "hit", 1, 0),
                  temporal_plus_minus_one_hit = ifelse((temporal_scored == "hitPlusOne" | temporal_scored == "hitMinusOne"), 1, 0),
                  temporal_miss = ifelse(temporal_scored == "incorrect", 1, 0),
                  question_exact_hit = ifelse(question_scored_exact == "hit", 1, 0),
                  question_same_category_hit = ifelse(question_scored_liberal == "hitSC", 1, 0),
                  question_exact_miss = ifelse(question_scored_exact == "incorrect", 1, 0),
                  question_liberal_miss = ifelse(question_scored_liberal == "incorrect", 1, 0),
                  exclude_behavioral = ifelse(item_recog_scored == "null" | exclude_trial_exceeds_trs == "exclude", 1, 0)) %>%
    # replace the NA values in the `temporal` and `question` columns w/ 0s
    # these seem to propogate from the `temporal_scored` and `question_scored_*` columns
    dplyr::mutate(temporal_exact_hit = ifelse(is.na(temporal_exact_hit), 0, temporal_exact_hit),
                  temporal_plus_minus_one_hit = ifelse(is.na(temporal_plus_minus_one_hit), 0, temporal_plus_minus_one_hit),
                  temporal_miss = ifelse(is.na(temporal_miss), 0, temporal_miss),
                  question_exact_hit = ifelse(is.na(question_exact_hit), 0, question_exact_hit),
                  question_same_category_hit = ifelse(is.na(question_same_category_hit), 0, question_same_category_hit),
                  question_exact_miss = ifelse(is.na(question_exact_miss), 0, question_exact_miss),
                  question_liberal_miss = ifelse(is.na(question_liberal_miss), 0, question_liberal_miss)) %>%
    # eliminate the rows that correspond to the fixation crosses from item_recog
    dplyr::filter(!is.na(trial_onset_time_from_start)) %>%
    # now that these garbage rows are removed, calculate the absoluate (ie, across run) trial ID from recog
    dplyr::mutate(recog_abs_trial_id = ((item_recog_run - 1) * 60) + recog_trial_id) %>%
    # select the columns that may be useful to filter on when selecting trial correlation pairs of interest
    # also include onset information so can logically check that have correct trials
    dplyr::select(subj, enc_trial_id, recog_trial_id, recog_abs_trial_id, objnum, item_recog_run, trial_onset_time_from_start, listnum, enctask, item_temporal_rhit:exclude_behavioral, exclude_trial_exceeds_trs) %>%
    # ensure the data is in the same order as the trial x trial correlation pairs will be
    dplyr::arrange(item_recog_run, recog_trial_id) %>%
    # duplicate the recog_abs_trial_id column for matching up in `rsa-load-data-btwn-runs.R`
    dplyr::mutate(row_recog_trial_id = recog_abs_trial_id,
                  col_recog_trial_id = recog_abs_trial_id)

  # --- add subject-specific trials to be excluded ---
  # exclude encoding list 8 for s8 because they may have seen items twice
  trial_info_mtx$exclude_behavioral[trial_info_mtx$subj == "s8" & trial_info_mtx$listnum == 8] <- 1
  # exclude trial 1 in run 1 for s8 b/c get a "file too small" error when read in single trial beta
  trial_info_mtx$exclude_behavioral[trial_info_mtx$subj == "s8" & trial_info_mtx$item_recog_run == 1 & trial_info_mtx$recog_trial_id == 1] <- 1
  # this beta also seems to be corrupted - exclude
  trial_info_mtx$exclude_behavioral[trial_info_mtx$subj == "s8" & trial_info_mtx$item_recog_run == 6 & trial_info_mtx$recog_trial_id == 32] <- 1
  # exclude s14 encoding list 8 items
  trial_info_mtx$exclude_behavioral[trial_info_mtx$subj == "s14" & trial_info_mtx$listnum == 8] <- 1
  # exclude s15 run 1 b/c subject re-saw first few items
  # since we don't know how many items were repeated, be conservative and drop entire run
  trial_info_mtx$exclude_behavioral[trial_info_mtx$subj == "s15" & trial_info_mtx$item_recog_run == 1] <- 1
  # exclude s25 encoding list 1 b/c re-started and may have seen items multiple times
  trial_info_mtx$exclude_behavioral[trial_info_mtx$subj == "s25" & trial_info_mtx$listnum == 1] <- 1

  head(trial_info_mtx)

  if(SAVE_FLAG == 1){
    trial_info_fname <- paste0(analyzed_behavioral_dir, halle::ensure_trailing_slash(isubj), "trial_information.csv")
    write.csv(trial_info_mtx, file = trial_info_fname, row.names=FALSE, quote=FALSE)
  }

} #for(isubj in subjects){

#' # Print out R session info
sessionInfo()
