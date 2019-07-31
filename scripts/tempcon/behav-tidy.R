#' ---
#' title: Load Tempcon behavioral data
#' author: Halle R. Dimsdale-Zucker
#' output:
#'  html_document:
#'    toc: true
#'    toc_depth: 5
#'    toc_float:
#'      collapsed: false
#'      smooth_scroll: false
#'    number_sections: true
#'    theme: spacelab
#'    code_folding: "hide"
#' ---

#' # Setup
#+ label="Setup stuff"
#+ initialize, warning = FALSE, message = FALSE
devtools::load_all()
# following instructions from http://kbroman.org/pkg_primer/pages/github.html
# this ensures that the most recent version of the `halle` package also gets added
devtools::install_github("hallez/halle", subdir="halle")

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
raw_behavioral_dir <- halle::ensure_trailing_slash(config$directories$raw_behavioral)
analyzed_behavioral_dir <- halle::ensure_trailing_slash(config$directories$analyzed_behavioral)
dropbox_dir <- halle::ensure_trailing_slash(config$directories$dropbox_tempcon)

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

# chuck out any subjects who need to be excluded
# at this point, only exclude subjects who we know have missing behavioral data
# (assume we're blind to knowing about whether or not certain subjects get included in the fMRI analyses)
exclude_subjects <- c(2, 3)
exclude_subjects_formatted <- NULL
for(i in 1:length(exclude_subjects)){
  exclude_subjects_formatted[i] <- halle::prepend_s_to_subject_id(exclude_subjects[i])
}
subjects <- subjects[!is.element(subjects, exclude_subjects_formatted)]

# save out excluded subjects so can easily pass to other scripts
save(exclude_subjects_formatted, file = file.path(analyzed_behavioral_dir, "behav-exclude-subjects.RData"))

#' ### List included subjects
print(subjects)
length(subjects)

#' ## Define other variables
SAVE_FLAG <- 1
GRAPH_FLAG <- 1

#' ## Setup blank variable for group data
all_subj <- NULL
all_subj <- data.frame()

#' # Load in the data
# these are the files that are created by `load_data.R`
# this means that we only need to deal w/ the Presentation annoyances of file naming, etc. once
for(isubj in subjects){
  print(sprintf("Current subject is %s", isubj))
  subj_analyzed_dir <- paste0(analyzed_behavioral_dir, halle::ensure_trailing_slash(isubj))

  # --- load in encoding, item recog, temporal source, and question source files ---
  # load doesn't like when you try to set it directly to a variable
  # so just rename once the data has been loaded
  enc_dat <- NULL
  item_recog_dat <- NULL
  temporal_source_dat <- NULL
  question_source_dat <- NULL

  load(paste0(subj_analyzed_dir,paste0(isubj,"_encoding.Rdata")))
  cur_enc_data <- enc_dat

  load(paste0(subj_analyzed_dir,paste0(isubj,"_item_recog.Rdata")))
  cur_item_recog_data <- item_recog_dat

  load(paste0(subj_analyzed_dir,paste0(isubj,"_temporal_source_mem.Rdata")))
  cur_temporal_source_data <- temporal_source_dat

  load(paste0(subj_analyzed_dir,paste0(isubj,"_question_source_mem.Rdata")))
  cur_question_source_data <- question_source_dat

  # --- clarify column labels that are different between files ---
  # item recog
  cur_item_recog_data_renamed <- NULL
  cur_item_recog_data_tidy <- NULL
  cur_item_recog_data_renamed <-
    cur_item_recog_data %>%
    dplyr::rename(item_recog_response = response,
                  item_recog_jitter = jitter,
                  item_recog_run = run,
                  item_recog_row_number = row_number,
                  item_scored_from_presentation = scored)
  head(cur_item_recog_data_renamed)
  # replace "other" label in the `listnum` column so can merge w/ other dataframes
  # these are the labels for the fixation trials that start each run
  class(cur_item_recog_data_renamed$listnum)
  cur_item_recog_data_renamed$listnum[cur_item_recog_data_renamed$listnum == "other"] <- 0
  cur_item_recog_data_renamed$listnum <- as.numeric(cur_item_recog_data_renamed$listnum)
  class(cur_item_recog_data_renamed$listnum)
  # replace blank values and `next` in `item_recog_response`
  class(cur_item_recog_data_renamed$item_recog_response)
  cur_item_recog_data_renamed$item_recog_response[cur_item_recog_data_renamed$item_recog_response == ""] <- NA
  cur_item_recog_data_renamed$item_recog_response[cur_item_recog_data_renamed$item_recog_response == "next"] <- NA
  cur_item_recog_data_renamed$item_recog_response <- as.numeric(cur_item_recog_data_renamed$item_recog_response)
  class(cur_item_recog_data_renamed$item_recog_response)
  # use consistent variable name to match other dataframes
  cur_item_recog_data_tidy <- cur_item_recog_data_renamed
  head(cur_item_recog_data_tidy)

  # temporal source
  cur_temporal_source_data_tidy <- NULL
  cur_temporal_source_data_tidy <-
    cur_temporal_source_data %>%
    dplyr::rename(temporal_source_response = response,
                  temporal_source_jitter = jitter,
                  temporal_source_run = run,
                  temporal_source_row_number = row_number,
                  temporal_scored_from_presentation = scored)
  head(cur_temporal_source_data_tidy)

  # question source
  cur_question_source_data_tidy <- NULL
  cur_question_source_data_tidy <-
    cur_question_source_data %>%
    dplyr::rename(question_source_response = response,
                  question_source_jitter = jitter,
                  question_source_run = run,
                  question_source_row_number = row_number,
                  question_scored_from_presentation = scored)
  head(cur_question_source_data_tidy)

  # --- merge together data from item recognition, temporal source memory, and question source memory ---
  enc_item <- NULL
  enc_item_temporal <- NULL
  cur_all_data <- NULL

  enc_item <- dplyr::full_join(cur_item_recog_data_tidy, cur_enc_data, by = intersect(names(cur_item_recog_data_tidy), names(cur_enc_data)))
  head(enc_item)
  dim(enc_item)

  enc_item_temporal <- dplyr::full_join(enc_item,cur_temporal_source_data_tidy, by = intersect(names(enc_item), names(cur_temporal_source_data_tidy)))
  head(enc_item_temporal)
  dim(enc_item_temporal)

  cur_all_data <- dplyr::full_join(enc_item_temporal,cur_question_source_data_tidy, by = intersect(names(enc_item_temporal), names(cur_question_source_data_tidy)))
  head(cur_all_data)
  dim(cur_all_data)

  # ---score the memory data---
  # we're at the mercy of MM's scripts in terms of figuring out what the trials and labels are representing
  # this is based on `onsets_rkn_hitmiss_0dur.m`
  cur_data_scored_for_onsets <- NULL
  cur_data_scored_for_onsets <-
    cur_all_data %>%
    # we're going to be relying on combining `mutate` and `ifelse`
    # based on: https://rstudio-pubs-static.s3.amazonaws.com/116317_e6922e81e72e4e3f83995485ce686c14.html#/6
    # ---temporal scoring---
    # if listnum == temporal_source_response, this is a temporal hit
    # temporal hits correspond to MM temporal(i,7) == 1
    # within +/- 1 correspond to MM temporal(i,8) == 1
    # based on `make_stimuli.m`, seems like `listnum == 9` is for new items
    # there are no notes about what it means when `temporal_source_response == 9`
    # (since they should have only had 8 options), so throwing out these responses
    # ---question scoring---
    # question hits correspond to MM question(i,7) == 1 (`hit`)
    # question same category hits correspond to MM question(i,8) == 1 (`hitSC`)
    # ---object recognition scoring---
    # for object recognition, response options were:
    # 1 = new, 2 = familiar, 3 = recollect non-temporal details, 4 = recollect temporal details
    # Presentation seems to leave column completely blank when there was no response therefore mistakenly gets read as next column over (and these are all very large values)
    dplyr::mutate(temporal_scored = as.factor(ifelse(100 > temporal_source_response & temporal_source_response > 8, "response exceeds scale",
                                           ifelse(temporal_source_response > 100, "no resp",
                                           ifelse(is.na(temporal_source_response), NA,
                                           ifelse(listnum == temporal_source_response, "hit",
                                           ifelse(listnum == temporal_source_response+1, "hitPlusOne",
                                           ifelse(listnum == temporal_source_response-1, "hitMinusOne", "incorrect"))))))),
                  temporal_scored_exact = as.factor(ifelse(100 > temporal_source_response & temporal_source_response > 8, "response exceeds scale",
                                                     ifelse(temporal_source_response > 100, "no resp",
                                                            ifelse(is.na(temporal_source_response), NA,
                                                                   ifelse(listnum == temporal_source_response, "correct", "incorrect"))))),
                  temporal_scored_list_firstsecond_half = as.factor(case_when(temporal_source_response %in% c(1:4) & listnum %in% c(1:4) ~ "first_half_hit",
                                                                              temporal_source_response %in% c(5:8) & listnum %in% c(5:8) ~ "second_half_hit",
                                                                              temporal_source_response %in% c(1:4) & listnum %in% c(5:8) ~ "wrong_half",
                                                                              temporal_source_response %in% c(5:8) & listnum %in% c(1:4) ~ "wrong_half",
                                                                              temporal_source_response > 8 ~ "response exceeds scale")),
                  # for doing list half bias analyses, also recode responses based on list half
                  # without regard for whether or not the response was correct
                  temporal_response_as_list_firstsecond_half = as.factor(case_when(temporal_source_response %in% c(1:4) ~ "first_half",
                                                                              temporal_source_response %in% c(5:8) ~ "second_half",
                                                                              temporal_source_response > 8 ~ "response exceeds scale")),
                  # "same_half" are all hits, and "diff_half" are misses
                  temporal_scored_list_samediff_half = as.factor(case_when(temporal_source_response %in% c(1:4) & listnum %in% c(1:4) ~ "same_half_hit",
                                                                           temporal_source_response %in% c(5:8) & listnum %in% c(5:8) ~ "same_half_hit",
                                                                           temporal_source_response %in% c(1:4) & listnum %in% c(5:8) ~ "wrong_half",
                                                                           temporal_source_response %in% c(5:8) & listnum %in% c(1:4) ~ "wrong_half",
                                                                           temporal_source_response > 8 ~ "response exceeds scale")),
                  question_scored_exact = as.factor(ifelse(100 > question_source_response & question_source_response > 4, "response exceeds scale",
                                                     ifelse(question_source_response > 100, "no resp",
                                                     ifelse(is.na(question_source_response), NA,
                                                     ifelse(enctask == question_source_response, "hit", "incorrect"))))),
                  # based on `make_stimuli.m` line 181, it seems like question 1 = fridge, 2 = bathtub, 3 = supermarket, 4 = convenience store,
                  question_scored_liberal = as.factor(ifelse(100 > question_source_response & question_source_response > 4, "response exceeds scale",
                                           ifelse(question_source_response > 100, "no resp",
                                           ifelse(is.na(question_source_response), NA,
                                           ifelse(enctask == question_source_response, "hitSC",
                                           ifelse(enctask == 1 & question_source_response == 2, "hitSC",
                                           ifelse(enctask == 2 & question_source_response == 1, "hitSC",
                                           ifelse(enctask == 3 & question_source_response == 4, "hitSC",
                                           ifelse(enctask == 4 & question_source_response == 3, "hitSC", "incorrect"))))))))),
                  item_recog_scored = as.factor(ifelse(100 > item_recog_response & item_recog_response > 4, "response exceeds scale",
                                             ifelse(item_recog_response > 100, "no resp",
                                             ifelse(is.na(item_recog_response), NA,
                                             ifelse(oldnew == 1 & item_recog_response == 1, "miss",
                                             ifelse(oldnew == 1 & item_recog_response == 2, "fhit",
                                             ifelse(oldnew == 1 & item_recog_response == 3, "rhitNT",
                                             ifelse(oldnew == 1 & item_recog_response == 4, "rhitT",
                                             ifelse(oldnew == 2 & item_recog_response == 1, "cr",
                                             ifelse(oldnew == 2 & item_recog_response > 1, "fa", "null")))))))))),
                  item_recog_scored_rhit_collapsed = as.factor(ifelse(100 > item_recog_response & item_recog_response > 4, "response exceeds scale",
                                                            ifelse(item_recog_response > 4, "no resp",
                                                            ifelse(is.na(item_recog_response), NA,
                                                            ifelse(oldnew == 1 & item_recog_response == 1, "miss",
                                                            ifelse(oldnew == 1 & item_recog_response == 2, "fhit",
                                                            ifelse(oldnew == 1 & item_recog_response == 3, "rhit",
                                                            ifelse(oldnew == 1 & item_recog_response == 4, "rhit",
                                                            ifelse(oldnew == 2 & item_recog_response == 1, "cr",
                                                            ifelse(oldnew == 2 & item_recog_response > 1, "fa", "null")))))))))),
                  item_recog_response_label = as.factor(ifelse(100 > item_recog_response & item_recog_response > 4, "response exceeds scale",
                                                     ifelse(item_recog_response > 4, "no resp",
                                                     ifelse(is.na(item_recog_response), NA,
                                                     ifelse(item_recog_response == 1, "new",
                                                     ifelse(item_recog_response == 2, "fam",
                                                     ifelse(item_recog_response == 3, "rem",
                                                     ifelse(item_recog_response == 4, "rem", "null")))))))),
                  item_recog_scored_corr_incorr = as.factor(ifelse(item_recog_scored_rhit_collapsed == "no resp", "no resp",
                                                         ifelse(item_recog_scored_rhit_collapsed == "rhit", "correct",
                                                         ifelse(item_recog_scored_rhit_collapsed == "fhit", "correct",
                                                         ifelse(item_recog_scored_rhit_collapsed == "cr", "correct",
                                                         ifelse(item_recog_scored_rhit_collapsed == "fa", "incorrect",
                                                         ifelse(item_recog_scored_rhit_collapsed == "miss", "incorrect", item_recog_scored_rhit_collapsed))))))))

  summary(cur_data_scored_for_onsets)

  # save out merged and minimally scored file for each subject to be used by `create-onset-files.R`
  if(SAVE_FLAG == 1){
    fname_out <- paste0(subj_analyzed_dir,paste0(isubj,"_concat_for_onsets.Rdata"))
    save(cur_data_scored_for_onsets, file = fname_out)
  }

  # --- moar data cleanup ---
  # now that we've saved out the information we'll need to generate the onsets,
  # let's clean up the dataframe
  # first, let's look at the summary to see why we've identified these columns for cleaning
  # essentially, looking for values where the `mean` or `max` don't fit with the rest of the data
  summary(cur_data_scored_for_onsets)

  # now, actually clean up
  cur_data_scored_clean <- NULL
  cur_data_scored_clean <-
    cur_data_scored_for_onsets %>%
    # remove the columns that just contain 0s and NaNs
    dplyr::select(-temporal_source_jitter, -question_source_jitter) %>%
    # remove rows that are the start of a run and not a trial of interest
    # conveniently, these super high values seem to get written into the `listnum` column
    # this must just have been based on where the presentation script was writing out to (since it's the first column)
    # because it doesn't make sense to put this data in this column
    # this should remove 6 rows of data (b/c 6 item recog runs)
    dplyr::filter(listnum <= 9) %>%
    # remove rows that correspond to fixation time
    # according to `onsets_rkn_hitmiss_0dur.m`,
    # column 4 (which now corresponds to `oldnew`)
    # is where these values are contained
    # again, there should be 6 rows removed here b/c of the 6 item recog runs
    dplyr::filter(oldnew <= 2)

  summary(cur_data_scored_clean)

  # let's sanity check that we now have 360 rows
  # this accounts for the 288 old items + 72 new items seen in the scanner
  if(dim(cur_data_scored_clean)[1] != 360){
    cat(print("WARNING: Incorrect number of rows of data."))
  }
  # another way to check is that should have 72 trials for which `enc_trial_id`==NA
  if(dim(cur_data_scored_clean[is.na(cur_data_scored_clean$enc_trial_id),])[1] != 72){
    cat(print("\nWARNING: Incorrect number of trial IDs not present at encoding"))
  }

  # let's also check that for other variables w/ extreme-ish values, we can see why they occur
  # item_recog_jitter
  # not because we care about jitter, but just to ensure that all the wonky rows have been removed
  # this should basically look like a left-skewed histogram
  cur_data_scored_clean %>%
    ggplot2::ggplot(ggplot2::aes(x = item_recog_jitter)) +
    ggplot2::geom_histogram()
  # it seems like `16000` was the jitter max
  # so report if any values exceed this
  cur_data_scored_clean %>%
    dplyr::filter(item_recog_jitter > 16000)

  # temporal_source_response
  # what do the rows look like where the response exceeds the max response?
  cur_data_scored_clean %>%
    dplyr::filter(temporal_source_response > 8)
  # were those values also in the original file?
  # do they return the same values as in the cleaned dataframe?
  # this should help figure out if the wonky values were introduced by somethign like misalignment
  cur_temporal_source_data %>%
    dplyr::filter(response > 8)
  # conclusion: these large values occur for trials in which subjects made an invalid response (9) or no repsonse (temporal_scored_from_presentation = miss)

  # question_source_response
  # let's take the same approach as w/ `temporal_source_response`
  # what are the other values in the row where there's an abnormally high value?
  cur_data_scored_clean %>%
    dplyr::filter(question_source_response > 4)
  # what does the original file look like?
  cur_question_source_data %>%
    dplyr::filter(response > 4)
  # let's again confirm that these trials are ones where the subject made an invalid response (>4) or just didn't respond ("miss")
  cur_data_scored_clean %>%
    dplyr::filter(question_scored_from_presentation == "incorrect",
                  question_source_response > 4)
  cur_data_scored_clean %>%
    dplyr::filter(question_scored_from_presentation == "miss")
  summary(cur_data_scored_clean$question_source_response)

  # score the lag between where an item was studied and where someone said it was studied
  # relabel old/new with verbal labels
  cur_data_scored_clean <- cur_data_scored_clean %>%
    dplyr::mutate(temporal_scored_lag = factor(x = temporal_source_response - listnum, levels = c(-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7), ordered = TRUE),
                  oldnew_label = (dplyr::recode_factor(oldnew, `1` = "old", `2` = "new")))

  # --- combine subjects into a group dataframe variable ---
  if(dim(all_subj)[1] == 0){
    all_subj <- cur_data_scored_clean
  } else {
    all_subj <- dplyr::full_join(all_subj, cur_data_scored_clean, by = intersect(names(all_subj), names(cur_data_scored_clean)))
  }

} #for(isubj in subjects

#' # Take an overall peek at `all_subj` dataframe
summary(all_subj)

#' # Save out group subject dataframe
dim(all_subj)
unique(all_subj$subj)
length(unique(all_subj$subj))

# check to see if there are 360 rows per subject
if(dim(all_subj)[1]/360 != length(unique(all_subj$subj))){
  warning("Not all subjects have 360 rows of data.")
}

save(all_subj, file = file.path(analyzed_behavioral_dir, sprintf("all_subj_behav_data_N%d.RData", length(unique(all_subj$subj)))))

#' # Print out R session info
sessionInfo()
