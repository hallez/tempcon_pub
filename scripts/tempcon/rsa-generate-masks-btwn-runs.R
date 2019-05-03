#' ---
#' title: Create TempCon masks &mdash; Between runs
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

#' ## Define functions
write_mask <- function(input_dat, out_path){
  # name of variable as string based on: https://stackoverflow.com/questions/14577412/how-to-convert-variable-object-name-into-string
  # saving out matrix to txt based on: https://stackoverflow.com/questions/12907198/how-to-export-r-matrix-object-to-txt-file
  var_name <- deparse(substitute(input_dat))
  stopifnot(isSymmetric(input_dat))
  write.table(input_dat,
              file = paste0(halle::ensure_trailing_slash(out_path),
                            sprintf("%s.txt", var_name)),
              row.names = FALSE, col.names = FALSE)
}

#+ label = "subjects"
#' ## Figure out subjects
# assume that subjects have folders in `raw_behavioral_dir` that start with `s` and are followed by one or two digits
subjects <- c(list.files(path=raw_behavioral_dir,pattern="^[s][(123456789)]$|^[s][(123456789)][(0123456789)]$"))
# chuck out any subjects who need to be excluded
exclude_subjects <- c("s1", "s2", "s3", "s11", "s13")
subjects <- subjects[!is.element(subjects, exclude_subjects)]
length(subjects)
subjects

#+ label = "setup vars"
#' ## Set other variables
runs <- c(1:6)
num_trials_per_run <- 60

#' ## Setup w/in run NaN matrix
zeros_trial_by_trial_mtx <- matrix(data = 0, nrow = ((length(runs)*num_trials_per_run)), ncol = ((length(runs)*num_trials_per_run)))
ones_trial_by_trial_mtx <- matrix(data = 1, nrow = ((length(runs)*num_trials_per_run)), ncol = ((length(runs)*num_trials_per_run)))
NaN_trial_by_trial_mtx <- matrix(data = NaN, nrow = ((length(runs)*num_trials_per_run)), ncol = ((length(runs)*num_trials_per_run)))

rownames(ones_trial_by_trial_mtx) <- c(1:360)
colnames(ones_trial_by_trial_mtx) <- c(1:360)
ones_trial_by_trial_mtx_w_NaN <- ones_trial_by_trial_mtx
ones_trial_by_trial_mtx_w_NaN[1:60, 1:60] <- NaN
# for runs 2:6, can figure out start and end values this way
for(irun in c(2:6)){
  start_value <- num_trials_per_run * (irun - 1) + 1
  end_value <- num_trials_per_run * (irun)
  ones_trial_by_trial_mtx_w_NaN[start_value:end_value, start_value:end_value] <- NaN
}

within_run_corr_as_NaN <- ones_trial_by_trial_mtx_w_NaN # this will later be used to mask out w/in run correlations
write_mask(within_run_corr_as_NaN, file.path(analyzed_behavioral_dir))
lower_tri_mtx <- upper.tri(ones_trial_by_trial_mtx, diag = FALSE)
lower_tri_mtx[lower_tri_mtx == FALSE] <- NaN # this will later be used to ensure just grab half of mtx (since it's symmetric)

# to check: ggplot2::ggplot(reshape2::melt(lower_tri_mtx), ggplot2::aes(Var2, Var1, fill = value)) + ggplot2::geom_tile()

#+ label="read in data"
#+ warning = FALSE
#' # Read in data: Loop across subjects, hemispheres, rois, and runs
for(isubj in subjects){
  time_subj_start <- Sys.time()
  cat(sprintf("\nWorking on %s\n", isubj))

  # read in trial information file
  # this gets created in `create-onset-files.R`
  cur_subj_info_file <- NULL
  cur_subj_info_file <- paste0(halle::ensure_trailing_slash(analyzed_behavioral_dir),
                               halle::ensure_trailing_slash(isubj),
                               "trial_information.csv")

  if(!file.exists(cur_subj_info_file)){
    print(sprintf("\nNo subject info file (%s) - skipping.\n", cur_subj_info_file))
    next
  } else {
    cur_subj_info <- NULL
    cur_subj_info <- read.csv(cur_subj_info_file)

    mask_out_path <- paste0(halle::ensure_trailing_slash(analyzed_mri_dir),
                            halle::ensure_trailing_slash("multivariate_itemHitBYType"),
                            halle::ensure_trailing_slash(isubj))

    dir.create(mask_out_path, recursive = TRUE)

    # --- setup masks ---
    # let's only analyze rhit trials
    # in a perfect world, would restrict to temporal rhit trials
    # but this won't be enough trial pairs w/in each run for each subject
    # so, include `temporal` and `non_temporal` rhits
    #
    # also, let's make sure we're only selecting include trials
    # ie, `exclude_behavioral == 0`
    #
    # to check: ggplot2::ggplot(reshape2::melt(rHits_list1_list2_sameQuestion_mask), ggplot2::aes(Var2, Var1, fill = value)) + ggplot2::geom_tile()
    # ggplot2::ggsave(filename = file.path(halle::ensure_trailing_slash(config$directories$dropbox_tempcon),
    #                                      halle::ensure_trailing_slash("presentations"), halle::ensure_trailing_slash("sfn-2017"), halle::ensure_trailing_slash("plots"), "list1_list2_mask.pdf"),
    #                 height = 10, width = 10)
    #
    # for matrices that are not pulling the same index for both dimension (e.g., the off-diagonal between trial matrices), check that they're symmetric:
    # ggplot2::ggplot(reshape2::melt(combo11), ggplot2::aes(Var2, Var1, fill = value)) + ggplot2::geom_tile()
    # NaNs <- matrix(data = NaN, nrow = 10, ncol = 10)
    # idx1 <- c(1,2,4,7,9)
    # idx2 <- c(2,3,6)
    # # this is not symmetric
    # combo12 <- NaNs
    # combo12[idx1, idx2] <- 1
    #
    # # this is symmetric
    # combo12c <- NaNs
    # combo12c[idx1,idx2] <- 1
    # combo12c[idx2,idx1] <- 1
    #
    # # use this function to verify that the matrix is symmetric
    # isSymmetric(combo12c)

    rhits_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$item_temporal_rhit == 1 | cur_subj_info$item_non_temporal_rhit == 1)

    fhits_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$item_fhit ==1)

    item_miss_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                             cur_subj_info$item_miss == 1)

    fhitsANDmisses_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                                  cur_subj_info$item_fhit ==1 | cur_subj_info$item_miss == 1)

    temporal_source_hits_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                                        cur_subj_info$temporal_exact_hit == 1 | cur_subj_info$temporal_plus_minus_one_hit == 1)

    temporal_source_miss_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                                        cur_subj_info$temporal_miss == 1)

    # Distinguish between when subjects got the exact question correct (out of 4 options) vs just the fit/find stem correct
    question_exact_source_hits_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                                        cur_subj_info$question_exact_hit == 1)

    # this variable gets created in `create-onset-files.R`, l.205
    # 1 = hit, same question stem; 0 = incorrect/all other responses
    question_liberal_source_hits_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                                              cur_subj_info$question_same_category_hit == 1)

    question_exact_source_miss_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                                        cur_subj_info$question_exact_miss == 1)

    question_liberal_source_miss_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                                              cur_subj_info$question_liberal_miss == 1)

    list1_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$listnum == 1)
    list2_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$listnum == 2)
    list3_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$listnum == 3)
    list4_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$listnum == 4)
    list5_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$listnum == 5)
    list6_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$listnum == 6)
    list7_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$listnum == 7)
    list8_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                         cur_subj_info$listnum == 8)

    list1_mask <- NaN_trial_by_trial_mtx
    list1_mask[list1_idx, list1_idx] <- 1
    write_mask(list1_mask, mask_out_path)
    list2_mask <- NaN_trial_by_trial_mtx
    list2_mask[list2_idx, list2_idx] <- 1
    write_mask(list2_mask, mask_out_path)
    list3_mask <- NaN_trial_by_trial_mtx
    list3_mask[list3_idx, list3_idx] <- 1
    write_mask(list3_mask, mask_out_path)
    list4_mask <- NaN_trial_by_trial_mtx
    list4_mask[list4_idx, list4_idx] <- 1
    write_mask(list4_mask, mask_out_path)
    list5_mask <- NaN_trial_by_trial_mtx
    list5_mask[list5_idx, list5_idx] <- 1
    write_mask(list5_mask, mask_out_path)
    list6_mask <- NaN_trial_by_trial_mtx
    list6_mask[list6_idx, list6_idx] <- 1
    write_mask(list6_mask, mask_out_path)
    list7_mask <- NaN_trial_by_trial_mtx
    list7_mask[list7_idx, list7_idx] <- 1
    write_mask(list7_mask, mask_out_path)
    list8_mask <- NaN_trial_by_trial_mtx
    list8_mask[list8_idx, list8_idx] <- 1
    write_mask(list8_mask, mask_out_path)

    # if we care about what question and item was associated with, can use fit vs find masks based on `enctask`
    # however, if what we care about is whether the subject got the question type source correct,
    # then it makes more sense to use the `question_same_category_hit` column
    question_fit_idx <- which(cur_subj_info$enctask == 1 | cur_subj_info$enctask == 2)
    question_find_idx <- which(cur_subj_info$enctask == 3 | cur_subj_info$enctask == 4)

    question1_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                             cur_subj_info$enctask == 1)
    question2_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                             cur_subj_info$enctask == 2)
    question3_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                             cur_subj_info$enctask == 3)
    question4_idx <- which(cur_subj_info$exclude_behavioral == 0 &
                             cur_subj_info$enctask == 4)

    question1_mask <- NaN_trial_by_trial_mtx
    question1_mask[question1_idx, question1_idx] <- 1
    write_mask(question1_mask, mask_out_path)
    question2_mask <- NaN_trial_by_trial_mtx
    question2_mask[question2_idx, question2_idx] <- 1
    write_mask(question2_mask, mask_out_path)
    question3_mask <- NaN_trial_by_trial_mtx
    question3_mask[question3_idx, question3_idx] <- 1
    write_mask(question3_mask, mask_out_path)
    question4_mask <- NaN_trial_by_trial_mtx
    question4_mask[question4_idx, question4_idx] <- 1
    write_mask(question4_mask, mask_out_path)

    rhits_mask <- NaN_trial_by_trial_mtx
    rhits_mask[rhits_idx, rhits_idx] <- 1
    # in order to save out so matlab can read mask it, can't have NA values
    # could use 1s/0s, but 0 could be a plausible value
    # so opting to use `NaN` which matlab can recognize
    write_mask(rhits_mask, mask_out_path)

    fhits_mask <- NaN_trial_by_trial_mtx
    fhits_mask[fhits_idx, fhits_idx] <- 1
    write_mask(fhits_mask, mask_out_path)

    fhitsANDmisses_mask <- NaN_trial_by_trial_mtx
    fhitsANDmisses_mask[fhitsANDmisses_idx, fhitsANDmisses_idx] <- 1
    write_mask(fhitsANDmisses_mask, mask_out_path)

    item_miss_mask <- NaN_trial_by_trial_mtx
    item_miss_mask[item_miss_idx, item_miss_idx] <- 1
    write_mask(item_miss_mask, mask_out_path)

    temporal_hits_mask <- NaN_trial_by_trial_mtx
    temporal_hits_mask[temporal_source_hits_idx, temporal_source_hits_idx] <- 1
    write_mask(temporal_hits_mask, mask_out_path)

    temporal_miss_mask <- NaN_trial_by_trial_mtx
    temporal_miss_mask[temporal_source_miss_idx, temporal_source_miss_idx] <- 1
    write_mask(temporal_miss_mask, mask_out_path)

    question_exact_hits_mask <- NaN_trial_by_trial_mtx
    question_exact_hits_mask[question_exact_source_hits_idx, question_exact_source_hits_idx] <- 1
    write_mask(question_exact_hits_mask, mask_out_path)

    question_liberal_hits_mask <- NaN_trial_by_trial_mtx
    question_liberal_hits_mask[question_liberal_source_hits_idx, question_liberal_source_hits_idx] <- 1
    write_mask(question_liberal_hits_mask, mask_out_path)

    question_exact_miss_mask <- NaN_trial_by_trial_mtx
    question_exact_miss_mask[question_exact_source_miss_idx, question_exact_source_miss_idx] <- 1
    write_mask(question_exact_miss_mask, mask_out_path)

    question_liberal_miss_mask <- NaN_trial_by_trial_mtx
    question_liberal_miss_mask[question_liberal_source_miss_idx, question_liberal_source_miss_idx] <- 1
    write_mask(question_liberal_miss_mask, mask_out_path)

    same_list_mask <- NaN_trial_by_trial_mtx
    same_list_mask[list1_idx, list1_idx] <- 1
    same_list_mask[list2_idx, list2_idx] <- 1
    same_list_mask[list3_idx, list3_idx] <- 1
    same_list_mask[list4_idx, list4_idx] <- 1
    same_list_mask[list5_idx, list5_idx] <- 1
    same_list_mask[list6_idx, list6_idx] <- 1
    same_list_mask[list7_idx, list7_idx] <- 1
    same_list_mask[list8_idx, list8_idx] <- 1
    write_mask(same_list_mask, mask_out_path)

    diff_list_mask <- NaN_trial_by_trial_mtx
    diff_list_mask[list1_idx, list2_idx] <- 1
    diff_list_mask[list1_idx, list3_idx] <- 1
    diff_list_mask[list1_idx, list4_idx] <- 1
    diff_list_mask[list1_idx, list5_idx] <- 1
    diff_list_mask[list1_idx, list6_idx] <- 1
    diff_list_mask[list1_idx, list7_idx] <- 1
    diff_list_mask[list1_idx, list8_idx] <- 1
    diff_list_mask[list2_idx, list1_idx] <- 1
    diff_list_mask[list2_idx, list3_idx] <- 1
    diff_list_mask[list2_idx, list4_idx] <- 1
    diff_list_mask[list2_idx, list5_idx] <- 1
    diff_list_mask[list2_idx, list6_idx] <- 1
    diff_list_mask[list2_idx, list7_idx] <- 1
    diff_list_mask[list2_idx, list8_idx] <- 1
    diff_list_mask[list3_idx, list1_idx] <- 1
    diff_list_mask[list3_idx, list2_idx] <- 1
    diff_list_mask[list3_idx, list3_idx] <- 1
    diff_list_mask[list3_idx, list4_idx] <- 1
    diff_list_mask[list3_idx, list5_idx] <- 1
    diff_list_mask[list3_idx, list6_idx] <- 1
    diff_list_mask[list3_idx, list7_idx] <- 1
    diff_list_mask[list3_idx, list8_idx] <- 1
    diff_list_mask[list4_idx, list1_idx] <- 1
    diff_list_mask[list4_idx, list2_idx] <- 1
    diff_list_mask[list4_idx, list3_idx] <- 1
    diff_list_mask[list4_idx, list5_idx] <- 1
    diff_list_mask[list4_idx, list6_idx] <- 1
    diff_list_mask[list4_idx, list7_idx] <- 1
    diff_list_mask[list4_idx, list8_idx] <- 1
    diff_list_mask[list5_idx, list1_idx] <- 1
    diff_list_mask[list5_idx, list2_idx] <- 1
    diff_list_mask[list5_idx, list3_idx] <- 1
    diff_list_mask[list5_idx, list4_idx] <- 1
    diff_list_mask[list5_idx, list6_idx] <- 1
    diff_list_mask[list5_idx, list7_idx] <- 1
    diff_list_mask[list5_idx, list8_idx] <- 1
    diff_list_mask[list6_idx, list1_idx] <- 1
    diff_list_mask[list6_idx, list2_idx] <- 1
    diff_list_mask[list6_idx, list3_idx] <- 1
    diff_list_mask[list6_idx, list4_idx] <- 1
    diff_list_mask[list6_idx, list5_idx] <- 1
    diff_list_mask[list6_idx, list7_idx] <- 1
    diff_list_mask[list6_idx, list8_idx] <- 1
    diff_list_mask[list7_idx, list1_idx] <- 1
    diff_list_mask[list7_idx, list2_idx] <- 1
    diff_list_mask[list7_idx, list3_idx] <- 1
    diff_list_mask[list7_idx, list4_idx] <- 1
    diff_list_mask[list7_idx, list5_idx] <- 1
    diff_list_mask[list7_idx, list6_idx] <- 1
    diff_list_mask[list7_idx, list8_idx] <- 1
    diff_list_mask[list8_idx, list1_idx] <- 1
    diff_list_mask[list8_idx, list2_idx] <- 1
    diff_list_mask[list8_idx, list3_idx] <- 1
    diff_list_mask[list8_idx, list4_idx] <- 1
    diff_list_mask[list8_idx, list5_idx] <- 1
    diff_list_mask[list8_idx, list6_idx] <- 1
    diff_list_mask[list8_idx, list7_idx] <- 1
    write_mask(diff_list_mask, mask_out_path)

    # to grab trial pairs from between lists, create separate masks
    # --- list 1 ---
    list1_list2_mask <- NaN_trial_by_trial_mtx
    list1_list2_mask[list1_idx, list2_idx] <- 1
    list1_list2_mask[list2_idx, list1_idx] <- 1
    write_mask(list1_list2_mask, mask_out_path)
    list1_list3_mask <- NaN_trial_by_trial_mtx
    list1_list3_mask[list1_idx, list3_idx] <- 1
    list1_list3_mask[list3_idx, list1_idx] <- 1
    write_mask(list1_list3_mask, mask_out_path)
    list1_list4_mask <- NaN_trial_by_trial_mtx
    list1_list4_mask[list1_idx, list4_idx] <- 1
    list1_list4_mask[list4_idx, list1_idx] <- 1
    write_mask(list1_list4_mask, mask_out_path)
    list1_list5_mask <- NaN_trial_by_trial_mtx
    list1_list5_mask[list1_idx, list5_idx] <- 1
    list1_list5_mask[list5_idx, list1_idx] <- 1
    write_mask(list1_list5_mask, mask_out_path)
    list1_list6_mask <- NaN_trial_by_trial_mtx
    list1_list6_mask[list1_idx, list6_idx] <- 1
    list1_list6_mask[list6_idx, list1_idx] <- 1
    write_mask(list1_list6_mask, mask_out_path)
    list1_list7_mask <- NaN_trial_by_trial_mtx
    list1_list7_mask[list1_idx, list7_idx] <- 1
    list1_list7_mask[list7_idx, list1_idx] <- 1
    write_mask(list1_list7_mask, mask_out_path)
    list1_list8_mask <- NaN_trial_by_trial_mtx
    list1_list8_mask[list1_idx, list8_idx] <- 1
    list1_list8_mask[list8_idx, list1_idx] <- 1
    write_mask(list1_list8_mask, mask_out_path)
    # --- list 2 ---
    list2_list1_mask <- NaN_trial_by_trial_mtx
    list2_list1_mask[list2_idx, list1_idx] <- 1
    list2_list1_mask[list1_idx, list2_idx] <- 1
    write_mask(list2_list1_mask, mask_out_path)
    list2_list3_mask <- NaN_trial_by_trial_mtx
    list2_list3_mask[list2_idx, list3_idx] <- 1
    list2_list3_mask[list3_idx, list2_idx] <- 1
    write_mask(list2_list3_mask, mask_out_path)
    list2_list4_mask <- NaN_trial_by_trial_mtx
    list2_list4_mask[list2_idx, list4_idx] <- 1
    list2_list4_mask[list4_idx, list2_idx] <- 1
    write_mask(list2_list4_mask, mask_out_path)
    list2_list5_mask <- NaN_trial_by_trial_mtx
    list2_list5_mask[list2_idx, list5_idx] <- 1
    list2_list5_mask[list5_idx, list2_idx] <- 1
    write_mask(list2_list5_mask, mask_out_path)
    list2_list6_mask <- NaN_trial_by_trial_mtx
    list2_list6_mask[list2_idx, list6_idx] <- 1
    list2_list6_mask[list6_idx, list2_idx] <- 1
    write_mask(list2_list6_mask, mask_out_path)
    list2_list7_mask <- NaN_trial_by_trial_mtx
    list2_list7_mask[list2_idx, list7_idx] <- 1
    list2_list7_mask[list7_idx, list2_idx] <- 1
    write_mask(list2_list7_mask, mask_out_path)
    list2_list8_mask <- NaN_trial_by_trial_mtx
    list2_list8_mask[list2_idx, list8_idx] <- 1
    list2_list8_mask[list8_idx, list2_idx] <- 1
    write_mask(list2_list8_mask, mask_out_path)
    # --- list 3 ---
    list3_list1_mask <- NaN_trial_by_trial_mtx
    list3_list1_mask[list3_idx, list1_idx] <- 1
    list3_list1_mask[list1_idx, list3_idx] <- 1
    write_mask(list3_list1_mask, mask_out_path)
    list3_list2_mask <- NaN_trial_by_trial_mtx
    list3_list2_mask[list3_idx, list2_idx] <- 1
    list3_list2_mask[list2_idx, list3_idx] <- 1
    write_mask(list3_list2_mask, mask_out_path)
    list3_list4_mask <- NaN_trial_by_trial_mtx
    list3_list4_mask[list3_idx, list4_idx] <- 1
    list3_list4_mask[list4_idx, list3_idx] <- 1
    write_mask(list3_list4_mask, mask_out_path)
    list3_list5_mask <- NaN_trial_by_trial_mtx
    list3_list5_mask[list3_idx, list5_idx] <- 1
    list3_list5_mask[list5_idx, list3_idx] <- 1
    write_mask(list3_list5_mask, mask_out_path)
    list3_list6_mask <- NaN_trial_by_trial_mtx
    list3_list6_mask[list3_idx, list6_idx] <- 1
    list3_list6_mask[list6_idx, list3_idx] <- 1
    write_mask(list3_list6_mask, mask_out_path)
    list3_list7_mask <- NaN_trial_by_trial_mtx
    list3_list7_mask[list3_idx, list7_idx] <- 1
    list3_list7_mask[list7_idx, list3_idx] <- 1
    write_mask(list3_list7_mask, mask_out_path)
    list3_list8_mask <- NaN_trial_by_trial_mtx
    list3_list8_mask[list3_idx, list8_idx] <- 1
    list3_list8_mask[list8_idx, list3_idx] <- 1
    write_mask(list3_list8_mask, mask_out_path)
    # --- list 4 ---
    list4_list1_mask <- NaN_trial_by_trial_mtx
    list4_list1_mask[list4_idx, list1_idx] <- 1
    list4_list1_mask[list1_idx, list4_idx] <- 1
    write_mask(list4_list1_mask, mask_out_path)
    list4_list2_mask <- NaN_trial_by_trial_mtx
    list4_list2_mask[list4_idx, list2_idx] <- 1
    list4_list2_mask[list2_idx, list4_idx] <- 1
    write_mask(list4_list2_mask, mask_out_path)
    list4_list3_mask <- NaN_trial_by_trial_mtx
    list4_list3_mask[list4_idx, list3_idx] <- 1
    list4_list3_mask[list3_idx, list4_idx] <- 1
    write_mask(list4_list3_mask, mask_out_path)
    list4_list5_mask <- NaN_trial_by_trial_mtx
    list4_list5_mask[list4_idx, list5_idx] <- 1
    list4_list5_mask[list5_idx, list4_idx] <- 1
    write_mask(list4_list5_mask, mask_out_path)
    list4_list6_mask <- NaN_trial_by_trial_mtx
    list4_list6_mask[list4_idx, list6_idx] <- 1
    list4_list6_mask[list6_idx, list4_idx] <- 1
    write_mask(list4_list6_mask, mask_out_path)
    list4_list7_mask <- NaN_trial_by_trial_mtx
    list4_list7_mask[list4_idx, list7_idx] <- 1
    list4_list7_mask[list7_idx, list4_idx] <- 1
    write_mask(list4_list7_mask, mask_out_path)
    list4_list8_mask <- NaN_trial_by_trial_mtx
    list4_list8_mask[list4_idx, list8_idx] <- 1
    list4_list8_mask[list8_idx, list4_idx] <- 1
    write_mask(list4_list8_mask, mask_out_path)
    # --- list 5 ---
    list5_list1_mask <- NaN_trial_by_trial_mtx
    list5_list1_mask[list5_idx, list1_idx] <- 1
    list5_list1_mask[list1_idx, list5_idx] <- 1
    write_mask(list5_list1_mask, mask_out_path)
    list5_list2_mask <- NaN_trial_by_trial_mtx
    list5_list2_mask[list5_idx, list2_idx] <- 1
    list5_list2_mask[list2_idx, list5_idx] <- 1
    write_mask(list5_list2_mask, mask_out_path)
    list5_list3_mask <- NaN_trial_by_trial_mtx
    list5_list3_mask[list5_idx, list3_idx] <- 1
    list5_list3_mask[list3_idx, list5_idx] <- 1
    write_mask(list5_list3_mask, mask_out_path)
    list5_list4_mask <- NaN_trial_by_trial_mtx
    list5_list4_mask[list5_idx, list4_idx] <- 1
    list5_list4_mask[list4_idx, list5_idx] <- 1
    write_mask(list5_list4_mask, mask_out_path)
    list5_list6_mask <- NaN_trial_by_trial_mtx
    list5_list6_mask[list5_idx, list6_idx] <- 1
    list5_list6_mask[list6_idx, list5_idx] <- 1
    write_mask(list5_list6_mask, mask_out_path)
    list5_list7_mask <- NaN_trial_by_trial_mtx
    list5_list7_mask[list5_idx, list7_idx] <- 1
    list5_list7_mask[list7_idx, list5_idx] <- 1
    write_mask(list5_list7_mask, mask_out_path)
    list5_list8_mask <- NaN_trial_by_trial_mtx
    list5_list8_mask[list5_idx, list8_idx] <- 1
    list5_list8_mask[list8_idx, list5_idx] <- 1
    write_mask(list5_list8_mask, mask_out_path)
    # --- list 6 ---
    list6_list1_mask <- NaN_trial_by_trial_mtx
    list6_list1_mask[list6_idx, list1_idx] <- 1
    list6_list1_mask[list1_idx, list6_idx] <- 1
    write_mask(list6_list1_mask, mask_out_path)
    list6_list2_mask <- NaN_trial_by_trial_mtx
    list6_list2_mask[list6_idx, list2_idx] <- 1
    list6_list2_mask[list2_idx, list6_idx] <- 1
    write_mask(list6_list2_mask, mask_out_path)
    list6_list3_mask <- NaN_trial_by_trial_mtx
    list6_list3_mask[list6_idx, list3_idx] <- 1
    list6_list3_mask[list3_idx, list6_idx] <- 1
    write_mask(list6_list3_mask, mask_out_path)
    list6_list4_mask <- NaN_trial_by_trial_mtx
    list6_list4_mask[list6_idx, list4_idx] <- 1
    list6_list4_mask[list4_idx, list6_idx] <- 1
    write_mask(list6_list4_mask, mask_out_path)
    list6_list5_mask <- NaN_trial_by_trial_mtx
    list6_list5_mask[list6_idx, list5_idx] <- 1
    list6_list5_mask[list5_idx, list6_idx] <- 1
    write_mask(list6_list5_mask, mask_out_path)
    list6_list7_mask <- NaN_trial_by_trial_mtx
    list6_list7_mask[list6_idx, list7_idx] <- 1
    list6_list7_mask[list7_idx, list6_idx] <- 1
    write_mask(list6_list7_mask, mask_out_path)
    list6_list8_mask <- NaN_trial_by_trial_mtx
    list6_list8_mask[list6_idx, list8_idx] <- 1
    list6_list8_mask[list8_idx, list6_idx] <- 1
    write_mask(list6_list8_mask, mask_out_path)
    # --- list 7 ---
    list7_list1_mask <- NaN_trial_by_trial_mtx
    list7_list1_mask[list7_idx, list1_idx] <- 1
    list7_list1_mask[list1_idx, list7_idx] <- 1
    write_mask(list7_list1_mask, mask_out_path)
    list7_list2_mask <- NaN_trial_by_trial_mtx
    list7_list2_mask[list7_idx, list2_idx] <- 1
    list7_list2_mask[list2_idx, list7_idx] <- 1
    write_mask(list7_list2_mask, mask_out_path)
    list7_list3_mask <- NaN_trial_by_trial_mtx
    list7_list3_mask[list7_idx, list3_idx] <- 1
    list7_list3_mask[list3_idx, list7_idx] <- 1
    write_mask(list7_list3_mask, mask_out_path)
    list7_list4_mask <- NaN_trial_by_trial_mtx
    list7_list4_mask[list7_idx, list4_idx] <- 1
    list7_list4_mask[list4_idx, list7_idx] <- 1
    write_mask(list7_list4_mask, mask_out_path)
    list7_list5_mask <- NaN_trial_by_trial_mtx
    list7_list5_mask[list7_idx, list5_idx] <- 1
    list7_list5_mask[list5_idx, list7_idx] <- 1
    write_mask(list7_list5_mask, mask_out_path)
    list7_list6_mask <- NaN_trial_by_trial_mtx
    list7_list6_mask[list7_idx, list6_idx] <- 1
    list7_list6_mask[list6_idx, list7_idx] <- 1
    write_mask(list7_list6_mask, mask_out_path)
    list7_list8_mask <- NaN_trial_by_trial_mtx
    list7_list8_mask[list7_idx, list8_idx] <- 1
    list7_list8_mask[list8_idx, list7_idx] <- 1
    write_mask(list7_list8_mask, mask_out_path)
    # --- list 8 ---
    list8_list1_mask <- NaN_trial_by_trial_mtx
    list8_list1_mask[list8_idx, list1_idx] <- 1
    list8_list1_mask[list1_idx, list8_idx] <- 1
    write_mask(list8_list1_mask, mask_out_path)
    list8_list2_mask <- NaN_trial_by_trial_mtx
    list8_list2_mask[list8_idx, list2_idx] <- 1
    list8_list2_mask[list2_idx, list8_idx] <- 1
    write_mask(list8_list2_mask, mask_out_path)
    list8_list3_mask <- NaN_trial_by_trial_mtx
    list8_list3_mask[list8_idx, list3_idx] <- 1
    list8_list3_mask[list3_idx, list8_idx] <- 1
    write_mask(list8_list3_mask, mask_out_path)
    list8_list4_mask <- NaN_trial_by_trial_mtx
    list8_list4_mask[list8_idx, list4_idx] <- 1
    list8_list4_mask[list4_idx, list8_idx] <- 1
    write_mask(list8_list4_mask, mask_out_path)
    list8_list5_mask <- NaN_trial_by_trial_mtx
    list8_list5_mask[list8_idx, list5_idx] <- 1
    list8_list5_mask[list5_idx, list8_idx] <- 1
    write_mask(list8_list5_mask, mask_out_path)
    list8_list6_mask <- NaN_trial_by_trial_mtx
    list8_list6_mask[list8_idx, list6_idx] <- 1
    list8_list6_mask[list6_idx, list8_idx] <- 1
    write_mask(list8_list6_mask, mask_out_path)
    list8_list7_mask <- NaN_trial_by_trial_mtx
    list8_list7_mask[list8_idx, list7_idx] <- 1
    list8_list7_mask[list7_idx, list8_idx] <- 1
    write_mask(list8_list7_mask, mask_out_path)

    same_question_mask <- NaN_trial_by_trial_mtx
    same_question_mask[question1_idx, question1_idx] <- 1
    same_question_mask[question2_idx, question2_idx] <- 1
    same_question_mask[question3_idx, question3_idx] <- 1
    same_question_mask[question4_idx, question4_idx] <- 1
    write_mask(same_question_mask, mask_out_path)

    diff_question_mask <- NaN_trial_by_trial_mtx
    diff_question_mask[question1_idx, question2_idx] <- 1
    diff_question_mask[question1_idx, question3_idx] <- 1
    diff_question_mask[question1_idx, question4_idx] <- 1
    diff_question_mask[question2_idx, question1_idx] <- 1
    diff_question_mask[question2_idx, question3_idx] <- 1
    diff_question_mask[question2_idx, question4_idx] <- 1
    diff_question_mask[question3_idx, question1_idx] <- 1
    diff_question_mask[question3_idx, question2_idx] <- 1
    diff_question_mask[question3_idx, question4_idx] <- 1
    diff_question_mask[question4_idx, question1_idx] <- 1
    diff_question_mask[question4_idx, question2_idx] <- 1
    diff_question_mask[question4_idx, question3_idx] <- 1
    write_mask(diff_question_mask, mask_out_path)

    # sanity checking that `*` is doing what I think
    # fake_0_mtx <- matrix(data = 0, nrow = 5, ncol = 5)
    # fake_0_mtx[1,3] <- 1
    # fake_0_mtx[3,5] <- 0.5
    #
    # fake_1_mtx <- matrix(data = 1, nrow = 5, ncol = 5)
    #
    # fake_mult_mtx <- fake_0_mtx * fake_1_mtx
    # fake_mult_mtx
    # fake_mtxmult_mtx <- fake_0_mtx %*% fake_1_mtx
    # fake_mtxmult_mtx

    ###################
    ### Rhits #########
    ###################
    # --- rhits: split up by each list ---
    rHits_list1_sameQuestion_mask <- rhits_mask * list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_sameQuestion_mask, mask_out_path)

    rHits_list2_sameQuestion_mask <- rhits_mask * list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_sameQuestion_mask, mask_out_path)

    rHits_list3_sameQuestion_mask <- rhits_mask * list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_sameQuestion_mask, mask_out_path)

    rHits_list4_sameQuestion_mask <- rhits_mask * list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_sameQuestion_mask, mask_out_path)

    rHits_list5_sameQuestion_mask <- rhits_mask * list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_sameQuestion_mask, mask_out_path)

    rHits_list6_sameQuestion_mask <- rhits_mask * list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_sameQuestion_mask, mask_out_path)

    rHits_list7_sameQuestion_mask <- rhits_mask * list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_sameQuestion_mask, mask_out_path)

    rHits_list8_sameQuestion_mask <- rhits_mask * list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_sameQuestion_mask, mask_out_path)

    rHits_list1_diffQuestion_mask <- rhits_mask * list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_diffQuestion_mask, mask_out_path)

    rHits_list2_diffQuestion_mask <- rhits_mask * list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_diffQuestion_mask, mask_out_path)

    rHits_list3_diffQuestion_mask <- rhits_mask * list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_diffQuestion_mask, mask_out_path)

    rHits_list4_diffQuestion_mask <- rhits_mask * list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_diffQuestion_mask, mask_out_path)

    rHits_list5_diffQuestion_mask <- rhits_mask * list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_diffQuestion_mask, mask_out_path)

    rHits_list6_diffQuestion_mask <- rhits_mask * list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_diffQuestion_mask, mask_out_path)

    rHits_list7_diffQuestion_mask <- rhits_mask * list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_diffQuestion_mask, mask_out_path)

    rHits_list8_diffQuestion_mask <- rhits_mask * list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_diffQuestion_mask, mask_out_path)

    # --- rhits: split up by each list (between list/off-diagonal elements) ---
    # list1
    rHits_list1_list2_sameQuestion_mask <- rhits_mask * list1_list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list2_sameQuestion_mask, mask_out_path)
    rHits_list1_list3_sameQuestion_mask <- rhits_mask * list1_list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list3_sameQuestion_mask, mask_out_path)
    rHits_list1_list4_sameQuestion_mask <- rhits_mask * list1_list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list4_sameQuestion_mask, mask_out_path)
    rHits_list1_list5_sameQuestion_mask <- rhits_mask * list1_list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list5_sameQuestion_mask, mask_out_path)
    rHits_list1_list6_sameQuestion_mask <- rhits_mask * list1_list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list6_sameQuestion_mask, mask_out_path)
    rHits_list1_list7_sameQuestion_mask <- rhits_mask * list1_list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list7_sameQuestion_mask, mask_out_path)
    rHits_list1_list8_sameQuestion_mask <- rhits_mask * list1_list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list8_sameQuestion_mask, mask_out_path)
    rHits_list1_list2_diffQuestion_mask <- rhits_mask * list1_list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list2_diffQuestion_mask, mask_out_path)
    rHits_list1_list3_diffQuestion_mask <- rhits_mask * list1_list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list3_diffQuestion_mask, mask_out_path)
    rHits_list1_list4_diffQuestion_mask <- rhits_mask * list1_list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list4_diffQuestion_mask, mask_out_path)
    rHits_list1_list5_diffQuestion_mask <- rhits_mask * list1_list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list5_diffQuestion_mask, mask_out_path)
    rHits_list1_list6_diffQuestion_mask <- rhits_mask * list1_list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list6_diffQuestion_mask, mask_out_path)
    rHits_list1_list7_diffQuestion_mask <- rhits_mask * list1_list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list7_diffQuestion_mask, mask_out_path)
    rHits_list1_list8_diffQuestion_mask <- rhits_mask * list1_list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list1_list8_diffQuestion_mask, mask_out_path)
    # list2
    rHits_list2_list1_sameQuestion_mask <- rhits_mask * list2_list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list1_sameQuestion_mask, mask_out_path)
    rHits_list2_list3_sameQuestion_mask <- rhits_mask * list2_list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list3_sameQuestion_mask, mask_out_path)
    rHits_list2_list4_sameQuestion_mask <- rhits_mask * list2_list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list4_sameQuestion_mask, mask_out_path)
    rHits_list2_list5_sameQuestion_mask <- rhits_mask * list2_list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list5_sameQuestion_mask, mask_out_path)
    rHits_list2_list6_sameQuestion_mask <- rhits_mask * list2_list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list6_sameQuestion_mask, mask_out_path)
    rHits_list2_list7_sameQuestion_mask <- rhits_mask * list2_list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list7_sameQuestion_mask, mask_out_path)
    rHits_list2_list8_sameQuestion_mask <- rhits_mask * list2_list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list8_sameQuestion_mask, mask_out_path)
    rHits_list2_list1_diffQuestion_mask <- rhits_mask * list2_list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list1_diffQuestion_mask, mask_out_path)
    rHits_list2_list3_diffQuestion_mask <- rhits_mask * list2_list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list3_diffQuestion_mask, mask_out_path)
    rHits_list2_list4_diffQuestion_mask <- rhits_mask * list2_list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list4_diffQuestion_mask, mask_out_path)
    rHits_list2_list5_diffQuestion_mask <- rhits_mask * list2_list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list5_diffQuestion_mask, mask_out_path)
    rHits_list2_list6_diffQuestion_mask <- rhits_mask * list2_list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list6_diffQuestion_mask, mask_out_path)
    rHits_list2_list7_diffQuestion_mask <- rhits_mask * list2_list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list7_diffQuestion_mask, mask_out_path)
    rHits_list2_list8_diffQuestion_mask <- rhits_mask * list2_list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list2_list8_diffQuestion_mask, mask_out_path)
    # list3
    rHits_list3_list1_sameQuestion_mask <- rhits_mask * list3_list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list1_sameQuestion_mask, mask_out_path)
    rHits_list3_list2_sameQuestion_mask <- rhits_mask * list3_list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list2_sameQuestion_mask, mask_out_path)
    rHits_list3_list4_sameQuestion_mask <- rhits_mask * list3_list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list4_sameQuestion_mask, mask_out_path)
    rHits_list3_list5_sameQuestion_mask <- rhits_mask * list3_list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list5_sameQuestion_mask, mask_out_path)
    rHits_list3_list6_sameQuestion_mask <- rhits_mask * list3_list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list6_sameQuestion_mask, mask_out_path)
    rHits_list3_list7_sameQuestion_mask <- rhits_mask * list3_list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list7_sameQuestion_mask, mask_out_path)
    rHits_list3_list8_sameQuestion_mask <- rhits_mask * list3_list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list8_sameQuestion_mask, mask_out_path)
    rHits_list3_list1_diffQuestion_mask <- rhits_mask * list3_list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list1_diffQuestion_mask, mask_out_path)
    rHits_list3_list2_diffQuestion_mask <- rhits_mask * list3_list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list2_diffQuestion_mask, mask_out_path)
    rHits_list3_list4_diffQuestion_mask <- rhits_mask * list3_list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list4_diffQuestion_mask, mask_out_path)
    rHits_list3_list5_diffQuestion_mask <- rhits_mask * list3_list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list5_diffQuestion_mask, mask_out_path)
    rHits_list3_list6_diffQuestion_mask <- rhits_mask * list3_list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list6_diffQuestion_mask, mask_out_path)
    rHits_list3_list7_diffQuestion_mask <- rhits_mask * list3_list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list7_diffQuestion_mask, mask_out_path)
    rHits_list3_list8_diffQuestion_mask <- rhits_mask * list3_list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list3_list8_diffQuestion_mask, mask_out_path)
    # list4
    rHits_list4_list1_sameQuestion_mask <- rhits_mask * list4_list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list1_sameQuestion_mask, mask_out_path)
    rHits_list4_list2_sameQuestion_mask <- rhits_mask * list4_list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list2_sameQuestion_mask, mask_out_path)
    rHits_list4_list3_sameQuestion_mask <- rhits_mask * list4_list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list3_sameQuestion_mask, mask_out_path)
    rHits_list4_list5_sameQuestion_mask <- rhits_mask * list4_list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list5_sameQuestion_mask, mask_out_path)
    rHits_list4_list6_sameQuestion_mask <- rhits_mask * list4_list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list6_sameQuestion_mask, mask_out_path)
    rHits_list4_list7_sameQuestion_mask <- rhits_mask * list4_list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list7_sameQuestion_mask, mask_out_path)
    rHits_list4_list8_sameQuestion_mask <- rhits_mask * list4_list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list8_sameQuestion_mask, mask_out_path)
    rHits_list4_list1_diffQuestion_mask <- rhits_mask * list4_list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list1_diffQuestion_mask, mask_out_path)
    rHits_list4_list2_diffQuestion_mask <- rhits_mask * list4_list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list2_diffQuestion_mask, mask_out_path)
    rHits_list4_list3_diffQuestion_mask <- rhits_mask * list4_list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list3_diffQuestion_mask, mask_out_path)
    rHits_list4_list5_diffQuestion_mask <- rhits_mask * list4_list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list5_diffQuestion_mask, mask_out_path)
    rHits_list4_list6_diffQuestion_mask <- rhits_mask * list4_list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list6_diffQuestion_mask, mask_out_path)
    rHits_list4_list7_diffQuestion_mask <- rhits_mask * list4_list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list7_diffQuestion_mask, mask_out_path)
    rHits_list4_list8_diffQuestion_mask <- rhits_mask * list4_list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list4_list8_diffQuestion_mask, mask_out_path)
    # list5
    rHits_list5_list1_sameQuestion_mask <- rhits_mask * list5_list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list1_sameQuestion_mask, mask_out_path)
    rHits_list5_list2_sameQuestion_mask <- rhits_mask * list5_list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list2_sameQuestion_mask, mask_out_path)
    rHits_list5_list3_sameQuestion_mask <- rhits_mask * list5_list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list3_sameQuestion_mask, mask_out_path)
    rHits_list5_list4_sameQuestion_mask <- rhits_mask * list5_list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list4_sameQuestion_mask, mask_out_path)
    rHits_list5_list6_sameQuestion_mask <- rhits_mask * list5_list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list6_sameQuestion_mask, mask_out_path)
    rHits_list5_list7_sameQuestion_mask <- rhits_mask * list5_list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list7_sameQuestion_mask, mask_out_path)
    rHits_list5_list8_sameQuestion_mask <- rhits_mask * list5_list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list8_sameQuestion_mask, mask_out_path)
    rHits_list5_list1_diffQuestion_mask <- rhits_mask * list5_list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list1_diffQuestion_mask, mask_out_path)
    rHits_list5_list2_diffQuestion_mask <- rhits_mask * list5_list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list2_diffQuestion_mask, mask_out_path)
    rHits_list5_list3_diffQuestion_mask <- rhits_mask * list5_list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list3_diffQuestion_mask, mask_out_path)
    rHits_list5_list4_diffQuestion_mask <- rhits_mask * list5_list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list4_diffQuestion_mask, mask_out_path)
    rHits_list5_list6_diffQuestion_mask <- rhits_mask * list5_list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list6_diffQuestion_mask, mask_out_path)
    rHits_list5_list7_diffQuestion_mask <- rhits_mask * list5_list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list7_diffQuestion_mask, mask_out_path)
    rHits_list5_list8_diffQuestion_mask <- rhits_mask * list5_list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list5_list8_diffQuestion_mask, mask_out_path)
    # list6
    rHits_list6_list1_sameQuestion_mask <- rhits_mask * list6_list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list1_sameQuestion_mask, mask_out_path)
    rHits_list6_list2_sameQuestion_mask <- rhits_mask * list6_list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list2_sameQuestion_mask, mask_out_path)
    rHits_list6_list3_sameQuestion_mask <- rhits_mask * list6_list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list3_sameQuestion_mask, mask_out_path)
    rHits_list6_list4_sameQuestion_mask <- rhits_mask * list6_list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list4_sameQuestion_mask, mask_out_path)
    rHits_list6_list5_sameQuestion_mask <- rhits_mask * list6_list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list5_sameQuestion_mask, mask_out_path)
    rHits_list6_list7_sameQuestion_mask <- rhits_mask * list6_list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list7_sameQuestion_mask, mask_out_path)
    rHits_list6_list8_sameQuestion_mask <- rhits_mask * list6_list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list8_sameQuestion_mask, mask_out_path)
    rHits_list6_list1_diffQuestion_mask <- rhits_mask * list6_list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list1_diffQuestion_mask, mask_out_path)
    rHits_list6_list2_diffQuestion_mask <- rhits_mask * list6_list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list2_diffQuestion_mask, mask_out_path)
    rHits_list6_list3_diffQuestion_mask <- rhits_mask * list6_list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list3_diffQuestion_mask, mask_out_path)
    rHits_list6_list4_diffQuestion_mask <- rhits_mask * list6_list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list4_diffQuestion_mask, mask_out_path)
    rHits_list6_list5_diffQuestion_mask <- rhits_mask * list6_list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list5_diffQuestion_mask, mask_out_path)
    rHits_list6_list7_diffQuestion_mask <- rhits_mask * list6_list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list7_diffQuestion_mask, mask_out_path)
    rHits_list6_list8_diffQuestion_mask <- rhits_mask * list6_list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list6_list8_diffQuestion_mask, mask_out_path)
    # list7
    rHits_list7_list1_sameQuestion_mask <- rhits_mask * list7_list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list1_sameQuestion_mask, mask_out_path)
    rHits_list7_list2_sameQuestion_mask <- rhits_mask * list7_list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list2_sameQuestion_mask, mask_out_path)
    rHits_list7_list3_sameQuestion_mask <- rhits_mask * list7_list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list3_sameQuestion_mask, mask_out_path)
    rHits_list7_list4_sameQuestion_mask <- rhits_mask * list7_list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list4_sameQuestion_mask, mask_out_path)
    rHits_list7_list5_sameQuestion_mask <- rhits_mask * list7_list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list5_sameQuestion_mask, mask_out_path)
    rHits_list7_list6_sameQuestion_mask <- rhits_mask * list7_list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list6_sameQuestion_mask, mask_out_path)
    rHits_list7_list8_sameQuestion_mask <- rhits_mask * list7_list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list8_sameQuestion_mask, mask_out_path)
    rHits_list7_list1_diffQuestion_mask <- rhits_mask * list7_list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list1_diffQuestion_mask, mask_out_path)
    rHits_list7_list2_diffQuestion_mask <- rhits_mask * list7_list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list2_diffQuestion_mask, mask_out_path)
    rHits_list7_list3_diffQuestion_mask <- rhits_mask * list7_list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list3_diffQuestion_mask, mask_out_path)
    rHits_list7_list4_diffQuestion_mask <- rhits_mask * list7_list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list4_diffQuestion_mask, mask_out_path)
    rHits_list7_list5_diffQuestion_mask <- rhits_mask * list7_list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list5_diffQuestion_mask, mask_out_path)
    rHits_list7_list6_diffQuestion_mask <- rhits_mask * list7_list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list6_diffQuestion_mask, mask_out_path)
    rHits_list7_list8_diffQuestion_mask <- rhits_mask * list7_list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list7_list8_diffQuestion_mask, mask_out_path)
    # list8
    rHits_list8_list1_sameQuestion_mask <- rhits_mask * list8_list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list1_sameQuestion_mask, mask_out_path)
    rHits_list8_list2_sameQuestion_mask <- rhits_mask * list8_list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list2_sameQuestion_mask, mask_out_path)
    rHits_list8_list3_sameQuestion_mask <- rhits_mask * list8_list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list3_sameQuestion_mask, mask_out_path)
    rHits_list8_list4_sameQuestion_mask <- rhits_mask * list8_list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list4_sameQuestion_mask, mask_out_path)
    rHits_list8_list5_sameQuestion_mask <- rhits_mask * list8_list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list5_sameQuestion_mask, mask_out_path)
    rHits_list8_list6_sameQuestion_mask <- rhits_mask * list8_list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list6_sameQuestion_mask, mask_out_path)
    rHits_list8_list7_sameQuestion_mask <- rhits_mask * list8_list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list7_sameQuestion_mask, mask_out_path)
    rHits_list8_list1_diffQuestion_mask <- rhits_mask * list8_list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list1_diffQuestion_mask, mask_out_path)
    rHits_list8_list2_diffQuestion_mask <- rhits_mask * list8_list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list2_diffQuestion_mask, mask_out_path)
    rHits_list8_list3_diffQuestion_mask <- rhits_mask * list8_list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list3_diffQuestion_mask, mask_out_path)
    rHits_list8_list4_diffQuestion_mask <- rhits_mask * list8_list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list4_diffQuestion_mask, mask_out_path)
    rHits_list8_list5_diffQuestion_mask <- rhits_mask * list8_list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list5_diffQuestion_mask, mask_out_path)
    rHits_list8_list6_diffQuestion_mask <- rhits_mask * list8_list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list6_diffQuestion_mask, mask_out_path)
    rHits_list8_list7_diffQuestion_mask <- rhits_mask * list8_list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_list8_list7_diffQuestion_mask, mask_out_path)

    # --- rhits and source hits: split up by each list ---
    rANDsourceHits_list1_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list1_sameQuestion_mask, mask_out_path)

    rANDsourceHits_list2_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list2_sameQuestion_mask, mask_out_path)

    rANDsourceHits_list3_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list3_sameQuestion_mask, mask_out_path)

    rANDsourceHits_list4_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list4_sameQuestion_mask, mask_out_path)

    rANDsourceHits_list5_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list5_sameQuestion_mask, mask_out_path)

    rANDsourceHits_list6_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list6_sameQuestion_mask, mask_out_path)

    rANDsourceHits_list7_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list7_sameQuestion_mask, mask_out_path)

    rANDsourceHits_list8_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list8_sameQuestion_mask, mask_out_path)

    rANDsourceHits_list1_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list1_diffQuestion_mask, mask_out_path)

    rANDsourceHits_list2_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list2_diffQuestion_mask, mask_out_path)

    rANDsourceHits_list3_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list3_diffQuestion_mask, mask_out_path)

    rANDsourceHits_list4_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list4_diffQuestion_mask, mask_out_path)

    rANDsourceHits_list5_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list5_diffQuestion_mask, mask_out_path)

    rANDsourceHits_list6_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list6_diffQuestion_mask, mask_out_path)

    rANDsourceHits_list7_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list7_diffQuestion_mask, mask_out_path)

    rANDsourceHits_list8_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_list8_diffQuestion_mask, mask_out_path)

    # --- rhits: split up by same/different list and question ---
    rHits_sameList_sameQuestion_mask <- rhits_mask * same_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_sameList_sameQuestion_mask, mask_out_path)

    rHits_diffList_sameQuestion_mask <- rhits_mask * diff_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rHits_diffList_sameQuestion_mask, mask_out_path)

    rHits_sameList_diffQuestion_mask <- rhits_mask * same_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_sameList_diffQuestion_mask, mask_out_path)

    rHits_diffList_diffQuestion_mask <- rhits_mask * diff_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rHits_diffList_diffQuestion_mask, mask_out_path)

    # --- rhits and source hits: split up by same/different list and question ---
    rANDsourceHits_sameList_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * same_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_sameList_sameQuestion_mask, mask_out_path)

    rANDsourceHits_diffList_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * diff_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_diffList_sameQuestion_mask, mask_out_path)

    rANDsourceHits_sameList_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * same_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_sameList_diffQuestion_mask, mask_out_path)

    rANDsourceHits_diffList_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * diff_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_diffList_diffQuestion_mask, mask_out_path)

    # --- rhits and source hits: split up by same/diff question (collapsing list) ---
    rANDsourceHits_sameQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_sameQuestion_mask, mask_out_path)

    rANDsourceHits_diffQuestion_mask <- rhits_mask * temporal_hits_mask * question_exact_hits_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(rANDsourceHits_diffQuestion_mask, mask_out_path)

    # --- rhits and source hits: do not require exact memory for source question - count as correct if get fit vs find category correct
    rHits_questionHitsLiberal_mask <- rhits_mask * question_liberal_hits_mask * within_run_corr_as_NaN
    write_mask(rHits_questionHitsLiberal_mask, mask_out_path)

    rHits_questionMissLiberal_mask <- rhits_mask * question_liberal_miss_mask * within_run_corr_as_NaN
    write_mask(rHits_questionMissLiberal_mask, mask_out_path)

    ###################
    ### ITEM MISSES ###
    ###################
    # --- item miss: split up by each list ---
    itemMiss_list1_sameQuestion_mask <- item_miss_mask * list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list1_sameQuestion_mask, mask_out_path)

    itemMiss_list2_sameQuestion_mask <- item_miss_mask * list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list2_sameQuestion_mask, mask_out_path)

    itemMiss_list3_sameQuestion_mask <- item_miss_mask * list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list3_sameQuestion_mask, mask_out_path)

    itemMiss_list4_sameQuestion_mask <- item_miss_mask * list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list4_sameQuestion_mask, mask_out_path)

    itemMiss_list5_sameQuestion_mask <- item_miss_mask * list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list5_sameQuestion_mask, mask_out_path)

    itemMiss_list6_sameQuestion_mask <- item_miss_mask * list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list6_sameQuestion_mask, mask_out_path)

    itemMiss_list7_sameQuestion_mask <- item_miss_mask * list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list7_sameQuestion_mask, mask_out_path)

    itemMiss_list8_sameQuestion_mask <- item_miss_mask * list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list8_sameQuestion_mask, mask_out_path)

    itemMiss_list1_diffQuestion_mask <- item_miss_mask * list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list1_diffQuestion_mask, mask_out_path)

    itemMiss_list2_diffQuestion_mask <- item_miss_mask * list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list2_diffQuestion_mask, mask_out_path)

    itemMiss_list3_diffQuestion_mask <- item_miss_mask * list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list3_diffQuestion_mask, mask_out_path)

    itemMiss_list4_diffQuestion_mask <- item_miss_mask * list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list4_diffQuestion_mask, mask_out_path)

    itemMiss_list5_diffQuestion_mask <- item_miss_mask * list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list5_diffQuestion_mask, mask_out_path)

    itemMiss_list6_diffQuestion_mask <- item_miss_mask * list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list6_diffQuestion_mask, mask_out_path)

    itemMiss_list7_diffQuestion_mask <- item_miss_mask * list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list7_diffQuestion_mask, mask_out_path)

    itemMiss_list8_diffQuestion_mask <- item_miss_mask * list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_list8_diffQuestion_mask, mask_out_path)

    # --- item miss and source miss: split up by each list ---
    itemMissANDsourceMiss_list1_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list1_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list2_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list2_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list3_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list3_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list4_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list4_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list5_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list5_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list6_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list6_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list7_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list7_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list8_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list8_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list1_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list1_diffQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list2_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list2_diffQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list3_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list3_diffQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list4_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list4_diffQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list5_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list5_diffQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list6_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list6_diffQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list7_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list7_diffQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_list8_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_list8_diffQuestion_mask, mask_out_path)

    # --- item miss: split up by same/different list and question ---
    itemMiss_sameList_sameQuestion_mask <- item_miss_mask * same_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_sameList_sameQuestion_mask, mask_out_path)

    itemMiss_diffList_sameQuestion_mask <- item_miss_mask * diff_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_diffList_sameQuestion_mask, mask_out_path)

    itemMiss_sameList_diffQuestion_mask <- item_miss_mask * same_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_sameList_diffQuestion_mask, mask_out_path)

    itemMiss_diffList_diffQuestion_mask <- item_miss_mask * diff_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMiss_diffList_diffQuestion_mask, mask_out_path)

    # --- item miss and source miss: split up by same/different list and question ---
    itemMissANDsourceMiss_sameList_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * same_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_sameList_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_diffList_sameQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * diff_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_diffList_sameQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_sameList_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * same_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_sameList_diffQuestion_mask, mask_out_path)

    itemMissANDsourceMiss_diffList_diffQuestion_mask <- item_miss_mask * temporal_miss_mask * question_exact_miss_mask * diff_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(itemMissANDsourceMiss_diffList_diffQuestion_mask, mask_out_path)

    ###########################
    ### FHITS & ITEM MISSES ###
    ###########################
    # --- fhits and misses: split up by each list ---
    fHitsANDmisses_list1_sameQuestion_mask <- fhitsANDmisses_mask * list1_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list1_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_list2_sameQuestion_mask <- fhitsANDmisses_mask * list2_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list2_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_list3_sameQuestion_mask <- fhitsANDmisses_mask * list3_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list3_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_list4_sameQuestion_mask <- fhitsANDmisses_mask * list4_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list4_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_list5_sameQuestion_mask <- fhitsANDmisses_mask * list5_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list5_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_list6_sameQuestion_mask <- fhitsANDmisses_mask * list6_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list6_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_list7_sameQuestion_mask <- fhitsANDmisses_mask * list7_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list7_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_list8_sameQuestion_mask <- fhitsANDmisses_mask * list8_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list8_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_list1_diffQuestion_mask <- fhitsANDmisses_mask * list1_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list1_diffQuestion_mask, mask_out_path)

    fHitsANDmisses_list2_diffQuestion_mask <- fhitsANDmisses_mask * list2_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list2_diffQuestion_mask, mask_out_path)

    fHitsANDmisses_list3_diffQuestion_mask <- fhitsANDmisses_mask * list3_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list3_diffQuestion_mask, mask_out_path)

    fHitsANDmisses_list4_diffQuestion_mask <- fhitsANDmisses_mask * list4_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list4_diffQuestion_mask, mask_out_path)

    fHitsANDmisses_list5_diffQuestion_mask <- fhitsANDmisses_mask * list5_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list5_diffQuestion_mask, mask_out_path)

    fHitsANDmisses_list6_diffQuestion_mask <- fhitsANDmisses_mask * list6_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list6_diffQuestion_mask, mask_out_path)

    fHitsANDmisses_list7_diffQuestion_mask <- fhitsANDmisses_mask * list7_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list7_diffQuestion_mask, mask_out_path)

    fHitsANDmisses_list8_diffQuestion_mask <- fhitsANDmisses_mask * list8_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_list8_diffQuestion_mask, mask_out_path)

    # --- fhits: split up by same/different list and question ---
    fHitsANDmisses_sameList_sameQuestion_mask <- fhitsANDmisses_mask * same_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_sameList_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_diffList_sameQuestion_mask <- fhitsANDmisses_mask * diff_list_mask * same_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_diffList_sameQuestion_mask, mask_out_path)

    fHitsANDmisses_sameList_diffQuestion_mask <- fhitsANDmisses_mask * same_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_sameList_diffQuestion_mask, mask_out_path)

    fHitsANDmisses_diffList_diffQuestion_mask <- fhitsANDmisses_mask * diff_list_mask * diff_question_mask * within_run_corr_as_NaN
    write_mask(fHitsANDmisses_diffList_diffQuestion_mask, mask_out_path)

  } #if(exists(cur_subj_info_file
  print(sprintf("%s took %.2f minutes", isubj, Sys.time() - time_subj_start))
} #for (isubj


#' # Print out R session info
sessionInfo()

#' # Print out script timer
sprintf("It took %.2f minutes to run all of this.", Sys.time() - time_script_start)
