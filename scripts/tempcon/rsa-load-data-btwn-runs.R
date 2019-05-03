#' ---
#' title: Load TempCon data &mdash; Between runs
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

#' ## Function definitions
gather_df <- function(input_df){
  input_df %>%
    tidyr::gather(col_name, r, -subj, -roi, -hemi, -question_id, -item_mem, -question_mem, -list_mem, -question, -list, -source_mem, -condition, -row_name) %>%
    dplyr::filter(!is.na(r))
}

filter_df <- function(input_df, info_df){
  input_df %>%
    dplyr::mutate(row_recog_trial_id = as.numeric(sub("x.", "", row_name)),
                  col_recog_trial_id = as.numeric(sub("x.", "", col_name))) %>%
    # figure out encoding trial ID for rows
    dplyr::left_join(info_df, by = "row_recog_trial_id") %>%
    dplyr::rename(row_enc_trial_id = enc_trial_id,
                  row_objnum = objnum,
                  row_enctask = enctask,
                  # b/c this column exists in both, remove the `x` that got added in the prior join
                  col_recog_trial_id = col_recog_trial_id.x,
                  subj = subj.x) %>%
    #figure out encoding trial ID for columns
    dplyr::left_join(info_df, by = "col_recog_trial_id") %>%
    dplyr::rename(col_enc_trial_id = enc_trial_id,
                  col_objnum = objnum,
                  col_enctask = enctask,
                  subj = subj.x,
                  row_recog_trial_id = row_recog_trial_id.x) %>%
    dplyr::select(-ends_with(".x"), -ends_with(".y"))
}

read_symmetric_mtx <- function(input_name, input_path){
  mask_file <- file.path(input_path,sprintf("%s.txt", input_name))
  if(!file.exists(mask_file)){
    sprintf("Error - mask file %s does not exist.", mask_file)
  }
  foo <- as.matrix(read.table(mask_file), as.is = TRUE)
  colnames(foo) <- NULL
  stopifnot(isSymmetric(foo))
  return(foo)
}

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

#+ label = "subjects"
#' ## Figure out subjects
# assume that subjects have folders in `raw_behavioral_dir` that start with `s` and are followed by one or two digits
subjects <- c(list.files(path=raw_behavioral_dir,pattern="^[s][(123456789)]$|^[s][(123456789)][(0123456789)]$"))
# chuck out any subjects who need to be excluded
exclude_subjects <- c("s1", "s2", "s3", "s7", "s11", "s13", "s14", "s24", "s26")
subjects <- subjects[!is.element(subjects, exclude_subjects)]
length(subjects)
subjects

cat(sprintf("\nPerforming initial setup\n"))

#+ label = "setup vars"
#' ## Set other variables
roi_dirs <- c("ashs_left","ashs_right", "manual_traced_rois")
rois <- c("CA1.nii","CA2.nii","CA3.nii","DG.nii","ERC.nii","subiculum.nii",
          "35_36.nii","CA2_3_DG.nii","CA3_DG.nii", "whole_hippo.nii",
          "CA1_head.nii", "CA1_body.nii", "CA1_tail.nii",
          "CA2_3_DG_head.nii", "CA2_3_DG_body.nii", "CA2_3_DG_tail.nii",
          "whole_hippo_head.nii", "whole_hippo_body.nii", "whole_hippo_tail.nii")

rois_of_most_interest <- c("subiculum_body.nii",
                           "CA1_body.nii", "CA2_3_DG_body.nii",
                           "PRC_L.nii", "PRC_R.nii",
                           "aPHC_L.nii", "aPHC_R.nii", "pPHC_L.nii", "pPHC_R.nii")

runs <- c(1:6)
num_trials_per_run <- 60
# read in mask for w/in run correlations
# this will later be used to mask out w/in run correlations (should already be masked out in trial masks, but can never be too safe!)
# assumes have 6 runs w/ 60 trials in each
within_run_corr_as_NaN <- read_symmetric_mtx("within_run_corr_as_NaN",analyzed_behavioral_dir)

#' ## Flags
REDUCED_ROIS_FLAG <- 1
MIXED_MODELS_FLAG <- 1
MAP_TYPE_FLAG <- "from_betas"

if(REDUCED_ROIS_FLAG==1){
  rois <- rois_of_most_interest
}

#+ label="read in data"
#+ warning = FALSE
#' # Read in data: Loop across subjects, hemispheres, rois, and runs
for(isubj in subjects){
  time_subj_start <- Sys.time()
  cat(sprintf("\nWorking on %s\n", isubj))
  # --- read in masks ---
  # these are created in `rsa-generate-masks-btwn-runs.R`
  # to check: library(superheat)
  # superheat::superheat(<MASK-NAME>)
  mask_path <- paste0(halle::ensure_trailing_slash(analyzed_mri_dir),
                          halle::ensure_trailing_slash("multivariate_itemHitBYType"),
                          halle::ensure_trailing_slash(isubj))

  #\\TODO: make this into a function. problem w/ dynamic variable name assignment
  # an alternative option would be to write out as `.RData` instead of `.txt`
  # --- encoding question regardless of memory ---
  question1_mask <- read_symmetric_mtx("question1_mask", mask_path)
  question2_mask <- read_symmetric_mtx("question2_mask", mask_path)
  question3_mask <- read_symmetric_mtx("question3_mask", mask_path)
  question4_mask <- read_symmetric_mtx("question4_mask", mask_path)
  # --- rhits: split up by each list ---
  rHits_list1_sameQuestion_mask <- read_symmetric_mtx("rHits_list1_sameQuestion_mask", mask_path)
  rHits_list2_sameQuestion_mask <- read_symmetric_mtx("rHits_list2_sameQuestion_mask", mask_path)
  rHits_list3_sameQuestion_mask <- read_symmetric_mtx("rHits_list3_sameQuestion_mask", mask_path)
  rHits_list4_sameQuestion_mask <- read_symmetric_mtx("rHits_list4_sameQuestion_mask", mask_path)
  rHits_list5_sameQuestion_mask <- read_symmetric_mtx("rHits_list5_sameQuestion_mask", mask_path)
  rHits_list6_sameQuestion_mask <- read_symmetric_mtx("rHits_list6_sameQuestion_mask", mask_path)
  rHits_list7_sameQuestion_mask <- read_symmetric_mtx("rHits_list7_sameQuestion_mask", mask_path)
  rHits_list8_sameQuestion_mask <- read_symmetric_mtx("rHits_list8_sameQuestion_mask", mask_path)
  rHits_list1_diffQuestion_mask <- read_symmetric_mtx("rHits_list1_diffQuestion_mask", mask_path)
  rHits_list2_diffQuestion_mask <- read_symmetric_mtx("rHits_list2_diffQuestion_mask", mask_path)
  rHits_list3_diffQuestion_mask <- read_symmetric_mtx("rHits_list3_diffQuestion_mask", mask_path)
  rHits_list4_diffQuestion_mask <- read_symmetric_mtx("rHits_list4_diffQuestion_mask", mask_path)
  rHits_list5_diffQuestion_mask <- read_symmetric_mtx("rHits_list5_diffQuestion_mask", mask_path)
  rHits_list6_diffQuestion_mask <- read_symmetric_mtx("rHits_list6_diffQuestion_mask", mask_path)
  rHits_list7_diffQuestion_mask <- read_symmetric_mtx("rHits_list7_diffQuestion_mask", mask_path)
  rHits_list8_diffQuestion_mask <- read_symmetric_mtx("rHits_list8_diffQuestion_mask", mask_path)
  # --- rhits: split up by each list, off-diagonal elements ---
  # list1
  rHits_list1_list2_sameQuestion_mask <- read_symmetric_mtx("rHits_list1_list2_sameQuestion_mask", mask_path)
  rHits_list1_list3_sameQuestion_mask <- read_symmetric_mtx("rHits_list1_list3_sameQuestion_mask", mask_path)
  rHits_list1_list4_sameQuestion_mask <- read_symmetric_mtx("rHits_list1_list4_sameQuestion_mask", mask_path)
  rHits_list1_list5_sameQuestion_mask <- read_symmetric_mtx("rHits_list1_list5_sameQuestion_mask", mask_path)
  rHits_list1_list6_sameQuestion_mask <- read_symmetric_mtx("rHits_list1_list6_sameQuestion_mask", mask_path)
  rHits_list1_list7_sameQuestion_mask <- read_symmetric_mtx("rHits_list1_list7_sameQuestion_mask", mask_path)
  rHits_list1_list8_sameQuestion_mask <- read_symmetric_mtx("rHits_list1_list8_sameQuestion_mask", mask_path)
  rHits_list1_list2_diffQuestion_mask <- read_symmetric_mtx("rHits_list1_list2_diffQuestion_mask", mask_path)
  rHits_list1_list3_diffQuestion_mask <- read_symmetric_mtx("rHits_list1_list3_diffQuestion_mask", mask_path)
  rHits_list1_list4_diffQuestion_mask <- read_symmetric_mtx("rHits_list1_list4_diffQuestion_mask", mask_path)
  rHits_list1_list5_diffQuestion_mask <- read_symmetric_mtx("rHits_list1_list5_diffQuestion_mask", mask_path)
  rHits_list1_list6_diffQuestion_mask <- read_symmetric_mtx("rHits_list1_list6_diffQuestion_mask", mask_path)
  rHits_list1_list7_diffQuestion_mask <- read_symmetric_mtx("rHits_list1_list7_diffQuestion_mask", mask_path)
  rHits_list1_list8_diffQuestion_mask <- read_symmetric_mtx("rHits_list1_list8_diffQuestion_mask", mask_path)
  # list2
  rHits_list2_list1_sameQuestion_mask <- read_symmetric_mtx("rHits_list2_list1_sameQuestion_mask", mask_path)
  rHits_list2_list3_sameQuestion_mask <- read_symmetric_mtx("rHits_list2_list3_sameQuestion_mask", mask_path)
  rHits_list2_list4_sameQuestion_mask <- read_symmetric_mtx("rHits_list2_list4_sameQuestion_mask", mask_path)
  rHits_list2_list5_sameQuestion_mask <- read_symmetric_mtx("rHits_list2_list5_sameQuestion_mask", mask_path)
  rHits_list2_list6_sameQuestion_mask <- read_symmetric_mtx("rHits_list2_list6_sameQuestion_mask", mask_path)
  rHits_list2_list7_sameQuestion_mask <- read_symmetric_mtx("rHits_list2_list7_sameQuestion_mask", mask_path)
  rHits_list2_list8_sameQuestion_mask <- read_symmetric_mtx("rHits_list2_list8_sameQuestion_mask", mask_path)
  rHits_list2_list1_diffQuestion_mask <- read_symmetric_mtx("rHits_list2_list1_diffQuestion_mask", mask_path)
  rHits_list2_list3_diffQuestion_mask <- read_symmetric_mtx("rHits_list2_list3_diffQuestion_mask", mask_path)
  rHits_list2_list4_diffQuestion_mask <- read_symmetric_mtx("rHits_list2_list4_diffQuestion_mask", mask_path)
  rHits_list2_list5_diffQuestion_mask <- read_symmetric_mtx("rHits_list2_list5_diffQuestion_mask", mask_path)
  rHits_list2_list6_diffQuestion_mask <- read_symmetric_mtx("rHits_list2_list6_diffQuestion_mask", mask_path)
  rHits_list2_list7_diffQuestion_mask <- read_symmetric_mtx("rHits_list2_list7_diffQuestion_mask", mask_path)
  rHits_list2_list8_diffQuestion_mask <- read_symmetric_mtx("rHits_list2_list8_diffQuestion_mask", mask_path)
  # list3
  rHits_list3_list1_sameQuestion_mask <- read_symmetric_mtx("rHits_list3_list1_sameQuestion_mask", mask_path)
  rHits_list3_list2_sameQuestion_mask <- read_symmetric_mtx("rHits_list3_list2_sameQuestion_mask", mask_path)
  rHits_list3_list4_sameQuestion_mask <- read_symmetric_mtx("rHits_list3_list4_sameQuestion_mask", mask_path)
  rHits_list3_list5_sameQuestion_mask <- read_symmetric_mtx("rHits_list3_list5_sameQuestion_mask", mask_path)
  rHits_list3_list6_sameQuestion_mask <- read_symmetric_mtx("rHits_list3_list6_sameQuestion_mask", mask_path)
  rHits_list3_list7_sameQuestion_mask <- read_symmetric_mtx("rHits_list3_list7_sameQuestion_mask", mask_path)
  rHits_list3_list8_sameQuestion_mask <- read_symmetric_mtx("rHits_list3_list8_sameQuestion_mask", mask_path)
  rHits_list3_list1_diffQuestion_mask <- read_symmetric_mtx("rHits_list3_list1_diffQuestion_mask", mask_path)
  rHits_list3_list2_diffQuestion_mask <- read_symmetric_mtx("rHits_list3_list2_diffQuestion_mask", mask_path)
  rHits_list3_list4_diffQuestion_mask <- read_symmetric_mtx("rHits_list3_list4_diffQuestion_mask", mask_path)
  rHits_list3_list5_diffQuestion_mask <- read_symmetric_mtx("rHits_list3_list5_diffQuestion_mask", mask_path)
  rHits_list3_list6_diffQuestion_mask <- read_symmetric_mtx("rHits_list3_list6_diffQuestion_mask", mask_path)
  rHits_list3_list7_diffQuestion_mask <- read_symmetric_mtx("rHits_list3_list7_diffQuestion_mask", mask_path)
  rHits_list3_list8_diffQuestion_mask <- read_symmetric_mtx("rHits_list3_list8_diffQuestion_mask", mask_path)
  # list4
  rHits_list4_list1_sameQuestion_mask <- read_symmetric_mtx("rHits_list4_list1_sameQuestion_mask", mask_path)
  rHits_list4_list2_sameQuestion_mask <- read_symmetric_mtx("rHits_list4_list2_sameQuestion_mask", mask_path)
  rHits_list4_list3_sameQuestion_mask <- read_symmetric_mtx("rHits_list4_list3_sameQuestion_mask", mask_path)
  rHits_list4_list5_sameQuestion_mask <- read_symmetric_mtx("rHits_list4_list5_sameQuestion_mask", mask_path)
  rHits_list4_list6_sameQuestion_mask <- read_symmetric_mtx("rHits_list4_list6_sameQuestion_mask", mask_path)
  rHits_list4_list7_sameQuestion_mask <- read_symmetric_mtx("rHits_list4_list7_sameQuestion_mask", mask_path)
  rHits_list4_list8_sameQuestion_mask <- read_symmetric_mtx("rHits_list4_list8_sameQuestion_mask", mask_path)
  rHits_list4_list1_diffQuestion_mask <- read_symmetric_mtx("rHits_list4_list1_diffQuestion_mask", mask_path)
  rHits_list4_list2_diffQuestion_mask <- read_symmetric_mtx("rHits_list4_list2_diffQuestion_mask", mask_path)
  rHits_list4_list3_diffQuestion_mask <- read_symmetric_mtx("rHits_list4_list3_diffQuestion_mask", mask_path)
  rHits_list4_list5_diffQuestion_mask <- read_symmetric_mtx("rHits_list4_list5_diffQuestion_mask", mask_path)
  rHits_list4_list6_diffQuestion_mask <- read_symmetric_mtx("rHits_list4_list6_diffQuestion_mask", mask_path)
  rHits_list4_list7_diffQuestion_mask <- read_symmetric_mtx("rHits_list4_list7_diffQuestion_mask", mask_path)
  rHits_list4_list8_diffQuestion_mask <- read_symmetric_mtx("rHits_list4_list8_diffQuestion_mask", mask_path)
  # list5
  rHits_list5_list1_sameQuestion_mask <- read_symmetric_mtx("rHits_list5_list1_sameQuestion_mask", mask_path)
  rHits_list5_list2_sameQuestion_mask <- read_symmetric_mtx("rHits_list5_list2_sameQuestion_mask", mask_path)
  rHits_list5_list3_sameQuestion_mask <- read_symmetric_mtx("rHits_list5_list3_sameQuestion_mask", mask_path)
  rHits_list5_list4_sameQuestion_mask <- read_symmetric_mtx("rHits_list5_list4_sameQuestion_mask", mask_path)
  rHits_list5_list6_sameQuestion_mask <- read_symmetric_mtx("rHits_list5_list6_sameQuestion_mask", mask_path)
  rHits_list5_list7_sameQuestion_mask <- read_symmetric_mtx("rHits_list5_list7_sameQuestion_mask", mask_path)
  rHits_list5_list8_sameQuestion_mask <- read_symmetric_mtx("rHits_list5_list8_sameQuestion_mask", mask_path)
  rHits_list5_list1_diffQuestion_mask <- read_symmetric_mtx("rHits_list5_list1_diffQuestion_mask", mask_path)
  rHits_list5_list2_diffQuestion_mask <- read_symmetric_mtx("rHits_list5_list2_diffQuestion_mask", mask_path)
  rHits_list5_list3_diffQuestion_mask <- read_symmetric_mtx("rHits_list5_list3_diffQuestion_mask", mask_path)
  rHits_list5_list4_diffQuestion_mask <- read_symmetric_mtx("rHits_list5_list4_diffQuestion_mask", mask_path)
  rHits_list5_list6_diffQuestion_mask <- read_symmetric_mtx("rHits_list5_list6_diffQuestion_mask", mask_path)
  rHits_list5_list7_diffQuestion_mask <- read_symmetric_mtx("rHits_list5_list7_diffQuestion_mask", mask_path)
  rHits_list5_list8_diffQuestion_mask <- read_symmetric_mtx("rHits_list5_list8_diffQuestion_mask", mask_path)
  # list6
  rHits_list6_list1_sameQuestion_mask <- read_symmetric_mtx("rHits_list6_list1_sameQuestion_mask", mask_path)
  rHits_list6_list2_sameQuestion_mask <- read_symmetric_mtx("rHits_list6_list2_sameQuestion_mask", mask_path)
  rHits_list6_list3_sameQuestion_mask <- read_symmetric_mtx("rHits_list6_list3_sameQuestion_mask", mask_path)
  rHits_list6_list4_sameQuestion_mask <- read_symmetric_mtx("rHits_list6_list4_sameQuestion_mask", mask_path)
  rHits_list6_list5_sameQuestion_mask <- read_symmetric_mtx("rHits_list6_list5_sameQuestion_mask", mask_path)
  rHits_list6_list7_sameQuestion_mask <- read_symmetric_mtx("rHits_list6_list7_sameQuestion_mask", mask_path)
  rHits_list6_list8_sameQuestion_mask <- read_symmetric_mtx("rHits_list6_list8_sameQuestion_mask", mask_path)
  rHits_list6_list1_diffQuestion_mask <- read_symmetric_mtx("rHits_list6_list1_diffQuestion_mask", mask_path)
  rHits_list6_list2_diffQuestion_mask <- read_symmetric_mtx("rHits_list6_list2_diffQuestion_mask", mask_path)
  rHits_list6_list3_diffQuestion_mask <- read_symmetric_mtx("rHits_list6_list3_diffQuestion_mask", mask_path)
  rHits_list6_list4_diffQuestion_mask <- read_symmetric_mtx("rHits_list6_list4_diffQuestion_mask", mask_path)
  rHits_list6_list5_diffQuestion_mask <- read_symmetric_mtx("rHits_list6_list5_diffQuestion_mask", mask_path)
  rHits_list6_list7_diffQuestion_mask <- read_symmetric_mtx("rHits_list6_list7_diffQuestion_mask", mask_path)
  rHits_list6_list8_diffQuestion_mask <- read_symmetric_mtx("rHits_list6_list8_diffQuestion_mask", mask_path)
  # list7
  rHits_list7_list1_sameQuestion_mask <- read_symmetric_mtx("rHits_list7_list1_sameQuestion_mask", mask_path)
  rHits_list7_list2_sameQuestion_mask <- read_symmetric_mtx("rHits_list7_list2_sameQuestion_mask", mask_path)
  rHits_list7_list3_sameQuestion_mask <- read_symmetric_mtx("rHits_list7_list3_sameQuestion_mask", mask_path)
  rHits_list7_list4_sameQuestion_mask <- read_symmetric_mtx("rHits_list7_list4_sameQuestion_mask", mask_path)
  rHits_list7_list5_sameQuestion_mask <- read_symmetric_mtx("rHits_list7_list5_sameQuestion_mask", mask_path)
  rHits_list7_list6_sameQuestion_mask <- read_symmetric_mtx("rHits_list7_list6_sameQuestion_mask", mask_path)
  rHits_list7_list8_sameQuestion_mask <- read_symmetric_mtx("rHits_list7_list8_sameQuestion_mask", mask_path)
  rHits_list7_list1_diffQuestion_mask <- read_symmetric_mtx("rHits_list7_list1_diffQuestion_mask", mask_path)
  rHits_list7_list2_diffQuestion_mask <- read_symmetric_mtx("rHits_list7_list2_diffQuestion_mask", mask_path)
  rHits_list7_list3_diffQuestion_mask <- read_symmetric_mtx("rHits_list7_list3_diffQuestion_mask", mask_path)
  rHits_list7_list4_diffQuestion_mask <- read_symmetric_mtx("rHits_list7_list4_diffQuestion_mask", mask_path)
  rHits_list7_list5_diffQuestion_mask <- read_symmetric_mtx("rHits_list7_list5_diffQuestion_mask", mask_path)
  rHits_list7_list6_diffQuestion_mask <- read_symmetric_mtx("rHits_list7_list6_diffQuestion_mask", mask_path)
  rHits_list7_list8_diffQuestion_mask <- read_symmetric_mtx("rHits_list7_list8_diffQuestion_mask", mask_path)
  # list8
  rHits_list8_list1_sameQuestion_mask <- read_symmetric_mtx("rHits_list8_list1_sameQuestion_mask", mask_path)
  rHits_list8_list2_sameQuestion_mask <- read_symmetric_mtx("rHits_list8_list2_sameQuestion_mask", mask_path)
  rHits_list8_list3_sameQuestion_mask <- read_symmetric_mtx("rHits_list8_list3_sameQuestion_mask", mask_path)
  rHits_list8_list4_sameQuestion_mask <- read_symmetric_mtx("rHits_list8_list4_sameQuestion_mask", mask_path)
  rHits_list8_list5_sameQuestion_mask <- read_symmetric_mtx("rHits_list8_list5_sameQuestion_mask", mask_path)
  rHits_list8_list6_sameQuestion_mask <- read_symmetric_mtx("rHits_list8_list6_sameQuestion_mask", mask_path)
  rHits_list8_list7_sameQuestion_mask <- read_symmetric_mtx("rHits_list8_list7_sameQuestion_mask", mask_path)
  rHits_list8_list1_diffQuestion_mask <- read_symmetric_mtx("rHits_list8_list1_diffQuestion_mask", mask_path)
  rHits_list8_list2_diffQuestion_mask <- read_symmetric_mtx("rHits_list8_list2_diffQuestion_mask", mask_path)
  rHits_list8_list3_diffQuestion_mask <- read_symmetric_mtx("rHits_list8_list3_diffQuestion_mask", mask_path)
  rHits_list8_list4_diffQuestion_mask <- read_symmetric_mtx("rHits_list8_list4_diffQuestion_mask", mask_path)
  rHits_list8_list5_diffQuestion_mask <- read_symmetric_mtx("rHits_list8_list5_diffQuestion_mask", mask_path)
  rHits_list8_list6_diffQuestion_mask <- read_symmetric_mtx("rHits_list8_list6_diffQuestion_mask", mask_path)
  rHits_list8_list7_diffQuestion_mask <- read_symmetric_mtx("rHits_list8_list7_diffQuestion_mask", mask_path)
  # --- rhits and source hits: split up by each list ---
  rANDsourceHits_list1_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list1_sameQuestion_mask", mask_path)
  rANDsourceHits_list2_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list2_sameQuestion_mask", mask_path)
  rANDsourceHits_list3_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list3_sameQuestion_mask", mask_path)
  rANDsourceHits_list4_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list4_sameQuestion_mask", mask_path)
  rANDsourceHits_list5_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list5_sameQuestion_mask", mask_path)
  rANDsourceHits_list6_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list6_sameQuestion_mask", mask_path)
  rANDsourceHits_list7_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list7_sameQuestion_mask", mask_path)
  rANDsourceHits_list8_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list8_sameQuestion_mask", mask_path)
  rANDsourceHits_list1_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list1_diffQuestion_mask", mask_path)
  rANDsourceHits_list2_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list2_diffQuestion_mask", mask_path)
  rANDsourceHits_list3_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list3_diffQuestion_mask", mask_path)
  rANDsourceHits_list4_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list4_diffQuestion_mask", mask_path)
  rANDsourceHits_list5_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list5_diffQuestion_mask", mask_path)
  rANDsourceHits_list6_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list6_diffQuestion_mask", mask_path)
  rANDsourceHits_list7_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list7_diffQuestion_mask", mask_path)
  rANDsourceHits_list8_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_list8_diffQuestion_mask", mask_path)
  # --- rhits: same/diff question and list ---
  rHits_sameList_sameQuestion_mask <- read_symmetric_mtx("rHits_sameList_sameQuestion_mask", mask_path)
  rHits_diffList_sameQuestion_mask <- read_symmetric_mtx("rHits_diffList_sameQuestion_mask", mask_path)
  rHits_sameList_diffQuestion_mask <- read_symmetric_mtx("rHits_sameList_diffQuestion_mask", mask_path)
  rHits_diffList_diffQuestion_mask <- read_symmetric_mtx("rHits_diffList_diffQuestion_mask", mask_path)
  # --- rhits and source hits: same/diff question and list ---
  rANDsourceHits_sameList_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_sameList_sameQuestion_mask", mask_path)
  rANDsourceHits_sameList_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_sameList_diffQuestion_mask", mask_path)
  rANDsourceHits_diffList_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_diffList_diffQuestion_mask", mask_path)
  rANDsourceHits_diffList_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_diffList_sameQuestion_mask", mask_path)
  rANDsourceHits_sameQuestion_mask <- read_symmetric_mtx("rANDsourceHits_sameQuestion_mask", mask_path)
  rANDsourceHits_diffQuestion_mask <- read_symmetric_mtx("rANDsourceHits_diffQuestion_mask", mask_path)
  # --- rhits and liberal question source scoring ---
  rHits_questionMissLiberal_mask <- read_symmetric_mtx("rHits_questionMissLiberal_mask", mask_path)
  rHits_questionHitsLiberal_mask <- read_symmetric_mtx("rHits_questionHitsLiberal_mask", mask_path)

  # --- item miss, question source miss, temporal source miss: split up by list ---
  # --- item miss, question source miss, temporal source miss ---
  # --- item miss: split up by list ---
  # --- item miss: split up by each list ---
  itemMiss_list1_sameQuestion_mask <- read_symmetric_mtx("itemMiss_list1_sameQuestion_mask", mask_path)
  itemMiss_list2_sameQuestion_mask <- read_symmetric_mtx("itemMiss_list2_sameQuestion_mask", mask_path)
  itemMiss_list3_sameQuestion_mask <- read_symmetric_mtx("itemMiss_list3_sameQuestion_mask", mask_path)
  itemMiss_list4_sameQuestion_mask <- read_symmetric_mtx("itemMiss_list4_sameQuestion_mask", mask_path)
  itemMiss_list5_sameQuestion_mask <- read_symmetric_mtx("itemMiss_list5_sameQuestion_mask", mask_path)
  itemMiss_list6_sameQuestion_mask <- read_symmetric_mtx("itemMiss_list6_sameQuestion_mask", mask_path)
  itemMiss_list7_sameQuestion_mask <- read_symmetric_mtx("itemMiss_list7_sameQuestion_mask", mask_path)
  itemMiss_list8_sameQuestion_mask <- read_symmetric_mtx("itemMiss_list8_sameQuestion_mask", mask_path)
  itemMiss_list1_diffQuestion_mask <- read_symmetric_mtx("itemMiss_list1_diffQuestion_mask", mask_path)
  itemMiss_list2_diffQuestion_mask <- read_symmetric_mtx("itemMiss_list2_diffQuestion_mask", mask_path)
  itemMiss_list3_diffQuestion_mask <- read_symmetric_mtx("itemMiss_list3_diffQuestion_mask", mask_path)
  itemMiss_list4_diffQuestion_mask <- read_symmetric_mtx("itemMiss_list4_diffQuestion_mask", mask_path)
  itemMiss_list5_diffQuestion_mask <- read_symmetric_mtx("itemMiss_list5_diffQuestion_mask", mask_path)
  itemMiss_list6_diffQuestion_mask <- read_symmetric_mtx("itemMiss_list6_diffQuestion_mask", mask_path)
  itemMiss_list7_diffQuestion_mask <- read_symmetric_mtx("itemMiss_list7_diffQuestion_mask", mask_path)
  itemMiss_list8_diffQuestion_mask <- read_symmetric_mtx("itemMiss_list8_diffQuestion_mask", mask_path)
  # --- item miss: same/diff question and list ---
  itemMiss_sameList_sameQuestion_mask <- read_symmetric_mtx("itemMiss_sameList_sameQuestion_mask", mask_path)
  itemMiss_diffList_sameQuestion_mask <- read_symmetric_mtx("itemMiss_diffList_sameQuestion_mask", mask_path)
  itemMiss_sameList_diffQuestion_mask <- read_symmetric_mtx("itemMiss_sameList_diffQuestion_mask", mask_path)
  itemMiss_diffList_diffQuestion_mask <- read_symmetric_mtx("itemMiss_diffList_diffQuestion_mask", mask_path)

  # --- fhits and misses: split up by each list ---
  fHitsANDmisses_list1_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list1_sameQuestion_mask", mask_path)
  fHitsANDmisses_list2_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list2_sameQuestion_mask", mask_path)
  fHitsANDmisses_list3_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list3_sameQuestion_mask", mask_path)
  fHitsANDmisses_list4_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list4_sameQuestion_mask", mask_path)
  fHitsANDmisses_list5_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list5_sameQuestion_mask", mask_path)
  fHitsANDmisses_list6_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list6_sameQuestion_mask", mask_path)
  fHitsANDmisses_list7_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list7_sameQuestion_mask", mask_path)
  fHitsANDmisses_list8_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list8_sameQuestion_mask", mask_path)
  fHitsANDmisses_list1_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list1_diffQuestion_mask", mask_path)
  fHitsANDmisses_list2_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list2_diffQuestion_mask", mask_path)
  fHitsANDmisses_list3_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list3_diffQuestion_mask", mask_path)
  fHitsANDmisses_list4_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list4_diffQuestion_mask", mask_path)
  fHitsANDmisses_list5_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list5_diffQuestion_mask", mask_path)
  fHitsANDmisses_list6_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list6_diffQuestion_mask", mask_path)
  fHitsANDmisses_list7_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list7_diffQuestion_mask", mask_path)
  fHitsANDmisses_list8_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_list8_diffQuestion_mask", mask_path)
  # --- fhits and misses: same/diff question and list ---
  fHitsANDmisses_sameList_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_sameList_sameQuestion_mask", mask_path)
  fHitsANDmisses_diffList_sameQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_diffList_sameQuestion_mask", mask_path)
  fHitsANDmisses_sameList_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_sameList_diffQuestion_mask", mask_path)
  fHitsANDmisses_diffList_diffQuestion_mask <- read_symmetric_mtx("fHitsANDmisses_diffList_diffQuestion_mask", mask_path)

  # read in trial information file
  # this gets created in `create-onset-files.R`
  cur_subj_info_file <- NULL
  cur_subj_info_file <- paste0(halle::ensure_trailing_slash(analyzed_behavioral_dir),
                               halle::ensure_trailing_slash(isubj),
                               "trial_information.csv")

  if(!file.exists(cur_subj_info_file)){
    print(sprintf("\nNo subject info file (%s) - skipping.\n", cur_subj_info_file))
    next
  } else{
    cur_subj_info <- NULL
    cur_subj_info <- read.csv(cur_subj_info_file)
  } # if(file.exists

  # loop across hemispheres
  for(idir in c(1:length(roi_dirs))){
    hemi_label <- NULL
    # figure out L vs R if using ASHS ROIs
    if(roi_dirs[idir]=="ashs_left"){
      hemi_label <- "left"
    } else if(roi_dirs[idir]=="ashs_right") {
      hemi_label <- "right"
    } else if(roi_dirs[idir] == "manual_traced_rois") {
      # we will actually set the hemisphere label later once we have the name of the ROI,
      # but setting some temporary value for now
      hemi_label <- "manual"
    }

    # loop across rois
    for(iroi in c(1:length(rois))){
      mm_trial_pairs <- data.frame()

      # define directories and filenames
      base_roi_dir <- paste(analyzed_mri_dir,halle::ensure_trailing_slash(isubj),halle::ensure_trailing_slash('ROIs'),sep="")
      cur_roi_dir <- paste(base_roi_dir,halle::ensure_trailing_slash(roi_dirs[idir]),sep="")
      cur_roi <- strsplit(rois[iroi],".nii")

      cat(sprintf("\nWorking on hemi: %s, cur_roi: %s\n",hemi_label,cur_roi))

      # --- read in the current subject's correlation matrix ---
      # this file gets created by `RSA_extract_betas_from_ROIs.m`
      cur_pattern_corr_file <- NULL
      cur_pattern_corr_file <- paste0(cur_roi_dir, paste0('br', cur_roi, '_pattern_corr_', MAP_TYPE_FLAG,'_all_runs.mat'))

      # check to make sure current run file exists before continuing on
      if(!file.exists(cur_pattern_corr_file)){
        cat(sprintf("\nPattern correlations file %s does not exist - skipping.\n", cur_pattern_corr_file))
        next
      } else{

        # for manually-traced ROIs, re-set the value for `hemi_label` to be meaningful
        if(roi_dirs[idir] == "manual_traced_rois") {
          # I think you need to double index b/c technically [1] isn't the string alone?
          # (this is a case of knowing that this works, but not fully knowing why - maybe it's a list?)
          roi_hemi_split <- strsplit(cur_roi[[1]], "_")
          hemi_only <- unlist(roi_hemi_split)[2]
          # now, make `hemi_label` match the formatting for the ashs rois
          if(hemi_only == "L"){
            hemi_label <- "left"
          } else if (hemi_only == "R") {
            hemi_label <- "right"
          }
        }

        cur_corr <- NULL
        cur_corr <- data.matrix(as.data.frame(R.matlab::readMat(cur_pattern_corr_file)))

        # add in rownames
        # this is easy to do b/c it's a symmetric matrix
        cur_corr_w_rownames <- cur_corr
        rownames(cur_corr_w_rownames) <- colnames(cur_corr)

        # --- notch out the within run correlations ---
        # because we know each run had 60 trials,
        # we've already setup a matrix of NaN values for w/in run correlations
        # and now we're just going to multiply that by our `cur_corr_w_rownames`
        cur_corr_no_autocorr <- NULL
        cur_corr_no_autocorr <- cur_corr_w_rownames * within_run_corr_as_NaN

        # if want to graphically check
        # ggplot2::ggplot(reshape2::melt(cur_corr_no_autocorr), ggplot2::aes(Var2, Var1, fill = value)) + ggplot2::geom_tile() + ggplot2::ggtitle(paste0(isubj, " ", hemi_label, " ", cur_roi))

        # --- select the trial pairs of interest ---
        # Split by question id, don't worry about memory
        question1 <- cur_corr_no_autocorr * question1_mask
        question1_mtx <-
          question1 %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "q1",
                        item_mem = "unknown",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "unknown",
                        list = "unknown",
                        source_mem = "unknown",
                        condition = "unknown",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        question1_trialnums <- length(question1_mtx$r)

        question2 <- cur_corr_no_autocorr * question2_mask
        question2_mtx <-
          question2 %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "q2",
                        item_mem = "unknown",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "unknown",
                        list = "unknown",
                        source_mem = "unknown",
                        condition = "unknown",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        question2_trialnums <- length(question2_mtx$r)

        question3 <- cur_corr_no_autocorr * question3_mask
        question3_mtx <-
          question3 %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "q3",
                        item_mem = "unknown",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "unknown",
                        list = "unknown",
                        source_mem = "unknown",
                        condition = "unknown",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        question3_trialnums <- length(question3_mtx$r)

        question4 <- cur_corr_no_autocorr * question4_mask
        question4_mtx <-
          question4 %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "q4",
                        item_mem = "unknown",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "unknown",
                        list = "unknown",
                        source_mem = "unknown",
                        condition = "unknown",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        question4_trialnums <- length(question4_mtx$r)

        # Full factorial: rhits, list source hits, question source hits
        # do this way so can plot and compare w/ mask
        rANDsourceHits_sameList_sameQuestion <- cur_corr_no_autocorr * rANDsourceHits_sameList_sameQuestion_mask
        rANDsourceHits_sameList_sameQuestion_mtx <-
          rANDsourceHits_sameList_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "correct",
                        list_mem = "correct",
                        question = "same",
                        list = "same",
                        source_mem = "rHit_listHit_questionHit",
                        condition = "sameList_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rANDsourceHits_sameList_sameQuestion_trialnums <- length(rANDsourceHits_sameList_sameQuestion_mtx$r)

        rANDsourceHits_diffList_sameQuestion <- cur_corr_no_autocorr * rANDsourceHits_diffList_sameQuestion_mask
        rANDsourceHits_diffList_sameQuestion_mtx <-
          rANDsourceHits_diffList_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "correct",
                        list_mem = "correct",
                        question = "same",
                        list = "diff",
                        source_mem = "rHit_listHit_questionHit",
                        condition = "diffList_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rANDsourceHits_diffList_sameQuestion_trialnums <- length(rANDsourceHits_diffList_sameQuestion_mtx$r)

        rANDsourceHits_sameList_diffQuestion <- cur_corr_no_autocorr * rANDsourceHits_sameList_diffQuestion_mask
        rANDsourceHits_sameList_diffQuestion_mtx <-
          rANDsourceHits_sameList_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "correct",
                        list_mem = "correct",
                        question = "diff",
                        list = "same",
                        source_mem = "rHit_listHit_questionHit",
                        condition = "sameList_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rANDsourceHits_sameList_diffQuestion_trialnums <- length(rANDsourceHits_sameList_diffQuestion_mtx$r)

        rANDsourceHits_diffList_diffQuestion <- cur_corr_no_autocorr * rANDsourceHits_diffList_diffQuestion_mask
        rANDsourceHits_diffList_diffQuestion_mtx <-
          rANDsourceHits_diffList_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "correct",
                        list_mem = "correct",
                        question = "diff",
                        list = "diff",
                        source_mem = "rHit_listHit_questionHit",
                        condition = "diffList_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rANDsourceHits_diffList_diffQuestion_trialnums <- length(rANDsourceHits_diffList_diffQuestion_mtx$r)

        # Full factorial: rhits
        rHits_sameList_sameQuestion <- cur_corr_no_autocorr * rHits_sameList_sameQuestion_mask
        rHits_sameList_sameQuestion_mtx <-
          rHits_sameList_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "same",
                        source_mem = "rHit",
                        condition = "sameList_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_sameList_sameQuestion_trialnums <- length(rHits_sameList_sameQuestion_mtx$r)

        rHits_diffList_sameQuestion <- cur_corr_no_autocorr * rHits_diffList_sameQuestion_mask
        rHits_diffList_sameQuestion_mtx <-
          rHits_diffList_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "diff",
                        source_mem = "rHit",
                        condition = "diffList_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_diffList_sameQuestion_trialnums <- length(rHits_diffList_sameQuestion_mtx$r)

        rHits_sameList_diffQuestion <- cur_corr_no_autocorr * rHits_sameList_diffQuestion_mask
        rHits_sameList_diffQuestion_mtx <-
          rHits_sameList_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "same",
                        source_mem = "rHit",
                        condition = "sameList_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_sameList_diffQuestion_trialnums <- length(rHits_sameList_diffQuestion_mtx$r)

        rHits_diffList_diffQuestion <- cur_corr_no_autocorr * rHits_diffList_diffQuestion_mask
        rHits_diffList_diffQuestion_mtx <-
          rHits_diffList_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "diff",
                        source_mem = "rHit",
                        condition = "diffList_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_diffList_diffQuestion_trialnums <- length(rHits_diffList_diffQuestion_mtx$r)

        # Same/diff question: rhits, question source hits
        rANDsourceHits_sameQuestion <- cur_corr_no_autocorr * rANDsourceHits_sameQuestion_mask
        rANDsourceHits_sameQuestion_mtx <-
          rANDsourceHits_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "correct",
                        list_mem = "correct",
                        question = "same",
                        list = "unknown",
                        source_mem = "rHit_questionHit",
                        condition = "sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rANDsourceHits_sameQuestion_trialnums <- length(rANDsourceHits_sameQuestion_mtx$r)

        rANDsourceHits_diffQuestion <- cur_corr_no_autocorr * rANDsourceHits_diffQuestion_mask
        rANDsourceHits_diffQuestion_mtx <-
          rANDsourceHits_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "correct",
                        list_mem = "correct",
                        question = "same",
                        list = "unknown",
                        source_mem = "rHit_questionHit",
                        condition = "diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rANDsourceHits_diffQuestion_trialnums <- length(rANDsourceHits_diffQuestion_mtx$r)

        # rhits and question source: liberal source scoring
        rHits_questionHitsLiberal <- cur_corr_no_autocorr * rHits_questionHitsLiberal_mask
        rHits_questionHitsLiberal_mtx <-
          rHits_questionHitsLiberal %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "correctQuestion",
                        list_mem = "unknown",
                        question = "unknown",
                        list = "unknown",
                        source_mem = "rHit_questionLiberalHit",
                        condition = "questionLiberalHit",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_questionHitsLiberal_trialnums <- length(rHits_questionHitsLiberal_mtx$r)

        rHits_questionMissLiberal <- cur_corr_no_autocorr * rHits_questionMissLiberal_mask
        rHits_questionMissLiberal_mtx <-
          rHits_questionMissLiberal %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "incorrectQuestion",
                        list_mem = "unknown",
                        question = "unknown",
                        list = "unknown",
                        source_mem = "rHit_questionLiberalMiss",
                        condition = "questionLiberalMiss",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_questionMissLiberal_trialnums <- length(rHits_questionMissLiberal_mtx$r)


        # Full factorial: item misses
        itemMiss_sameList_sameQuestion <- cur_corr_no_autocorr * itemMiss_sameList_sameQuestion_mask
        itemMiss_sameList_sameQuestion_mtx <-
          itemMiss_sameList_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "miss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "same",
                        source_mem = "itemMiss",
                        condition = "sameList_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        itemMiss_sameList_sameQuestion_trialnums <- length(itemMiss_sameList_sameQuestion_mtx$r)

        itemMiss_diffList_sameQuestion <- cur_corr_no_autocorr * itemMiss_diffList_sameQuestion_mask
        itemMiss_diffList_sameQuestion_mtx <-
          itemMiss_diffList_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "miss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "diff",
                        source_mem = "itemMiss",
                        condition = "diffList_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        itemMiss_diffList_sameQuestion_trialnums <- length(itemMiss_diffList_sameQuestion_mtx$r)

        itemMiss_sameList_diffQuestion <- cur_corr_no_autocorr * itemMiss_sameList_diffQuestion_mask
        itemMiss_sameList_diffQuestion_mtx <-
          itemMiss_sameList_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "miss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "same",
                        source_mem = "itemMiss",
                        condition = "sameList_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        itemMiss_sameList_diffQuestion_trialnums <- length(itemMiss_sameList_diffQuestion_mtx$r)

        itemMiss_diffList_diffQuestion <- cur_corr_no_autocorr * itemMiss_diffList_diffQuestion_mask
        itemMiss_diffList_diffQuestion_mtx <-
          itemMiss_diffList_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "miss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "diff",
                        source_mem = "itemMiss",
                        condition = "diffList_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        itemMiss_diffList_diffQuestion_trialnums <- length(itemMiss_diffList_diffQuestion_mtx$r)

        # --- Full factorial: rhits, list source hits, question source hits split by list ---
        #\\TODO: add these in later (not starting here b/c probably won't have very many trials)
        # do this way so can plot and compare w/ mask
        rANDsourceHits_list1_sameQuestion <- cur_corr_no_autocorr * rANDsourceHits_list1_sameQuestion_mask
        rANDsourceHits_list1_sameQuestion_mtx <-
          rANDsourceHits_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "correct",
                        list_mem = "correct",
                        question = "same",
                        list = "list1",
                        source_mem = "rHit_listHit_questionHit",
                        condition = "list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rANDsourceHits_list1_sameQuestion_trialnums <- length(rANDsourceHits_list1_sameQuestion_mtx$r)

        # Full factorial: rhits split by list
        rHits_list1_sameQuestion <- cur_corr_no_autocorr * rHits_list1_sameQuestion_mask
        rHits_list1_sameQuestion_mtx <-
          rHits_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1",
                        source_mem = "rHit",
                        condition = "list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_sameQuestion_trialnums <- length(rHits_list1_sameQuestion_mtx$r)

        rHits_list2_sameQuestion <- cur_corr_no_autocorr * rHits_list2_sameQuestion_mask
        rHits_list2_sameQuestion_mtx <-
          rHits_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2",
                        source_mem = "rHit",
                        condition = "list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_sameQuestion_trialnums <- length(rHits_list2_sameQuestion_mtx$r)

        rHits_list3_sameQuestion <- cur_corr_no_autocorr * rHits_list3_sameQuestion_mask
        rHits_list3_sameQuestion_mtx <-
          rHits_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3",
                        source_mem = "rHit",
                        condition = "list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_sameQuestion_trialnums <- length(rHits_list3_sameQuestion_mtx$r)

        rHits_list4_sameQuestion <- cur_corr_no_autocorr * rHits_list4_sameQuestion_mask
        rHits_list4_sameQuestion_mtx <-
          rHits_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4",
                        source_mem = "rHit",
                        condition = "list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_sameQuestion_trialnums <- length(rHits_list4_sameQuestion_mtx$r)

        rHits_list5_sameQuestion <- cur_corr_no_autocorr * rHits_list5_sameQuestion_mask
        rHits_list5_sameQuestion_mtx <-
          rHits_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5",
                        source_mem = "rHit",
                        condition = "list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_sameQuestion_trialnums <- length(rHits_list5_sameQuestion_mtx$r)

        rHits_list6_sameQuestion <- cur_corr_no_autocorr * rHits_list6_sameQuestion_mask
        rHits_list6_sameQuestion_mtx <-
          rHits_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6",
                        source_mem = "rHit",
                        condition = "list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_sameQuestion_trialnums <- length(rHits_list6_sameQuestion_mtx$r)

        rHits_list7_sameQuestion <- cur_corr_no_autocorr * rHits_list7_sameQuestion_mask
        rHits_list7_sameQuestion_mtx <-
          rHits_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7",
                        source_mem = "rHit",
                        condition = "list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_sameQuestion_trialnums <- length(rHits_list7_sameQuestion_mtx$r)

        rHits_list8_sameQuestion <- cur_corr_no_autocorr * rHits_list8_sameQuestion_mask
        rHits_list8_sameQuestion_mtx <-
          rHits_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8",
                        source_mem = "rHit",
                        condition = "list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_sameQuestion_trialnums <- length(rHits_list8_sameQuestion_mtx$r)

        rHits_list1_diffQuestion <- cur_corr_no_autocorr * rHits_list1_diffQuestion_mask
        rHits_list1_diffQuestion_mtx <-
          rHits_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1",
                        source_mem = "rHit",
                        condition = "list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_diffQuestion_trialnums <- length(rHits_list1_diffQuestion_mtx$r)

        rHits_list2_diffQuestion <- cur_corr_no_autocorr * rHits_list2_diffQuestion_mask
        rHits_list2_diffQuestion_mtx <-
          rHits_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2",
                        source_mem = "rHit",
                        condition = "list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_diffQuestion_trialnums <- length(rHits_list2_diffQuestion_mtx$r)

        rHits_list3_diffQuestion <- cur_corr_no_autocorr * rHits_list3_diffQuestion_mask
        rHits_list3_diffQuestion_mtx <-
          rHits_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3",
                        source_mem = "rHit",
                        condition = "list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_diffQuestion_trialnums <- length(rHits_list3_diffQuestion_mtx$r)

        rHits_list4_diffQuestion <- cur_corr_no_autocorr * rHits_list4_diffQuestion_mask
        rHits_list4_diffQuestion_mtx <-
          rHits_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4",
                        source_mem = "rHit",
                        condition = "list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_diffQuestion_trialnums <- length(rHits_list4_diffQuestion_mtx$r)

        rHits_list5_diffQuestion <- cur_corr_no_autocorr * rHits_list5_diffQuestion_mask
        rHits_list5_diffQuestion_mtx <-
          rHits_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5",
                        source_mem = "rHit",
                        condition = "list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_diffQuestion_trialnums <- length(rHits_list5_diffQuestion_mtx$r)

        rHits_list6_diffQuestion <- cur_corr_no_autocorr * rHits_list6_diffQuestion_mask
        rHits_list6_diffQuestion_mtx <-
          rHits_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6",
                        source_mem = "rHit",
                        condition = "list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_diffQuestion_trialnums <- length(rHits_list6_diffQuestion_mtx$r)

        rHits_list7_diffQuestion <- cur_corr_no_autocorr * rHits_list7_diffQuestion_mask
        rHits_list7_diffQuestion_mtx <-
          rHits_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7",
                        source_mem = "rHit",
                        condition = "list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_diffQuestion_trialnums <- length(rHits_list7_diffQuestion_mtx$r)

        rHits_list8_diffQuestion <- cur_corr_no_autocorr * rHits_list8_diffQuestion_mask
        rHits_list8_diffQuestion_mtx <-
          rHits_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8",
                        source_mem = "rHit",
                        condition = "list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_diffQuestion_trialnums <- length(rHits_list8_diffQuestion_mtx$r)

        # --- rHits by list, off-diagonal elements ---
        # list1
        rHits_list1_list2_sameQuestion <- cur_corr_no_autocorr * rHits_list1_list2_sameQuestion_mask
        rHits_list1_list2_sameQuestion_mtx <-
          rHits_list1_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1_list2",
                        source_mem = "rHit",
                        condition = "list1_list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list3_sameQuestion <- cur_corr_no_autocorr * rHits_list1_list3_sameQuestion_mask
        rHits_list1_list3_sameQuestion_mtx <-
          rHits_list1_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1_list3",
                        source_mem = "rHit",
                        condition = "list1_list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list4_sameQuestion <- cur_corr_no_autocorr * rHits_list1_list4_sameQuestion_mask
        rHits_list1_list4_sameQuestion_mtx <-
          rHits_list1_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1_list4",
                        source_mem = "rHit",
                        condition = "list1_list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list5_sameQuestion <- cur_corr_no_autocorr * rHits_list1_list5_sameQuestion_mask
        rHits_list1_list5_sameQuestion_mtx <-
          rHits_list1_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1_list5",
                        source_mem = "rHit",
                        condition = "list1_list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list6_sameQuestion <- cur_corr_no_autocorr * rHits_list1_list6_sameQuestion_mask
        rHits_list1_list6_sameQuestion_mtx <-
          rHits_list1_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1_list6",
                        source_mem = "rHit",
                        condition = "list1_list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list7_sameQuestion <- cur_corr_no_autocorr * rHits_list1_list7_sameQuestion_mask
        rHits_list1_list7_sameQuestion_mtx <-
          rHits_list1_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1_list7",
                        source_mem = "rHit",
                        condition = "list1_list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list8_sameQuestion <- cur_corr_no_autocorr * rHits_list1_list8_sameQuestion_mask
        rHits_list1_list8_sameQuestion_mtx <-
          rHits_list1_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1_list8",
                        source_mem = "rHit",
                        condition = "list1_list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list2_diffQuestion <- cur_corr_no_autocorr * rHits_list1_list2_diffQuestion_mask
        rHits_list1_list2_diffQuestion_mtx <-
          rHits_list1_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1_list2",
                        source_mem = "rHit",
                        condition = "list1_list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list3_diffQuestion <- cur_corr_no_autocorr * rHits_list1_list3_diffQuestion_mask
        rHits_list1_list3_diffQuestion_mtx <-
          rHits_list1_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1_list3",
                        source_mem = "rHit",
                        condition = "list1_list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list4_diffQuestion <- cur_corr_no_autocorr * rHits_list1_list4_diffQuestion_mask
        rHits_list1_list4_diffQuestion_mtx <-
          rHits_list1_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1_list4",
                        source_mem = "rHit",
                        condition = "list1_list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list5_diffQuestion <- cur_corr_no_autocorr * rHits_list1_list5_diffQuestion_mask
        rHits_list1_list5_diffQuestion_mtx <-
          rHits_list1_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1_list5",
                        source_mem = "rHit",
                        condition = "list1_list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list6_diffQuestion <- cur_corr_no_autocorr * rHits_list1_list6_diffQuestion_mask
        rHits_list1_list6_diffQuestion_mtx <-
          rHits_list1_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1_list6",
                        source_mem = "rHit",
                        condition = "list1_list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list7_diffQuestion <- cur_corr_no_autocorr * rHits_list1_list7_diffQuestion_mask
        rHits_list1_list7_diffQuestion_mtx <-
          rHits_list1_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1_list7",
                        source_mem = "rHit",
                        condition = "list1_list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list1_list8_diffQuestion <- cur_corr_no_autocorr * rHits_list1_list8_diffQuestion_mask
        rHits_list1_list8_diffQuestion_mtx <-
          rHits_list1_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1_list8",
                        source_mem = "rHit",
                        condition = "list1_list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        # list2
        rHits_list2_list1_sameQuestion <- cur_corr_no_autocorr * rHits_list2_list1_sameQuestion_mask
        rHits_list2_list1_sameQuestion_mtx <-
          rHits_list2_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2_list1",
                        source_mem = "rHit",
                        condition = "list2_list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list3_sameQuestion <- cur_corr_no_autocorr * rHits_list2_list3_sameQuestion_mask
        rHits_list2_list3_sameQuestion_mtx <-
          rHits_list2_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2_list3",
                        source_mem = "rHit",
                        condition = "list2_list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list4_sameQuestion <- cur_corr_no_autocorr * rHits_list2_list4_sameQuestion_mask
        rHits_list2_list4_sameQuestion_mtx <-
          rHits_list2_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2_list4",
                        source_mem = "rHit",
                        condition = "list2_list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list5_sameQuestion <- cur_corr_no_autocorr * rHits_list2_list5_sameQuestion_mask
        rHits_list2_list5_sameQuestion_mtx <-
          rHits_list2_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2_list5",
                        source_mem = "rHit",
                        condition = "list2_list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list6_sameQuestion <- cur_corr_no_autocorr * rHits_list2_list6_sameQuestion_mask
        rHits_list2_list6_sameQuestion_mtx <-
          rHits_list2_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2_list6",
                        source_mem = "rHit",
                        condition = "list2_list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list7_sameQuestion <- cur_corr_no_autocorr * rHits_list2_list7_sameQuestion_mask
        rHits_list2_list7_sameQuestion_mtx <-
          rHits_list2_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2_list7",
                        source_mem = "rHit",
                        condition = "list2_list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list8_sameQuestion <- cur_corr_no_autocorr * rHits_list2_list8_sameQuestion_mask
        rHits_list2_list8_sameQuestion_mtx <-
          rHits_list2_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2_list8",
                        source_mem = "rHit",
                        condition = "list2_list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list1_diffQuestion <- cur_corr_no_autocorr * rHits_list2_list1_diffQuestion_mask
        rHits_list2_list1_diffQuestion_mtx <-
          rHits_list2_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2_list1",
                        source_mem = "rHit",
                        condition = "list2_list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list3_diffQuestion <- cur_corr_no_autocorr * rHits_list2_list3_diffQuestion_mask
        rHits_list2_list3_diffQuestion_mtx <-
          rHits_list2_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2_list3",
                        source_mem = "rHit",
                        condition = "list2_list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list4_diffQuestion <- cur_corr_no_autocorr * rHits_list2_list4_diffQuestion_mask
        rHits_list2_list4_diffQuestion_mtx <-
          rHits_list2_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2_list4",
                        source_mem = "rHit",
                        condition = "list2_list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list5_diffQuestion <- cur_corr_no_autocorr * rHits_list2_list5_diffQuestion_mask
        rHits_list2_list5_diffQuestion_mtx <-
          rHits_list2_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2_list5",
                        source_mem = "rHit",
                        condition = "list2_list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list6_diffQuestion <- cur_corr_no_autocorr * rHits_list2_list6_diffQuestion_mask
        rHits_list2_list6_diffQuestion_mtx <-
          rHits_list2_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2_list6",
                        source_mem = "rHit",
                        condition = "list2_list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list7_diffQuestion <- cur_corr_no_autocorr * rHits_list2_list7_diffQuestion_mask
        rHits_list2_list7_diffQuestion_mtx <-
          rHits_list2_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2_list7",
                        source_mem = "rHit",
                        condition = "list2_list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list2_list8_diffQuestion <- cur_corr_no_autocorr * rHits_list2_list8_diffQuestion_mask
        rHits_list2_list8_diffQuestion_mtx <-
          rHits_list2_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2_list8",
                        source_mem = "rHit",
                        condition = "list2_list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        # list3
        rHits_list3_list1_sameQuestion <- cur_corr_no_autocorr * rHits_list3_list1_sameQuestion_mask
        rHits_list3_list1_sameQuestion_mtx <-
          rHits_list3_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3_list1",
                        source_mem = "rHit",
                        condition = "list3_list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list2_sameQuestion <- cur_corr_no_autocorr * rHits_list3_list2_sameQuestion_mask
        rHits_list3_list2_sameQuestion_mtx <-
          rHits_list3_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3_list2",
                        source_mem = "rHit",
                        condition = "list3_list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list4_sameQuestion <- cur_corr_no_autocorr * rHits_list3_list4_sameQuestion_mask
        rHits_list3_list4_sameQuestion_mtx <-
          rHits_list3_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3_list4",
                        source_mem = "rHit",
                        condition = "list3_list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list5_sameQuestion <- cur_corr_no_autocorr * rHits_list3_list5_sameQuestion_mask
        rHits_list3_list5_sameQuestion_mtx <-
          rHits_list3_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3_list5",
                        source_mem = "rHit",
                        condition = "list3_list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list6_sameQuestion <- cur_corr_no_autocorr * rHits_list3_list6_sameQuestion_mask
        rHits_list3_list6_sameQuestion_mtx <-
          rHits_list3_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3_list6",
                        source_mem = "rHit",
                        condition = "list3_list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list7_sameQuestion <- cur_corr_no_autocorr * rHits_list3_list7_sameQuestion_mask
        rHits_list3_list7_sameQuestion_mtx <-
          rHits_list3_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3_list7",
                        source_mem = "rHit",
                        condition = "list3_list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list8_sameQuestion <- cur_corr_no_autocorr * rHits_list3_list8_sameQuestion_mask
        rHits_list3_list8_sameQuestion_mtx <-
          rHits_list3_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3_list8",
                        source_mem = "rHit",
                        condition = "list3_list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list1_diffQuestion <- cur_corr_no_autocorr * rHits_list3_list1_diffQuestion_mask
        rHits_list3_list1_diffQuestion_mtx <-
          rHits_list3_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3_list1",
                        source_mem = "rHit",
                        condition = "list3_list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list2_diffQuestion <- cur_corr_no_autocorr * rHits_list3_list2_diffQuestion_mask
        rHits_list3_list2_diffQuestion_mtx <-
          rHits_list3_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3_list2",
                        source_mem = "rHit",
                        condition = "list3_list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list4_diffQuestion <- cur_corr_no_autocorr * rHits_list3_list4_diffQuestion_mask
        rHits_list3_list4_diffQuestion_mtx <-
          rHits_list3_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3_list4",
                        source_mem = "rHit",
                        condition = "list3_list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list5_diffQuestion <- cur_corr_no_autocorr * rHits_list3_list5_diffQuestion_mask
        rHits_list3_list5_diffQuestion_mtx <-
          rHits_list3_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3_list5",
                        source_mem = "rHit",
                        condition = "list3_list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list6_diffQuestion <- cur_corr_no_autocorr * rHits_list3_list6_diffQuestion_mask
        rHits_list3_list6_diffQuestion_mtx <-
          rHits_list3_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3_list6",
                        source_mem = "rHit",
                        condition = "list3_list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list7_diffQuestion <- cur_corr_no_autocorr * rHits_list3_list7_diffQuestion_mask
        rHits_list3_list7_diffQuestion_mtx <-
          rHits_list3_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3_list7",
                        source_mem = "rHit",
                        condition = "list3_list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list3_list8_diffQuestion <- cur_corr_no_autocorr * rHits_list3_list8_diffQuestion_mask
        rHits_list3_list8_diffQuestion_mtx <-
          rHits_list3_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3_list8",
                        source_mem = "rHit",
                        condition = "list3_list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        # list4
        rHits_list4_list1_sameQuestion <- cur_corr_no_autocorr * rHits_list4_list1_sameQuestion_mask
        rHits_list4_list1_sameQuestion_mtx <-
          rHits_list4_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4_list1",
                        source_mem = "rHit",
                        condition = "list4_list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list2_sameQuestion <- cur_corr_no_autocorr * rHits_list4_list2_sameQuestion_mask
        rHits_list4_list2_sameQuestion_mtx <-
          rHits_list4_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4_list2",
                        source_mem = "rHit",
                        condition = "list4_list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list3_sameQuestion <- cur_corr_no_autocorr * rHits_list4_list3_sameQuestion_mask
        rHits_list4_list3_sameQuestion_mtx <-
          rHits_list4_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4_list3",
                        source_mem = "rHit",
                        condition = "list4_list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list5_sameQuestion <- cur_corr_no_autocorr * rHits_list4_list5_sameQuestion_mask
        rHits_list4_list5_sameQuestion_mtx <-
          rHits_list4_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4_list5",
                        source_mem = "rHit",
                        condition = "list4_list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list6_sameQuestion <- cur_corr_no_autocorr * rHits_list4_list6_sameQuestion_mask
        rHits_list4_list6_sameQuestion_mtx <-
          rHits_list4_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4_list6",
                        source_mem = "rHit",
                        condition = "list4_list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list7_sameQuestion <- cur_corr_no_autocorr * rHits_list4_list7_sameQuestion_mask
        rHits_list4_list7_sameQuestion_mtx <-
          rHits_list4_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4_list7",
                        source_mem = "rHit",
                        condition = "list4_list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list8_sameQuestion <- cur_corr_no_autocorr * rHits_list4_list8_sameQuestion_mask
        rHits_list4_list8_sameQuestion_mtx <-
          rHits_list4_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4_list8",
                        source_mem = "rHit",
                        condition = "list4_list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list1_diffQuestion <- cur_corr_no_autocorr * rHits_list4_list1_diffQuestion_mask
        rHits_list4_list1_diffQuestion_mtx <-
          rHits_list4_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4_list1",
                        source_mem = "rHit",
                        condition = "list4_list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list2_diffQuestion <- cur_corr_no_autocorr * rHits_list4_list2_diffQuestion_mask
        rHits_list4_list2_diffQuestion_mtx <-
          rHits_list4_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4_list2",
                        source_mem = "rHit",
                        condition = "list4_list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list3_diffQuestion <- cur_corr_no_autocorr * rHits_list4_list3_diffQuestion_mask
        rHits_list4_list3_diffQuestion_mtx <-
          rHits_list4_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4_list3",
                        source_mem = "rHit",
                        condition = "list4_list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list5_diffQuestion <- cur_corr_no_autocorr * rHits_list4_list5_diffQuestion_mask
        rHits_list4_list5_diffQuestion_mtx <-
          rHits_list4_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4_list5",
                        source_mem = "rHit",
                        condition = "list4_list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list6_diffQuestion <- cur_corr_no_autocorr * rHits_list4_list6_diffQuestion_mask
        rHits_list4_list6_diffQuestion_mtx <-
          rHits_list4_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4_list6",
                        source_mem = "rHit",
                        condition = "list4_list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list7_diffQuestion <- cur_corr_no_autocorr * rHits_list4_list7_diffQuestion_mask
        rHits_list4_list7_diffQuestion_mtx <-
          rHits_list4_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4_list7",
                        source_mem = "rHit",
                        condition = "list4_list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list4_list8_diffQuestion <- cur_corr_no_autocorr * rHits_list4_list8_diffQuestion_mask
        rHits_list4_list8_diffQuestion_mtx <-
          rHits_list4_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4_list8",
                        source_mem = "rHit",
                        condition = "list4_list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        # list5
        rHits_list5_list1_sameQuestion <- cur_corr_no_autocorr * rHits_list5_list1_sameQuestion_mask
        rHits_list5_list1_sameQuestion_mtx <-
          rHits_list5_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5_list1",
                        source_mem = "rHit",
                        condition = "list5_list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list2_sameQuestion <- cur_corr_no_autocorr * rHits_list5_list2_sameQuestion_mask
        rHits_list5_list2_sameQuestion_mtx <-
          rHits_list5_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5_list2",
                        source_mem = "rHit",
                        condition = "list5_list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list3_sameQuestion <- cur_corr_no_autocorr * rHits_list5_list3_sameQuestion_mask
        rHits_list5_list3_sameQuestion_mtx <-
          rHits_list5_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5_list3",
                        source_mem = "rHit",
                        condition = "list5_list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list4_sameQuestion <- cur_corr_no_autocorr * rHits_list5_list4_sameQuestion_mask
        rHits_list5_list4_sameQuestion_mtx <-
          rHits_list5_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5_list4",
                        source_mem = "rHit",
                        condition = "list5_list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list6_sameQuestion <- cur_corr_no_autocorr * rHits_list5_list6_sameQuestion_mask
        rHits_list5_list6_sameQuestion_mtx <-
          rHits_list5_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5_list6",
                        source_mem = "rHit",
                        condition = "list5_list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list7_sameQuestion <- cur_corr_no_autocorr * rHits_list5_list7_sameQuestion_mask
        rHits_list5_list7_sameQuestion_mtx <-
          rHits_list5_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5_list7",
                        source_mem = "rHit",
                        condition = "list5_list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list8_sameQuestion <- cur_corr_no_autocorr * rHits_list5_list8_sameQuestion_mask
        rHits_list5_list8_sameQuestion_mtx <-
          rHits_list5_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5_list8",
                        source_mem = "rHit",
                        condition = "list5_list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list1_diffQuestion <- cur_corr_no_autocorr * rHits_list5_list1_diffQuestion_mask
        rHits_list5_list1_diffQuestion_mtx <-
          rHits_list5_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5_list1",
                        source_mem = "rHit",
                        condition = "list5_list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list2_diffQuestion <- cur_corr_no_autocorr * rHits_list5_list2_diffQuestion_mask
        rHits_list5_list2_diffQuestion_mtx <-
          rHits_list5_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5_list2",
                        source_mem = "rHit",
                        condition = "list5_list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list3_diffQuestion <- cur_corr_no_autocorr * rHits_list5_list3_diffQuestion_mask
        rHits_list5_list3_diffQuestion_mtx <-
          rHits_list5_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5_list3",
                        source_mem = "rHit",
                        condition = "list5_list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list4_diffQuestion <- cur_corr_no_autocorr * rHits_list5_list4_diffQuestion_mask
        rHits_list5_list4_diffQuestion_mtx <-
          rHits_list5_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5_list4",
                        source_mem = "rHit",
                        condition = "list5_list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list6_diffQuestion <- cur_corr_no_autocorr * rHits_list5_list6_diffQuestion_mask
        rHits_list5_list6_diffQuestion_mtx <-
          rHits_list5_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5_list6",
                        source_mem = "rHit",
                        condition = "list5_list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list7_diffQuestion <- cur_corr_no_autocorr * rHits_list5_list7_diffQuestion_mask
        rHits_list5_list7_diffQuestion_mtx <-
          rHits_list5_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5_list7",
                        source_mem = "rHit",
                        condition = "list5_list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list5_list8_diffQuestion <- cur_corr_no_autocorr * rHits_list5_list8_diffQuestion_mask
        rHits_list5_list8_diffQuestion_mtx <-
          rHits_list5_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5_list8",
                        source_mem = "rHit",
                        condition = "list5_list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        # list6
        rHits_list6_list1_sameQuestion <- cur_corr_no_autocorr * rHits_list6_list1_sameQuestion_mask
        rHits_list6_list1_sameQuestion_mtx <-
          rHits_list6_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6_list1",
                        source_mem = "rHit",
                        condition = "list6_list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list2_sameQuestion <- cur_corr_no_autocorr * rHits_list6_list2_sameQuestion_mask
        rHits_list6_list2_sameQuestion_mtx <-
          rHits_list6_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6_list2",
                        source_mem = "rHit",
                        condition = "list6_list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list3_sameQuestion <- cur_corr_no_autocorr * rHits_list6_list3_sameQuestion_mask
        rHits_list6_list3_sameQuestion_mtx <-
          rHits_list6_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6_list3",
                        source_mem = "rHit",
                        condition = "list6_list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list4_sameQuestion <- cur_corr_no_autocorr * rHits_list6_list4_sameQuestion_mask
        rHits_list6_list4_sameQuestion_mtx <-
          rHits_list6_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6_list4",
                        source_mem = "rHit",
                        condition = "list6_list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list5_sameQuestion <- cur_corr_no_autocorr * rHits_list6_list5_sameQuestion_mask
        rHits_list6_list5_sameQuestion_mtx <-
          rHits_list6_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6_list5",
                        source_mem = "rHit",
                        condition = "list6_list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list7_sameQuestion <- cur_corr_no_autocorr * rHits_list6_list7_sameQuestion_mask
        rHits_list6_list7_sameQuestion_mtx <-
          rHits_list6_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6_list7",
                        source_mem = "rHit",
                        condition = "list6_list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list8_sameQuestion <- cur_corr_no_autocorr * rHits_list6_list8_sameQuestion_mask
        rHits_list6_list8_sameQuestion_mtx <-
          rHits_list6_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6_list8",
                        source_mem = "rHit",
                        condition = "list6_list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list1_diffQuestion <- cur_corr_no_autocorr * rHits_list6_list1_diffQuestion_mask
        rHits_list6_list1_diffQuestion_mtx <-
          rHits_list6_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6_list1",
                        source_mem = "rHit",
                        condition = "list6_list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list2_diffQuestion <- cur_corr_no_autocorr * rHits_list6_list2_diffQuestion_mask
        rHits_list6_list2_diffQuestion_mtx <-
          rHits_list6_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6_list2",
                        source_mem = "rHit",
                        condition = "list6_list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list3_diffQuestion <- cur_corr_no_autocorr * rHits_list6_list3_diffQuestion_mask
        rHits_list6_list3_diffQuestion_mtx <-
          rHits_list6_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6_list3",
                        source_mem = "rHit",
                        condition = "list6_list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list4_diffQuestion <- cur_corr_no_autocorr * rHits_list6_list4_diffQuestion_mask
        rHits_list6_list4_diffQuestion_mtx <-
          rHits_list6_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6_list4",
                        source_mem = "rHit",
                        condition = "list6_list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list5_diffQuestion <- cur_corr_no_autocorr * rHits_list6_list5_diffQuestion_mask
        rHits_list6_list5_diffQuestion_mtx <-
          rHits_list6_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6_list5",
                        source_mem = "rHit",
                        condition = "list6_list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list7_diffQuestion <- cur_corr_no_autocorr * rHits_list6_list7_diffQuestion_mask
        rHits_list6_list7_diffQuestion_mtx <-
          rHits_list6_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6_list7",
                        source_mem = "rHit",
                        condition = "list6_list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list6_list8_diffQuestion <- cur_corr_no_autocorr * rHits_list6_list8_diffQuestion_mask
        rHits_list6_list8_diffQuestion_mtx <-
          rHits_list6_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6_list8",
                        source_mem = "rHit",
                        condition = "list6_list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        # list7
        rHits_list7_list1_sameQuestion <- cur_corr_no_autocorr * rHits_list7_list1_sameQuestion_mask
        rHits_list7_list1_sameQuestion_mtx <-
          rHits_list7_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7_list1",
                        source_mem = "rHit",
                        condition = "list7_list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list2_sameQuestion <- cur_corr_no_autocorr * rHits_list7_list2_sameQuestion_mask
        rHits_list7_list2_sameQuestion_mtx <-
          rHits_list7_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7_list2",
                        source_mem = "rHit",
                        condition = "list7_list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list3_sameQuestion <- cur_corr_no_autocorr * rHits_list7_list3_sameQuestion_mask
        rHits_list7_list3_sameQuestion_mtx <-
          rHits_list7_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7_list3",
                        source_mem = "rHit",
                        condition = "list7_list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list4_sameQuestion <- cur_corr_no_autocorr * rHits_list7_list4_sameQuestion_mask
        rHits_list7_list4_sameQuestion_mtx <-
          rHits_list7_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7_list4",
                        source_mem = "rHit",
                        condition = "list7_list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list5_sameQuestion <- cur_corr_no_autocorr * rHits_list7_list5_sameQuestion_mask
        rHits_list7_list5_sameQuestion_mtx <-
          rHits_list7_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7_list5",
                        source_mem = "rHit",
                        condition = "list7_list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list6_sameQuestion <- cur_corr_no_autocorr * rHits_list7_list6_sameQuestion_mask
        rHits_list7_list6_sameQuestion_mtx <-
          rHits_list7_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7_list6",
                        source_mem = "rHit",
                        condition = "list7_list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list8_sameQuestion <- cur_corr_no_autocorr * rHits_list7_list8_sameQuestion_mask
        rHits_list7_list8_sameQuestion_mtx <-
          rHits_list7_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7_list8",
                        source_mem = "rHit",
                        condition = "list7_list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list1_diffQuestion <- cur_corr_no_autocorr * rHits_list7_list1_diffQuestion_mask
        rHits_list7_list1_diffQuestion_mtx <-
          rHits_list7_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7_list1",
                        source_mem = "rHit",
                        condition = "list7_list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list2_diffQuestion <- cur_corr_no_autocorr * rHits_list7_list2_diffQuestion_mask
        rHits_list7_list2_diffQuestion_mtx <-
          rHits_list7_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7_list2",
                        source_mem = "rHit",
                        condition = "list7_list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list3_diffQuestion <- cur_corr_no_autocorr * rHits_list7_list3_diffQuestion_mask
        rHits_list7_list3_diffQuestion_mtx <-
          rHits_list7_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7_list3",
                        source_mem = "rHit",
                        condition = "list7_list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list4_diffQuestion <- cur_corr_no_autocorr * rHits_list7_list4_diffQuestion_mask
        rHits_list7_list4_diffQuestion_mtx <-
          rHits_list7_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7_list4",
                        source_mem = "rHit",
                        condition = "list7_list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list5_diffQuestion <- cur_corr_no_autocorr * rHits_list7_list5_diffQuestion_mask
        rHits_list7_list5_diffQuestion_mtx <-
          rHits_list7_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7_list5",
                        source_mem = "rHit",
                        condition = "list7_list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list6_diffQuestion <- cur_corr_no_autocorr * rHits_list7_list6_diffQuestion_mask
        rHits_list7_list6_diffQuestion_mtx <-
          rHits_list7_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7_list6",
                        source_mem = "rHit",
                        condition = "list7_list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list7_list8_diffQuestion <- cur_corr_no_autocorr * rHits_list7_list8_diffQuestion_mask
        rHits_list7_list8_diffQuestion_mtx <-
          rHits_list7_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7_list8",
                        source_mem = "rHit",
                        condition = "list7_list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        # list8
        rHits_list8_list1_sameQuestion <- cur_corr_no_autocorr * rHits_list8_list1_sameQuestion_mask
        rHits_list8_list1_sameQuestion_mtx <-
          rHits_list8_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8_list1",
                        source_mem = "rHit",
                        condition = "list8_list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list2_sameQuestion <- cur_corr_no_autocorr * rHits_list8_list2_sameQuestion_mask
        rHits_list8_list2_sameQuestion_mtx <-
          rHits_list8_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8_list2",
                        source_mem = "rHit",
                        condition = "list8_list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list3_sameQuestion <- cur_corr_no_autocorr * rHits_list8_list3_sameQuestion_mask
        rHits_list8_list3_sameQuestion_mtx <-
          rHits_list8_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8_list3",
                        source_mem = "rHit",
                        condition = "list8_list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list4_sameQuestion <- cur_corr_no_autocorr * rHits_list8_list4_sameQuestion_mask
        rHits_list8_list4_sameQuestion_mtx <-
          rHits_list8_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8_list4",
                        source_mem = "rHit",
                        condition = "list8_list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list5_sameQuestion <- cur_corr_no_autocorr * rHits_list8_list5_sameQuestion_mask
        rHits_list8_list5_sameQuestion_mtx <-
          rHits_list8_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8_list5",
                        source_mem = "rHit",
                        condition = "list8_list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list6_sameQuestion <- cur_corr_no_autocorr * rHits_list8_list6_sameQuestion_mask
        rHits_list8_list6_sameQuestion_mtx <-
          rHits_list8_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8_list6",
                        source_mem = "rHit",
                        condition = "list8_list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list7_sameQuestion <- cur_corr_no_autocorr * rHits_list8_list7_sameQuestion_mask
        rHits_list8_list7_sameQuestion_mtx <-
          rHits_list8_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8_list7",
                        source_mem = "rHit",
                        condition = "list8_list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list1_diffQuestion <- cur_corr_no_autocorr * rHits_list8_list1_diffQuestion_mask
        rHits_list8_list1_diffQuestion_mtx <-
          rHits_list8_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8_list1",
                        source_mem = "rHit",
                        condition = "list8_list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list2_diffQuestion <- cur_corr_no_autocorr * rHits_list8_list2_diffQuestion_mask
        rHits_list8_list2_diffQuestion_mtx <-
          rHits_list8_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8_list2",
                        source_mem = "rHit",
                        condition = "list8_list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list3_diffQuestion <- cur_corr_no_autocorr * rHits_list8_list3_diffQuestion_mask
        rHits_list8_list3_diffQuestion_mtx <-
          rHits_list8_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8_list3",
                        source_mem = "rHit",
                        condition = "list8_list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list4_diffQuestion <- cur_corr_no_autocorr * rHits_list8_list4_diffQuestion_mask
        rHits_list8_list4_diffQuestion_mtx <-
          rHits_list8_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8_list4",
                        source_mem = "rHit",
                        condition = "list8_list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list5_diffQuestion <- cur_corr_no_autocorr * rHits_list8_list5_diffQuestion_mask
        rHits_list8_list5_diffQuestion_mtx <-
          rHits_list8_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8_list5",
                        source_mem = "rHit",
                        condition = "list8_list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list6_diffQuestion <- cur_corr_no_autocorr * rHits_list8_list6_diffQuestion_mask
        rHits_list8_list6_diffQuestion_mtx <-
          rHits_list8_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8_list6",
                        source_mem = "rHit",
                        condition = "list8_list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        rHits_list8_list7_diffQuestion <- cur_corr_no_autocorr * rHits_list8_list7_diffQuestion_mask
        rHits_list8_list7_diffQuestion_mtx <-
          rHits_list8_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "rHit",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8_list7",
                        source_mem = "rHit",
                        condition = "list8_list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        #############################
        ### FHITS and ITEM MISSES ###
        #############################
        # fhits and misses: same/diff question and list
        fHitsANDmisses_sameList_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_sameList_sameQuestion_mask
        fHitsANDmisses_sameList_sameQuestion_mtx <-
          fHitsANDmisses_sameList_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "same",
                        source_mem = "fHitANDmiss",
                        condition = "sameList_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_sameList_sameQuestion_trialnums <- length(fHitsANDmisses_sameList_sameQuestion_mtx$r)

        fHitsANDmisses_diffList_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_diffList_sameQuestion_mask
        fHitsANDmisses_diffList_sameQuestion_mtx <-
          fHitsANDmisses_diffList_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "diff",
                        source_mem = "fHitANDmiss",
                        condition = "diffList_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_diffList_sameQuestion_trialnums <- length(fHitsANDmisses_diffList_sameQuestion_mtx$r)

        fHitsANDmisses_sameList_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_sameList_diffQuestion_mask
        fHitsANDmisses_sameList_diffQuestion_mtx <-
          fHitsANDmisses_sameList_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "same",
                        source_mem = "fHitANDmiss",
                        condition = "sameList_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_sameList_diffQuestion_trialnums <- length(fHitsANDmisses_sameList_diffQuestion_mtx$r)

        fHitsANDmisses_diffList_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_diffList_diffQuestion_mask
        fHitsANDmisses_diffList_diffQuestion_mtx <-
          fHitsANDmisses_diffList_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "diff",
                        source_mem = "fHitANDmiss",
                        condition = "diffList_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_diffList_diffQuestion_trialnums <- length(fHitsANDmisses_diffList_diffQuestion_mtx$r)

        # fhits and misses, split by list: same/diff question
        fHitsANDmisses_list1_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list1_sameQuestion_mask
        fHitsANDmisses_list1_sameQuestion_mtx <-
          fHitsANDmisses_list1_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list1",
                        source_mem = "fHitANDmiss",
                        condition = "list1_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list1_sameQuestion_trialnums <- length(fHitsANDmisses_list1_sameQuestion_mtx$r)

        fHitsANDmisses_list2_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list2_sameQuestion_mask
        fHitsANDmisses_list2_sameQuestion_mtx <-
          fHitsANDmisses_list2_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list2",
                        source_mem = "fHitANDmiss",
                        condition = "list2_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list2_sameQuestion_trialnums <- length(fHitsANDmisses_list2_sameQuestion_mtx$r)

        fHitsANDmisses_list3_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list3_sameQuestion_mask
        fHitsANDmisses_list3_sameQuestion_mtx <-
          fHitsANDmisses_list3_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list3",
                        source_mem = "fHitANDmiss",
                        condition = "list3_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list3_sameQuestion_trialnums <- length(fHitsANDmisses_list3_sameQuestion_mtx$r)

        fHitsANDmisses_list4_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list4_sameQuestion_mask
        fHitsANDmisses_list4_sameQuestion_mtx <-
          fHitsANDmisses_list4_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list4",
                        source_mem = "fHitANDmiss",
                        condition = "list4_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list4_sameQuestion_trialnums <- length(fHitsANDmisses_list4_sameQuestion_mtx$r)

        fHitsANDmisses_list5_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list5_sameQuestion_mask
        fHitsANDmisses_list5_sameQuestion_mtx <-
          fHitsANDmisses_list5_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list5",
                        source_mem = "fHitANDmiss",
                        condition = "list5_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list5_sameQuestion_trialnums <- length(fHitsANDmisses_list5_sameQuestion_mtx$r)

        fHitsANDmisses_list6_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list6_sameQuestion_mask
        fHitsANDmisses_list6_sameQuestion_mtx <-
          fHitsANDmisses_list6_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list6",
                        source_mem = "fHitANDmiss",
                        condition = "list6_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list6_sameQuestion_trialnums <- length(fHitsANDmisses_list6_sameQuestion_mtx$r)

        fHitsANDmisses_list7_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list7_sameQuestion_mask
        fHitsANDmisses_list7_sameQuestion_mtx <-
          fHitsANDmisses_list7_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list7",
                        source_mem = "fHitANDmiss",
                        condition = "list7_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list7_sameQuestion_trialnums <- length(fHitsANDmisses_list7_sameQuestion_mtx$r)

        fHitsANDmisses_list8_sameQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list8_sameQuestion_mask
        fHitsANDmisses_list8_sameQuestion_mtx <-
          fHitsANDmisses_list8_sameQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "same",
                        list = "list8",
                        source_mem = "fHitANDmiss",
                        condition = "list8_sameQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list8_sameQuestion_trialnums <- length(fHitsANDmisses_list8_sameQuestion_mtx$r)

        fHitsANDmisses_list1_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list1_diffQuestion_mask
        fHitsANDmisses_list1_diffQuestion_mtx <-
          fHitsANDmisses_list1_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list1",
                        source_mem = "fHitANDmiss",
                        condition = "list1_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list1_diffQuestion_trialnums <- length(fHitsANDmisses_list1_diffQuestion_mtx$r)

        fHitsANDmisses_list2_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list2_diffQuestion_mask
        fHitsANDmisses_list2_diffQuestion_mtx <-
          fHitsANDmisses_list2_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list2",
                        source_mem = "fHitANDmiss",
                        condition = "list2_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list2_diffQuestion_trialnums <- length(fHitsANDmisses_list2_diffQuestion_mtx$r)

        fHitsANDmisses_list3_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list3_diffQuestion_mask
        fHitsANDmisses_list3_diffQuestion_mtx <-
          fHitsANDmisses_list3_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list3",
                        source_mem = "fHitANDmiss",
                        condition = "list3_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list3_diffQuestion_trialnums <- length(fHitsANDmisses_list3_diffQuestion_mtx$r)

        fHitsANDmisses_list4_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list4_diffQuestion_mask
        fHitsANDmisses_list4_diffQuestion_mtx <-
          fHitsANDmisses_list4_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list4",
                        source_mem = "fHitANDmiss",
                        condition = "list4_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list4_diffQuestion_trialnums <- length(fHitsANDmisses_list4_diffQuestion_mtx$r)

        fHitsANDmisses_list5_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list5_diffQuestion_mask
        fHitsANDmisses_list5_diffQuestion_mtx <-
          fHitsANDmisses_list5_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list5",
                        source_mem = "fHitANDmiss",
                        condition = "list5_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list5_diffQuestion_trialnums <- length(fHitsANDmisses_list5_diffQuestion_mtx$r)

        fHitsANDmisses_list6_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list6_diffQuestion_mask
        fHitsANDmisses_list6_diffQuestion_mtx <-
          fHitsANDmisses_list6_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list6",
                        source_mem = "fHitANDmiss",
                        condition = "list6_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list6_diffQuestion_trialnums <- length(fHitsANDmisses_list6_diffQuestion_mtx$r)

        fHitsANDmisses_list7_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list7_diffQuestion_mask
        fHitsANDmisses_list7_diffQuestion_mtx <-
          fHitsANDmisses_list7_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list7",
                        source_mem = "fHitANDmiss",
                        condition = "list7_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list7_diffQuestion_trialnums <- length(fHitsANDmisses_list7_diffQuestion_mtx$r)

        fHitsANDmisses_list8_diffQuestion <- cur_corr_no_autocorr * fHitsANDmisses_list8_diffQuestion_mask
        fHitsANDmisses_list8_diffQuestion_mtx <-
          fHitsANDmisses_list8_diffQuestion %>%
          as.data.frame() %>%
          dplyr::mutate(subj = isubj,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        question_id = "unknown",
                        item_mem = "fHitANDmiss",
                        question_mem = "unknown",
                        list_mem = "unknown",
                        question = "diff",
                        list = "list8",
                        source_mem = "fHitANDmiss",
                        condition = "list8_diffQuestion",
                        row_name = rownames(.)) %>%
          gather_df(.) %>%
          filter_df(., cur_subj_info)

        fHitsANDmisses_list8_diffQuestion_trialnums <- length(fHitsANDmisses_list8_diffQuestion_mtx$r)

        ##########################
        ## Summarize for output ##
        ##########################
        # --- save out as trial pairs ---
        if(MIXED_MODELS_FLAG == 1){
          cat(sprintf("\nCreating mixed model dataframe for %s\n", isubj))

          join_col_names <- names(rANDsourceHits_sameList_sameQuestion_mtx)

          if(dim(mm_trial_pairs)[1] == 0){
            # nested join based on: http://stackoverflow.com/questions/32066402/how-to-perform-multiple-left-joins-using-dplyr-in-r
            mm_trial_pairs <- dplyr::full_join(rHits_sameList_sameQuestion_mtx, rHits_diffList_sameQuestion_mtx,
                                               by = join_col_names) %>%
              dplyr::full_join(., rHits_sameList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_diffList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_sameList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_diffList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_sameList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_diffList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., itemMiss_sameList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., itemMiss_diffList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., itemMiss_sameList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., itemMiss_diffList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_sameList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_diffList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_sameList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_diffList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              # off-diagonal elements
              dplyr::full_join(., rHits_list1_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              # include liberal scoring of question source memory
              dplyr::full_join(., rHits_questionMissLiberal_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_questionHitsLiberal_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., question1_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., question2_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., question3_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., question4_mtx,
                               by = join_col_names)
          } else {
            mm_trial_pairs <- dplyr::full_join(mm_trial_pairs, rHits_sameList_sameQuestion_mtx,
                                               by = join_col_names) %>%
              dplyr::full_join(., rHits_diffList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_sameList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_diffList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_sameList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_diffList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_sameList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_diffList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rANDsourceHits_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(mm_trial_pairs, itemMiss_sameList_sameQuestion_mtx,
                                 by = join_col_names) %>%
              dplyr::full_join(., itemMiss_diffList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., itemMiss_sameList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., itemMiss_diffList_diffQuestion_mtx,
                               by = join_col_names)  %>%
              dplyr::full_join(., rHits_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_sameList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_diffList_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_sameList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_diffList_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., fHitsANDmisses_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              # off-diagonal elements
              dplyr::full_join(., rHits_list1_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list1_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list2_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list3_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list4_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list5_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list6_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list8_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list7_list8_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list1_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list2_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list3_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list4_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list5_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list6_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list7_sameQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list1_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list2_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list3_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list4_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list5_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list6_diffQuestion_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_list8_list7_diffQuestion_mtx,
                               by = join_col_names) %>%
              # include liberal scoring of question source memory
              dplyr::full_join(., rHits_questionMissLiberal_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., rHits_questionHitsLiberal_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., question1_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., question2_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., question3_mtx,
                               by = join_col_names) %>%
              dplyr::full_join(., question4_mtx,
                               by = join_col_names)
          } #if(dim(mm_trial_pairs)[1] == 0

          mm_filepath <- paste0(analyzed_mri_dir,
                                halle::ensure_trailing_slash(isubj),
                                sprintf("mm_trial_pairs_%s_%s.RData",hemi_label,cur_roi))
          cat(sprintf("\nWriting single trial data to file: %s\n",mm_filepath))
          save(mm_trial_pairs,file=mm_filepath)

          mm_trial_pairs <- NULL
          print(Sys.time() - time_subj_start)

        } #if(MIXED_MODELS_FLAG == 1
      } #if(exist(cur_pattern_corr_file
    } #for (iroi
  } #for (idir

  #' # Subject timer
  print(sprintf("%s complete.", isubj))
  print(Sys.time() - time_subj_start)
} #for (isubj

#' # Print out R session info
sessionInfo()

#' # Print out script timer
print(sprintf("Complete."))
print(Sys.time() - time_script_start)
