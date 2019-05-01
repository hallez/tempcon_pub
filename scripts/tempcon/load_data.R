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
exclude_subjects <- c()
exclude_subjects_formatted <- NULL
for(i in 1:length(exclude_subjects)){
  exclude_subjects_formatted[i] <- halle::prepend_s_to_subject_id(exclude_subjects[i])
}
subjects <- subjects[!is.element(subjects, exclude_subjects_formatted)]

#' ### List included subjects
print(subjects)
length(subjects)

#' ## Set flags and other variables
SAVE_FLAG <- 1
num_enc_blocks <- 8

#' ## Setup columns of interest and colnames
# names are based on what MM wrote out as the headers for the group file in her `load_data_matlab.m` script and comparing against what's written out in `test_run*_combined.txt` files
# column numbers are based on what she subset
# (I feel better that she pulled the same columns from each logfile; makes it seem like Presentation output has a similar format. But, admittedly, this is still a gamble)
# I added in reading column 8 (`scored`)
colnames_from_MM <- c("listnum","oldnew","enctask","jitter","objnum","scored","response")
colnums_from_MM <- c(3:9)
# columns 3:9 names are based on what's in `load_data_matlab.m`
# column 12 column name is inferred from `onsets_rkn_hitmiss_0dur.m`
item_recog_colnames_from_MM <- c("listnum","oldnew","enctask","jitter","objnum","scored","response","trial_onset_time")
item_recog_colnums_from_MM <- c(3:9,12)

#+ label="Read in subject data"
#' # Read in subject data
#+ warning = FALSE
for(isubj in subjects){
  cat(sprintf("\nWorking on %s\n", isubj))

  # reset run numbers for wonky subjects
  if(isubj == "s15" | isubj == "s19"){
    item_recog_runs <- 2:7
  } else if(isubj == "s23"){
    item_recog_runs <- c(2,3,6:9)
  } else if(isubj %in% c("s29", "s4")){
    enc_blocks <- 2:9
  } else {
    enc_blocks <- 1:8
    item_recog_runs <- 1:6
    temporal_source_runs <- 1:4
    question_source_runs <- 1:4
  }

  # setup directories
  # deal w/ weird naming for backroom subjects' raw data directories
  if(isubj == "s11" | isubj == "s12"){
    subj_raw_dir <- paste0(raw_behavioral_dir,halle::ensure_trailing_slash(paste0(isubj, "_backroom")))
  } else {
    subj_raw_dir <- paste0(raw_behavioral_dir,halle::ensure_trailing_slash(isubj))
  }

  # create output directory (if doesn't already exist)
  subj_analyzed_dir <- paste0(analyzed_behavioral_dir, halle::ensure_trailing_slash(isubj))
  if(!file.exists(subj_analyzed_dir)){
    cat(sprintf("Creating %s\n", subj_analyzed_dir))
    dir.create(subj_analyzed_dir, recursive = TRUE)
  }

  #' ## --- Encoding ---
  enc_block_counter <- 0

  for(ienc in 1:num_enc_blocks){
    enc_block_counter <- enc_block_counter + 1

    cur_block <- enc_blocks[ienc]

    if(isubj == "s14" & cur_block == 8){
      enc_fpath <- file.path(subj_raw_dir, sprintf("%s-enc_%d_backroom3.log", isubj, cur_block))
    } else if(isubj == "s17" & cur_block == 1){
      enc_fpath <- file.path(subj_raw_dir, sprintf("%s-enc_%d_backroom1.log", isubj, cur_block))
    } else if(isubj == "s22" & cur_block == 3){
      enc_fpath <- file.path(subj_raw_dir, sprintf("%s-enc_%d_backroom..log", isubj, cur_block))
    } else if(isubj == "s22" & cur_block == 5){
      enc_fpath <- file.path(subj_raw_dir, sprintf("%s-enc_%d_backroom...log", isubj, cur_block))
    } else if(isubj == "s22" & cur_block == 6){
      enc_fpath <- file.path(subj_raw_dir, sprintf("%s-enc_%d_backroom..log", isubj, cur_block))
    } else if(isubj == "s25" & cur_block == 1){
      enc_fpath <- file.path(subj_raw_dir, sprintf("%s-enc_%d..log", isubj, cur_block))
    } else if(is.element(isubj, backroom_subjects_formatted) == TRUE) {
      enc_fpath <- file.path(subj_raw_dir, sprintf("%s-enc_%d_backroom.log", isubj, cur_block))
    } else {
      enc_fpath <- file.path(subj_raw_dir, sprintf("%s-enc_%d.log", isubj, cur_block))
    }

    if(file.exists(enc_fpath)){
      enc_start_line <- grep("^Event",readLines(enc_fpath))+1

      cur_enc <- read.table(enc_fpath,
                            header = FALSE,
                            skip = enc_start_line,
                            fill = TRUE,
                            stringsAsFactors = FALSE)

      cur_enc_tidy <- cur_enc %>%
        # these are based on looking at `s1-enc_1.log`
        # start by just duplicating the column names in the file
        dplyr::rename("type" = "V3",
                      "response" = "V4") %>%
        # add trial order (this is based on the order of the rows)
        # based on: https://stackoverflow.com/questions/36893957/dplyr-generate-row-number-row-position-in-group-by
        dplyr::mutate(enc_trial_id = 1:n(),
                      subj = isubj,
                      listnum = ienc) %>%
        dplyr::select(-starts_with("V")) %>%
        # these are based on comparing with `enc1.txt` which was created by `make_stimuli.m`
        dplyr::rename("enctask" = "type",
                      "objnum" = "response")

        if(enc_block_counter==1){
          # create this variable if this is the first iteration through the loop
          enc_dat <- cur_enc_tidy
        } else {
          # merge together
          enc_dat <- dplyr::full_join(enc_dat, cur_enc_tidy, by = intersect(names(enc_dat), names(cur_enc_tidy)))
        }
      }else {
        print(sprintf("Encoding file %s not found.", enc_fpath))
      }# file.exists
  } #ienc

  #' ## Check the datafile
  if(!exists("enc_dat")){
    message(sprintf("WARNING: %s HAS NO ENCODING DATA", isubj))
  } else if (is.null(enc_dat)){
    message(sprintf("WARNING: %s HAS NO ENCODING DATA", isubj))
  } else if(dim(enc_dat)[1] != 8*36){
    print(sprintf("WARNING: %s does not have the correct number of encoding trials.", isubj))
  }

  #' ## Save out encoding data
  if(SAVE_FLAG == 1){
    fname_out <- paste0(subj_analyzed_dir,paste0(isubj,"_encoding.Rdata"))
    save(enc_dat, file = fname_out)
  }

  #' ## ---Item recognition---
  recog_run_counter <- 0
  for(irun in item_recog_runs){

    recog_run_counter <- recog_run_counter + 1

    if(isubj == "s11" | isubj == "s12"){
      # since these files were generated in the scanner, `backroom` won't be added onto the filenames themselves
      fname <- paste0(subj_raw_dir,paste0(isubj,"-test_",irun,".log"))
    } else {
      fname <- paste0(subj_raw_dir,paste0(isubj,"-test_",irun,".log"))
    }

    if(file.exists(fname)){
      # read in the data
      # we know the data we want is always 2 lines below the header
      start_line <- grep("^Event",readLines(fname))+1
      cur_recog_data <- read.table(fname,
                             header = FALSE,
                             skip = start_line,
                             fill = TRUE,
                             stringsAsFactors = FALSE)

      head(cur_recog_data)

      # subset the data based on how MM does it in `load_data_matlab.m`
      cur_recog_data_subset <- cur_recog_data[,item_recog_colnums_from_MM]
      head(cur_recog_data_subset)

      # add the column labels to the columns of interest
      # again, we're totally trusting that MM has pulled the right columns and labelled them correctly
      colnames(cur_recog_data_subset) <- unlist(item_recog_colnames_from_MM)

      # add in columns for run number, subject ID, and row number (this will be how we determine trial number)
      cur_recog_data_subset$run <- recog_run_counter
      cur_recog_data_subset$subj <- isubj
      cur_recog_data_subset$row_number <- as.numeric(rownames(cur_recog_data_subset))

      head(cur_recog_data_subset)

      if(recog_run_counter == 1){
        # create this variable if this is the first iteration through the loop
        item_recog_dat <- NULL
        item_recog_dat <- cur_recog_data_subset
      } else {
        # merge together
        # NB: this will throw a warning about joining factors w different levels
        item_recog_dat <- dplyr::full_join(item_recog_dat, cur_recog_data_subset, by = intersect(names(item_recog_dat), names(cur_recog_data_subset)))
      }

    } else {
      print(sprintf("Item recog file %s not found.", fname))
    }#if(file.exists)
  } #irun

  #' ## Check the datafile
  if(!exists("item_recog_dat")){
    message(sprintf("WARNING: %s HAS NO ITEM RECOG DATA", isubj))
  } else if (is.null(item_recog_dat)){
    message(sprintf("WARNING: %s HAS NO ITEM RECOG DATA", isubj))
  } else if(dim(item_recog_dat)[1] != 6*62){
    print(sprintf("WARNING: %s does not have the correct number of item recognition trials.", isubj))
  }

  #' ## Save out item recog data
  if(SAVE_FLAG == 1){
    fname_out <- paste0(subj_analyzed_dir,paste0(isubj,"_item_recog.Rdata"))
    save(item_recog_dat, file = fname_out)
  }

  #' ## ---Temporal source recognition---
  temp_run_counter <- 0
  for(irun in temporal_source_runs){
    temp_run_counter <- temp_run_counter + 1

    if(isubj == "s11" | isubj == "s12") {
      fname <- paste0(subj_raw_dir,paste0(isubj,"-tempcon",irun,"_backroom.log"))
    } else if(is.element(isubj, backroom_subjects_formatted) == TRUE) {
      fname <- paste0(subj_raw_dir,paste0(isubj,"-tempcon",irun,"_backroom.log"))
    } else {
      fname <- paste0(subj_raw_dir,paste0(isubj,"-tempcon",irun,".log"))
    }

    if(file.exists(fname)){
      # read in the data
      # we know the data we want is always 2 lines below the header
      start_line <- grep("^Event",readLines(fname))+1
      cur_temp_data <- read.table(fname,
                             header = FALSE,
                             skip = start_line,
                             fill = TRUE,
                             stringsAsFactors = FALSE)

      head(cur_temp_data)

      # subset the data
      cur_temp_data_subset <- cur_temp_data[,colnums_from_MM]
      head(cur_temp_data_subset)

      # add column names
      colnames(cur_temp_data_subset) <- unlist(colnames_from_MM)

      # add in columns for run number, subject ID, and row number (this will be how we determine trial number)
      cur_temp_data_subset$run <- temp_run_counter
      cur_temp_data_subset$subj <- isubj
      cur_temp_data_subset$row_number <- as.numeric(rownames(cur_temp_data_subset))

      head(cur_temp_data_subset)

      if(temp_run_counter == 1){
        # create this variable if this is the first iteration through the loop
        temporal_source_dat <- NULL
        temporal_source_dat <- cur_temp_data_subset
      } else {
        # merge together
        # NB: this will throw a warning about joining factors w different levels
        temporal_source_dat <- dplyr::full_join(temporal_source_dat, cur_temp_data_subset, by = intersect(names(temporal_source_dat), names(cur_temp_data_subset)))
      }
    } else {
      print(sprintf("Temporal source file %s not found.", fname))
    }#if(file.exists)
  } #irun

  #' ## Check the datafile
  if(!exists("temporal_source_dat")){
    message(sprintf("WARNING: %s HAS NO TEMPORAL SOURCE MEMORY DATA", isubj))
  } else if (is.null(temporal_source_dat)){
    message(sprintf("WARNING: %s HAS NO TEMPORAL SOURCE MEMORY DATA", isubj))
  } else if(dim(temporal_source_dat)[1] != 4*72){
    print(sprintf("WARNING: %s does not have the correct number of temporal source memory trials.", isubj))
  }

  #' ## Save out temporal source memory data
  if(SAVE_FLAG == 1){
    fname_out <- paste0(subj_analyzed_dir,paste0(isubj,"_temporal_source_mem.Rdata"))
    save(temporal_source_dat, file = fname_out)
  }

  #' ## ---Question source recognition---
  question_run_counter <- 0
  for(irun in question_source_runs){
    question_run_counter <- question_run_counter + 1

    if(isubj == "s11" | isubj == "s12") {
      fname <- paste0(subj_raw_dir,paste0(isubj,"-assoc",irun,"_backroom.log"))
    } else if(is.element(isubj, backroom_subjects_formatted) == TRUE) {
      fname <- paste0(subj_raw_dir,paste0(isubj,"-assoc",irun,"_backroom.log"))
    } else {
      fname <- paste0(subj_raw_dir,paste0(isubj,"-assoc",irun,".log"))
    }

    if(file.exists(fname)){
      # read in the data
      # we know the data we want is always 2 lines below the header
      start_line <- grep("^Event",readLines(fname))+1
      cur_question_data <- read.table(fname,
                             header = FALSE,
                             skip = start_line,
                             fill = TRUE,
                             stringsAsFactors = FALSE)

      head(cur_question_data)

      # subset the data
      cur_question_data_subset <- cur_question_data[,colnums_from_MM]
      head(cur_question_data_subset)

      # add column names
      colnames(cur_question_data_subset) <- unlist(colnames_from_MM)

      # add in columns for run number, subject ID, and row number (this will be how we determine trial number)
      cur_question_data_subset$run <- question_run_counter
      cur_question_data_subset$subj <- isubj
      cur_question_data_subset$row_number <- as.numeric(rownames(cur_question_data_subset))

      head(cur_question_data_subset)

      if(question_run_counter == 1){
        # create this variable if this is the first iteration through the loop
        question_source_dat <- NULL
        question_source_dat <- cur_question_data_subset
      } else {
        # merge together
        # NB: this will throw a warning about joining factors w different levels
        question_source_dat <- dplyr::full_join(question_source_dat, cur_question_data_subset, by = intersect(names(question_source_dat), names(cur_question_data_subset)))
      }

    } else {
      print(sprintf("Question source file %s not found.", fname))
    }#if(file.exists)
  } #irun

  #' ## Check the datafile
  if(!exists("question_source_dat")){
    message(sprintf("WARNING: %s HAS NO QUESTION SOURCE MEMORY DATA", isubj))
  } else if (is.null(question_source_dat)){
    message(sprintf("WARNING: %s HAS NO QUESTION SOURCE MEMORY DATA", isubj))
  } else if(dim(question_source_dat)[1] != 4*72){
    print(sprintf("WARNING: %s does not have the correct number of question source memory trials.", isubj))
  }

  #' ## Save out question source memory
  if(SAVE_FLAG == 1){
    fname_out <- paste0(subj_analyzed_dir,paste0(isubj,"_question_source_mem.Rdata"))
    save(question_source_dat, file = fname_out)
  }

  # clear out between participants so don't accidentally save out previous subjects' data if one subject is missing data
  enc_dat <- NULL
  item_recog_dat <- NULL
  temporal_source_dat <- NULL
  question_source_dat <- NULL

} #isubj

#' # Print out R session info
sessionInfo()
