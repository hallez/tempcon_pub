#' ---
#' title: Stats for Tempcon behavioral data
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
analyzed_mri_dir <- halle::ensure_trailing_slash(config$directories$analyzed_mri)
raw_behavioral_dir <- halle::ensure_trailing_slash(config$directories$raw_behavioral)
analyzed_behavioral_dir <- halle::ensure_trailing_slash(config$directories$analyzed_behavioral)
dropbox_dir <- halle::ensure_trailing_slash(config$directories$dropbox_tempcon)
graphs_dir <- paste0(halle::ensure_trailing_slash(analyzed_behavioral_dir),
                     halle::ensure_trailing_slash("plots"))

if(!file.exists(graphs_dir)){
  cat(sprintf("Creating graphs_dir\n"))
  dir.create(graphs_dir, recursive = TRUE)
}

#' ## Define other variables
SAVE_FLAG <- 1
GRAPH_FLAG <- 1
num_subj <- 30

#' # Load data
# this file is created by `behav-tidy.R`
load(file.path(analyzed_behavioral_dir, sprintf("all_subj_behav_data_N%d.RData", num_subj)))

#' ## Load excluded subjects
# this file gets created by `rsa-tidy-data-btwn-runs.R`
rsa_exclude_subjects_fpath <- file.path(analyzed_mri_dir, "excluded_subjects_variables.RData")
if(file.exists(rsa_exclude_subjects_fpath)) {
  load(file = rsa_exclude_subjects_fpath)
} else {
  exclude_subjects <- c("s1", "s2", "s3", "s11", "s13", "s15", "s24", "s26")
}

message(sprintf("Excluded subjects: %s", paste0(exclude_subjects, collapse = ', ')))

#' ## Read in subject information files
# These are needed to drop behavioral trials for which we drop the imaging data
#+ warning = FALSE
all_subj_info <- data.frame()

# assume that subjects have folders that start with `s` and are followed by one or two digits
subjects <- c(list.files(path = analyzed_behavioral_dir, pattern = "^[s][(123456789)]$|^[s][(123456789)][(0123456789)]$"))

for(isubj in subjects){
  cur_fpath <- file.path(analyzed_behavioral_dir, isubj, "trial_information.csv")

  if(!file.exists(cur_fpath)){
    message(sprintf("trial information file %s does not exist - skipping\n", cur_fpath))
    next
  }

  cur_info_file <- read.csv(cur_fpath)

  if(dim(all_subj_info)[1] == 0){
    all_subj_info <- cur_info_file
  } else {
    # NB: will get error about merging `subj` column (vector vs factor) - can ignore
    all_subj_info <- dplyr::full_join(all_subj_info, cur_info_file, by = names(all_subj_info))
  }

} #for(isubj

#' ## Merge behavioral data with trial information files
all_subj_with_info <- all_subj %>%
  dplyr::left_join(., all_subj_info, by = intersect(names(all_subj), names(all_subj_info)))

#' ## Filter out excluded subjects and excluded trials
all_subj_filt <- all_subj_with_info %>%
  dplyr::filter(!subj %in% exclude_subjects) %>%
  dplyr::filter(exclude_behavioral == 0)

#' ## Take a quick look at the data
summary(all_subj_filt)

#' ## Report current subjects
unique(all_subj_filt$subj)
length(unique(all_subj_filt$subj))

#' # Response rates
#' ## Item recognition
itemrec_denom <- all_subj_filt %>%
  dplyr::group_by(subj, oldnew_label) %>%
  dplyr::summarise(numresp = length(item_recog_scored_rhit_collapsed))

itemrec_rates <- all_subj_filt %>%
  dplyr::group_by(subj, oldnew_label, item_recog_scored_rhit_collapsed) %>%
  dplyr::summarise(respcount = length(item_recog_scored_rhit_collapsed)) %>%
  dplyr::left_join(itemrec_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

itemrec_rates %>%
  dplyr::group_by(oldnew_label, item_recog_scored_rhit_collapsed) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

save("itemrec_rates", file = file.path(analyzed_behavioral_dir, "itemrec_rates.RData"))

#' ## Question source memory
questmem_denom <- all_subj_filt %>%
  dplyr::group_by(subj) %>%
  # filter out non responses so that rates will add up to 1
  dplyr::filter(question_scored_exact %in% c("hit", "incorrect")) %>%
  dplyr::summarise(numresp = length(question_scored_exact))

questmem_rates <- all_subj_filt %>%
  dplyr::group_by(subj, question_scored_exact) %>%
  dplyr::filter(question_scored_exact %in% c("hit", "incorrect")) %>%
  dplyr::summarise(respcount = length(question_scored_exact)) %>%
  dplyr::left_join(questmem_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

questmem_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(question_scored_exact) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

#' ## Question source memory by item recognition
questmem_denom <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed) %>%
  # filter out non responses so that rates will add up to 1
  dplyr::filter(question_scored_exact %in% c("hit", "incorrect")) %>%
  dplyr::summarise(numresp = length(question_scored_exact))

questmem_rates <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, question_scored_exact) %>%
  dplyr::filter(question_scored_exact %in% c("hit", "incorrect")) %>%
  dplyr::summarise(respcount = length(question_scored_exact)) %>%
  dplyr::left_join(questmem_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

questmem_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(item_recog_scored_rhit_collapsed, question_scored_exact) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

save("questmem_rates", file = file.path(analyzed_behavioral_dir, "questmem_by_itemmem.RData"))

#' ## List source memory (exact list)
listmem_exact_denom <- all_subj_filt %>%
  dplyr::group_by(subj) %>%
  # filter out non responses so that rates will add up to 1
  dplyr::filter(temporal_scored_exact %in% c("correct", "incorrect")) %>%
  dplyr::summarise(numresp = length(temporal_scored_exact))

listmem_exact_rates <- all_subj_filt %>%
  dplyr::group_by(subj, temporal_scored_exact) %>%
  dplyr::filter(temporal_scored_exact %in% c("correct", "incorrect")) %>%
  dplyr::summarise(respcount = length(temporal_scored_exact)) %>%
  dplyr::left_join(listmem_exact_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

listmem_exact_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(temporal_scored_exact) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

#' ## List source memory (exact list) by item recognition
listmem_exact_denom <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed) %>%
  # filter out non responses so that rates will add up to 1
  dplyr::filter(temporal_scored_exact %in% c("correct", "incorrect")) %>%
  dplyr::summarise(numresp = length(temporal_scored_exact))

listmem_exact_rates <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_exact) %>%
  dplyr::filter(temporal_scored_exact %in% c("correct", "incorrect")) %>%
  dplyr::summarise(respcount = length(temporal_scored_exact)) %>%
  dplyr::left_join(listmem_exact_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

listmem_exact_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(item_recog_scored_rhit_collapsed, temporal_scored_exact) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

save("listmem_exact_rates", file = file.path(analyzed_behavioral_dir, "listmem_exact_by_itemmem.RData"))

#' ## List source memory (exact list) by item recognition and question memory
questmemXlistmem_exact_denom <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, question_scored_exact) %>%
  # filter out non responses so that rates will add up to 1
  dplyr::filter(temporal_scored_exact %in% c("correct", "incorrect"),
                question_scored_exact %in% c("hit", "incorrect")) %>%
  dplyr::summarise(numresp = length(temporal_scored_exact))

questmemXlistmem_exact_rates <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_exact, question_scored_exact) %>%
  dplyr::filter(temporal_scored_exact %in% c("correct", "incorrect"),
                question_scored_exact %in% c("hit", "incorrect")) %>%
  dplyr::summarise(respcount = length(temporal_scored_exact)) %>%
  dplyr::left_join(questmemXlistmem_exact_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

questmemXlistmem_exact_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(item_recog_scored_rhit_collapsed, temporal_scored_exact, question_scored_exact) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

save("questmemXlistmem_exact_rates", file = file.path(analyzed_behavioral_dir, "questmemXlistmem_exact_by_itemmem.RData"))

#' ## List source memory (list half)
listmem_denom <- all_subj_filt %>%
  dplyr::group_by(subj) %>%
  # filter out non responses so that rates will add up to 1
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half")) %>%
  dplyr::summarise(numresp = length(temporal_scored_list_samediff_half))

listmem_rates <- all_subj_filt %>%
  dplyr::group_by(subj, temporal_scored_list_samediff_half) %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half")) %>%
  dplyr::summarise(respcount = length(temporal_scored_list_samediff_half)) %>%
  dplyr::left_join(listmem_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

listmem_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(temporal_scored_list_samediff_half) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

#' ## List source memory (list half) by item recognition
listmem_denom <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed) %>%
  # filter out non responses so that rates will add up to 1
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half")) %>%
  dplyr::summarise(numresp = length(temporal_scored_list_samediff_half))

listmem_rates <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_list_samediff_half) %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half")) %>%
  dplyr::summarise(respcount = length(temporal_scored_list_samediff_half)) %>%
  dplyr::left_join(listmem_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

listmem_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(item_recog_scored_rhit_collapsed, temporal_scored_list_samediff_half) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

save("listmem_rates", file = file.path(analyzed_behavioral_dir, "listmem_half_by_itemmem.RData"))

#' ## List source memory (list half) by item recognition and question memory
questmemXlistmem_denom <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, question_scored_exact) %>%
  # filter out non responses so that rates will add up to 1
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half"),
                question_scored_exact %in% c("hit", "incorrect")) %>%
  dplyr::summarise(numresp = length(temporal_scored_list_samediff_half))

questmemXlistmem_rates <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_list_samediff_half, question_scored_exact) %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half"),
                question_scored_exact %in% c("hit", "incorrect")) %>%
  dplyr::summarise(respcount = length(temporal_scored_list_samediff_half)) %>%
  dplyr::left_join(questmemXlistmem_denom) %>%
  dplyr::mutate(resprate = respcount / numresp)

questmemXlistmem_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(item_recog_scored_rhit_collapsed, temporal_scored_list_samediff_half, question_scored_exact) %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate)) %>%
  kableExtra::kable()

save("questmemXlistmem_rates", file = file.path(analyzed_behavioral_dir, "questmemXlistmem_half_by_itemmem.RData"))

#' # Item recollection by list half (first/second)
#' ## Graph counts
all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_list_firstsecond_half) %>%
  dplyr::summarise(resp_count = n()) %>%
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit", "wrong_half"),
                !is.na(item_recog_scored_rhit_collapsed)) %>%
  ggplot2::ggplot(ggplot2::aes(x = item_recog_scored_rhit_collapsed, y = resp_count, fill = temporal_scored_list_firstsecond_half)) +
  ggplot2::geom_boxplot()

#' ## Stats (counts)
all_subj_filt %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_list_firstsecond_half) %>%
  dplyr::summarise(resp_count = n()) %>%
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit"),
                !is.na(item_recog_scored_rhit_collapsed),
                # since trying to get this analysis to parallel the RSA, actually only care about rhits
                item_recog_scored_rhit_collapsed == "rhit") %>%
  ez::ezANOVA(data = .,
              dv = .(resp_count),
              wid = .(subj),
              within = .(temporal_scored_list_firstsecond_half),
              type = 2)

#' ## Graph proportions
itemrec_firstsecond_respcounts <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed) %>%
  dplyr::filter(item_recog_scored_rhit_collapsed == "rhit") %>%
  # filter out nonsensical temporal responses so that proportions sum to 1
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit", "wrong_half")) %>%
  dplyr::summarise(rhit_resp = length(item_recog_scored_rhit_collapsed))

itemrec_firstsecond_resprate <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_list_firstsecond_half) %>%
  dplyr::filter(item_recog_scored_rhit_collapsed == "rhit") %>%
  dplyr::summarise(list_resp = length(temporal_scored_list_firstsecond_half)) %>%
  dplyr::left_join(itemrec_firstsecond_respcounts) %>%
  dplyr::group_by(subj) %>%
  dplyr::mutate(resprate = list_resp / rhit_resp)

itemrec_firstsecond_resprate %>%
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit", "wrong_half")) %>%
  ggplot2::ggplot(ggplot2::aes(x = item_recog_scored_rhit_collapsed, y = resprate, fill = temporal_scored_list_firstsecond_half)) +
  ggplot2::geom_boxplot()

ggplot2::ggsave(filename = file.path(graphs_dir,
                                     sprintf("list-firstsecond-half-mem_for-rhits_N%d.pdf", length(unique(all_subj_filt$subj)))),
                width = 6, height = 6)

#' ### Sanity check ourselves
# are there subjects w/ more than 360 responses?
itemrec_firstsecond_respcounts %>%
  dplyr::filter(rhit_resp > 360)

# do the means add up to 1?
itemrec_firstsecond_resprate %>%
  dplyr::group_by(subj) %>%
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit", "wrong_half")) %>%
  dplyr::summarise(total_resprate = sum(resprate))

#' ## Stats (proportions)
# is there a difference between 1st/2nd half memory
# no need to test for hits vs misses b/c that's accounted for in the item recollection by list half (same/diff) analyses
itemrec_firstsecond_resprate %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit"),
                !is.na(item_recog_scored_rhit_collapsed),
                # since trying to get this analysis to parallel the RSA, actually only care about rhits
                item_recog_scored_rhit_collapsed == "rhit") %>%
  ez::ezANOVA(data = .,
              dv = .(resprate),
              wid = .(subj),
              within = .(temporal_scored_list_firstsecond_half),
              type = 2)

# does it make more sense to use a paired t-test?
itemrec_firstsecond_tt <- itemrec_firstsecond_resprate %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit"),
                !is.na(item_recog_scored_rhit_collapsed),
                # since trying to get this analysis to parallel the RSA, actually only care about rhits
                item_recog_scored_rhit_collapsed == "rhit") %>%
  # drop non-unique values so that when spread don't get NAs
  dplyr::select(-list_resp) %>%
  tidyr::spread(key = temporal_scored_list_firstsecond_half, value = resprate)

t.test(itemrec_firstsecond_tt$first_half_hit, itemrec_firstsecond_tt$second_half_hit, paired = TRUE)

#' ## Half bias (proportion)
itemrec_firstsecond_bias_respcounts <- all_subj_filt %>%
  dplyr::filter(item_recog_scored_rhit_collapsed == "rhit") %>%
  dplyr::group_by(subj, temporal_response_as_list_firstsecond_half) %>%
  dplyr::summarise(resp_count = n()) %>%
  # filter out nonsensical temporal responses so that proportions sum to 1
  dplyr::filter(temporal_response_as_list_firstsecond_half %in% c("first_half", "second_half"))

itemrec_firstsecond_bias_resprate <- all_subj_filt %>%
  dplyr::filter(item_recog_scored_rhit_collapsed == "rhit") %>%
  # filter out nonsensical temporal responses so that proportions sum to 1
  dplyr::filter(temporal_response_as_list_firstsecond_half %in% c("first_half", "second_half")) %>%
  dplyr::group_by(subj, temporal_scored_list_firstsecond_half) %>%
  dplyr::summarise(list_resp = length(temporal_scored_list_firstsecond_half)) %>%
  dplyr::left_join(itemrec_firstsecond_bias_respcounts) %>%
  dplyr::group_by(subj) %>%
  dplyr::mutate(resprate = list_resp / resp_count)

itemrec_firstsecond_bias_resprate %>%
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit")) %>%
  ggplot2::ggplot(ggplot2::aes(x = temporal_scored_list_firstsecond_half, y = resprate, fill = temporal_scored_list_firstsecond_half)) +
  ggplot2::geom_boxplot() +
  ggplot2::ylab("temporal half hits / all responses by half") +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(), legend.position = "none")

ggplot2::ggsave(filename = file.path(graphs_dir,
                                     sprintf("list-firstsecond-bias-half-mem_for-rhits_N%d.pdf", length(unique(all_subj_filt$subj)))),
                width = 6, height = 6)

#' # Item recollection by list half (same/diff)
#' ## Graph (counts)
all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_list_samediff_half) %>%
  dplyr::summarise(resp_count = n()) %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half"),
                !is.na(item_recog_scored_rhit_collapsed)) %>%
  ggplot2::ggplot(ggplot2::aes(x = item_recog_scored_rhit_collapsed, y = resp_count, fill = temporal_scored_list_samediff_half)) +
  ggplot2::geom_boxplot()

#' ## Stats (counts)
all_subj_filt %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_list_samediff_half) %>%
  dplyr::summarise(resp_count = n()) %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half"),
                !is.na(item_recog_scored_rhit_collapsed),
                # since trying to get this analysis to parallel the RSA, actually only care about rhits
                item_recog_scored_rhit_collapsed == "rhit") %>%
  ez::ezANOVA(data = .,
              dv = .(resp_count),
              wid = .(subj),
              within = .(temporal_scored_list_samediff_half),
              type = 2)

#' ## Graph proportions
itemrec_samediff_respcounts <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed) %>%
  dplyr::filter(item_recog_scored_rhit_collapsed == "rhit") %>%
  # filter out nonsensical temporal responses so that proportions sum to 1
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half")) %>%
  dplyr::summarise(rhit_resp = length(item_recog_scored_rhit_collapsed))

itemrec_samediff_resprate <- all_subj_filt %>%
  dplyr::group_by(subj, item_recog_scored_rhit_collapsed, temporal_scored_list_samediff_half) %>%
  dplyr::filter(item_recog_scored_rhit_collapsed == "rhit") %>%
  dplyr::summarise(list_resp = length(temporal_scored_list_samediff_half)) %>%
  dplyr::left_join(itemrec_samediff_respcounts) %>%
  dplyr::group_by(subj) %>%
  dplyr::mutate(resprate = list_resp / rhit_resp)

itemrec_samediff_resprate %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half")) %>%
  ggplot2::ggplot(ggplot2::aes(x = item_recog_scored_rhit_collapsed, y = resprate, fill = temporal_scored_list_samediff_half)) +
  ggplot2::geom_boxplot()

ggplot2::ggsave(filename = file.path(graphs_dir,
                                     sprintf("list-samediff-half-mem_for-rhits_N%d.pdf", length(unique(all_subj_filt$subj)))),
                width = 6, height = 6)

#' ### Sanity check ourselves
# are there subjects w/ more than 360 responses?
itemrec_samediff_respcounts %>%
  dplyr::filter(rhit_resp > 360)

# do the means add up to 1?
itemrec_samediff_resprate %>%
  dplyr::group_by(subj) %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half")) %>%
  dplyr::summarise(total_resprate = sum(resprate))

#' ## Stats (proportions)
itemrec_samediff_resprate %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half"),
                !is.na(item_recog_scored_rhit_collapsed),
                # since trying to get this analysis to parallel the RSA, actually only care about rhits
                item_recog_scored_rhit_collapsed == "rhit") %>%
  ez::ezANOVA(data = .,
              dv = .(resprate),
              wid = .(subj),
              within = .(temporal_scored_list_samediff_half),
              type = 2)

# does it make more sense to use a paired t-test?
itemrec_samediff_tt <- itemrec_samediff_resprate %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::filter(temporal_scored_list_samediff_half %in% c("same_half_hit", "wrong_half"),
                !is.na(item_recog_scored_rhit_collapsed),
                # since trying to get this analysis to parallel the RSA, actually only care about rhits
                item_recog_scored_rhit_collapsed == "rhit") %>%
  # drop non-unique values so that when spread don't get NAs
  dplyr::select(-list_resp) %>%
  tidyr::spread(key = temporal_scored_list_samediff_half, value = resprate)

t.test(itemrec_samediff_tt$same_half_hit, itemrec_samediff_tt$wrong_half, paired = TRUE)

#' # Temporal memory accuracy by list half (first/second)
#' ## Graph
# first, figure out how many responses each person made (to get their denominator)
firstsecond_all_resp <- all_subj_filt %>%
  dplyr::filter(temporal_source_response %in% c(1:8)) %>%
  dplyr::group_by(subj) %>%
  dplyr::summarise(allresp = length(temporal_scored_list_firstsecond_half))

firstsecond_rates <- all_subj_filt %>%
  dplyr::filter(temporal_source_response %in% c(1:8)) %>%
  dplyr::group_by(subj, temporal_scored_list_firstsecond_half) %>%
  dplyr::summarise(respcount = length(temporal_scored_list_firstsecond_half)) %>%
  dplyr::left_join(firstsecond_all_resp) %>%
  dplyr::mutate(resprate = respcount / allresp)

firstsecond_rates %>%
  ggplot2::ggplot(ggplot2::aes(x = temporal_scored_list_firstsecond_half, y = resprate)) +
  ggplot2::geom_boxplot()

#' ## Stats
# somehow this makes more sense to me as an anova than as a t-test - does it matter?
firstsecond_rates %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::filter(temporal_scored_list_firstsecond_half %in% c("first_half_hit", "second_half_hit")) %>%
  ez::ezANOVA(data = .,
            dv = .(resprate),
            wid = .(subj),
            within = .(temporal_scored_list_firstsecond_half),
            type = 2)

#' # Temporal memory accuracy by list half (same/diff)
#' ## Graph
samediff_all_resp <- all_subj_filt %>%
  dplyr::filter(temporal_source_response %in% c(1:8)) %>%
  dplyr::group_by(subj) %>%
  dplyr::summarise(allresp = length(temporal_scored_list_samediff_half))

samediff_rates <- all_subj_filt %>%
  dplyr::filter(temporal_source_response %in% c(1:8)) %>%
  dplyr::group_by(subj, temporal_scored_list_samediff_half) %>%
  dplyr::summarise(respcount = length(temporal_scored_list_samediff_half)) %>%
  dplyr::left_join(samediff_all_resp) %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  dplyr::mutate(resprate = respcount / allresp)

samediff_rates %>%
  ggplot2::ggplot(ggplot2::aes(x = temporal_scored_list_samediff_half, y = resprate)) +
  ggplot2::geom_boxplot()

#' ## Stats (differ vs. chance)
samediff_rates %>%
  dplyr::filter(temporal_scored_list_samediff_half == "same_half_hit") %>%
  t.test(data = ., x = .$resprate, mu = 0.5)

#' # Question source memory, split by question
count_questmem_by_enctask <- all_subj_filt %>%
  # drop s14 because they made insufficient responses for each list half
  # NB: they do not get dropped from the rsa analyses or most other behavioral analyses
  dplyr::filter(subj != "s14") %>%
  # restrict to items that were correctly remembered
  # and for which they made a correct question source judgment
  dplyr::filter(item_recog_scored_rhit_collapsed == "rhit") %>%
  dplyr::filter(question_scored_exact == "hit") %>%
  # do not include unstudied items
  # (these should have already been dropped since they cannot be rhits, but just being safe)
  dplyr::filter(enctask %in% c(1, 2, 3, 4)) %>%
  dplyr::group_by(subj, enctask) %>%
  dplyr::summarise(numresp = length(question_scored_exact))

#' ## Graph
count_questmem_by_enctask %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(enctask), y = numresp)) +
  ggplot2::geom_boxplot()

#' ## Stats
count_questmem_by_enctask %>%
  ez::ezANOVA(data = .,
              dv = .(numresp),
              wid = .(subj),
              within = .(enctask),
              type = 2)
