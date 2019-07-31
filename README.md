# Public repository for scripts related to:
Dimsdale-Zucker, H.R., Montchal, M.E., Reagh, Z.M., Wang, S.-F., Libby, L.A., & Ranganath, C. (in prep)

## Data access:
1. Raw behavioral data, de-faced raw MRI images, single-trial betas, and resliced & binaraized ROIs can be found on [osf.io](https://osf.io/qfcjg/)
1. A reproducible code capsule is available from [Code Ocean](https://codeocean.com/capsule/0129473)
  * this capsule contains the data and code to re-generate behavioral plots and stats as well as RSA plots and mixed models (but not full permutation results, due to computational limitations)

## Environment setup and requirements:
1. If you would simply like to re-generate analyses, that can be done from [Code Ocean](https://codeocean.com/capsule/0129473)
  * If you sign up with your own account, I think you should also be able to modify the code and re-generate the results to see how this changes the output
1. If you would like to run on your own compute environment, here are the basic requirements:
  * MATLAB r2014b (note: scripts may work in newer versions, but in the limited testing I've done in 2016 and 2017 versions there seem to be some SPM compatibility issues)
  * [SPM 8](https://www.fil.ion.ucl.ac.uk/spm/software/download/)
  * some [custom MATLAB functions](https://github.com/hallez/MatlabFunctions)
  * Python 2.7.15
    - dependency packages can be installed via `pip install -r requirements.txt`
  * R >= 3.4.4
    - see the DESCRIPTION file for required packages
  * Copy config.yml.example and update the paths `cp config.yml.example config.yml`
    - NB: I saved some results to a dropbox folder; when setting this path (`dropbox_tempcon` and `dropbox_base_dir`), you can specify anywhere you would like the results to be saved
  * Update the paths for `project_dir` in `scripts/mri_analyses/initialize_temp_con.m`
    - Some of these analyses were run on our HPC. That is what `SLURM_FLAG` is referring to.
1. Code from other sources has been included in the `vendor` directory. This is not code that I wrote and therefore is subject to the usage license and instructions of those authors. It is included here out of convenience for those trying to run the scripts in the current repository.
2. R scripts are setup as an [R package](http://r-pkgs.had.co.nz/) - this is probably not quite the right time to use a package, but I did this as a convenient way to load the dependencies. What this means is that all of the R scripts assume that you are first opening the `tempcon.Rproj` file and then running the R scripts from within RStudio.
3. Matlab scripts assume you are running from `<repo-dir>/scripts/mri_analyses` with SPM8, `vendor` and my [custom MATLAB functions](https://github.com/hallez/MatlabFunctions) as part of your path.
  - If they aren't on your path, a good place to add them would be `scripts/mri_analyses/initialize_temp_con.m`
1. Python scripts assume you are running them from terminal within the `<repo-dir`.

## Excluded subjects
R-formatted list: `exclude_subjects <- c("s1", "s2", "s3", "s11", "s13", "s15", "s24", "s26")`
Matlab-formatted list: `config.exclude_subjects = [1, 2, 3, 11, 13, 15, 24, 26];`
To print out information about excluded runs and subjects, `python scripts/tabulate_exclusions.py`
* s1: problem at data collection, unsure what happened but missing one trial in each run; excluding entirely just to be safe
* s2: data were incorrectly oriented at acquisition, and don't have item recog data
* s3: no item recog data
* s7: exclude runs 4 and 5 because incorrect logfiles were used
* s8: exclude encoding list 8, had to restart this list at encoding and subject may have double-seen items
* s11: excessive motion in T2; unable to reliably trace subfields
* s13: experimenter error in which files were run; subject did not see proper counterbalancing of stimuli
* s14: exclude encoding list 8 items b/c had to restart list and subject may have seen some items twice
   * missing some temporal and question source memory data; excluded from relevant behavior analyses, but not rsa (since this only depends on item recognition data)
* s15: exclude run 1 because subject saw first few stimuli twice when re-started this run at the scanner
  * excluding run 1 in addition to excluding runs 3-6 due to motion means only one run (run 2) remains; therefore, exclude this subject entirely.
* s24: too much motion to trace T2
* s25: exclude encoding list 1, had to restart this list and subject saw multiple times
* s26: missing temporal source memory data

## Notes about number of TRs
trials that do not have at least a 16 second buffer from the end of the run get dropped in `create-onset-files.R`
* s1-s10 have a different number of TRs than s12-s32
* s12-s32 have 268 TRs in each run
* s1:s5, s10 have 262 TRs in each run
* s6:s9 have 267 TRs in each run

## Behavioral analyses:
1. Load in data from Presentation and re-format into R: `load_data.R`
2. Tidy up and create some columns that are helpful for scoring: `behav-tidy.R`
3. Generate onset and trial information files: `create-onset-files.R`
4. Plot behavioral data and generate some plots: `behav-stats.R`

## Prepping fMRI data:
1. I like to rename subject folders so everyone has the same directory structure (i.e., instead of dealing with individualized paths in every script). The standard directory names are contained in `scripts/mri_analyses/default_dirs.yml`. Subject-specific directories are saved in `raw_behavioral/<subject-dir>/<subject-id>.yml`
  - Since `raw_behavioral` is in the `.gitignore` file, these subject-specific directory files can be downloaded from [osf.io](https://osf.io/qfcjg/)
1. If starting from raw scanner files, unzip (`python scripts/mri_analyses/unzip_raw_mri.py`), do dicom conversion (`scripts/mri_analyses/dicom_conversion.m`), and rename subject directories (`python scripts/mri_analyses/rename_unpacked_folders.py`)
  - NB: This is not necessary since I am not supplying zipped dicoms. Scripts are included for reference.
1. Preprocess the data (this will generate `rf*.nii` files): `scripts/mri_analyses/preproc_moco_quickcoreg.m`
1. Run QA. I use code that is based on [Maureen Ritchey's code](https://github.com/memobc/memolab-fmri-qa)
  - The key thing is that the `spike_regs_rp.txt` motion files get created so they can be read into the models.
1. Based on the QA, decide what runs/subjects need to be dropped. Manually update `subject_exclusion_info.csv` with any runs that need to be dropped based on QA. Run `python scripts/create_subject_exclude_run_files.py` to create exclusion files for each subject.  
  - This can be used in combination with `run_exceptions_TempCon.m` (although currently I haven't updated the subject-specific run information in this file - instead I deal with subject-specific things in the tidying scripts)

## ROI segmentation: ASHS (automated)
This is how I generate the hippocampal subfield ROIs.

NB: I've already uploaded ROI masks to [osf.io](https://osf.io/qfcjg/), but this is here for reference if you wanted to re-generate the ROIs.

1. Install ashs. Check on [NITRC](https://www.nitrc.org/projects/ashs/) for latest install. Git code below for version downloaded October 2016.
```
git clone https://github.com/pyushkevich/ashs/tree/d03ddc919e1d38777736dec8b95286cf6afa786d
```
1. Test install using [their instructions](https://sites.google.com/site/hipposubfields/installation)
1. Download atlas of choice. I use: `ashs_atlas_upennpmc_20140416.tar`
1. Update the paths in `scripts/mri_analyses/run-ashs.sh`
1. Run the following command for whichever subjects currently working on.
```
for SUBJ in `seq 1 32`; do export SUBJ="s$SUBJ"; scripts/mri_analyses/run-ashs.sh; done
```

## ROI segmentation: Manual
This is how I generate the MTLc ROIs.

NB: I've already uploaded ROI masks to [osf.io](https://osf.io/qfcjg/). Original tracings are available upon request.

## ROI prep
Even though the ROIs have been traced, we need to extract them from the tracing files and ensure they'll play nicely when we use them to mask our beta images. Again, all of these steps have already been done if you're downloading them from [osf.io](https://osf.io/qfcjg/). Scripts are provided here as a reference.

1. Divide the ASHS and manually traced ROIs into single files: `scripts/mri_analyses/RSA_extractROIs.m`. In order to run for both sets of ROIs, will need to reset the value of `ASHS_FLAG` and `MANUAL_FLAG`
1. Combine ROIs of interest (e.g., CA23DG) - only relevant for ASHS ROIs: `scripts/mri_analyses/RSA_combine_ROIs_batch.m`
1. Split ROIs into head/body tail based on these boundaries - only relevant for ASHS ROIs: `RSA_split_HBT.m`
  * This uses the subject-specific transition files (`<subject-id>_hc_transitions.yml`) that are included with the raw data [osf.io](https://osf.io/qfcjg/).
1. Reslice the ROIs into EPI space: `scripts/mri_analyses/RSA_reslice_t2_and_ROIs_batch.m`
1. Binairize the ROIs so can use as masks: `scripts/mri_analyses/RSA_binarize_ROIs_batch.m`

## Multivariate fMRI analyses
Many of these scripts are very time-intensive to run. If possible, I recommend running on an HPC cluster. However, I have included both the single trial beta images as well as the extracted pattern correlations on [osf.io](https://osf.io/qfcjg/). You can also re-generate the output from the mixed models in the [Code Ocean capsule](https://codeocean.com/capsule/0129473). If you want to run all 1000 iterations of the permutations, you can either modify the code on Code Ocean or run on your own HPC. Even my very powerful iMac with 32GB of RAM struggled with these.

1. Generate regressors for single trial models: `scripts/mri_analyses/RSA_generate_single_trial_regressors.m`
1. Compute single trial models. This takes ~7.5 hours/subject to run. NB: This assumes you have preprocessed to generate the `rf` files:  `scripts/mri_analyses/RSA_single_trial_models_batch.m`
1. Extract the beta values from within each ROI of interest and compute trial-by-trial correlation matrices: `scripts/mri_analyses/RSA_extract_betas_from_ROIs.m`

***Move from matlab into R (just my personal preference)***

1. Create masks based on trials pairs of interest: `scripts/tempcon/rsa-generate-masks-btwn-runs.R`
1. Read in the trial-by-trial pattern correlation matrices from matlab into R and extract the trial pairs of interest: `scripts/tempcon/rsa-load-data-btwn-runs.R`
1. Tidy up these trial pair correlations so easier to work with: `scripts/tempcon/rsa-tidy-data-btwn-runs.R`
2. Run mixed models: `scripts/tempcon/mixed-models-btwn-runs.R`
3. To evaluate the significance of these mixed models, run permutations (with 1000 iterations). These are split into four separate files to speed things up since on most clusters they can be run in parallel. If you just want to see how the scripts work, you can check them out on [Code Ocean](https://codeocean.com/capsule/0129473) where they run with a reduced number of permutation iterations: `scripts/tempcon/perms_*.R`
4. To generate the figures as shown in the paper that combine means, difference scores, and permutations use `scripts/tempcon/graphs-of-mixed-models-data.R`

## Searchlight
**NB: these are just here as example scripts, these are not analyses included in this project.**

1. These scripts assume that trial-by-trial masks have been created (`rsa-generate-masks-btwn-runs.R`), the gray matter masks exists (`run1/c1meanf*.nii`; these can be generated with SPM's segment and deformations procedures which generate probabilistic brain masks - an example file for s10 can be found on [osf.io](https://osf.io/qfcjg/)), and single trial betas already exist (`RSA_single_trial_models_batch.m`).
2. Run searchlight (this takes about 45 minutes per subject): `searchlight_btwn_runs.m`
3. To evaluate these searchlights, you would need to take a subtraction between the maps from each condition. Then, to make inferences about group-level results, the maps would need to be normalized into some common space (e.g., MNI). Since I don't actually include these analyses, I am not posting scripts for these steps.

## Univariate
(These scripts/analyses are still in progress)
