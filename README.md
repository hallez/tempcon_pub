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
* s1: only has 59 trials in run 3; currently excluding until figure out a subject-specific way to handle different trial numbers
* s2: data were incorrectly oriented at acquisition, and don't have item recog data
* s3: no item recog data
* s7: exclude runs 4 and 5 because incorrect logfiles were used
* s8: exclude encoding list 8, had to restart this list at encoding and subject may have double-seen items
* s11: excessive motion in T2; unable to reliably trace subfields
* s13: experimenter error in which files were run; subject did not see proper counterbalancing of stimuli
* s15: exclude run 1 because subject saw first few stimuli twice when re-started this run at the scanner
* s24: too much motion to trace T2
* s25: exclude encoding list 1, had to restart this list and subject saw multiple times
* s26: missing temporal source memory data

## Behavioral analyses:
1. Load in data from Presentation and re-format into R: `load_data.R`
2. Tidy up and create some columns that are helpful for scoring: `behav-tidy.R`
3. Generate onset and trial information files: `create-onset-files.R`
4. Plot behavioral data and generate some plots: `behav-stats.R`

## Prepping fMRI data:
1. I like to rename subject folders so everyone has the same directory structure (i.e., instead of dealing with individualized paths in every script). The standard directory names are contained in `scripts/mri_analyses/default_dirs.yml`. Subject-specific directories are saved in `raw_behavioral/<subject-dir>/<subject-id>.yml`
  - Since `raw_behavioral` is in the `.gitignore` file, these subject-specific directory files can be downloaded from [osf.io](https://osf.io/qfcjg/)
1. If starting from raw scanner files, unzip (`python scripts/mri_analyses/unzip_raw_mri.py`), do dicom conversion (`dicom_conversion.m`), and rename subject directories (`python scripts/mri_analyses/rename_unpacked_folders.py`)
  - NB: This is not necessary since I am not supplying zipped dicoms. Scripts are included for reference.
1. Preprocess the data (this will generate `rf*.nii` files): `preproc_moco_quickcoreg.m`
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

1. Divide the ASHS and manually traced ROIs into single files `RSA_extractROIs.m`. In order to run for both sets of ROIs, will need to reset the value of `ASHS_FLAG` and `MANUAL_FLAG`
