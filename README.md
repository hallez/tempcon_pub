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
