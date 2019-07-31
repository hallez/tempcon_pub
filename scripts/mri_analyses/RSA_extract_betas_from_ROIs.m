function [] = RSA_extract_betas_from_ROIs
% Loop across ROIs, extract (mask) single-trial betas within these ROIs,
% and write out voxel x trial pattern matrices and trial pattern
% correlation matrices
%
% Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;


script_start = tic; 

% set current directory - this is a SPM oddity
curdir = pwd;
SLURM_FLAG = 0;
file_type = 'beta';

if SLURM_FLAG
    pc = parcluster('local');
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
else
    num_cores = feature('numcores');
    parpool(num_cores);
end
myPool = gcp;

parfor isub=1:32
    subj_start = tic;
    
    config = initialize_temp_con;
    
    REDUCED_ROIS_FLAG = 1;
    
    b.cur_subj = ['s' num2str(isub,'%01d')];
    b.analyzed_mri_dir = fullfile(config.analyzed_mri_dir, b.cur_subj, filesep);
    b.raw_behavioral_dir = fullfile(config.raw_behavioral_dir, b.cur_subj, filesep);
    b.analyzed_behavioral_dir = fullfile(config.analyzed_behavioral_dir, b.cur_subj);
    b.roi_dirs = {'ashs_left', 'ashs_right', 'manual_traced_rois'};
    
    % hardcoding this is necessary for setting up `pattern_corr_all_runs`
    b.num_trials_per_run = 60;

    b.model_dir = fullfile(config.analyzed_mri_dir,'multivariate_itemHitBYType');
    b.subj_model_dir = fullfile(b.model_dir, b.cur_subj);
    subj_model_dir_runs = dir(fullfile(b.subj_model_dir,'run*'));
    b.runs = {subj_model_dir_runs.name};
    b.num_runs = length(b.runs);
    
    if ismember(isub, config.exclude_subjects)
        fprintf('Current subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end

    fprintf('\n----Working on extracting %ss for %s----\n',file_type, b.cur_subj)
    
    b = run_exceptions_TempCon(b);
    
    exclude_runs_fname = fullfile(b.raw_behavioral_dir, 'exclude_runs.txt');
    b.exclude_runs = importdata(exclude_runs_fname);

    trial_lbls_fname = fullfile(b.analyzed_behavioral_dir, 'trial_information.csv');
    if(~exist(trial_lbls_fname, 'file'))
        fprintf('%s: no trial labels file.\n', b.cur_subj)
        continue;
    end
    
    ids = tdfread(trial_lbls_fname,',');

    for idir=1:length(b.roi_dirs)

        fprintf('%s: working on ROI dir: %s.\n', b.cur_subj, b.roi_dirs{idir})

        b.cur_ROI_dir = fullfile(b.analyzed_mri_dir, 'ROIs', b.roi_dirs{idir});
        b.all_ROIs = dir(b.cur_ROI_dir);
        b.all_ROIs = {b.all_ROIs.name};
        % make sure only getting binarized ROIs
        idx = regexp(b.all_ROIs,'^b*');
        idx = ~cellfun(@isempty,idx);
        b.all_ROIs = {b.all_ROIs{idx}};
        % make sure only getting ROIs (this is critical if re-running the
        % script since you'll have br*.mat files)
        idx2 = regexp(b.all_ROIs,'\>.nii');
        idx2 = ~cellfun(@isempty,idx2);
        b.all_ROIs = {b.all_ROIs{idx2}};

        % remove other ROIs we're never going to analyze to speed up the script
        rois_to_include = 1:length(b.all_ROIs);
        for iroi = 1:length(b.all_ROIs)
           if contains_str(b.all_ROIs{iroi}, 'zeros') ||  contains_str(b.all_ROIs{iroi}, 'Clear_Label') || contains_str(b.all_ROIs{iroi}, 'MISC') || contains_str(b.all_ROIs{iroi}, 'CS')
              rois_to_include(iroi) = 0;
           end
        end %iroi

        rois_to_include = logical(rois_to_include);
        b.all_ROIs = {b.all_ROIs{rois_to_include}};

        if REDUCED_ROIS_FLAG
            b.all_ROIs = {'brCA1_body.nii', 'brCA2_3_DG_body.nii', 'brsubiculum_body.nii',...
                'brwhole_hippo.nii',...
                'brPRC_L.nii', 'brPRC_R.nii',...
                'braPHC_L.nii', 'braPHC_R.nii', 'brpPHC_L.nii', 'brpPHC_R.nii',...
                'bralEC_L.nii', 'bralEC_R.nii', 'brpmEC_L.nii', 'brpmEC_R.nii'}; 
        end

        for iroi=1:length(b.all_ROIs)

            sprintf('Working on ROI %d of %d ROI: %s\n',iroi,length(b.all_ROIs),b.all_ROIs{iroi})

            cur_roi_fpath =fullfile(b.cur_ROI_dir,b.all_ROIs{iroi});
            
            if ~exist(cur_roi_fpath, 'file')
                fprintf('%s: ROI %s does not exist. Continuing.\n', b.cur_subj, cur_roi_fpath)
                continue;
            end
            cur_roi = spm_vol(cur_roi_fpath);
            cur_roi.img = spm_read_vols(cur_roi);

            if nansum(unique(cur_roi.img))> 0
                col_counter = 0;
                
                pattern_all_runs = nan(size(find(cur_roi.img > 0),1),b.num_trials_per_run*length(b.runs));

                for irun=1:length(b.runs)
                    
                    b.cur_run = b.runs{irun};
                    b.cur_run_num = str2double(strtok(b.cur_run,'run'));

                    b.current_run_dir = fullfile(b.subj_model_dir,b.cur_run);
                    current_trial_dirs = dir(b.current_run_dir);
                    trialCounter = 0;

                    % figure out how many trials were in the current run
                    for j = 1:size(current_trial_dirs,1)
                        s=current_trial_dirs(j).name;
                        if strfind(s,'trial_')
                            trialCounter = trialCounter + 1;
                        end 
                    end 

                    % subjects should have had 6 runs of 60 trials for a total of 360 trials
                    if trialCounter <= (b.num_trials_per_run * b.num_runs / b.num_runs)
                        numCurTrials = 1:trialCounter;
                    else
                        sprintf('Implausible number of trials for current run: %s.\n',b.cur_run)
                        continue;
                    end %if trialCounter

                    if(ismember(b.cur_run_num,b.exclude_runs))
                        % fill in a run of NaN values and iterate the
                        % column counter if the run is excluded so that
                        % the structure is preserved
                        fprintf('%s: excluding run %d - filling with NaNs.\n', b.cur_subj, b.cur_run_num)
                        pattern_all_runs(:,col_counter + 1:col_counter + b.num_trials_per_run) = nan(size(find(cur_roi.img > 0),1),b.num_trials_per_run);
                    else
                        fprintf('%s: working on: %s\n', b.cur_subj, b.cur_run);

                        % subset structure, based on:
                        % http://www.mathworks.com/matlabcentral/newsreader/view_thread/236371
                        cur_run_idx = ids.item_recog_run == b.cur_run_num;
                        b.cur_ids = structfun(@(v) v(cur_run_idx), ids, 'Uniform',0);

                        voxel_by_trial_patterns = nan(size(find(cur_roi.img > 0),1),trialCounter);

                        for itrial=1:trialCounter

                            b.current_trial = ['trial_',num2str(numCurTrials(itrial),'%03d')];

                            fprintf('%s: current trial is %s\n', b.cur_subj, b.current_trial)

                            % these beta files are created by
                            % 'RSA_single_trial_models_batch.m'
                            b.current_trial_dir = fullfile(b.model_dir,b.cur_subj,b.cur_run,b.current_trial);
                            if(~exist(fullfile(b.current_trial_dir,sprintf('%s_0001.img', file_type)), 'file'))
                                fprintf('\t %s does not exist - continuing.\n', sprintf('%s_0001.img', file_type))
                                continue;
                            end
                            b.cur_beta=spm_vol(fullfile(b.current_trial_dir,sprintf('%s_0001.img', file_type)));
                            b.cur_beta.img = spm_read_vols(b.cur_beta);

                            voxel_by_trial_patterns(:,itrial) = b.cur_beta.img(cur_roi.img > 0);
                            
                        end %itrial=
                        
                        % in order to collapse across trials, save out
                        % a version that includes NaN values
                        if size(voxel_by_trial_patterns,1) > 0                         
                          % put the current pattern into a
                          % variable that will store the pattern across
                          % all runs so that the ROI dimensions are
                          % preserved
                          pattern_all_runs(:,col_counter + 1:col_counter + b.num_trials_per_run) = voxel_by_trial_patterns;
                            
                          % first, let's get rid of rows w/ bad voxels and then NaN
                          % out trials that need to be removed (either on the basis
                          % of behavior or being a "bad" beta -- this is the same
                          % procedure for removing trials that's used in searchlight)
                          % it's okay to do this here b/c not removing any
                          % trial information, just bad **voxels** (ie, can
                          % still index based on trial ID matrix)
                          %
                          % NB: this elimination of bad voxels (rows) and bad trials
                          % (columns) is done in a 2-step process so don't end up in
                          % a situation where have NaN values in rows AND columns
                          % because then when go to compute correlations, this would
                          % result in a matrix of NaNs
                          %
                          % based on: http://www.mathworks.com/matlabcentral/answers/68510-remove-rows-or-cols-whose-elements-are-all-nan
                          pattern_without_bad_voxels = voxel_by_trial_patterns(all(~isnan(voxel_by_trial_patterns),2),:);

                          % now, NaN out bad trials
                          % we do NOT want to eliminate these columns b/c then this
                          % will get tricky when want to match up with trial IDs
                          pattern = pattern_without_bad_voxels;
                          pattern(:,logical(b.cur_ids.exclude_behavioral == 1)) = NaN;

                          % TODO: add in option to NaN out trials based on
                          % outlier betas 

                          % save out pattern matrix
                          [~,roi_name,~] = fileparts(cur_roi.fname);
                          b.fname_pattern_out = cellstr(fullfile(b.cur_ROI_dir,sprintf('%s_pattern_mtx_from_%ss_run_%s.mat', roi_name, file_type,num2str(b.cur_run_num,'%02d'))));
                          parsave(b.fname_pattern_out{:},pattern)
                        end %if size(pattern,1) > 0

                    end %(if(ismember(b.cur_run_num, b.exclude_runs 
                    
                    % iterate counter for pattern_all_runs
                    col_counter = col_counter + b.num_trials_per_run;

                end %irun=

                % save out pattern matrix
                % deal w/ when have no patterns
                % will happen for zeros ROIs
                if size(pattern_all_runs,1) > 0
                    % NaN out bad trials 
                    % we do NOT want to eliminate these columns b/c then this
                    % will get tricky when want to match up with trial IDs
                    pattern_all_runs_no_bad_trials = pattern_all_runs;
                    pattern_all_runs_no_bad_trials(:,logical(ids.exclude_behavioral == 99)) = NaN;

                    % NOW, we can finally save some stuff out
                    fname_pattern_out = cellstr(fullfile(b.cur_ROI_dir,sprintf('%s_pattern_mtx_from_%ss_all_runs.mat', roi_name, file_type)));
                    parsave(fname_pattern_out{:},pattern_all_runs_no_bad_trials)

                    % take correlation across all voxels and all trials
                    % and save this out too 
                    % test to see if matrix ONLY has NaN values - this will
                    % happen if have excluded full runs
                    % based on: https://stackoverflow.com/questions/38071444/check-if-matrix-is-not-nan-in-matlab#
                    if(~any(~isnan(pattern_all_runs_no_bad_trials(:))))
                        fprintf('%s: can''t run correlation - ''pattern_all_runs_no_bad_trials'' is all NaNs so just saving out.\n', b.cur_subj)
                        pattern_corr_all_runs = pattern_all_runs_no_bad_trials;
                    else
                        % use the 'rows', 'pairwise' option to only compute
                        % correlations for rows that contain non-NaN
                        % values. this is necessary for when there are both
                        % bad voxels (NaNs in rows) and bad trials/excluded
                        % runs (NaNs in columns)
                        pattern_corr_all_runs = corrcoef(pattern_all_runs_no_bad_trials, 'rows', 'pairwise');
                    end
                    fname_corr_out = cellstr(fullfile(b.cur_ROI_dir,sprintf('%s_pattern_corr_from_%ss_all_runs.mat', roi_name, file_type)));
                    parsave(fname_corr_out{:},pattern_corr_all_runs)
                end %if(size(pattern_all_runs

            end %if(sum(unique(cur_roi.img
                            
        end %iroi
                    
    end %idir=
    
    tEnd_subj =  toc(subj_start);
    sprintf('Extracting %ss for %s took %d minutes and %f seconds.\n',file_type, b.cur_subj,floor(tEnd_subj/60),rem(tEnd_subj,60))
    
end %isub

delete(myPool);

time_end_script =  toc(script_start);
sprintf('Extracting %ss for all subjects took %d minutes and %f seconds.\n', file_type, floor(time_end_script/60), rem(time_end_script,60))

end

% based on: http://www.mathworks.com/matlabcentral/answers/135285-how-do-i-use-save-with-a-parfor-loop-using-parallel-computing-toolbox
function parsave(fname, x)
    save(fname, 'x')
end
