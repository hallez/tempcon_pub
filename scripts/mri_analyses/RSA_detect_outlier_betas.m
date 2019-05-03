% Read in betas and compute the mean absolute beta weight to detect "outlier" 
% betas. 
%
% Halle R. Dimsdale-Zucker

% \\TODO: write out diaries/logs for future reference (problem: this won't
% play nicely w/ parfor)
all_subjects_timer = tic; % time script so know how long things take for future use

outlier_threshold = 3; % mark a beta for exclusion if the mean absolute z-score exceeds 3 SD from the mean

mask_type = 'whole_brain'; % can use this to set ROIs or whole brain in future

% if want to test out z-scoring procedure, use this test code:
% for tmp_trial = 1:4
%     if tmp_trial == 1 % this should probably actually be a switch/case loop instead
%         tmp_cur_beta = [-5.0131, 17.9966, -15.2866, NaN; 21.2882, -8.3428, 23.7833, -8.4788; NaN, -2.8608, 5.6037, 5.0752];
%     elseif tmp_trial == 2
%         tmp_cur_beta = [NaN, 4.8909, NaN, 17.0109; -11.8961, -11.8837, -4.3014, -35.1398; -30.7396, -2.6523, NaN, 3.8816];
%     elseif tmp_trial == 3
%         tmp_cur_beta = [6.2717, 16.2217, -7.4785, NaN; -9.8486, 2.0320, NaN, 25.8610; -8.2944, 4.8142, NaN, -12.4082];
%     elseif tmp_trial == 4
%         tmp_cur_beta = [NaN, NaN, -2.3983, -4.5782; 3.2712, 0.2500, -4.4423, -1.2393; -3.1110, -4.1794, NaN, 1.2477];
%     end
%     tmp_mask = [1,1,0,0;1,1,1,1;1,0,1,1]; % eventually this will be a brain mask or an ROI mask
%     tmp_all_betas(:, tmp_trial) = tmp_cur_beta(tmp_mask > 0);
% end
% 
% % we want to z-score the values of a voxel across trials so
% % take nanmean and nanstd across trials dimension (ie,
% % dimension 2)
% nanmean_tmp_all_betas = nanmean(tmp_all_betas,2); % should yield a column w/ as many rows as 1s in tmp_all_betas
% 
% % remember that for nanstd to take SD across a different
% % dimension, need to specify a flag (1 or 0) for how to compute
% % this SD. opting to use flag = 1 which means that normalizing
% % by n (instead of n-1)
% nanstd_tmp_all_betas = nanstd(tmp_all_betas, 1, 2);
% 
% % use repmat so that the summary (ie, mean, SD) matrices are
% % the same size as the matrix of all betas
% nanmean_mtx = repmat(nanmean_tmp_all_betas, [1,size(tmp_all_betas,2)]); 
% nanstd_mtx = repmat(nanstd_tmp_all_betas, [1, size(tmp_all_betas,2)]);
% 
% % compute z-score ((value - mean) / SD)
% tmp_numerator = tmp_all_betas - nanmean_mtx; % only breaking this out so can check division
% z_tmp_all_betas = tmp_numerator ./ nanstd_mtx;
% 
% % what we actually care about is the overall (ie, mean absolute
% % value) of these z-scores for each trial (column)
% abs_z_tmp_all_betas = abs(z_tmp_all_betas);
% mean_abs_z_tmp_all_betas = nanmean(abs_z_tmp_all_betas);


for isub = 1:32
    % let's also time the length for each subject 
    cur_subj_timer = tic;
    
    config = initialize_temp_con; % do w/in isub loop in case want to use parfor

    b.cur_subj = ['s' num2str(isub,'%01d')];
    fprintf('------Working on subject %s------\n',b.cur_subj)

    if ismember(isub, config.exclude_subjects)
        fprintf('\tCurrent subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end
        
    % setup subject-specific directories, etc.
    b.cur_subj_fpath = fullfile(config.analyzed_mri_dir,...
        'multivariate_itemHitBYType',...
        b.cur_subj);
    run_dirs = dir(fullfile(b.cur_subj_fpath, 'run*')); % define number of runs based on directories that are present
    num_runs = length(run_dirs); % is `size(run_dirs, 1)` more appropriate here?
    b.runs = {run_dirs.name};
    
    % run exceptions for subject-specific naming conventions
    % don't actually think there are any exceptions for tempcon subjects,
    % but keeping this here in case need to add any
    b = run_exceptions_TempCon(b);
        
    % read in subject-specific brain mask. these brain masks are created by
    % the QA script and are always stored in the first run folder (w/ the
    % assumption that the subject should always be in the same location
    % across runs - which is true since the files that are used to compute
    % betas are realigned EPIs)
    if strcmp(mask_type, 'whole_brain')
        subj_brain_mask_fpath = fullfile(config.analyzed_mri_dir, b.cur_subj, 'run1', 'brainmask.nii');
    else 
        % in future, can use this for different ROIs
    end
    
    if ~exist(subj_brain_mask_fpath, 'file')
        fprintf('\tCurrent brain mask %s does not exist. Skipping.\n', subj_brain_mask_fpath)
        continue;
    end
    subj_brain_mask = spm_read_vols(spm_vol(subj_brain_mask_fpath));
    
    % read in beta maps across trials & runs
    for irun = 1:num_runs
        b.cur_run = b.runs{irun};
        b.cur_run_fpath = fullfile(b.cur_subj_fpath,...
                b.cur_run);
        trial_dirs = dir(fullfile(b.cur_run_fpath, 'trial*'));
        num_trials = length(trial_dirs);
        
        for itrial = 1:num_trials
            b.cur_trial = ['trial_', num2str(itrial, '%03d')];
            b.cur_trial_fpath = fullfile(b.cur_run_fpath,...
                b.cur_trial);
            % only the first beta image is ever used for RSA b/c of how I 
            % build the single trial models
            b.cur_beta_fpath = fullfile(b.cur_trial_fpath, 'beta_0001.img'); 
            
            if ~exist(b.cur_beta_fpath, 'file')
                fprintf('\tCurrent beta file %s does not exist. Skipping.\n', b.cur_beta_fpath)
                continue;
            end %if exist
            b.cur_beta = spm_read_vols(spm_vol(b.cur_beta_fpath)); 
            
            % only grab values that are w/in the brain mask (ie, don't 
            % include all of the voxels that are outside of the brain b/c
            % these are NaN values and make things more challenging to work
            % with. the brain mask only contains 1s/0s, but, technically, the
            % 1s are 1.000 (or something like this). grabbing voxels for which the
            % mask is > 0 is essentially a way of binarizing it. after
            % extract voxels that are w/in the brain mask, there is no
            % longer any sense of x/y/z dimensions. this is ok b/c the
            % voxels retain the same index position (which is what we will
            % need to z-score a voxel across trials) even if we no longer
            % know where it was in x/y/z coordiante space.
            % //TODO can also adapt this mask to be an ROI of
            % interest if using regional variance to identify bad betas
            % instead of the whole brain
            masked_beta_vals(:, itrial) = b.cur_beta(subj_brain_mask > 0);
            
        end %itrial
        
        % if for some reason there are no trial for this run,
        % `masked_beta_vals` shouldn't exist so skip trying to z-score.
        if ~exist('masked_beta_vals') 
            continue;
        end
        % z-scoring: this is done for each voxel across all trials. we
        % can't z-score w/in a beta b/c then this wouldn't allow us to
        % detect if the beta itself is an outlier. things get a little 
        % tricky though because there are sometimes still NaN values. 
        % why? could be slight motion/misalignment
        % of betas w/ brain mask (although we hope not!), could be due
        % to problem when converting multiband,
        % could be CSF (since the brain masks are created by combining
        % GM, WM, CSF). but, what matters, is that we can disregard
        % these NaN voxels for the purposes of detecting outlier betas.
        nanmean_masked_beta_vals = nanmean(masked_beta_vals,2);
        % opting to use flag = 1 which means that normalizing
        % by n (instead of n-1)
        nanstd_masked_beta_vals = nanstd(masked_beta_vals, 1, 2);
        
        % compute z-score ((value - mean) / SD)
        z_masked_beta_vals = (masked_beta_vals - repmat(nanmean_masked_beta_vals, [1,size(masked_beta_vals,2)])) ./ repmat(nanstd_masked_beta_vals, [1, size(masked_beta_vals,2)]);
        mean_abs_betas(irun,:) = nanmean(abs(z_masked_beta_vals)); % aggregate across runs for plotting; this will fail for s1 b/c run 3 has a different number of trials

        % save out bad beta information
        if sum(mean_abs_betas(irun,:) > outlier_threshold) > 0
            fprintf('\tWriting out bad timepoints.\n')
            fname_bad_timepoints = sprintf('%s_bad_timepoints_%s.dat', b.cur_subj, b.cur_run);
            cur_run_outliers = find(mean_abs_betas(irun,:) > outlier_threshold);
            csvwrite(fullfile(b.cur_subj_fpath, fname_bad_timepoints), cur_run_outliers)
        end 
        
        clear masked_beta_vals % do this so can re-set between runs incase one run gets skipped
        
    end %irun
    
    % plot betas across runs 
    % //TODO make all of the legend labels visible (gets truncated b/c
    % > 50 characters)
    figure;
    hist(mean_abs_betas');
    title(sprintf('%s - betas from %s', b.cur_subj, mask_type),'interpreter','none'); % tell matlab underscores aren't subscripts
    legend(sprintf('%s\t', b.runs{:})) % this dynamically can change if number of runs changes
    box off; % hide the outside axis lines
    ylabel('Number of trials')
    xlabel('Mean abs z-score')
    graph_fname_out = fullfile(b.cur_subj_fpath, sprintf('%s_mean-abs-z-betas_distribution.png', b.cur_subj));
    saveas(gcf, graph_fname_out);
    
    tEnd_subj =  toc(cur_subj_timer);
    fprintf('\tProcessing %s took %d minutes and %f seconds.\n',b.cur_subj,floor(tEnd_subj/60),rem(tEnd_subj,60))
end %isub

tEnd_all_subjects = toc(all_subjects_timer);
fprintf('\nRunning all those subjects took %d minutes and %f seconds. \n',floor(tEnd_all_subjects/60),rem(tEnd_all_subjects,60))

