function [] = searchlight_btwn_runs()

% This function conducts a searchligh pattern similarity analysis. The
% basic approach taken is to iterate through all voxels in the brain and to
% extract the pattern from the given volume (cube, sphere, weird 3D diamond
% thing).
% 
% Author: Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;

script_start = tic; 

allModels(1).name = 'itemHitBYType_model';
allModels(1).type = 'searchlight';

modelSelect = 0;
for imodel=1:size(allModels,2)
    curModel = allModels(imodel).name;
    curModelBase = strtok(curModel,'_');
    curModelType = allModels(imodel).type;
    modelSelect = input(['Do you want to analyze ' curModel '? (Y=1,N=0): ']);
    if modelSelect
        break;
    end %if modelSelect
end %imodel=
fprintf('Selected to analyze %s\n',curModel)

num_cores = feature('numcores');
parpool(num_cores);
myPool = gcp;

parfor isub=1:32
    subj_start = tic;
    
    config = initialize_temp_con;
    
    b.cur_subj = ['s' num2str(isub,'%01d')];
    b.analyzed_mri_dir = fullfile(config.analyzed_mri_dir, b.cur_subj, filesep);
    b.analyzed_behavioral_dir = fullfile(config.analyzed_behavioral_dir, b.cur_subj, filesep);
    b.cur_single_trials_dir = fullfile(config.analyzed_mri_dir, ['multivariate_', curModelBase], b.cur_subj);
    searchlight_dirname = [curModelType, '_', curModelBase];
    searchlight_dir = fullfile(config.analyzed_mri_dir, searchlight_dirname);
    % create output directory if it doesn't already exist
    if(~isdir(searchlight_dir))
       mkdir(searchlight_dir); 
    end
    b.cur_searchlight_dir = fullfile(config.analyzed_mri_dir, searchlight_dirname, b.cur_subj);
    b.runs = {'run1' 'run2' 'run3' 'run4' 'run5' 'run6'};
    b.numruns = length(b.runs);
    b.numtrials = 60;
    
    if ismember(isub, config.exclude_subjects)
        fprintf('Current subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end

    fprintf('\n----Working on searchlight for %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);
    
    ZSCORE_FLAG = 0;
    
    % this is used to preallocate the size of `mean_PS_by_trial_type`
    trial_types = {'rANDsourceHits_sameQuestion', 'rANDsourceHits_diffQuestion'};
    num_trial_types = size(trial_types,2);
    
    %%% read in mask files
    % these are created by `rsa-load-data-btwn-runs.R`
    % set loaded file equal to a variable (this is necessary for parfor)
    rANDsourceHits_sameQuestion_mask = load(fullfile(b.cur_single_trials_dir, 'rANDsourceHits_sameQuestion_mask.txt'));
    rANDsourceHits_diffQuestion_mask = load(fullfile(b.cur_single_trials_dir, 'rANDsourceHits_diffQuestion_mask.txt'));
    
    
    %%% setup stuff about the searchlight itself
    % let's define the voxel extent off of the center cube
    % MCI and LL call this `slrad` or `ceil_radius` 
    % (even though it's not technically the radius)
    %
    % LL and MCI use a 5-voxel "radius" w/ standard res 
    % thus, they are getting a 5*3 = 15mm diadmeter
    % one approach would be to match the diameter (in mm), but, since we went
    % to the trouble of getting the high res data, seems like we want to
    % leverage this. thus, sticking with the total number of voxels in the
    % searchlight being around 20-25 seems like a good plan
    searchlight_arm_length = 2; 

    % setup distance metric
    % 1 = correlation
    % 2 = Euclidian distances
    % 3 = Mahalanobis 
    % 4 = cosine
    % ADD SPEARMANS
    distance_metric = 1;
    
    if(~isdir(b.cur_searchlight_dir))
       mkdir(b.cur_searchlight_dir); 
    end
    
    % limit to GM-only voxels
    % this file is created during the `segmentation` step
    % by script `univariate_normalize_to_MNI.m`
    % 
    % remember that even though this file doesn't have an `r`, it's in
    % register w/ the `rf` images b/c the `rf` images are created when you
    % align to the `meanf`
    run1_files = dir(fullfile(b.analyzed_mri_dir, 'run1', 'c1meanf*.nii'));
    b.maskname = fullfile(b.analyzed_mri_dir,'run1',run1_files.name);
    if(~exist(b.maskname, 'file'))
        fprintf('Missing brain mask for %s.\n', b.cur_subj)
        continue; 
    end
    mask = spm_read_vols(spm_vol(b.maskname));
    % set to NaN where there's no brain
    mask(mask==0) = NaN;
    
    mean_PS_by_trial_type = zeros(size(mask));
    
    % read in betas across trials (so not doing in the larger forloop
    % later)
    all_betas = zeros(size(mask, 1), size(mask,2), size(mask,3), b.numtrials, b.numruns);
    for irun = 1:b.numruns
        b.cur_run = b.runs{irun};
        b.cur_run_num = str2double(strtok(b.cur_run, 'run'));

        % load beta patterns for each trial
        b.cur_betas_dir = fullfile(b.cur_single_trials_dir, b.cur_run);  
        % first, we need to generate a list of filenames
        % [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
        % 'ExtFPListRec' means search recursively within subdirectories 
        % if omit `frames`, will default to 1
        [files,~] = spm_select('ExtFPListRec',b.cur_betas_dir,'beta_0001.img'); 

        % read in the image volume information for each beta
        V = spm_vol(files);

        for itrial = 1:size(V,1)
           all_betas(:,:,:,itrial,irun) = spm_read_vols(V{itrial}) .* mask; 
        end

    end %irun =
        
    % and now we loop across all the voxels in each trial's beta
    % MR uses adjblock to eliminate looping across dimensions where 
    % half of the searchlight is hanging off into non-brain space
    %
    % MCI and LL use a different way of buffering by actually changing
    % the coordinates the the spotlight is covering
    % also, the size of the MCI/LL spotlight can change (although this
    % only occurs near the edges of the brain)
    for x = 1:size(mask,1)
        for y = 1:size(mask,2)
            for z = 1:size(mask,3)
                % check to make sure there's actually brain at the current
                % x,y,z coordinate
                if ~isempty(find(mask(x,y,z) > 0, 1))
                    
                    % set the size/shape of the searchlight
                    % (spotlight)
                    spotlight = [];
                    for xdim = -1*searchlight_arm_length:searchlight_arm_length
                        for ydim = -1*searchlight_arm_length:searchlight_arm_length
                            for zdim = -1*searchlight_arm_length:searchlight_arm_length
                                % check the make sure the
                                % spotlight won't be butting up
                                % against the edge of the brain
                                %
                                % check that each coordinate is >0
                                % (x+xdim > 0) but does not
                                % exceed the maximum mask dimension
                                % (x+xdim <
                                % size(mask(relevant_dimension,1))
                                % 
                                % also check that current
                                % spotlight does not exceed
                                % the dimensions that the 
                                % spotlight would take up if
                                % it were an inflated cube
                                % volume of a cube = edge^3
                                % but, since xdim, ydim,
                                % zdim may be be equal
                                % (although when would they
                                % not be???), compute it as
                                % if there are 3 squares and
                                % take the area of each
                                % separately then sum
                                % together
                                if x + xdim > 0 && ...
                                   x + xdim < size(mask,1) && ...
                                   y + ydim > 0 && ...
                                   y + ydim < size(mask,2) && ...
                                   z + zdim > 0 && ...
                                   z + zdim < size(mask,3) && ...
                                   sqrt(xdim^2 + ydim^2 + zdim^2) <= searchlight_arm_length

                                   % also check that the
                                   % current x, y, z voxel is
                                   % w/in the brain mask
                                   % previously, this would
                                   % compare mask(x + xdim, y +
                                   % ydim, z + zdim) == 1, but
                                   % this doesn't work (as per
                                   % https://www.mathworks.com/matlabcentral/newsreader/view_thread/66553)
                                   % instead. additionally, 
                                   % cannot just testing if
                                   % mask(x + xdim, y + ydim, z
                                   % + zdim) is true because
                                   % this runs into issues with
                                   % this value is NaN. thus,
                                   % entering loop when the
                                   % current mask coordinate is
                                   % NON a NaN value. so
                                   % convoluted.
                                   if ~isnan(mask(x + xdim, y + ydim, z + zdim))
                                       % append the current
                                       % dimensions to
                                       % spotlight for each of
                                       % the x, y, and z (1:3) 
                                       % dims
                                       spotlight(end+1,1:3) = [x + xdim, y + ydim, z + zdim];
                                   end % if mask(x+xdim...

                                end % if x+xdim...
                            end %zdim
                        end % ydim
                    end %xdim, spotlight definition
                    
                    % checks that the the pattern is at
                    % least some minimum size
                    if size(spotlight,1) <= 9
                        fprintf('Spotlight < 9 voxels for x = %d, y = %d, z = %d.\n', x, y, z)
                        continue;
                    end
                        
                    % initialize a column counter 
                    % this will be used to make sure that the current run's pattern
                    % gets put into the correct place in the across run matrix
                    col_counter = 0;

                    % initialize a matrix of NaNs where can put pattern
                    % matrices from each run
                    pattern_all_runs = NaN(size(spotlight,1), b.numtrials*b.numruns);

                    % extract the pattern at each trial across
                    % each run
                    for irun = 1:b.numruns
                        b.cur_run = b.runs{irun};
                        b.cur_run_num = str2double(strtok(b.cur_run, 'run'));

                        % create a temporary variable the size 
                        % of the searchlight by the size of the 
                        % number of trials
                        patt_spotlight_by_trials = NaN((size(spotlight,1)), b.numtrials); 

                        % loop through trials extracting the beta
                        % values for the current searchlight
                        % coordinates
                        % to get the indexing to work properly,
                        % you need to loop across trials and
                        % then select the appropriate row
                        % dimension for the current coordinates
                        for itrial = 1:b.numtrials
                            for idim = 1:size(spotlight,1)
                                patt_spotlight_by_trials(idim,itrial) = all_betas(spotlight(idim,1),spotlight(idim,2),spotlight(idim,3), itrial, irun); 
                            end
                        end

                        % put the current pattern (and ids) into a variable
                        % across runs
                        % this method is nice b/c we're left with a 2D
                        % variable which is easier to use in `corrcoef`
                        pattern_all_runs(:,col_counter + 1:col_counter + b.numtrials) = patt_spotlight_by_trials;
                        col_counter = col_counter + b.numtrials;

                    end %irun =

                    % check that patt_spotlight_by_trials has at least 2
                    % rows with data and 2 columns
                    % because, if not, when try to take
                    % correlations (ie, pattsim) will fail
                    % logic based on: http://www.mathworks.com/matlabcentral/answers/31104-deleting-nan-rows-in-a-matrix
                    if sum(any(~isnan(pattern_all_runs),1))>=2 && sum(any(~isnan(pattern_all_runs),2))>=2

                        % eliminate rows that contain NaNs because
                        % this will mess up `corrcoeff`
                        pattern_all_runs_no_nans = pattern_all_runs(~any(isnan(pattern_all_runs),2),:);

                        % z_score (if required)
                        % this can't have a different name b/c then
                        % later references will break (or would need if
                        % statements and that's a overkill)
                        if ZSCORE_FLAG 
                            pattern_all_runs_no_nans = bsxfun(@rdivide,bsxfun(@minus,pattern_all_runs_no_nans,nanmean(pattern_all_runs_no_nans)),nanstd(pattern_all_runs_no_nans));
                        end %ZSCORE_FLAG

                        % compute pairwise correlations for voxels
                        % w/in current searchlight
                        if distance_metric == 1 % correlationa
                            pattsim = corrcoef(pattern_all_runs_no_nans);
                        elseif distance_metric == 2 % Euclidian distance (aka, pairwise distance) 
                                                    % this is the default behavior of pdist
                            pattsim = pdist(pattern_all_runs_no_nans);
                        elseif distance_metric == 3
                            pattsim = pdist(pattern_all_runs_no_nans, 'mahalanobis');
                        elseif distance_metric == 4
                            pattsim = pdist(pattern_all_runs_no_nans, 'cosine');
                        else
                            error('No distance metric selected.\n');
                        end

                        % append trial patterns to a variable to be
                        % written out later
                        % take the mean across rows and across columns so have one summary value
                        mean_PS_by_trial_type(x, y, z, 1) = nanmean(nanmean(rANDsourceHits_sameQuestion_mask .* pattsim));
                        mean_PS_by_trial_type(x, y, z, 2) = nanmean(nanmean(rANDsourceHits_diffQuestion_mask .* pattsim));
                        
                        
                    end %if sum(any(~isnan(patt_spotlight_by_trials),1))>=2 && sum(any(~isnan(patt_spotlight_by_trials),2))>=2
                end %if ~isempty(find(mask(x,y,z...
            end %z
        end %y
    end %x
    
    % save out data
    for itrial_type=1:num_trial_types
        fname_out = fullfile(b.cur_searchlight_dir, sprintf('%s_searchlight_%s.nii', b.cur_subj, trial_types{itrial_type}));

        % just pick any of the files in V to get the formatting of the
        % structure in a way that spm_write_vol likes it
        T = V{1};
        T.fname = fname_out;
        T.descrip = trial_types{itrial_type};
        spm_write_vol(T, mean_PS_by_trial_type(:,:,:,itrial_type));
    end
    
    % fisher transform data (but only if using correlations as the distance
    % metric and save this version out as well
    if distance_metric==1
        zmean_PS_by_trial_type = fisherz(mean_PS_by_trial_type);
        
        % save out data
        for itrial_type=1:num_trial_types
            fname_out = fullfile(b.cur_searchlight_dir, sprintf('z_%s_searchlight_%s.nii', b.cur_subj, trial_types{itrial_type}));

            % just pick any of the files in V to get the formatting of the
            % structure in a way that spm_write_vol likes it
            T = V{1};
            T.fname = fname_out;
            T.descrip = trial_types{itrial_type};
            spm_write_vol(T,zmean_PS_by_trial_type(:,:,:,itrial_type));
        end
    end

    % end timer for current subject
    tEnd_subj =  toc(subj_start);
    fprintf('Running the searchlight for %s took %d minutes and %f seconds.\n\n', b.cur_subj, floor(tEnd_subj/60), rem(tEnd_subj,60))
    
end

delete(myPool);

tEnd_all_searchlights = toc(script_start);
fprintf('\nRunning searchlights took %d minutes and %f seconds.\n', floor(tEnd_all_searchlights/60),rem(tEnd_all_searchlights,60))

end %function