function []= RSA_binarize_ROIs_batch()

% Read in coregistered ROIs and binarize
%
% Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;

% set current directory - this is a SPM oddity
curdir = pwd;

script_start = tic; 

num_cores = feature('numcores');
parpool(num_cores);
myPool = gcp;

%% Loop across subjects
parfor isub=1:32

    subj_start = tic;
    
    config = initialize_temp_con;
    
    b.cur_subj = ['s' num2str(isub,'%01d')];
    b.analyzed_mri_dir = fullfile(config.analyzed_mri_dir, b.cur_subj, filesep);
    b.subj_ROI_base_dir = fullfile(b.analyzed_mri_dir, 'ROIs');

    b.curModel = 'itemHitBYType_model';
    b.curModelBase = strtok(b.curModel,'_');
    b.curModelType = 'multivariate';
    b.model_dir = fullfile(config.analyzed_mri_dir, [b.curModelType, '_', b.curModelBase]);
    b.runs = {'run1','run2','run3','run4','run5','run6'};
    
    if ismember(isub, config.exclude_subjects)
        fprintf('Current subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end

    fprintf('\n----Working on binarizing ROIs for %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);
    
    if ~exist(b.subj_ROI_base_dir,'dir');
       sprintf('ROI directory does not exist.\n')
       continue;
    end % if ~exist

    roi_dirs = {'ashs_left','ashs_right','manual_traced_rois'};
    for idir = 1:size(roi_dirs,2)
        fprintf('\nCurrent ROI directory: %s\n',roi_dirs{idir})
        b.cur_path = fullfile(b.subj_ROI_base_dir,roi_dirs{idir});
        all_ROIs = dir([b.cur_path,filesep,'r*.nii']);
        all_ROIs = {all_ROIs.name};
        num_rois = length(all_ROIs);

        for iroi=1:num_rois
            cur_roi_name = fullfile(b.cur_path,all_ROIs{iroi});
            [~,NAME,EXT] = fileparts(cur_roi_name);
            b.fname_in = [NAME,EXT];
            b.fname_out = ['b',NAME,EXT];

            matlabbatch = batch_job(b);

            cd(b.analyzed_mri_dir);
            try
                spm_jobman('initcfg')
                spm('defaults', 'FMRI');
                spm_jobman('serial', matlabbatch);
            catch
                cd(curdir);
                continue;
            end %try

            cd(curdir);
        end
    end %idir=
    
    tEnd_subj =  toc(subj_start);
    sprintf('Binarizing ROIs for %s took %d minutes and %f seconds.\n', b.cur_subj, floor(tEnd_subj/60), rem(tEnd_subj,60))

end %for isub=

delete(myPool);

time_end_script =  toc(script_start);
sprintf('Binarizing ROIs for all subjects took %d minutes and %f seconds.\n', floor(time_end_script/60), rem(time_end_script,60))


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
function [matlabbatch]=batch_job(b)
    matlabbatch{1}.spm.util.imcalc.input = {fullfile(b.cur_path, b.fname_in)};
    matlabbatch{1}.spm.util.imcalc.output = b.fname_out;
    matlabbatch{1}.spm.util.imcalc.outdir = {b.cur_path};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1>0';
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    % uint8=2; int16=4; int32=8; float32=16; float64=64
    matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
end
