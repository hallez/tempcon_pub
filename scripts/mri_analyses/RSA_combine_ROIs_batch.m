function []= RSA_combine_ROIs_batch()

% Create combined ROIs (eg, CA23DG) and split into head/body/tail ROIs.
% This should run relatively quickly, so no need to use `parfor`. 
%
% Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;

% set current directory - this is a SPM oddity
curdir = pwd;

roi_dirs = {'ashs_left','ashs_right','manual'};
% create a variable that stores the name of the combined ROI and the ROIs that are needed to make it up
% is there a better way of doing this???
rois(1).combined = 'CA2_3_DG.nii';
rois(1).parts = {'CA2.nii', 'CA3.nii', 'DG.nii'};
rois(1).expression = 'i1+i2+i3';

rois(2).combined = 'CA3_DG.nii';
rois(2).parts = {'CA3.nii', 'DG.nii'};
rois(2).expression = 'i1+i2';

rois(3).combined = 'whole_hippo.nii';
rois(3).parts = {'CA1.nii', 'CA2.nii', 'CA3.nii', 'DG.nii', 'subiculum.nii'};
rois(3).expression = 'i1+i2+i3+i4+i5';

num_rois = length(rois);

%% Loop across subjects
for isub=1:32
    subj_start = tic;
    
    config = initialize_temp_con;
        
    b.cur_subj = ['s' num2str(isub,'%01d')];
    b.unzipped_mri_dir = config.unzipped_mri_dir;
    b.analyzed_mri_dir = fullfile(config.analyzed_mri_dir, b.cur_subj, filesep);
    b.subj_ROI_base_dir = fullfile(b.analyzed_mri_dir,'ROIs');
    
    if ismember(isub, config.exclude_subjects)
        fprintf('\tCurrent subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end

    fprintf('\n----Working on combining ROIs for %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);
        
    if ~exist(b.subj_ROI_base_dir,'dir');
       error('ROI directory does not exist.\n')
    end 

    for idir = 1:size(roi_dirs,2)
        fprintf('\nCurrent ROI directory: %s\n', roi_dirs{idir})
        b.cur_path = fullfile(b.subj_ROI_base_dir, roi_dirs{idir});
        
        if exist(b.cur_path,'dir')
            
            for iroi=1:num_rois
                % generate the input files based on which combined ROI we're working on
                for ipart=1:length(rois(iroi).parts)
                    b.fname_in(ipart,1) = cellstr(strcat(fullfile(b.cur_path, rois(iroi).parts(ipart)), ',1'));
                end

                b.fname_out = rois(iroi).combined;

                b.expression = rois(iroi).expression;

                % specify matlabbatch variable with subject-specific inputs
                matlabbatch = batch_job(b);

                % run matlabbatch job
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
            end %for iroi
            
        end %if exist(b.cur_path

    end %idir=
    
    tEnd_subj =  toc(subj_start);
    fprintf('Combining ROIs for %s took %d minutes and %f seconds.\n', b.cur_subj, floor(tEnd_subj/60), rem(tEnd_subj,60))

end %for isub=

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
function [matlabbatch]=batch_job(b)
    matlabbatch{1}.spm.util.imcalc.input = b.fname_in;
    matlabbatch{1}.spm.util.imcalc.output = b.fname_out;
    matlabbatch{1}.spm.util.imcalc.outdir = {b.cur_path};
    matlabbatch{1}.spm.util.imcalc.expression = b.expression;
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    % uint8=2; int16=4; int32=8; float32=16; float64=64
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
end
