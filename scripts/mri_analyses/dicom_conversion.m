function [] = dicom_conversion()

% Convert .dcm to .nii
%
% Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;

% set current directory - this is a SPM oddity
curdir = pwd;
num_cores = feature('numcores');
parpool(num_cores);
myPool = gcp;

parfor isub = 1:32 
    subj_start = tic;
    
    config = initialize_temp_con;
        
    b.cur_subj = ['s' num2str(isub,'%01d')];
    b.unzipped_mri_dir = config.unzipped_mri_dir;
    b.analyzed_mri_dir = fullfile(config.analyzed_mri_dir, b.cur_subj, filesep);
    
    if ismember(isub, config.exclude_subjects)
        fprintf('\tCurrent subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end

    fprintf('\n----Working on dicom conversion for %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);
    
    corrupt_zip_subjects = 2:9;
    if ismember(b.cur_subj, corrupt_zip_subjects)
        fprintf('\tCurrent subject has a corrupt zip file. Skipping dicom conversion.\n')
        continue;
    end
    
    ds_dir_info = dir(fullfile(b.unzipped_mri_dir, b.cur_subj, 'ds*/'));
    % this handles the case where there is no 'ds*' directory
    % when there is no directory, this creates a 0x1 struct array hence why
    % checking for size(<var>,1) handles this case
    if ~size(ds_dir_info,1)
        warning('\tRaw ds* directory does not exist. Skipping.\n')
        continue;
    end
    
    if strcmp(b.cur_subj, 's1') && (size(ds_dir_info,1) == 2)
        ds_dir_info = ds_dir_info(1);
    end
    
    b.unzipped_mri_dir = fullfile(b.unzipped_mri_dir, b.cur_subj, ds_dir_info.name, filesep);
    
    if ~isdir(b.unzipped_mri_dir)
       warning('\tUnzipped mri directory does not exist. Skipping.\n')
       continue;
    end
    
    if ~isdir(b.analyzed_mri_dir)
        fprintf('\tCreating analyzed mri dir: %s\n', b.analyzed_mri_dir)
        mkdir(b.analyzed_mri_dir)
    end

    cur_dirs = dir(b.unzipped_mri_dir);
    for idir = 1:size(cur_dirs,1)
        % check that the current value is a directory
        % and also check that it's not just a "dot" directory
        % based on: https://www.mathworks.com/matlabcentral/newsreader/view_thread/258220
        if ~cur_dirs(idir).isdir || strcmpi(cur_dirs(idir).name, '.') || strcmpi(cur_dirs(idir).name, '..')
            continue;
        end
        
        b.cur_dir_name = cur_dirs(idir).name;

        matlabbatch = batch_job(b);

        cd(b.analyzed_mri_dir);
        try
            spm_jobman('initcfg')
            spm('defaults', 'FMRI');
            spm_jobman('serial', matlabbatch);
        catch
            cd(curdir);
            continue;
        end 
    end %idir

    cd(curdir);
            
    tEnd_subj =  toc(subj_start);
    fprintf('Dicom conversion for %s took %d minutes and %f seconds.\n', b.cur_subj, floor(tEnd_subj/60), rem(tEnd_subj,60))
    
end %for isub=

delete(myPool);

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions

function [matlabbatch]=batch_job(b)
    fprintf('\tFiguring out files in %s.\n', b.cur_dir_name)
    matlabbatch{1}.cfg_basicio.file_fplist.dir = {fullfile(b.unzipped_mri_dir, b.cur_dir_name)};
    matlabbatch{1}.cfg_basicio.file_fplist.filter = 'dcm';
    matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPList';

    %% Dicom conversion 
    fprintf('\tStarting dicom conversion.\n')
    %% Specfiy files
    matlabbatch{2}.spm.util.dicom.data(1) = cfg_dep;
    matlabbatch{2}.spm.util.dicom.data(1).tname = 'DICOM files';
    matlabbatch{2}.spm.util.dicom.data(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{2}.spm.util.dicom.data(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{2}.spm.util.dicom.data(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{2}.spm.util.dicom.data(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{2}.spm.util.dicom.data(1).sname = 'File Selector (Batch Mode): Selected Files (dcm)';
    matlabbatch{2}.spm.util.dicom.data(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.spm.util.dicom.data(1).src_output = substruct('.','files');
    matlabbatch{2}.spm.util.dicom.root = 'series'; %other option is 'flat'
    matlabbatch{2}.spm.util.dicom.outdir = {b.analyzed_mri_dir};
    matlabbatch{2}.spm.util.dicom.convopts.format = 'nii';
    matlabbatch{2}.spm.util.dicom.convopts.icedims = 0;
end
