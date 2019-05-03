function [] = RSA_reslice_t2_and_ROIs_batch()

% Apply parameters needed to coregister and reslice T2 into mean EPI space
% to all ROIs.
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

parfor isub=1:32
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

    fprintf('\n----Working on reslicing T2 and ROIs for %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);
    
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
    
    tEnd_subj =  toc(subj_start);
    fprintf('Reslicing for %s took %d minutes and %f seconds.\n', b.cur_subj, floor(tEnd_subj/60), rem(tEnd_subj,60))

end %isub

delete(myPool);

time_end_script =  toc(script_start);
sprintf('Reslicing for all subjects took %d minutes and %f seconds.\n', floor(time_end_script/60), rem(time_end_script,60))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions

function [matlabbatch]=batch_job(b)

%% Specify files
matlabbatch{1}.cfg_basicio.file_fplist.dir = {fullfile(b.analyzed_mri_dir,'run1/')};
matlabbatch{1}.cfg_basicio.file_fplist.filter = '^meanf';
matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPList';
matlabbatch{2}.cfg_basicio.file_fplist.dir = {fullfile(b.analyzed_mri_dir,'t2_19/')};
matlabbatch{2}.cfg_basicio.file_fplist.filter = '^s';
matlabbatch{2}.cfg_basicio.file_fplist.rec = 'FPList';
matlabbatch{3}.cfg_basicio.file_fplist.dir = {fullfile(b.analyzed_mri_dir,'ROIs/')};
% based on http://stackoverflow.com/questions/2116328/regexp-matching-string-not-starting-with-my
% and http://www.regextester.com/15
% skip over ROIs that start with 'r' or 'br' (ie, don't re-reslice already
% resliced OR already binarized and resliced files)
matlabbatch{3}.cfg_basicio.file_fplist.filter = '^((?![br]).|(?![r]).)*$';
matlabbatch{3}.cfg_basicio.file_fplist.rec = 'FPListRec';

%% Coregister and reslice
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep;
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tname = 'Reference Image';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).name = 'class';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).sname = 'File Selector (Batch Mode): Selected Files (^meanf)';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).src_output = substruct('.','files');
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1) = cfg_dep;
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tname = 'Source Image';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).name = 'class';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).sname = 'File Selector (Batch Mode): Selected Files (^s)';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).src_output = substruct('.','files');
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1) = cfg_dep;
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tname = 'Other Images';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tgt_spec{1}(1).name = 'class';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).sname = 'File Selector (Batch Mode): Selected Files (^((?![br]).|(?![r]).)*$)';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).src_output = substruct('.','files');
matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

end
