function [] = preproc_moco_quickcoreg()

% Runs basic preprocessing steps including motion correction (realign &
% reslice) and quick coreg (ie, not good enough if using high res, but
% good enough to display results on anatomicals for visualization purposes)
%
% Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;

% set current directory - this is a SPM oddity
curdir = pwd;
SLURM_FLAG = 0;

if SLURM_FLAG
    pc = parcluster('local');
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
else
    num_cores = feature('numcores');
    parpool(num_cores);
end
myPool = gcp;

parfor isub = 1:32
    subj_start = tic;
    
    config = initialize_temp_con;
    
    b.cur_subj = ['s' num2str(isub,'%01d')];
    b.analyzed_mri_dir = fullfile(config.analyzed_mri_dir, b.cur_subj, filesep);

    if ismember(isub, config.exclude_subjects)
        fprintf('\tCurrent subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end

    fprintf('\n----Working on preprocessing %s----\n', b.cur_subj)

    b = run_exceptions_TempCon(b);

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

    fprintf('Done with preprocessing and motion correction.')
    time_end_subj =  toc(subj_start);
    sprintf('Running preproc for %s took %d minutes and %f seconds.\n', b.cur_subj, floor(time_end_subj/60), rem(time_end_subj,60))
    
end %isub

delete(myPool);

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions

function [matlabbatch]=batch_job(b)
%This function generates the matlabbatch variable: can be copied in
%directly from the batch job output, then modify lines as necessary to
%generalize the paths, etc, using b variables

    %% Specfiy files
    fprintf('\tFiguring out files for subject %s.\n',b.cur_subj)
    matlabbatch{1}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep,'mprage/']};
    matlabbatch{1}.cfg_basicio.file_fplist.filter = '^s';
    matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{2}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep, 't2_19/']};
    matlabbatch{2}.cfg_basicio.file_fplist.filter = '^s';
    matlabbatch{2}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{3}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep,'run1/']};
    matlabbatch{3}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{3}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{4}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep,'run2/']};
    matlabbatch{4}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{4}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{5}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep,'run3/']};
    matlabbatch{5}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{5}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{6}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep,'run4/']};
    matlabbatch{6}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{6}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{7}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep,'run4/']};
    matlabbatch{7}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{7}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{8}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep,'run5/']};
    matlabbatch{8}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{8}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{9}.cfg_basicio.file_fplist.dir = {[b.analyzed_mri_dir,filesep,'run6/']};
    matlabbatch{9}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{9}.cfg_basicio.file_fplist.rec = 'FPList';

    %% Motion correction (realignment + reslice)
    fprintf('\tDoing motion correction for subject %s.\n',b.cur_subj)
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep;
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1).tname = 'Session';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
    matlabbatch{10}.spm.spatial.realign.estwrite.data{1}(1).src_output = substruct('.','files');
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep;
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1).tname = 'Session';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1});
    matlabbatch{10}.spm.spatial.realign.estwrite.data{2}(1).src_output = substruct('.','files');
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1) = cfg_dep;
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1).tname = 'Session';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1});
    matlabbatch{10}.spm.spatial.realign.estwrite.data{3}(1).src_output = substruct('.','files');
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1) = cfg_dep;
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1).tname = 'Session';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1});
    matlabbatch{10}.spm.spatial.realign.estwrite.data{4}(1).src_output = substruct('.','files');
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1) = cfg_dep;
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1).tname = 'Session';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1});
    matlabbatch{10}.spm.spatial.realign.estwrite.data{5}(1).src_output = substruct('.','files');
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1) = cfg_dep;
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1).tname = 'Session';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1).src_exbranch = substruct('.','val', '{}',{8}, '.','val', '{}',{1});
    matlabbatch{10}.spm.spatial.realign.estwrite.data{6}(1).src_output = substruct('.','files');
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1) = cfg_dep;
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1).tname = 'Session';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1).src_exbranch = substruct('.','val', '{}',{9}, '.','val', '{}',{1});
    matlabbatch{10}.spm.spatial.realign.estwrite.data{7}(1).src_output = substruct('.','files');
    matlabbatch{10}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{10}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{10}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{10}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{10}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{10}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{10}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{10}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{10}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{10}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{10}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{10}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

    %% Quick n dirty coreg
    fprintf('\tDoing coreg with MPRAGE to allow for displaying results for subject %s.\n',b.cur_subj)
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep;
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1).tname = 'Reference Image';
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1).sname = 'Realign: Estimate & Reslice: Mean Image';
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1).src_exbranch = substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{11}.spm.spatial.coreg.estwrite.ref(1).src_output = substruct('.','rmean');
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1) = cfg_dep;
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1).tname = 'Source Image';
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1).sname = 'File Selector (Batch Mode): Selected Files (^s)';
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{11}.spm.spatial.coreg.estwrite.source(1).src_output = substruct('.','files');
    matlabbatch{11}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{11}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{11}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{11}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{11}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{11}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{11}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{11}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{11}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

    fprintf('\tDoing coreg with T2 to allow for displaying results for subject %s.\n',b.cur_subj)
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep;
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1).tname = 'Reference Image';
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1).sname = 'Realign: Estimate & Reslice: Mean Image';
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1).src_exbranch = substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{12}.spm.spatial.coreg.estwrite.ref(1).src_output = substruct('.','rmean');
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1) = cfg_dep;
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1).tname = 'Source Image';
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1).sname = 'File Selector (Batch Mode): Selected Files (^s)';
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
    matlabbatch{12}.spm.spatial.coreg.estwrite.source(1).src_output = substruct('.','files');
    matlabbatch{12}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{12}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{12}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{12}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{12}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{12}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{12}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{12}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{12}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

end
