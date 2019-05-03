function [] = RSA_single_trial_models_batch()

% Loop across subjects, runs, and trials to estimate single trial betas for
% multivariate (pattern similarity) analyses.
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
    parpool(pc, 8);
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
    b.curModel = 'itemHitBYType_model';
    b.curModelBase = strtok(b.curModel,'_');
    b.curModelType = 'multivariate';
    b.model_dir = fullfile(config.analyzed_mri_dir, [b.curModelType, '_', b.curModelBase]);
    b.runs = {'run1','run2','run3','run4','run5','run6'};
    
    if ismember(isub, config.exclude_subjects)
        fprintf('Current subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end

    fprintf('\n----Working on single trial modelling for %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);
    
    subj_dir_out = fullfile(b.model_dir, b.cur_subj);
    if ~isdir(subj_dir_out)
        mkdir(subj_dir_out);
    end
    
    for irun=1:length(b.runs)
        b.cur_run = b.runs{irun};

        curRunDir = fullfile(b.model_dir, b.cur_subj, b.cur_run);
        curDirs = dir(curRunDir);
        trialCounter = 0;

        for j = 1:size(curDirs,1)
            s=curDirs(j).name;
            if strfind(s,'trial_')
                trialCounter = trialCounter + 1;
            end %if strfind
        end %j=

        % 360 = total number of item recog trials per run
        % 6 = number of item recog runs
        if trialCounter <= (360/6)
            numCurTrials = 1:trialCounter;
        else 
            error('Implausible number of trials for current run: %s.\n',b.cur_run)
        end %if trialCounter

        for itrial=1:trialCounter

            b.cur_trial = ['trial_',num2str(numCurTrials(itrial),'%03d')];

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

        end %itrial=  

    end %irun=
    
    cd(curdir);
    
    tEnd_subj =  toc(subj_start);
    sprintf('Running the models for %s took %d minutes and %f seconds.\n', b.cur_subj, floor(tEnd_subj/60), rem(tEnd_subj,60))

end %parfor isub

delete(myPool);

end %function single_trial_parallel_trials
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
function [matlabbatch]=batch_job(b)
%This function generates the matlabbatch variable: can be copied in
%directly from the batch job output, then modify lines as necessary to
%generalize the paths, etc, using b variables

    %% Specfiy files
    fprintf('Figuring out files for subject %s.\n', b.cur_subj)
    matlabbatch{1}.cfg_basicio.file_fplist.dir = {fullfile(b.analyzed_mri_dir, b.cur_run, filesep)};
    matlabbatch{1}.cfg_basicio.file_fplist.filter = '^rf';
    matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{2}.cfg_basicio.file_fplist.dir = {fullfile(b.model_dir, b.cur_subj, b.cur_run, b.cur_trial, filesep)};
    matlabbatch{2}.cfg_basicio.file_fplist.filter = 'regs.mat';
    matlabbatch{2}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{3}.cfg_basicio.file_fplist.dir = {fullfile(b.analyzed_mri_dir ,b.cur_run, filesep)};
    matlabbatch{3}.cfg_basicio.file_fplist.filter = 'spike_regs_rp.txt';
    matlabbatch{3}.cfg_basicio.file_fplist.rec = 'FPList';
    
    %% Run model
    matlabbatch{4}.spm.stats.fmri_spec.dir = {fullfile(b.model_dir, b.cur_subj, b.cur_run, b.cur_trial, filesep)};
    matlabbatch{4}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{4}.spm.stats.fmri_spec.timing.RT = 2.01;
    matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep;
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tname = 'Scans';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).sname = 'File Selector (Batch Mode): Selected Files (^rf)';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).src_output = substruct('.','files');
    matlabbatch{4}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1) = cfg_dep;
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tname = 'Multiple conditions';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).sname = 'File Selector (Batch Mode): Selected Files (regs.mat)';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).src_output = substruct('.','files');
    matlabbatch{4}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1) = cfg_dep;
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tname = 'Multiple regressors';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).sname = 'File Selector (Batch Mode): Selected Files (spike_regs_rp.txt)';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).src_output = substruct('.','files');
    matlabbatch{4}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{4}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{4}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{4}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{4}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{4}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{4}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% Estimate model
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
    matlabbatch{5}.spm.stats.fmri_est.method.Classical = 1;
    
end %batch_job
