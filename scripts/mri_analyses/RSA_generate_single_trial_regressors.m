% This script generates regressors for running single trial event models. 
%
% Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;

config = initialize_temp_con;

%% MODEL INFO
curModel = 'itemHitBYType_model';
curModelBase = strtok(curModel,'_');
curModelType = 'multivariate';

outputDir = [config.analyzed_mri_dir, filesep, curModelType '_' curModelBase];
if ~isdir(outputDir)
    fprintf('Creating model directory: %s\n',curModel)
    mkdir(outputDir)
end %if isdir(

%% VARIABLES
excludereg = 99;

%% OPTIONS
write_reg = 1; %do you want to write out regressor mats?
include_dur = 0; %do you want to include durations? 

%% END OF MODIFICATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(strcat(outputDir,filesep,'model_info'));

for isub = 1:32
    b.cur_subj = ['s' num2str(isub,'%01d')];  
    b.raw_mri_dir = fullfile(config.raw_mri_dir, b.cur_subj, filesep);
    b.runs = {'run1','run2','run3','run4','run5','run6'};
    
    if ismember(isub, config.exclude_subjects)
        fprintf('Current subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end
    
    fprintf('\n----Generating single trial regressors for subject %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);

    subj_dir_out = [outputDir,filesep, b.cur_subj];
    if ~isdir(subj_dir_out)
        mkdir(subj_dir_out);
    end
    
    diaryname = fullfile(subj_dir_out,filesep,'RSA_single_trial_regressors_diary.txt');
    diary(diaryname);
            
    datafile = fullfile(config.analyzed_behavioral_dir, b.cur_subj,'onsets_file.csv');
    if(exist(datafile,'file'))
        conditions = tdfread(datafile,','); 
        
        if ~isfield(conditions,curModel) 
            error('There is no column for the current model (%s) in onsets file (%s).\n',curModel,datafile);
            continue;
        end %if 
        
        %%% Generate regressor files
        for irun = 1:length(b.runs)

            cur_run = b.runs{irun};
            cur_run_num = str2double(strtok(cur_run, 'run'));

            runDir = fullfile(outputDir,b.cur_subj,cur_run);
            if ~isdir(runDir)
                    fprintf('Creating %s directory.\n',cur_run)
                    mkdir(runDir)
            end %if isdir(

            % generate variable to represent where current functional run's
            % data is stored
            curfuncrun = fullfile(config.analyzed_mri_dir,b.cur_subj,cur_run);

            % create a logical index for the current run (using
            % conditions.item_recog_run)
            curRunIdx = conditions.item_recog_run == cur_run_num;
            
            if(sum(curRunIdx) ~= 60)
                warning('Incorrect number of trials for %s. Should be 60 and instead is %d.\n', cur_run, sum(curRunIdx))
            end

            % subset information w/in conditions struct for the current run
            % using this logical index
            % NB: a more standard Matlab way to do this would be to have a
            % matrix where columns represent each of these variables (would
            % then also need a key to know which columns represent which
            % variables)
            curConditions.onset = conditions.onset(curRunIdx);
            curConditions.item_recog_row_number = conditions.item_recog_row_number(curRunIdx);
            curConditions.temporal_scored = cellstr(conditions.temporal_scored(curRunIdx,:));
            curConditions.question_scored_exact = cellstr(conditions.question_scored_exact(curRunIdx,:));
            curConditions.itemHitBYType_model = conditions.itemHitBYType_model(curRunIdx);
            curConditions.itemHitBYType_key = cellstr(conditions.itemHitBYType_key(curRunIdx,:));
            curConditions.duration = conditions.duration(curRunIdx);
            numCurTrials = length(unique(curConditions.item_recog_row_number));

            for itrial = 1:numCurTrials

                % create an index for the current trial
                % remember, rows in onset_file.csv are NOT (necessarily) ordered
                % by trial number which is why need an index
                % we're using itrial + 1 because the fixation cross got
                % labelled as trial 1 in the output file
                curTrialIdx = curConditions.item_recog_row_number == itrial+1;

                if curConditions.itemHitBYType_model(curTrialIdx) <= excludereg

                    trialDir = [runDir,filesep,'trial_',num2str(itrial,'%03d')];
                    if ~isdir(trialDir)
                        mkdir(trialDir)
                    end %if isdir(

                    filename = [trialDir,filesep,'regs'];

                    onsets = {}; names = {}; durations = {};

                    if include_dur 
                        add_dur_cur_trial = curConditions.duration(curTrialIdx);
                        add_dur_all_other_trials = curConditions.duration(~curTrialIdx)';
                    else
                        add_dur_cur_trial = 0;
                        add_dur_all_other_trials = zeros(size(curConditions.duration(~curTrialIdx)))';
                    end

                    % figure out onsets
                    % set the first value to be a cell and then will group
                    % values for each regressor together (this is a hack)
                    onsets = [num2cell(curConditions.onset(curTrialIdx))];
                    onsets = [onsets curConditions.onset(~curTrialIdx)'];

                    % figure out durations
                    % use the same trick as w/ onsets to get information about
                    % regressors to stay together
                    durations = [num2cell(add_dur_cur_trial)];
                    durations = [durations add_dur_all_other_trials]; 

                    % append name of current regressor 
                    cur_trial_name = ['run',num2str(cur_run_num,'%02d'),...
                                      '_trial',num2str(itrial,'%03d'),...
                                      '_',curConditions.itemHitBYType_key{curTrialIdx},...
                                      '_',curConditions.temporal_scored{curTrialIdx},...
                                      '_',curConditions.question_scored_exact{curTrialIdx}];
                    names = [cellstr(cur_trial_name) cellstr(['NOT_' cur_trial_name])]; 

                    % save out the regressor file
                    if write_reg 
                        if exist(filename,'file')
                            regsOverwrite = input('Regs file already exists. Overwrite? (Y=1,N=0):');
                            if regsOverwrite
                                save(filename,'names','onsets','durations');
                            else
                                save([filename '+'],'names','onsets','durations'); 
                            end %if regsOverwrite
                        else
                            save(filename,'names','onsets','durations');
                        end %if exist(filename
                    end %if write_reg

                end %if curConditions.itemHitBYType_model(itrial)

            end %itrial

        end %irun
    end % if exist(datafile
    
    diary off
     
end %isub