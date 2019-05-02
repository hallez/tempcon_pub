function [config] = initialize_temp_con()
    % Initialize script for TempCon directories
    % Author: Halle R. Dimsdale-Zucker
    
    SLURM_FLAG = 0;
    if SLURM_FLAG 
        project_dir = <absolute-path-as-string>;
    else
        project_dir = <absolute-path-as-string>;
    end
    addpath(genpath(fullfile(project_dir,'vendor','yamlmatlab','0.4.3')));
    addpath(genpath(fullfile(project_dir, 'vendor')))
    config_file = ReadYaml(fullfile(project_dir,'config.yml'));
    config.vendor_dir = config_file.directories.vendor_dir; 

    % data directories
    config.base_dir = config_file.directories.base_dir;
    config.raw_behavioral_dir = fullfile(config_file.directories.raw_behavioral);
    config.raw_mri_dir = config_file.directories.raw_mri;
    config.analyzed_behavioral_dir = fullfile(config_file.directories.analyzed_behavioral);
    config.unzipped_mri_dir = fullfile(config_file.directories.local_unzipped_dir);
    config.analyzed_mri_dir = config_file.directories.analyzed_mri;

    % script directories
    config.scripts_mri_analyses = fullfile(config_file.directories.mri_analysis_scripts);
    config.scripts_behavioral_analyses = fullfile(config_file.directories.behavioral_analysis_scripts);
    
    % notes and diaries directories
    config.diaries_dir = fullfile(project_dir, 'diary_files');
    
    %%% Subjects
    config.exclude_subjects = [1, 2, 3, 11, 13, 14, 24, 26];  
    config.backroom_subjects = [11, 12, 14:23, 26:28, 31];
end 
