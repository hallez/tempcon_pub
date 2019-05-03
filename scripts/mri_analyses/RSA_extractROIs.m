function [] = RSA_extractROIs()

% Code to read in traced ROIs (manual or ASHS) and extract ROIs (based on
% values in labels file) and save out as individual files.
%
% Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;

num_cores = feature('numcores');
parpool(num_cores);
myPool = gcp;

%% Loop across subjects
parfor isub=1:32
    
    subj_start = tic;
    
    config = initialize_temp_con;
        
    b.cur_subj = ['s' num2str(isub,'%01d')];
    b.unzipped_mri_dir = config.unzipped_mri_dir;
    b.analyzed_mri_dir = fullfile(config.analyzed_mri_dir, b.cur_subj, filesep);
    b.subj_ROI_base_dir = fullfile(b.analyzed_mri_dir, 'ROIs');

    if ismember(isub, config.exclude_subjects)
        fprintf('\tCurrent subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end

    fprintf('\n----Working on extracting ROIs for %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);
    
    if ~exist(b.subj_ROI_base_dir,'dir');
       fprintf('Making base ROIs directory.\n')
       mkdir(b.subj_ROI_base_dir);
    end % if ~exist
    
    ASHS_FLAG = 0;
    MANUAL_FLAG = 1;

    T2_19_FLAG = 1;

    hemis = {'left', 'right'};    

    if ASHS_FLAG
      for ihemi = 1:length(hemis)
        % this will loop over left and right hemis
        cur_hemi = hemis{ihemi};

        % define where split ROIs will be written out to
        b.cur_ROI_dir = fullfile(b.subj_ROI_base_dir,['ashs_',cur_hemi]);

        % define where to-be-split file lives
        b.ashs_rois_dir = fullfile(config.analyzed_mri_dir,b.cur_subj,'ashs','final');

        % define name of the to-be split file
        b.nii_file = fullfile(b.ashs_rois_dir,[b.cur_subj, '_', cur_hemi,'_lfseg_corr_usegray.nii']);

        fprintf('b.nii_file is %s\n', b.nii_file)

        % read in the labels file
        b.lbls=tdfread(fullfile(config.base_dir,'snaplabels_forMatlab.txt'),',');

        % Actually run the extraction
        extract(b);
      end %for ihemi

    elseif MANUAL_FLAG
        % define where split ROIs will be written out to
        b.cur_ROI_dir = fullfile(b.subj_ROI_base_dir,'manual_traced_rois');
        
        % generate a list of to-be-split files since manual tracings were
        % done with one roi in each file
        manual_roi_files = dir(fullfile(b.cur_ROI_dir, 's*.nii*'));
        manual_roi_fnames = {manual_roi_files.name};

        for ifile = 1:length(manual_roi_fnames)
 
            % define name of the to-be split file
            cur_fname = manual_roi_fnames{ifile};
            
            % extract function is expecting file names that end in .nii so
            % make sure that is what gets passed
            [~, b.cur_fname_base, b.cur_fname_ext] = fileparts(cur_fname);
            if(strcmp('.gz', b.cur_fname_ext))
                cur_fname = b.cur_fname_base;
            end
            
            b.nii_file=fullfile(b.cur_ROI_dir, cur_fname);
            
            % read in the labels file
            % NB: in the alEC/pmEC tracings, PRC is also present. however,
            % it is tagged with different labels than the PRC tracings. in
            % case need for future reference, here are the label id values:
            % IDX,R,G,B,A,VIS,MSH,LABEL
            % 1,46,139,87,1,1,1,PRC_L
            % 2,255,228,225,1,1,1,PRC_R
    
            b.lbls=tdfread(fullfile(config.base_dir,'snaplabels_manualtracings_formatlab.txt'),',');

            % Actually run the extraction
            extract_manual_traced(b);
        end

    end %if ASHS_FLAG
    
    tEnd_subj =  toc(subj_start);
    fprintf('\nExtracting ROI for %s took %d minutes and %f seconds.\n', b.cur_subj, floor(tEnd_subj/60), rem(tEnd_subj,60))

end %for isub=
delete(myPool);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
function []=extract(b)

    if ~exist(b.cur_ROI_dir,'dir');
       fprintf('Making current ROI directory.\n')
       mkdir(b.cur_ROI_dir);
    end % if ~exist

    % Unzip the .nii.gz file
    % Should create a .nii file with the same name in the pwd
    if exist(b.nii_file,'file')==0
        fprintf('Unzipping %s\n',b.nii_file)
        gunzip([b.nii_file,'.gz']);
    end %if ~exist

    % Read in the .nii file
    b.all_ROIs=spm_vol(b.nii_file);
    b.all_ROIs.img = spm_read_vols(b.all_ROIs);

    % Get info about the labels
    % 'stable' will return them in the order they occur in the file; for
    % ITK-SNAP output, not an issue b/c ordered in descending numerical
    % order (but just a precaution)
    b.lbl_ids=unique(b.lbls.IDX,'stable');
    b.nlbls=length(b.lbl_ids);

    % Do a sanity check that the image contains all of the values from the ROI
    % labels
    b.lbl_idx=ismember(b.lbl_ids,unique(b.all_ROIs.img,'stable'));

    for ilbl=1:b.nlbls

        % skip over labels that don't exist
        if sum(strcmp(cellstr(b.lbls.LABEL(ilbl,:)),cellstr(b.lbls.LABEL(~b.lbl_idx,:))))>0
            fprintf('Skipping %s. Does not exist in current file.\n',(b.lbls.LABEL(ilbl,:)))
        else
            % set the current label (name, index value)
            cur_lbl_name=b.lbls.LABEL(ilbl,:);
            cur_lbl_idx=b.lbls.IDX(ilbl);

            fprintf('Working on %s.\n',cur_lbl_name);

            %deal with whitespace in cur_lbl_name
            for i=1:length(cur_lbl_name)
                if cur_lbl_name(i)==' '
                    cur_lbl_name=cur_lbl_name(1:i-1);
                    break;
                end %if cur_lbl_name(i)==' '
            end %i=1:length(cur_lbl_name)

            % extract the voxels from b.all_ROIs that match the current label's index
            % value
            cur_ROI=b.all_ROIs.img==cur_lbl_idx;

            % save out new ROI but make sure don't overwrite if already exists
            % If need to check:
            % imagesc(cur_ROI_nii(:,:,i) where i iterates through the slices in the
            % z dimension
            if exist([b.cur_ROI_dir,filesep,cur_lbl_name,'.nii'],'file')==0
                fname = fullfile(b.cur_ROI_dir,[cur_lbl_name,'.nii']);
                tmp_b_all_rois = b.all_ROIs;
                tmp_b_all_rois.fname = fname;
                spm_write_vol(tmp_b_all_rois,cur_ROI);
            end %if exist
        end %sum(strcmp(

    end %ilbl=

end

function []=extract_manual_traced(b)

    if ~exist(b.cur_ROI_dir,'dir');
       fprintf('Making current ROI directory.\n')
       mkdir(b.cur_ROI_dir);
    end % if ~exist
    
    fprintf('\nCurrent ROI: %s.\n', b.nii_file);

    % Unzip the .nii.gz file
    % Should create a .nii file with the same name in the pwd
    if exist(b.nii_file,'file')==0
        fprintf('Unzipping %s\n',b.nii_file)
        gunzip([b.nii_file,'.gz']);
    end %if ~exist

    % Read in the .nii file
    b.all_ROIs=spm_vol(b.nii_file);
    b.all_ROIs.img = spm_read_vols(b.all_ROIs);

    % Get info about the labels
    % 'stable' will return them in the order they occur in the file; for
    % ITK-SNAP output, not an issue b/c ordered in descending numerical
    % order (but just a precaution)
    b.lbl_ids=unique(b.lbls.IDX,'stable');

    % for manual traced rois, most of them have a single label in each
    % tracing file. only alEC/pmEC have multiple labels in each tracing
    all_labels = unique(b.all_ROIs.img);
    b.nlbls=length(all_labels);

    for ilbl=1:b.nlbls
        
        b.cur_lbl_id = all_labels(ilbl);
        b.lbl_idx=ismember(b.lbl_ids, b.cur_lbl_id);

        % skip over labels that don't exist
        if sum(b.lbl_idx) <= 0
            cur_lbl_name=b.lbls.LABEL(b.lbl_idx,:);
            if isempty(cur_lbl_name)
                fprintf('Skipping label %d. Does not exist as a possible tracing label.\n', b.cur_lbl_id)
            else
                fprintf('Skipping %s. Does not exist in current file.\n', cur_lbl_name)
            end
            
        else
            % set the current label (name, index value)
            cur_lbl_name=b.lbls.LABEL(b.lbl_idx,:);
            cur_lbl_idx=b.lbls.IDX(b.lbl_idx);

            fprintf('Working on %s.\n',cur_lbl_name);

            %deal with whitespace in cur_lbl_name
            for i=1:length(cur_lbl_name)
                if cur_lbl_name(i)==' '
                    cur_lbl_name=cur_lbl_name(1:i-1);
                    break;
                end %if cur_lbl_name(i)==' '
            end %i=1:length(cur_lbl_name)

            % extract the voxels from b.all_ROIs that match the current label's index
            % value
            cur_ROI=b.all_ROIs.img==cur_lbl_idx;

            % save out new ROI but make sure don't overwrite if already exists
            % If need to check:
            % imagesc(cur_ROI_nii(:,:,i) where i iterates through the slices in the
            % z dimension
            if exist([b.cur_ROI_dir,filesep,cur_lbl_name,'.nii'],'file')==0
                fname = fullfile(b.cur_ROI_dir,[cur_lbl_name,'.nii']);
                tmp_b_all_rois = b.all_ROIs;
                tmp_b_all_rois.fname = fname;
                spm_write_vol(tmp_b_all_rois,cur_ROI);
            end %if exist
        end %sum(strcmp(

    end %ilbl=

end
