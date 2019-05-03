function [] = RSA_split_HBT()
% split ROIs of interest into head, body, and tail segments based on the (manually-determined) tranistion slice values
%
% Halle R. Dimsdale-Zucker

close all;
fclose('all');  
clc;

script_start = tic; 

for isub = 1:32
    
    subj_start = tic;
    
    config = initialize_temp_con;
    
    b.cur_subj = ['s' num2str(isub,'%01d')];
    b.analyzed_mri_dir = fullfile(config.analyzed_mri_dir, b.cur_subj, filesep);
    b.raw_behavioral_dir = fullfile(config.raw_behavioral_dir, b.cur_subj); 
    b.rois_to_split = {'CA1', 'CA3', 'CA2_3_DG', 'CA3_DG', 'subiculum', 'whole_hippo'};
    b.roi_dirs = {'ashs_left', 'ashs_right'};
    
    if ismember(isub, config.exclude_subjects)
        fprintf('\tCurrent subject %s is marked for exclusion. Skipping.\n', b.cur_subj)
        continue;
    end
    
    fprintf('\n----Working on splitting ROIs into head, body, tail for %s----\n',b.cur_subj)
    
    b = run_exceptions_TempCon(b);
    
    transitions_fname = fullfile(b.raw_behavioral_dir, sprintf('%s_hc_transitions.yml',b.cur_subj));
    
    if exist(transitions_fname, 'file')
        transitions = ReadYaml(transitions_fname);

        b.left_head = transitions.left.head_slice;
        b.right_head = transitions.right.head_slice;
        b.left_body_tail = transitions.left.body_tail_transition;
        b.right_body_tail = transitions.right.body_tail_transition;
    else
       warning('No transitions file for %s. Moving onto next subject', b.cur_subj)
    end

    for idir = 1:length(b.roi_dirs)
        b.cur_roi_dir = b.roi_dirs{idir};
        % this is a bit of a hack since `strtok` won't look for a delimiter
        % **string**
        b.cur_hemi = b.cur_roi_dir(strfind(b.cur_roi_dir,'_')+1:end);

        for iroi = 1:length(b.rois_to_split)
            b.cur_roi = b.rois_to_split{iroi};
            b.cur_roi_fpath = fullfile(b.analyzed_mri_dir, 'ROIs', b.cur_roi_dir);
            b.cur_roi_fname = fullfile(b.cur_roi_fpath, sprintf('%s.nii', b.cur_roi));

            if ~exist(b.cur_roi_fname, 'file')
                sprintf('%s does not exist. Skipping.\n', b.cur_roi_fname)
                continue;
            end
            
            Q = spm_vol(b.cur_roi_fname);
            QV = spm_read_vols(Q);

            if strcmp(b.cur_hemi, 'left')
                QV_head = zeros(size(QV));
                QV_head(:,:,1:b.left_head) = QV(:,:,1:b.left_head);

                QV_body = zeros(size(QV));
                QV_body(:,:,b.left_head+1:b.left_body_tail) = QV(:,:,b.left_head+1:b.left_body_tail);

                QV_tail = zeros(size(QV));
                QV_tail(:,:,b.left_body_tail+1:end) = QV(:,:,b.left_body_tail+1:end);

                QV_body_tail = zeros(size(QV));
                QV_body_tail(:,:,b.left_head+1:end) = QV(:,:,b.left_head+1:end);
            elseif strcmp(b.cur_hemi, 'right')
                QV_head = zeros(size(QV));
                QV_head(:,:,1:b.right_head) = QV(:,:,1:b.right_head);

                QV_body = zeros(size(QV));
                QV_body(:,:,b.right_head+1:b.right_body_tail) = QV(:,:,b.right_head+1:b.right_body_tail);

                QV_tail = zeros(size(QV));
                QV_tail(:,:,b.right_body_tail+1:end) = QV(:,:,b.right_body_tail+1:end);

                QV_body_tail = zeros(size(QV));
                QV_body_tail(:,:,b.right_head+1:end) = QV(:,:,b.right_head+1:end);
            end %if strcmp(b.cur_hemi

            % write out new ROIs
            % use 'Q' to get the structure in the way
            % spm_write_vol likes
            Q_head = Q;
            Q_body = Q;
            Q_tail = Q;
            Q_body_tail = Q;

            Q_head.fname = fullfile(b.cur_roi_fpath, sprintf('%s_head.nii', b.cur_roi));
            Q_body.fname = fullfile(b.cur_roi_fpath, sprintf('%s_body.nii', b.cur_roi));
            Q_tail.fname = fullfile(b.cur_roi_fpath, sprintf('%s_tail.nii', b.cur_roi));
            Q_body_tail.fname = fullfile(b.cur_roi_fpath, sprintf('%s_body_tail.nii', b.cur_roi));

            spm_write_vol(Q_head, QV_head);
            spm_write_vol(Q_body, QV_body);
            spm_write_vol(Q_tail, QV_tail);
            spm_write_vol(Q_body_tail, QV_body_tail);

            clear Q QV Q_head QV_head Q_body QV_body Q_tail QV_tail Q_body_tail QV_body_tail
        end %iroi

    end % idir
    
    tEnd_subj =  toc(subj_start);
    fprintf('Splitting ROIs into HBT for %s took %d minutes and %f seconds.\n', b.cur_subj, floor(tEnd_subj/60), rem(tEnd_subj,60))

end %isub

time_end_script =  toc(script_start);
sprintf('Splitting HBT for all subjects took %d minutes and %f seconds.\n', floor(time_end_script/60), rem(time_end_script,60))

end