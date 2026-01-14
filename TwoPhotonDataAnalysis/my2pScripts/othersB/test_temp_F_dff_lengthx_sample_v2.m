% F_dff_length1_period1_raw;
% length1_period1_interval;
temp0 = length1_period1_interval(3:4);
temp0(end) = temp0(end)-1;
temp_range = temp0(1):temp0(end);
length1_sample_interval = temp0-temp0(1)+1;
F_dff_length1_sample_raw = F_dff_length1_period1_raw(:,:,temp_range);

% F_dff_length2_period1_raw;
% length2_period1_interval;
temp0 = length2_period1_interval(3:5);
temp0(end) = temp0(end)-1;
temp_range = temp0(1):temp0(end);
length2_sample_interval = temp0-temp0(1)+1;
F_dff_length2_sample_raw = F_dff_length2_period1_raw(:,:,temp_range);


% F_dff_length3_period1_raw;
% length3_period1_interval;
temp0 = length3_period1_interval(3:6);
temp0(end) = temp0(end)-1;
temp_range = temp0(1):temp0(end);
length3_sample_interval = temp0-temp0(1)+1;
F_dff_length3_sample_raw = F_dff_length3_period1_raw(:,:,temp_range);

F_dff_length1_sample_raw;
F_dff_length2_sample_raw;
F_dff_length3_sample_raw;

F_dff_length1_sample_raw_allTrial = nan(size(F_dff_length1_sample_raw,1),trial_para.trial_count,size(F_dff_length1_sample_raw,3));
F_dff_length2_sample_raw_allTrial = nan(size(F_dff_length2_sample_raw,1),trial_para.trial_count,size(F_dff_length2_sample_raw,3));
F_dff_length3_sample_raw_allTrial = nan(size(F_dff_length3_sample_raw,1),trial_para.trial_count,size(F_dff_length3_sample_raw,3));

F_dff_length1_sample_raw_allTrial(:,trialIndex_length1_bool,:) = F_dff_length1_sample_raw;
F_dff_length2_sample_raw_allTrial(:,trialIndex_length2_bool,:) = F_dff_length2_sample_raw;
F_dff_length3_sample_raw_allTrial(:,trialIndex_length3_bool,:) = F_dff_length3_sample_raw;


a1 = squeeze(F_dff_length1_sample_raw_allTrial(:,1,:));
a2 = squeeze(F_dff_length1_sample_raw(:,1,:));
a11 = squeeze(F_dff_length1_sample_raw_allTrial(:,3,:));

F_dff_length1_sample_raw_allTrial;
F_dff_length2_sample_raw_allTrial;
F_dff_length3_sample_raw_allTrial;
length1_sample_interval;
length2_sample_interval;
length3_sample_interval;

%% End