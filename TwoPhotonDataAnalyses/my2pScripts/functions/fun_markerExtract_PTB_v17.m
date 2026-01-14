function output = fun_markerExtract_PTB_v17(fullMarkerName,fullPTBName,edf)

% fullMarkerName = 'D:\twoPhotonRawData\113Recording_20230109A_Ding_Site09\2023-01-09-PTB-marker-001-fromMarkerCount3.MKdata';
% load D:\twoPhotonRawData\113Recording_20230109A_Ding_Site09\2023-01-09trainDing_Step15_eye_for113Record_niMarker_15-1.mat

load(fullPTBName,'basic_para','trial_para')

marker_PTB_fullName = fullMarkerName;
[~,temp_fname,temp_ext] = fileparts(marker_PTB_fullName);
marker_PTB_fileName = [temp_fname,temp_ext];

%% Load marker file
fid_PTB = fopen(marker_PTB_fullName,'r');
current_seek = ftell(fid_PTB);
fseek(fid_PTB,0,'eof');
file_length_PTB_raw = ftell(fid_PTB);
fseek(fid_PTB,current_seek,'bof');
[marker_data_PTB, file_length_PTB] = fread(fid_PTB,Inf,'*uint8');
fprintf("Load %s ....\n",marker_PTB_fileName);
% fprintf("------  Load %s ....\n",marker_PTB_fileName);
if file_length_PTB == file_length_PTB_raw
    clear file_length_PTB_raw
else
    fprintf('file_length_PTB_raw is different from file_length_PTB!\n');
end
fclose(fid_PTB);

a = 1;

temp_str_start = strfind(fullMarkerName,'fromMarkerCount')+length('fromMarkerCount');
temp_ok_flag = 0;
temp_index = temp_str_start;
while temp_ok_flag == 0
    temp_char = fullMarkerName(temp_index);
    if temp_char >= uint8('0') && temp_char <= uint8('9')
        temp_index = temp_index + 1;        
    else
        temp_ok_flag = 1;
        temp_index = temp_index - 1;
    end    
end
temp_str_end = temp_index;

% fromMarkerCount = str2double(fullMarkerName(strfind(fullMarkerName,'fromMarkerCount')+length('fromMarkerCount')));
fromMarkerCount = str2double(fullMarkerName(temp_str_start:temp_str_end));

trial_marker_count_ground = 0;
firstTRIALID_markerIndex = 0;
for tempi=(fromMarkerCount+1):trial_para.marker.count
    k = strfind(trial_para.marker.content{tempi}, 'TRIALID ');% if find, k shouldn't be empty
    if k == 1
        if firstTRIALID_markerIndex == 0
            firstTRIALID_markerIndex = tempi;
        end        
        trial_marker_count_ground = trial_marker_count_ground + 1;
    end
end
fprintf("PTB theoretically recorded total_marker_count (from %d) = %d, ",trial_para.marker.count,trial_para.marker.count-firstTRIALID_markerIndex+1);
fprintf("trial_marker_count = %d.\n", trial_marker_count_ground);

%% Parse marker file
temp_ok_flag = 0;
validStartSampleIndex_fromENDMarker = 1;
while temp_ok_flag == 0
    high_count_continous_all = zeros(1,ceil(file_length_PTB/10));
    high_count_continous = 0;
    high_count = 0;
    sample_count = 0;
    trial_marker_count = 0;
    innerTrial_marker_count = 0;
    end_marker_count = 0;
    total_marker_count = 0; %#ok<*NASGU>
    PTB_marker_boolIndex = false(1,file_length_PTB*8);
    PTB_marker_sampleIndex = zeros(1,ceil(file_length_PTB/10));
    PTB_marker_TRIALID_sampleIndex = zeros(1,ceil(file_length_PTB/10));
    PTB_marker_END_sampleIndex = [];
    temp_count = 0;
    
    for tempi=1:file_length_PTB
        for tempj=1:8            
            sample_count = sample_count + 1;
            if sample_count < validStartSampleIndex_fromENDMarker
                continue
            end
            
            if bitget(marker_data_PTB(tempi),tempj) == 1
                high_count_continous = high_count_continous + 1;
                high_count = high_count + 1;
            else
                if high_count_continous > 0
                    %fprintf("There are %d high samples! (sample_count %d to sample_count %d).\n",...
                    %    high_count_continous, sample_count-high_count_continous, sample_count-1);
                    temp_count = temp_count + 1;
                    high_count_continous_all(temp_count) = high_count_continous;
                    current_sampleRange = sample_count-high_count_continous+1:sample_count;
                    
                    valid_markerFlag = 0;
                    %if high_count_continous >= 9 && high_count_continous <= 15
                    %if high_count_continous >= 11 && high_count_continous <= 13
                    if high_count_continous >= 100 && high_count_continous <= 140
                        valid_markerFlag = 1;
                        trial_marker_count = trial_marker_count + 1;
                        PTB_marker_TRIALID_sampleIndex(trial_marker_count) = current_sampleRange(1);
                    end
                    %if high_count_continous >= 3 && high_count_continous <= 8
                    %if high_count_continous >= 4 && high_count_continous <= 8
                    if high_count_continous >= 40 && high_count_continous <= 80
                        if trial_marker_count > 0
                            valid_markerFlag = 1;
                            innerTrial_marker_count = innerTrial_marker_count + 1;
                        end
                    end
                    %if high_count_continous > 150
                    if high_count_continous > 350
                        valid_markerFlag = 1;
                        end_marker_count = end_marker_count + 1;
                        PTB_marker_END_sampleIndex = [PTB_marker_END_sampleIndex current_sampleRange(end)];
                    end
                    if valid_markerFlag == 1
                        total_marker_count = total_marker_count + 1;
                        PTB_marker_boolIndex(current_sampleRange) = true;
                        PTB_marker_sampleIndex(total_marker_count) = current_sampleRange(1);
                    end
                    high_count_continous = 0;
                    
                end
            end
            
        end
    end
    
    if length(PTB_marker_END_sampleIndex) == 1
        temp_ok_flag = 1;
    else
        validStartSampleIndex_fromENDMarker = PTB_marker_END_sampleIndex(end-1) + 1;
    end
    
end


PTB_marker_sampleIndex(PTB_marker_sampleIndex==0) = [];
PTB_marker_TRIALID_sampleIndex(PTB_marker_TRIALID_sampleIndex==0) = [];
high_count_continous_all(high_count_continous_all==0) = [];

fprintf("total_marker_count = %d, ", total_marker_count);
fprintf("trial_marker_count = %d, ", trial_marker_count);
fprintf("innerTrial_marker_count = %d, ", innerTrial_marker_count);
fprintf("end_marker_count = %d.\n", end_marker_count);



if false
    PTB_marker_content_trial_para = trial_para.marker.content; %#ok<*UNRCH>
    
    % PTB_marker_content_edf = edf.Messages.info(11:end);
    PTB_marker_sampleIndex_edf = edf.Messages.time(11:end);
    
    PTB_marker_sampleIndex_edf4K_raw = PTB_marker_sampleIndex_edf*4;
    PTB_marker_sampleIndex_edf4K = PTB_marker_sampleIndex_edf4K_raw-PTB_marker_sampleIndex_edf4K_raw(1)+PTB_marker_sampleIndex(1);
    
    temp_diff1 = diff(PTB_marker_sampleIndex);
    temp_diff2 = diff(PTB_marker_sampleIndex_edf4K);
    temp_diff12 = [temp_diff1 temp_diff1(end)] - temp_diff2;
    
    PTB_marker_sampleIndex; %#ok<*VUNUS>
    PTB_marker_content_trial_para;
    
    compensate_sampleIndex = 137;
    compensate_interval = 44;
    
    PTB_marker_sampleIndex_compensate = nan(1,length(PTB_marker_content_trial_para));
    PTB_marker_sampleIndex_compensate(1:(compensate_sampleIndex-1)) = PTB_marker_sampleIndex(1:(compensate_sampleIndex-1));
    
    PTB_marker_sampleIndex_compensate(compensate_sampleIndex) = PTB_marker_sampleIndex(compensate_sampleIndex-1) + compensate_interval;
    
    PTB_marker_sampleIndex_compensate((compensate_sampleIndex+1):end) = PTB_marker_sampleIndex(compensate_sampleIndex:end);
    
    temp_diff1_compensate = diff(PTB_marker_sampleIndex_compensate);
    
    
    PTB_marker_sampleIndex;
    PTB_marker_sampleIndex_compensate;
    
    PTB_marker_sampleIndex_raw = PTB_marker_sampleIndex;
    PTB_marker_sampleIndex = PTB_marker_sampleIndex_compensate;    
    
    
    total_marker_count_raw = total_marker_count;
    innerTrial_marker_count_raw = innerTrial_marker_count;
    
    total_marker_count = length(PTB_marker_sampleIndex);
    innerTrial_marker_count = total_marker_count - trial_marker_count;   
    
    fprintf("Compensate now!\n");

    fprintf("total_marker_count = %d, ", total_marker_count);
    fprintf("trial_marker_count = %d, ", trial_marker_count);
    fprintf("innerTrial_marker_count = %d, ", innerTrial_marker_count);
    fprintf("end_marker_count = %d.\n", end_marker_count);
    
end

a = 1;
% If lose some markers at begining due to not open giveMarker program in time, then correct it here.
if trial_marker_count < trial_marker_count_ground    
    trial_marker_count_ground_raw = trial_marker_count_ground;
    firstTRIALID_markerIndex_raw = firstTRIALID_markerIndex;    
    redundant_trial_marker_num = trial_marker_count_ground_raw - trial_marker_count;
            
    temp_TRIALID_count = 0;
    for tempi=(firstTRIALID_markerIndex_raw+1):trial_para.marker.count
        k = strfind(trial_para.marker.content{tempi}, 'TRIALID ');% if find, k shouldn't be empty
        if k == 1
            temp_TRIALID_count = temp_TRIALID_count + 1;
            if temp_TRIALID_count == redundant_trial_marker_num
                firstTRIALID_markerIndex = tempi;
                break
            end
        end
    end
    trial_marker_count_ground = trial_marker_count;
    fprintf("(Corrected) PTB theoretically recorded total_marker_count (from %d) = %d, ",trial_para.marker.count,trial_para.marker.count-firstTRIALID_markerIndex+1);
    fprintf("trial_marker_count = %d.\n", trial_marker_count_ground);     
end
    
% if total_marker_count ~= (trial_para.marker.count-fromMarkerCount)
if total_marker_count ~= (trial_para.marker.count-firstTRIALID_markerIndex+1)
    fprintf('Warning! Lose some PTB marker!\n');
end
a = 1;


%% Double check for TRIALID sample index
temp_1_markerCount = 0;
temp_1_markerIndex = [];
% for tempi=(fromMarkerCount+1):trial_para.marker.count
for tempi=firstTRIALID_markerIndex:trial_para.marker.count
    k = strfind(trial_para.marker.content{tempi}, 'TRIALID ');% if find, k shouldn't be empty
    if k == 1
        %temp_1_markerIndex = [temp_1_markerIndex tempi-fromMarkerCount]; %#ok<*AGROW>
        temp_1_markerIndex = [temp_1_markerIndex tempi-firstTRIALID_markerIndex+1]; %#ok<*AGROW>        
        temp_1_markerCount = temp_1_markerCount + 1;
    end
end
PTB_marker_1_sampleIndex = PTB_marker_sampleIndex(temp_1_markerIndex);
if sum(ismember(PTB_marker_TRIALID_sampleIndex,PTB_marker_1_sampleIndex)) ~= trial_marker_count
    fprintf('Warning! Trial marker inconsistent!\n');
end
a = 1;

%% Get REWARD sample index
temp_SUBMIT_markerCount = 0;
temp_SUBMIT_markerIndex = [];
temp_SUBMIT_markerBoolIndex = false(1,length(trial_para.isCorrect));
temp_SUBMIT_groundTruthCount = 0;
% for tempi=(fromMarkerCount+1):trial_para.marker.count
for tempi=1:trial_para.marker.count
    k = strfind(trial_para.marker.content{tempi}, 'SUBMIT');% if find, k shouldn't be empty
    if k == 1
        temp_SUBMIT_groundTruthCount = temp_SUBMIT_groundTruthCount + 1;
        %if tempi >= (fromMarkerCount+1)
        if tempi >= firstTRIALID_markerIndex
            %temp_SUBMIT_markerIndex = [temp_SUBMIT_markerIndex tempi-fromMarkerCount]; %#ok<*AGROW>
            temp_SUBMIT_markerIndex = [temp_SUBMIT_markerIndex tempi-firstTRIALID_markerIndex+1]; %#ok<*AGROW>
            temp_SUBMIT_markerCount = temp_SUBMIT_markerCount + 1;
            temp_SUBMIT_markerBoolIndex(temp_SUBMIT_groundTruthCount) = true;
        end        
    end
end

temp_REWARD_markerIndex = temp_SUBMIT_markerIndex(trial_para.isCorrect(temp_SUBMIT_markerBoolIndex)==1);
PTB_marker_REWARD_sampleIndex = PTB_marker_sampleIndex(temp_REWARD_markerIndex);



% PTB_marker_REWARD_sampleIndex;
PTB_marker_correctTRIALID_sampleIndex = zeros(1,length(PTB_marker_REWARD_sampleIndex));
for tempi=1:length(PTB_marker_REWARD_sampleIndex)
    tempIndex = find(PTB_marker_TRIALID_sampleIndex<PTB_marker_REWARD_sampleIndex(tempi),1,'last');
    PTB_marker_correctTRIALID_sampleIndex(tempi) = PTB_marker_TRIALID_sampleIndex(tempIndex);
end

PTB_marker_SUBMIT_sampleIndex = PTB_marker_sampleIndex(temp_SUBMIT_markerIndex);
PTB_marker_validTRIALID_sampleIndex = zeros(1,length(PTB_marker_SUBMIT_sampleIndex));
for tempi=1:length(PTB_marker_SUBMIT_sampleIndex)
    tempIndex = find(PTB_marker_TRIALID_sampleIndex<PTB_marker_SUBMIT_sampleIndex(tempi),1,'last');
    PTB_marker_validTRIALID_sampleIndex(tempi) = PTB_marker_TRIALID_sampleIndex(tempIndex);
end

% PTB_marker_REWARD_sampleIndex;
% PTB_marker_SUBMIT_sampleIndex;
% PTB_marker_correctTRIALID_sampleIndex;
% PTB_marker_validTRIALID_sampleIndex;
firstValidTrialIndex = find(temp_SUBMIT_markerBoolIndex==true,1);
validTrial_marker_count = length(PTB_marker_SUBMIT_sampleIndex);
correctTrial_marker_count = length(PTB_marker_REWARD_sampleIndex);



output = struct;
output.basic_para = basic_para;
output.trial_para = trial_para;
output.file_length = file_length_PTB;
output.high_count = high_count;
output.total_marker_count = total_marker_count;
output.trial_marker_count = trial_marker_count;
output.innerTrial_marker_count = innerTrial_marker_count;
output.end_marker_count = end_marker_count;
output.firstValidTrialIndex = firstValidTrialIndex;
output.validTrial_marker_count = validTrial_marker_count;
output.correctTrial_marker_count = correctTrial_marker_count;
output.firstTRIALID_markerIndex = firstTRIALID_markerIndex;
output.marker_content = trial_para.marker.content(firstTRIALID_markerIndex:end);
output.marker_data = marker_data_PTB;
output.marker_boolIndex = PTB_marker_boolIndex;
output.marker_sampleIndex = PTB_marker_sampleIndex;
output.marker_TRIALID_sampleIndex = PTB_marker_TRIALID_sampleIndex;
output.marker_validTRIALID_sampleIndex = PTB_marker_validTRIALID_sampleIndex;
output.marker_correctTRIALID_sampleIndex = PTB_marker_correctTRIALID_sampleIndex;
output.marker_SUBMIT_sampleIndex = PTB_marker_SUBMIT_sampleIndex;
output.marker_REWARD_sampleIndex = PTB_marker_REWARD_sampleIndex;


a = 1;
%% End