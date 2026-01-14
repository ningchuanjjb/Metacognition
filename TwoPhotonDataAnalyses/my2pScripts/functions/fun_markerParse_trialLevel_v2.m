function markerParse_trialLevel = fun_markerParse_trialLevel_v2(extractOutput_PTB)
% clear
% close all
% 
% load extractOutput_PTB
% load extractOutput_2p

trial_para = extractOutput_PTB.trial_para;

firstValidTrialIndex = extractOutput_PTB.firstValidTrialIndex;
validTrial_marker_count = extractOutput_PTB.validTrial_marker_count; %#ok<*NASGU>

marker_content = extractOutput_PTB.marker_content;
marker_sampleIndex = extractOutput_PTB.marker_sampleIndex;
marker_SUBMIT_sampleIndex = extractOutput_PTB.marker_SUBMIT_sampleIndex;

TOTAL_downSampling = extractOutput_PTB.TOTAL_downSampling;
validTRIALID_downSampling = extractOutput_PTB.validTRIALID_downSampling;

firstValidTrialIndex_2p = find(validTRIALID_downSampling.valid_marker_X_trialIndex==true,1);
validTrial_marker_count_2p = sum(validTRIALID_downSampling.valid_marker_X_trialIndex);

markerParse_trialLevel = cell(1,validTrial_marker_count_2p);

fprintf('validTrial_marker_count_2p = %d.\n',validTrial_marker_count_2p);
a = 1;
% tempi = 1;
for tempi=1:validTrial_marker_count_2p
    
    
    trial_index = (firstValidTrialIndex_2p+firstValidTrialIndex-1) + tempi - 1;
    
    currentSequence = trial_para.currentSequence{trial_index};
    
    current_TRIALID_sampleIndex = validTRIALID_downSampling.valid_marker_X_sampleIndex(tempi);
    
    temp_index = find(marker_SUBMIT_sampleIndex>current_TRIALID_sampleIndex,1);
    current_SUBMIT_sampleIndex = marker_SUBMIT_sampleIndex(temp_index);
    
    temp1 = find(TOTAL_downSampling.valid_marker_X_sampleIndex==current_TRIALID_sampleIndex,1);
    temp2 = find(TOTAL_downSampling.valid_marker_X_sampleIndex==current_SUBMIT_sampleIndex,1);
    currentTrialTotalMarkers_rangeForFrame = temp1:temp2;
    currentTrialTotalMarkers_frameIndex = TOTAL_downSampling.valid_marker_X_frameIndex(currentTrialTotalMarkers_rangeForFrame);
    
    temp1 = find(marker_sampleIndex==current_TRIALID_sampleIndex,1);
    temp2 = find(marker_sampleIndex==current_SUBMIT_sampleIndex,1);
    currentTrialTotalMarkers_rangeForSample = temp1:temp2;
    currentTrialTotalMarkers_sampleIndex = marker_sampleIndex(currentTrialTotalMarkers_rangeForSample);
    
    currentTrialTotalMarkers = marker_content(currentTrialTotalMarkers_rangeForSample);
    
    temp_order = [];
    temp_location = [];
    for tempj=1:length(currentTrialTotalMarkers)
        k = strfind(currentTrialTotalMarkers{tempj}, 'TARGET ');% if find, k shouldn't be empty
        if k == 1
            %currentTrialTotalMarkers{tempj};
            k_order = 1;
            k_location = strfind(currentTrialTotalMarkers{tempj}, 'ITEM ');
            temp_order = [temp_order str2double(currentTrialTotalMarkers{tempj}(length('TARGET ')+k_order))]; %#ok<*AGROW>
            temp_location = [temp_location str2double(currentTrialTotalMarkers{tempj}(length('ITEM ')+k_location))];            
        end                
    end
    if sum(ismember(temp_location,currentSequence)) ~= length(currentSequence)
        fprintf('Warning! Sequence info is inconsistent between trial_para and marker_content.\n');
    end
    
    
    temp_struct = struct;
    
    temp_struct.trial_index = trial_index;
    temp_struct.currentSequence = currentSequence;
    temp_struct.current_TRIALID_sampleIndex = current_TRIALID_sampleIndex;
    temp_struct.current_SUBMIT_sampleIndex = current_SUBMIT_sampleIndex;
    temp_struct.currentTrialTotalMarkers_rangeForFrame = currentTrialTotalMarkers_rangeForFrame;
    temp_struct.currentTrialTotalMarkers_frameIndex = currentTrialTotalMarkers_frameIndex;
    temp_struct.currentTrialTotalMarkers_rangeForSample = currentTrialTotalMarkers_rangeForSample;
    temp_struct.currentTrialTotalMarkers_sampleIndex = currentTrialTotalMarkers_sampleIndex;
    temp_struct.currentTrialTotalMarkers = currentTrialTotalMarkers;
    
    markerParse_trialLevel{tempi} = temp_struct;
    
    a = 1;
    
    %if isempty(currentTrialTotalMarkers_rangeForFrame) == 1
    %    continue
    %end
    %home
    %fprintf('Current trial index is %d.\n',trial_index);
    %fprintf('CurrentSequence = %s.\n',num2str(currentSequence));
    %fprintf('Current trial marker sample count is %d.\n',length(currentTrialTotalMarkers_rangeForSample));
    %fprintf('Current trial marker frame count is %d.\n',length(currentTrialTotalMarkers_rangeForFrame));
    %fprintf('Current trial sample index range from %d to %d.\n',current_TRIALID_sampleIndex,current_SUBMIT_sampleIndex);
    %fprintf('Current trial frame index range from %d to %d.\n',currentTrialTotalMarkers_frameIndex(1),currentTrialTotalMarkers_frameIndex(end));
end

a = 1;
% markerParse_trialLevel;





