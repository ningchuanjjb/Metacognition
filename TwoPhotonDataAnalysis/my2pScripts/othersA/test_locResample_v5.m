
if_resample = 1;

num_roi = size(F_dff,1);


%% resample
if if_resample == 1
    %% seqCount1
    seqCount1 = zeros(numFrames,2);
    for temploc=1:numFrames
        currentLabel = temp_svm_Y_valid(temploc,:);
        uniqueSeq = unique(currentLabel);
        for tempi=1:length(uniqueSeq)
            seqCount1(temploc,tempi) = sum(currentLabel==uniqueSeq(tempi),'all');
        end
    end
    
    %% resample_trialIndex
    resampleCount1 = min(seqCount1(:,1));
    resampleCount2 = min(seqCount1(:,2));
    resampleCount = min(resampleCount1,resampleCount2);
    
    resample_targetLoc_trialIndex = zeros(numFrames,resampleCount);
    resample_nonTargetLoc_trialIndex = zeros(numFrames,resampleCount);
    for tempi=1:numFrames
        currentLabel = temp_svm_Y_valid(tempi,:);
        
        temp_targetLoc_trialIndex = find(currentLabel==true);
        temp_targetLoc_index_resampled = sort(randperm(length(temp_targetLoc_trialIndex),resampleCount));
        temp_targetLoc_trialIndex_resampled = temp_targetLoc_trialIndex(temp_targetLoc_index_resampled);
        resample_targetLoc_trialIndex(tempi,:) = temp_targetLoc_trialIndex_resampled;
        
        temp_nonTargetLoc_trialIndex = find(currentLabel==false);
        temp_nonTargetLoc_index_resampled = sort(randperm(length(temp_nonTargetLoc_trialIndex),resampleCount));
        temp_nonTargetLoc_trialIndex_resampled = temp_nonTargetLoc_trialIndex(temp_nonTargetLoc_index_resampled);
        resample_nonTargetLoc_trialIndex(tempi,:) = temp_nonTargetLoc_trialIndex_resampled;
    end
    resample_trialIndex = [resample_targetLoc_trialIndex,resample_nonTargetLoc_trialIndex];
    
    %% F_dff_loc & temp_svm_Y_valid_loc
    F_dff_loc = cell(numFrames,1);
    temp_svm_Y_valid_loc = cell(numFrames,1);
    for tempi=1:numFrames
        F_dff_loc{tempi} = F_dff(:,resample_trialIndex(tempi,:));
        temp_svm_Y_valid_loc{tempi} = temp_svm_Y_valid(:,resample_trialIndex(tempi,:));
    end
        
    F_dff_loc;
    temp_svm_Y_valid_loc;
end


