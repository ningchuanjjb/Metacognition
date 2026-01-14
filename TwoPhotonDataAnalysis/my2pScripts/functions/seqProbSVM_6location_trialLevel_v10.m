function [svm_posterior_seqLevel,validBoolIndex_seqLevel] = ...
    seqProbSVM_6location_trialLevel_v10(F_dff_raw,temp_svm_Y,boolIndex_location_seq,t_SVM,numFrames) %#ok<*INUSD>

%% Initialization

if_resample = 0;


F_dff = F_dff_raw; %#ok<*NASGU>

sum_F_dff = sum(F_dff,1);
validSeqBoolIndex = ~isnan(sum_F_dff);

temp_svm_Y_valid = temp_svm_Y(:,validSeqBoolIndex);
temp_svm_Y_valid_T = temp_svm_Y_valid';
F_dff = F_dff(:,validSeqBoolIndex);

F_dff_T = F_dff';

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
        
    F_dff_loc; %#ok<*VUNUS>
    temp_svm_Y_valid_loc;
end

a = 1;

%% temp_Mdl_CV_binary
KFold_num = 10;%10
temp_Mdl_CV_binary = cell(numFrames,1);
c_all = cell(1,numFrames);
for tempi=1:numFrames
    if if_resample == 0
        temp_svm_X = F_dff_T;        
        currentLabel = temp_svm_Y_valid(tempi,:);
    elseif if_resample == 1
        temp_svm_X = F_dff_loc{tempi}';
        currentLabel = temp_svm_Y_valid_loc{tempi}(tempi,:);        
    end
    
    c = cvpartition(currentLabel,'KFold',KFold_num);
    c_all{tempi} = c;
    
    temp_Mdl_binary = fitcecoc(temp_svm_X,currentLabel,'Learners',t_SVM,'FitPosterior',true); %#ok<*PFOUS>
    temp_Mdl_CV_binary{tempi} = crossval(temp_Mdl_binary,'CVPartition',c); % Very time-consuming!!!    
    a = 1;
end

%% multi_Posterior_cell
multi_Posterior_cell = cell(1, KFold_num);
for temploc=1:numFrames
    for tempk=1:KFold_num
        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{temploc}.ModelParameters.Generator.UseObsForIter(:,tempk);
        if if_resample == 0
            temp_F_dff_T = F_dff_T(~tempTrialBoolIndex_fold,:); %#ok<*PFBNS>
        elseif if_resample == 1
            temp_F_dff_T = F_dff_loc{temploc}(:,~tempTrialBoolIndex_fold)'; %#ok<*PFBNS>            
        end
        temp_F_dff_T_2d = temp_F_dff_T;
        [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},temp_F_dff_T_2d);
        tempPosterior_2 = tempPosterior(:,2);
        multi_Posterior_cell{tempk}(:,temploc) = tempPosterior_2;
    end
end

%% Posterior_2d
if if_resample == 0
    Posterior_2d = zeros(sum(validSeqBoolIndex),numFrames);
elseif if_resample == 1
    Posterior_2d = zeros(resampleCount*2,numFrames);
end
for tempk=1:KFold_num
    temp_Posterior = multi_Posterior_cell{tempk};
    for temploc=1:numFrames
        tempTrialBoolIndex_fold = temp_Mdl_CV_binary{temploc}.ModelParameters.Generator.UseObsForIter(:,tempk);
        Posterior_2d(~tempTrialBoolIndex_fold,temploc) = temp_Posterior(:,temploc);
    end
end

predictLabel_boolIndex = Posterior_2d > 0.5;%0.5

if if_resample == 0
    tempBoolIndex = predictLabel_boolIndex == temp_svm_Y_valid_T;
elseif if_resample == 1
    
end
predict_boolIndex = sum(tempBoolIndex,2) == numFrames;

%% svm_posterior_seqLevel
svm_posterior_seqLevel = zeros(1,size(boolIndex_location_seq,2));
validBoolIndex_seqLevel = true(1,size(boolIndex_location_seq,2));
a = 1;
for tempi=1:size(boolIndex_location_seq,2)
    tempSeq_bool = boolIndex_location_seq(:,tempi)';

    temp1_boolIndex = sum(temp_svm_Y_valid_T == tempSeq_bool,2) == numFrames;
        
    if sum(temp1_boolIndex) == 0
        validBoolIndex_seqLevel(tempi) = false;
        continue
    end
    
    temp_svm_posterior_seqLevel = sum(predict_boolIndex(temp1_boolIndex))/sum(temp1_boolIndex);    
    svm_posterior_seqLevel(tempi) = temp_svm_posterior_seqLevel;        
end
a = 1;



% prob_mean = mean(svm_posterior_seqLevel)
% prob_median = median(svm_posterior_seqLevel)
% a = 1;
% figure(2)
% plot(svm_posterior_seqLevel)


%% End
