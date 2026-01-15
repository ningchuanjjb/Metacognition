% Chuan's 7th script (20251214)
% This script: 
% (1) Compute single-neuron selectivity in one FOV.
% (2) Compute single-trial memory strength, related to figure2.
% (3) Conduct within-sequence analysis of memory strength, related to figure2.

% please run 'stepD3_train123_test123.m' first

%% Initialization
% close all

% if_dff_singleSession1_twoSession2 = 2;

% if if_monkey_D0_Z1 == 0
%     currentSessionIndex_A = 19;
% elseif if_monkey_D0_Z1 == 1
%     currentSessionIndex_A = 7;
% end

% currentSessionIndex = 20;%20-->10-->28-->5-->16-->21-->19

% locCorr_trialLevel
% memoryPrecision_trialLevel


if_loadDemoData = if_loadDemoData;

if_dff_singleSession1_twoSession2 = if_dff_singleSession1_twoSession2_allSession3;

if_compute_step2 = 1;%1

if_locTuningFit_lasso0_fitglm1 = 1;

if_tuning_location0_precision1 = 0;

if_memoryPrecision_accuracy0_sigma1 = if_memoryPrecision_accuracy0_sigma1; %#ok<*ASGSL>

if_entropy = if_entropy;

if_pProd_trialLevel_n11n = 1;


if_twoThreshold_median0_optimal1 = 0;%0,1

if if_twoThreshold_median0_optimal1 == 0
    if_memoryPrecisionThreshold_median0_optimal1 = 0;
elseif if_twoThreshold_median0_optimal1 == 1
    if_memoryPrecisionThreshold_median0_optimal1 = 1;
end


if_ePAIRS = 0;

if_trialEvidenceExampleSeq_histogram0_violinPlot1 = 1;
if_trialEvidenceAllSeq_violinplot0_pairline1 = 1;

if_memeoryPrecision_stimuli0_response1 = 0;%1, 0(entropy)
if_plot_memoryPrecision_trialLevelEvidence = 0;%0
if_plot_memoryPrecision_trialLevelEvidence_timeCourse = 0;%0
if_colormap_loadEnhanced = 1;

if exist('if_compute_summary','var') == 1
    if_memeoryPrecision_stimuli0_response1 = if_memeoryPrecision_stimuli0_response1_summary;
    if_plot_memoryPrecision_trialLevelEvidence = 1; % only temporary ban it !!!!
end

% if_plot_memoryPrecision_trialLevelEvidence = 0;

if_plot_additionalSmooth = 1;


if_memoryPrecisionTuning_decoder0_behavior1 = 1;%0

% if_memoryPrecision_accuracy0_sigma1 = 1;%1

if_exampleNeuron_trialTuning0_seqTuning1 = 1;%1

if_plotBeta_trialTuning0_seqTuning1_mix2 = 2;%2

if_beta_n11n = 0;%1

if_inferOffloadResponse = 0;

if_plot = 1;

if_plot_3d_beta = 0;
if_plot_2d_3sub_beta = 1;

color_choiceMemory = [1,133+20,113]/255;
color_choiceOffload = [166+40,97,26]/255;


%lowThreshold_memoryPrecision

fprintf('if_memeoryPrecision_stimuli0_response1=%d.\n',if_memeoryPrecision_stimuli0_response1);

%% Load trial data
trialIndex_bool_memoryCorrect;

% temp_load = load([output_path,'\','trial_para.mat']);
% trial_para = temp_load.trial_para;
% trial_para_choiceCondition_flag = trial_para.choiceCondition_flag;
trial_para_choiceCondition_flag;

temp_trial_para_ifSelectOffloading = trial_para_ifSelectOffloading==1;

choiceBoolIndex = trial_para_choiceCondition_flag == 2;
allMemoryBoolIndex = ~temp_trial_para_ifSelectOffloading;
choiceMemoryBoolIndex = choiceBoolIndex & allMemoryBoolIndex;
choiceOffloadBoolIndex = choiceBoolIndex & temp_trial_para_ifSelectOffloading;

allMemoryCorrectBoolIndex = trialIndex_bool_memoryCorrect;
allMemoryErrorBoolIndex = allMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;

choiceMemoryCorrectBoolIndex = choiceMemoryBoolIndex & trialIndex_bool_memoryCorrect;
choiceMemoryErrorBoolIndex = choiceMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;


% length error
allMemoryErrorBoolIndex;
allMemoryLengthErrorBoolIndex = false(size(allMemoryErrorBoolIndex));
allMemoryNonLengthErrorBoolIndex = false(size(allMemoryErrorBoolIndex));

seqIndex;
seqIndex_response;
for tempi=1:length(allMemoryErrorBoolIndex)
    if allMemoryErrorBoolIndex(tempi) == false
        continue
    end
    temp_seqLength = sum(boolIndex_location_seq_T(seqIndex(tempi),:));
    if ~isnan(seqIndex_response(tempi))
        temp_seqLength_response = sum(boolIndex_location_seq_T(seqIndex_response(tempi),:));
    else
        temp_seqLength_response = 0;
    end
    
    if temp_seqLength ~= temp_seqLength_response
        allMemoryLengthErrorBoolIndex(tempi) = true;
    else
        allMemoryNonLengthErrorBoolIndex(tempi) = true;
    end
end

% choiceMemoryBoolIndex_mildSeq = false(size(choiceBoolIndex));
% choiceOffloadBoolIndex_mildSeq = false(size(choiceBoolIndex));


offloadingProb_inOne;
seqIndex;

rProbLowThreshold_mildSeq = 0.40;%0.4-->0.45-->0.47
rProbHighThreshold_mildSeq = 1-rProbLowThreshold_mildSeq;

% exampleSeq = 13;%17
exampleSeq = 16;% Accuracy: 17,13,16
% exampleSeq = 10;% Offloading rate: 17,10,21

% mildSeqBoolIndex_seqLevel = (offloadingProb_inOne>=rProbLowThreshold_mildSeq) & (offloadingProb_inOne<=rProbHighThreshold_mildSeq);
% a1 = find(mildSeqBoolIndex_seqLevel==true);

mildSeqBoolIndex_seqLevel = false(1,length(offloadingProb_inOne));
mildSeqBoolIndex_seqLevel(exampleSeq) = true;

mildSeqBoolIndex = false(size(choiceBoolIndex));
for tempi=1:length(mildSeqBoolIndex)
    if mildSeqBoolIndex_seqLevel(seqIndex(tempi)) == true
        mildSeqBoolIndex(tempi) = true;
    end
end
mildSeqBoolIndex;
a1 = sum(mildSeqBoolIndex);

choiceMemoryBoolIndex_mildSeq = choiceMemoryBoolIndex & mildSeqBoolIndex;
choiceOffloadBoolIndex_mildSeq = choiceOffloadBoolIndex & mildSeqBoolIndex;
a2 = sum(choiceMemoryBoolIndex_mildSeq);
a3 = sum(choiceOffloadBoolIndex_mildSeq);

choiceMemoryCorrectBoolIndex_mildSeq = choiceMemoryCorrectBoolIndex & mildSeqBoolIndex;
choiceMemoryErrorBoolIndex_mildSeq = choiceMemoryErrorBoolIndex & mildSeqBoolIndex;

allMemoryCorrectBoolIndex_mildSeq = allMemoryCorrectBoolIndex & mildSeqBoolIndex;
allMemoryErrorBoolIndex_mildSeq = allMemoryErrorBoolIndex & mildSeqBoolIndex;


if if_loadDemoData == 0
    if if_dff_singleSession1_twoSession2 == 1
        FOVName_currentFOV = selevtivity_multiFOV.FOVName_multiFOV(currentSessionIndex,:);
        FOVName_currentFOV2 = FOVName_currentFOV(8:(8+4));
    elseif if_dff_singleSession1_twoSession2 == 2
        FOVName_currentFOV = currentSession;
        FOVName_currentFOV2 = [FOVName_currentFOV(5:9),'+',FOVName_currentFOV(19:23)];
    end
elseif if_loadDemoData == 1
    FOVName_currentFOV2 = 'demo (monkey D, FOV 8)';
end

% if if_dff_singleSession1_twoSession2 == 1
%     temp_range = FOVAllCellRange_multiFOV(currentSessionIndex,1):FOVAllCellRange_multiFOV(currentSessionIndex,2);
% elseif if_dff_singleSession1_twoSession2 == 2
%     temp_range_A = FOVAllCellRange_multiFOV(currentSessionIndex_A,1):FOVAllCellRange_multiFOV(currentSessionIndex_A,2);
%     temp_range_AB = decodingDataSimplified.extraForMerged.tempMappingCellIndex(:,1)' + temp_range_A(1) - 1;
%     temp_range = temp_range_AB;
% end
% cellIndex_suite2p = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range);
%
% if false
%     cellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1);
% end

if if_loadDemoData == 0
    if if_dff_singleSession1_twoSession2 == 1
        temp_range = FOVAllCellRange_multiFOV(currentSessionIndex,1):FOVAllCellRange_multiFOV(currentSessionIndex,2);
        cellIndex_suite2p = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range);
        
    elseif if_dff_singleSession1_twoSession2 == 2
        %temp_range_A = FOVAllCellRange_multiFOV(currentSessionIndex_A,1):FOVAllCellRange_multiFOV(currentSessionIndex_A,2);
        %temp_range_AB = decodingDataSimplified.extraForMerged.tempMappingCellIndex(:,1)' + temp_range_A(1) - 1;
        %temp_range = temp_range_AB;
        %cellIndex_suite2p = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range);
        
        cellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1);
        
        
        temp_range_B = FOVAllCellRange_multiFOV(currentSessionIndex_B,1):FOVAllCellRange_multiFOV(currentSessionIndex_B,2);
        cellIndex_suite2p_B_raw = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range_B);
        cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
        temp_range_AB = nan(1,length(cellIndex_suite2p_B));
        for tempi=1:length(cellIndex_suite2p_B)
            temp_range_AB(tempi) = find(cellIndex_suite2p_B_raw==cellIndex_suite2p_B(tempi));
        end
        temp_range_AB = temp_range_AB + temp_range_B(1) - 1;
        
        cellIndex_suite2p_B2 = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range_AB);
        
        r2_6loc = glm_r2_6locMean_delay1Bin_multiFOV(temp_range_AB);
        beta_6loc = glm_beta_6locMean_delay1Bin_multiFOV(temp_range_AB,:);
        std_beta_6loc = std(beta_6loc,1,2);
        
        
        boolIndex_6loc = tempBoolIndex123_or(temp_range_AB);
        p_6loc = nan(size(r2_6loc));
        p_6loc(boolIndex_6loc) = 0.001;
        p_6loc(~boolIndex_6loc) = 0.1;
        
        a = 1;
        
        
        temp1 = max(beta_6loc,[],2);
        temp2 = min(beta_6loc,[],2);
        temp3 = nan(length(beta_6loc),1);
        for tempi=1:length(beta_6loc)
            if abs(temp1(tempi)) > abs(temp2(tempi))
                temp3(tempi) = temp1(tempi);
            else
                temp3(tempi) = temp2(tempi);
            end
        end
        beta_6loc_peak = temp3;
        
    end
    
elseif if_loadDemoData == 1
    %cellIndex_suite2p = 1:size(F_dff_decisionBin_demo,1);    
    cellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1);
    cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
end


a = 1;

%% Step 1: Tuning for seq-level offloading rate

% selectiveCellBoolIndex_rProb_glm_delayx_currentFOV = selectiveCellBoolIndex_rProb_glm_delayx_multiFOV(temp_range);
% r2_rProb_glm_delayx_currentFOV = r2_rProb_glm_delayx_multiFOV(temp_range);
% p_rProb_glm_delayx_currentFOV = p_rProb_glm_delayx_multiFOV(temp_range);

% [M,I] = sort(r2_rProb_glm_delayx_currentFOV,'descend');
% I_suite2p = cellIndex_suite2p(I);
% r2_rProb_glm_delayx_currentFOV;
% r2_rProb = r2_rProb_glm_delayx_currentFOV;
% p_rProb = p_rProb_glm_delayx_currentFOV;

if if_loadDemoData == 0
    F_dff_decisionBin1 = mean(F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2)),3);
elseif if_loadDemoData == 1
    F_dff_decisionBin1 = F_dff_decisionBin_demo;
end

F_dff_decisionBin1 = F_dff_decisionBin1(neuronBoolIndex,:);


F_dff_decisionBin1_seqMerged = nan(size(F_dff_decisionBin1,1),sum(numSeq));
for tempi=1:sum(numSeq)
    temp_dff = F_dff_decisionBin1(:,seqIndex==tempi);
    F_dff_decisionBin1_seqMerged(:,tempi) = mean(temp_dff,2);
end

F_dff_decisionBin1_zScore = nan(size(F_dff_decisionBin1));
for tempi=1:size(F_dff_decisionBin1_zScore,1)
    temp1 = F_dff_decisionBin1(tempi,:);
    temp2 =(temp1-mean(temp1))/std(temp1); % z-score
    F_dff_decisionBin1_zScore(tempi,:) = temp2;
end
F_dff_decisionBin1_zScore_seqMerged = nan(size(F_dff_decisionBin1_zScore,1),sum(numSeq));
for tempi=1:sum(numSeq)
    temp_dff = F_dff_decisionBin1_zScore(:,seqIndex==tempi);
    F_dff_decisionBin1_zScore_seqMerged(:,tempi) = mean(temp_dff,2);
end



if if_locTuningFit_lasso0_fitglm1 == 1
    
    %% r2_6loc_allMemoryCorrect
    x_raw = boolIndex_location_seq_T(seqIndex,:);
    x_raw2 = x_raw(allMemoryCorrectBoolIndex,:);
    x = (x_raw2-mean(x_raw2,1))./std(x_raw2,0,1);
    
    y_raw = F_dff_decisionBin1_zScore';
    y_raw2 = y_raw(allMemoryCorrectBoolIndex,:);
    y = (y_raw2-mean(y_raw2,1))./std(y_raw2,0,1);
    
    r2 = zeros(roiNum,1);
    p = zeros(roiNum,1);
    beta = zeros(roiNum,numFrames);
    parfor tempi=1:roiNum
        %for tempi=1:roiNum
        warning('off');
        temp_mdl = fitglm(x,y(:,tempi));
        r2(tempi) = temp_mdl.Rsquared.Adjusted;
        p(tempi) = coefTest(temp_mdl);
        beta(tempi,:) = temp_mdl.Coefficients.Estimate(2:end);
        warning('on');
    end
    r2_6loc_allMemoryCorrect = r2;
    beta_6loc_allMemoryCorrect = beta;
    p_6loc_allMemoryCorrect = p;
    %boolIndex_6loc_allMemoryCorrect = p_6loc < 0.01;
    boolIndex_6loc_allMemoryCorrect = p_6loc_allMemoryCorrect < 0.01;
    
    
    %% r2_6loc_noChoiceCorrect
    x_raw = boolIndex_location_seq_T(seqIndex,:);
    x_raw2 = x_raw(allMemoryCorrectBoolIndex&(~choiceBoolIndex),:);
    x = (x_raw2-mean(x_raw2,1))./std(x_raw2,0,1);
    
    y_raw = F_dff_decisionBin1_zScore';
    y_raw2 = y_raw(allMemoryCorrectBoolIndex&(~choiceBoolIndex),:);
    y = (y_raw2-mean(y_raw2,1))./std(y_raw2,0,1);
    
    r2 = zeros(roiNum,1);
    p = zeros(roiNum,1);
    beta = zeros(roiNum,numFrames);
    parfor tempi=1:roiNum
        %for tempi=1:roiNum
        warning('off');
        temp_mdl = fitglm(x,y(:,tempi));
        r2(tempi) = temp_mdl.Rsquared.Adjusted;
        p(tempi) = coefTest(temp_mdl);
        beta(tempi,:) = temp_mdl.Coefficients.Estimate(2:end);
        warning('on');
    end
    r2_6loc_noChoiceCorrect = r2;
    beta_6loc_noChoiceCorrect = beta;
    p_6loc_noChoiceCorrect = p;
    boolIndex_6loc_noChoiceCorrect = p_6loc_noChoiceCorrect < 0.01;
    
    
    %%
    r2_6loc = r2_6loc_noChoiceCorrect;
    beta_6loc = beta_6loc_noChoiceCorrect;
    p_6loc = p_6loc_noChoiceCorrect;
    boolIndex_6loc = boolIndex_6loc_noChoiceCorrect;
    
    
    
    std_beta_6loc = std(beta_6loc,1,2);
    
    temp1 = max(beta_6loc,[],2);
    temp2 = min(beta_6loc,[],2);
    temp3 = nan(length(beta_6loc),1);
    for tempi=1:length(beta_6loc)
        if abs(temp1(tempi)) > abs(temp2(tempi))
            temp3(tempi) = temp1(tempi);
        else
            temp3(tempi) = temp2(tempi);
        end
    end
    beta_6loc_peak = temp3;
end



if false
    if false
        
        temp1 = selectiveCellBoolIndex_rProb_glm_delay1_multiFOV(temp_range_AB);
        temptempBoolIndex_rProb = temp1;
        
        dff_seqMean_delay1 = dff_seqMean_delay1_multiFOV(temp_range_AB,:);
        
        %dff_seqMean_delay1 = nan(size(F_dff_decisionBin1,1),sum(numSeq));
        %for tempi=1:sum(numSeq)
        %    temp_dff = F_dff_decisionBin1(:,seqIndex==tempi & trialIndex_bool_memoryCorrect);
        %    dff_seqMean_delay1(:,tempi) = mean(temp_dff,2);
        %end
        
        %% Seq-level glm for location tuning
        x = boolIndex_location_seq_T(1:41,:);
        y = dff_seqMean_delay1(:,1:41);
        
        temp_glm_beta = zeros(size(y,1),numFrames);
        temp_glm_r2 = zeros(size(y,1),1);
        parfor tempi=1:size(y,1)
            warning('off');
            temp_x = x;
            %temp_y = y(tempi,:)';
            
            temp_y0 = y(tempi,:)';
            temp_y = (temp_y0-mean(temp_y0))/std(temp_y0); % z-score
            
            temp_mdl= fitglm(temp_x,temp_y,'linear','Distribution','normal','Intercept',false);
            temp_glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
            temp_glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
            warning('on');
        end
        
        x;
        y = dff_seqMean_delay1(:,1:41);
        
        roiNum = size(y,1);
        
        temp_shuffleNum = 16*7;
        temp_glm_r2_shuffled = zeros(size(y,1),temp_shuffleNum);
        parfor tempShuffleIndex=1:temp_shuffleNum
            %for tempShuffleIndex=1:10
            temp1 = 1:41;
            temp1_shuffled = temp1(randperm(length(temp1)));
            x_shuffled = x(temp1_shuffled,:);
            
            for tempi=1:roiNum
                temp_x = x_shuffled;
                %temp_y = y(tempi,:)';
                
                temp_y0 = y(tempi,:)';
                temp_y = (temp_y0-mean(temp_y0))/std(temp_y0); % z-score
                
                temp_mdl= fitglm(temp_x,temp_y,'linear','Distribution','normal','Intercept',false);
                temp_glm_r2_shuffled(tempi,tempShuffleIndex) = temp_mdl.Rsquared.Adjusted;
            end
        end
        
        temp1 = mean(temp_glm_r2_shuffled,2);
        temp2 = prctile(temp_glm_r2_shuffled,99.9,2);
        
        sum(temp_glm_r2>temp2)
        
        temptempBoolIndex_6loc = temp_glm_r2>temp2;
        
    end
    
    
    if false
        %% Seq-level anova for location tuning
        x = [1:6,1:6,1:6]';
        y = nan(size(dff_seqMean_delay1,1),length(x));
        
        temp1 = boolIndex_location_seq_T(1:41,:);
        for tempj=1:3
            for tempi=1:6
                temp_range = sum(numSeq(1:tempj-1))+1:sum(numSeq(1:tempj));
                
                temp0 = false(41,1);
                temp0(temp_range) = true;
                
                temp1_locx = temp1(:,tempi)==true & temp0;
                temp2_locx = dff_seqMean_delay1(:,temp1_locx);
                temp2_locx_mean = mean(temp2_locx,2);
                
                tempij = (tempj-1)*6 + tempi;
                y(:,tempij) = temp2_locx_mean;
            end
        end
        
        temp_p = zeros(size(y,1),1);
        parfor tempi=1:size(y,1)
            %for tempi=1:1
            temp_y = y(tempi,:)';
            [temp_p(tempi),~,~] = anova1(temp_y, x,'off');
        end
        p_anova_6loc = temp_p;
        
    end
    
    
    if false
        %% Trial-level glm for location tuning
        temptempBoolIndex = trialIndex_bool_memoryCorrect & seqIndex<=sum(numSeq(1:3));
        x = boolIndex_location_allTrial(:,temptempBoolIndex)';
        y = F_dff_decisionBin1(:,temptempBoolIndex);
        
        temp_glm_beta = zeros(size(y,1),numFrames);
        temp_glm_p = zeros(size(y,1),numFrames);
        temp_glm_r2 = zeros(size(y,1),1);
        parfor tempi=1:size(y,1)
            temp_x = x;
            
            %temp_y = y(tempi,:)';
            
            temp_y0 = y(tempi,:)';
            temp_y = (temp_y0-mean(temp_y0))/std(temp_y0); % z-score
            
            temp_mdl= fitglm(temp_x,temp_y,'linear','Distribution','normal','Intercept',false);
            temp_glm_beta(tempi,:) = temp_mdl.Coefficients.Estimate;
            temp_glm_p(tempi,:) = temp_mdl.Coefficients.pValue;
            temp_glm_r2(tempi) = temp_mdl.Rsquared.Adjusted;
        end
        
        temp1 = sum(temp_glm_p<0.01,2);
        temp2 = temp1>=4;
        sum(temp2)
        temptempBoolIndex_6loc = temp2;
        
        x = boolIndex_location_allTrial(:,temptempBoolIndex)';
        y = F_dff_decisionBin1(:,temptempBoolIndex);
        
        roiNum = size(y,1);
        
        temp_shuffleNum = 16*7;
        temp_glm_r2_shuffled = zeros(size(y,1),temp_shuffleNum);
        parfor tempShuffleIndex=1:temp_shuffleNum
            %for tempShuffleIndex=1:temp_shuffleNum
            temp1 = 1:size(x,1);
            temp1_shuffled = temp1(randperm(length(temp1)));
            x_shuffled = x(temp1_shuffled,:);
            
            for tempi=1:roiNum
                temp_x = x_shuffled;
                
                %temp_y = y(tempi,:)';
                
                temp_y0 = y(tempi,:)';
                temp_y = (temp_y0-mean(temp_y0))/std(temp_y0); % z-score
                
                temp_mdl= fitglm(temp_x,temp_y,'linear','Distribution','normal','Intercept',false);
                temp_glm_r2_shuffled(tempi,tempShuffleIndex) = temp_mdl.Rsquared.Adjusted;
            end
        end
        
        temp1 = mean(temp_glm_r2_shuffled,2);
        temp2 = prctile(temp_glm_r2_shuffled,99.9,2);
        
        sum(temp_glm_r2>temp2)
        
        temptempBoolIndex_6loc = temp_glm_r2>temp2;
        
    end
    
    if false
        temptempBoolIndex_A = trialIndex_bool_memoryCorrect & seqIndex<=sum(numSeq(1:3));
        temptempBoolIndex_B = seqIndex<=6;
        temptempBoolIndex = temptempBoolIndex_A & temptempBoolIndex_B;
        x = seqIndex(temptempBoolIndex);
        y = F_dff_decisionBin1(:,temptempBoolIndex);
        
        p_anova = zeros(size(y,1),1);
        parfor tempi=1:size(y,1)
            %for tempi=1:1
            temp_x = x';
            temp_y = y(tempi,:)';
            
            [p_anova(tempi),~,~] = anova1(temp_y, temp_x,'off');
        end
        sum(p_anova<0.01)
        
    end
    
end


% figglm
r2_rProb = nan(size(F_dff_decisionBin1,1),1);
beta_rProb = nan(size(F_dff_decisionBin1,1),1);
p_rProb = nan(size(F_dff_decisionBin1,1),1);
parfor tempi=1:size(F_dff_decisionBin1,1)
    %x = F_dff_decisionBin1_seqMerged(tempi,:)';
    %y = offloadingProb_inOne';
    %temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
    
    x = offloadingProb_inOne';
    %y = F_dff_decisionBin1_seqMerged(tempi,:)';
    %y0 = F_dff_decisionBin1_seqMerged(tempi,:)';
    %y = (y0-mean(y0))/std(y0); % z-score
    y = F_dff_decisionBin1_zScore_seqMerged(tempi,:)';
    temp_mdl = fitglm(x,y);
    
    r2_rProb(tempi) = temp_mdl.Rsquared.Adjusted;
    beta_rProb(tempi) = temp_mdl.Coefficients.Estimate(2);
    p_rProb(tempi) = temp_mdl.Coefficients.pValue(2);
end
beta_rProb_raw = beta_rProb;
if if_beta_n11n == 1
    beta_rProb = rescale(beta_rProb_raw,-1,1);
end

r2_gProb = nan(size(F_dff_decisionBin1,1),1);
beta_gProb = nan(size(F_dff_decisionBin1,1),1);
p_gProb = nan(size(F_dff_decisionBin1,1),1);
parfor tempi=1:size(F_dff_decisionBin1,1)
    %x = F_dff_decisionBin1_seqMerged(tempi,:)';
    %y = (1-offloadingProb_inOne)';
    %temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
    
    x = (1-offloadingProb_inOne)';
    y = F_dff_decisionBin1_zScore_seqMerged(tempi,:)';
    temp_mdl = fitglm(x,y);
    
    r2_gProb(tempi) = temp_mdl.Rsquared.Adjusted;
    beta_gProb(tempi) = temp_mdl.Coefficients.Estimate(2);
    p_gProb(tempi) = temp_mdl.Coefficients.pValue(2);
end
beta_gProb_raw = beta_gProb;
if if_beta_n11n == 1
    beta_gProb = rescale(beta_gProb_raw,-1,1);
end

r2_seqPrecision = nan(size(F_dff_decisionBin1,1),1);
beta_seqPrecision = nan(size(F_dff_decisionBin1,1),1);
p_seqPrecision = nan(size(F_dff_decisionBin1,1),1);
parfor tempi=1:size(F_dff_decisionBin1,1)
    %y = seqPrecision_behavior;
    %y = seqPrecision_behavior_56;
    %x = F_dff_decisionBin1_seqMerged(tempi,1:length(y))';
    %temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
    
    x = seqPrecision_behavior_56;
    y = F_dff_decisionBin1_zScore_seqMerged(tempi,:)';
    temp_mdl = fitglm(x,y);
    
    r2_seqPrecision(tempi) = temp_mdl.Rsquared.Adjusted;
    beta_seqPrecision(tempi) = temp_mdl.Coefficients.Estimate(2);
    p_seqPrecision(tempi) = temp_mdl.Coefficients.pValue(2);
end
beta_seqPrecision_raw = beta_seqPrecision;
if if_beta_n11n == 1
    beta_seqPrecision = rescale(beta_seqPrecision_raw,-1,1);
end

[M,I] = max(r2_seqPrecision);


temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
% temptemp_range = (baselinePeriod_interval(1)):baselinePeriod_interval(3);
if if_loadDemoData == 0
    F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);
elseif if_loadDemoData == 1
    F_dff_baselineBin = F_dff_baselineBin_demo;
end
F_dff_baselineBin_seqMerged = nan(size(F_dff_baselineBin,1),sum(numSeq));
for tempi=1:sum(numSeq)
    temp_dff = F_dff_baselineBin(:,seqIndex==tempi);
    F_dff_baselineBin_seqMerged(:,tempi) = mean(temp_dff,2);
end

F_dff_baselineBin_zScore = nan(size(F_dff_baselineBin));
for tempi=1:size(F_dff_baselineBin_zScore,1)
    temp1 = F_dff_baselineBin(tempi,:);
    temp2 =(temp1-mean(temp1))/std(temp1); % z-score
    F_dff_baselineBin_zScore(tempi,:) = temp2;
end
F_dff_baselineBin_zScore_seqMerged = nan(size(F_dff_baselineBin_zScore,1),sum(numSeq));
for tempi=1:sum(numSeq)
    temp_dff = F_dff_baselineBin_zScore(:,seqIndex==tempi);
    F_dff_baselineBin_zScore_seqMerged(:,tempi) = mean(temp_dff,2);
end

r2_rProb_baseline = nan(size(F_dff_baselineBin,1),1);
beta_rProb_baseline = nan(size(F_dff_decisionBin1,1),1);
p_rProb_baseline = nan(size(F_dff_baselineBin,1),1);
parfor tempi=1:size(F_dff_baselineBin,1)
    %x = F_dff_baselineBin_seqMerged(tempi,:)';
    %y = offloadingProb_inOne';
    %temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
    
    x = offloadingProb_inOne';
    y = F_dff_baselineBin_zScore_seqMerged(tempi,:)';
    temp_mdl = fitglm(x,y);
    
    r2_rProb_baseline(tempi) = temp_mdl.Rsquared.Adjusted;
    beta_rProb_baseline(tempi) = temp_mdl.Rsquared.Adjusted;
    p_rProb_baseline(tempi) = temp_mdl.Coefficients.pValue(2);
end
beta_rProb_baseline_raw = beta_rProb_baseline;
if if_beta_n11n == 1
    beta_rProb_baseline = rescale(beta_rProb_baseline_raw,-1,1);
end

r2_gProb_baseline = nan(size(F_dff_baselineBin,1),1);
beta_gProb_baseline = nan(size(F_dff_decisionBin1,1),1);
p_gProb_baseline = nan(size(F_dff_baselineBin,1),1);
parfor tempi=1:size(F_dff_baselineBin,1)
    %x = F_dff_baselineBin_seqMerged(tempi,:)';
    %y = (1-offloadingProb_inOne)';
    %temp_mdl = fitglm(x(~isnan(x)),y(~isnan(x)));
    
    x = (1-offloadingProb_inOne)';
    y = F_dff_baselineBin_zScore_seqMerged(tempi,:)';
    temp_mdl = fitglm(x,y);
    
    r2_gProb_baseline(tempi) = temp_mdl.Rsquared.Adjusted;
    beta_gProb_baseline(tempi) = temp_mdl.Rsquared.Adjusted;
    p_gProb_baseline(tempi) = temp_mdl.Coefficients.pValue(2);
end
beta_gProb_baseline_raw = beta_gProb_baseline;
if if_beta_n11n == 1
    beta_gProb_baseline = rescale(beta_gProb_baseline_raw,-1,1);
end



threshold_p = 0.01;%0.05,0.01



temp_r2_rProb = r2_rProb;
temp_r2_rProb(p_rProb>=threshold_p) = -1;
threshold_r2_rProb = min(temp_r2_rProb(temp_r2_rProb>-1));

temp_r2_seqPrecision = r2_seqPrecision;
temp_r2_seqPrecision(p_seqPrecision>=threshold_p) = -1;
threshold_r2_seqPrecision = min(temp_r2_seqPrecision(temp_r2_seqPrecision>-1));

threshold_r2_rProb_baseline = (threshold_r2_rProb+threshold_r2_seqPrecision)/2;


if false
    %     fun_sigmaBased_singleTrialPrecision_Name_v = autoGetFunName_myScripts('fun_sigmaBased_singleTrialPrecision', [targetPATH '\functions']);
    %     fun_sigmaBased_singleTrialPrecision = str2func(fun_sigmaBased_singleTrialPrecision_Name_v);
    %
    %     temp_options = struct;
    %     temp_options.tempSeqIndex = tempSeqIndex;
    %     temp_options.boolIndex_location_seq_T = boolIndex_location_seq_T;
    %     temp_options.Posterior_2d_n11n_mean = Posterior_2d_n11n_mean;
    %     temp_options.numSeq = numSeq;
    %     temp_options.valid_length = valid_length;
    %     temp_options.score_stimuli_to_response = score_stimuli_to_response;
    %     temp_options.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
    %
    %     temp_precision = fun_sigmaBased_singleTrialPrecision(temp_options);
end

%% Real-time r2
if false
    % if if_compute_step2 == 1
    if_delay1_forward0_backward1 = 0; % new
    
    %% Preparation for F_dff_baseline
    F_dff_baseline = F_dff_baselinePeriod;
    F_dff_baseline_seqMerged = nan(size(F_dff_baseline,1),sum(numSeq),size(F_dff_baseline,3));
    for tempi=1:sum(numSeq)
        temp_dff = F_dff_baseline(:,seqIndex==tempi,:);
        F_dff_baseline_seqMerged(:,tempi,:) = mean(temp_dff,2);
    end
    
    
    %% Preparation for F_dff_lastT_sample
    for tempi=1:3
        temp_range = sum(numSeq(1:tempi-1))+1:sum(numSeq(1:tempi));
        
        temptempBoolIndex = ismember(seqIndex,temp_range)';
        if tempi==1
            trialBoolIndex_length1 = temptempBoolIndex;
        elseif tempi==2
            trialBoolIndex_length2 = temptempBoolIndex;
        elseif tempi==3
            trialBoolIndex_length3 = temptempBoolIndex;
        end
    end
    
    lengthx_sample_interval = length1_sample_interval;
    temp_trialBoolIndex_lengthx = trialBoolIndex_length1 | trialBoolIndex_length2 | trialBoolIndex_length3;
    
    
    F_dff_lastT_sample = nan(size(F_dff_length1_sample));
    for tempi=1:size(F_dff_lastT_sample,2)
        if trialBoolIndex_length1(tempi) == true
            temp1 = squeeze(F_dff_length1_sample(:,tempi,:));
        end
        if trialBoolIndex_length2(tempi) == true
            temp1 = squeeze(F_dff_length2_sample(:,tempi,:));
        end
        if trialBoolIndex_length3(tempi) == true
            temp1 = squeeze(F_dff_length3_sample(:,tempi,:));
        end
        if isnan(temp1(1,1)) == true
            continue
        end
        temp1 = temp1(:,end-lengthx_sample_interval(end)+1:end);
        F_dff_lastT_sample(:,tempi,:) = temp1;
    end
    
    
    F_dff_lastT_sample_seqMerged = nan(size(F_dff_lastT_sample,1),sum(numSeq),size(F_dff_lastT_sample,3));
    for tempi=1:sum(numSeq)
        temp_dff = F_dff_lastT_sample(:,seqIndex==tempi,:);
        F_dff_lastT_sample_seqMerged(:,tempi,:) = mean(temp_dff,2);
    end
    
    
    %% Preparation for F_dff_decision
    if if_delay1_forward0_backward1 == 1
        F_dff_decision = F_dff_decisionPeriodA(:,:,decisionPeriodA_interval(1):decisionPeriodA_interval(2));
    elseif if_delay1_forward0_backward1 == 0
        F_dff_decision = F_dff_decisionPeriodB(:,:,decisionPeriodB_interval(1):decisionPeriodB_interval(2));
    end
    F_dff_decision_seqMerged = nan(size(F_dff_decision,1),sum(numSeq),size(F_dff_decision,3));
    for tempi=1:sum(numSeq)
        temp_dff = F_dff_decision(:,seqIndex==tempi,:);
        F_dff_decision_seqMerged(:,tempi,:) = mean(temp_dff,2);
    end
    
    %% Preparation end
    F_dff_baseline;
    F_dff_baseline_seqMerged;
    F_dff_lastT_sample;
    F_dff_lastT_sample_seqMerged;
    F_dff_decision;
    F_dff_decision_seqMerged;
    
    
    %% Real-time r2 of meta-memory
    r2_rProb_timeCourse_baseline = nan(size(F_dff_baseline,1),size(F_dff_baseline,3));
    temp_timeSize = size(F_dff_baseline,3);
    parfor tempi=1:size(F_dff_baseline,1)
        for tempj=1:temp_timeSize
            x = offloadingProb_inOne';
            y = F_dff_baseline_seqMerged(tempi,:,tempj)';
            temp_mdl = fitglm(x,y);
            r2_rProb_timeCourse_baseline(tempi,tempj) = temp_mdl.Rsquared.Adjusted;
        end
    end
    %r2_rProb_timeCourse_baseline_mean = mean(r2_rProb_timeCourse_baseline(p_rProb<0.01,:),1);
    
    
    r2_rProb_timeCourse_lastT = nan(size(F_dff_lastT_sample,1),size(F_dff_lastT_sample,3));
    temp_timeSize = size(F_dff_lastT_sample,3);
    parfor tempi=1:size(F_dff_lastT_sample,1)
        for tempj=1:temp_timeSize
            x = offloadingProb_inOne';
            y = F_dff_lastT_sample_seqMerged(tempi,:,tempj)';
            temp_mdl = fitglm(x,y);
            r2_rProb_timeCourse_lastT(tempi,tempj) = temp_mdl.Rsquared.Adjusted;
        end
    end
    %r2_rProb_timeCourse_lastT_mean = mean(r2_rProb_timeCourse_lastT(p_rProb<0.01,:),1);
    
    
    r2_rProb_timeCourse_delay1 = nan(size(F_dff_decision,1),size(F_dff_decision,3));
    temp_timeSize = size(F_dff_decision,3);
    parfor tempi=1:size(F_dff_decision,1)
        for tempj=1:temp_timeSize
            x = offloadingProb_inOne';
            y = F_dff_decision_seqMerged(tempi,:,tempj)';
            temp_mdl = fitglm(x,y);
            r2_rProb_timeCourse_delay1(tempi,tempj) = temp_mdl.Rsquared.Adjusted;
        end
    end
    %r2_rProb_timeCourse_delay1_mean = mean(r2_rProb_timeCourse_delay1(p_rProb<0.01,:),1);
    
    
    
    
    %% Real-time r2 of location, trial-level
    temptempBoolIndex = trialIndex_bool_memoryCorrect & seqIndex<=sum(numSeq(1:3));
    
    
    temp_timeSize = size(F_dff_baseline,3);
    x = boolIndex_location_allTrial(:,temptempBoolIndex)';
    y = F_dff_baseline(:,temptempBoolIndex,:);
    r2_6loc_baseline_timeCourse = nan(size(F_dff_baseline,1),size(F_dff_baseline,3));
    parfor tempi=1:size(F_dff_baseline,1)
        for tempj=1:temp_timeSize
            temp_x = x;
            temp_y = y(tempi,:,tempj)'; %#ok<*PFTUSW>
            
            temp_mdl= fitglm(temp_x,temp_y,'linear','Distribution','normal','Intercept',false);
            r2_6loc_baseline_timeCourse(tempi,tempj) = temp_mdl.Rsquared.Adjusted;
        end
    end
    %r2_6loc_baseline_timeCourse_mean = mean(r2_6loc_baseline_timeCourse(r2_6loc>0.05,:),1);%0.05
    
    
    temp_timeSize = size(F_dff_lastT_sample,3);
    x = boolIndex_location_allTrial(:,temptempBoolIndex)';
    y = F_dff_lastT_sample(:,temptempBoolIndex,:);
    r2_6loc_lastT_timeCourse = nan(size(F_dff_lastT_sample,1),size(F_dff_lastT_sample,3));
    parfor tempi=1:size(F_dff_lastT_sample,1)
        for tempj=1:temp_timeSize
            temp_x = x;
            temp_y = y(tempi,:,tempj)'; %#ok<*PFTUSW>
            
            temp_mdl= fitglm(temp_x,temp_y,'linear','Distribution','normal','Intercept',false);
            r2_6loc_lastT_timeCourse(tempi,tempj) = temp_mdl.Rsquared.Adjusted;
        end
    end
    %r2_6loc_lastT_timeCourse_mean = mean(r2_6loc_lastT_timeCourse(r2_6loc>0.05,:),1);%0.05
    
    
    temp_timeSize = size(F_dff_decision,3);
    x = boolIndex_location_allTrial(:,temptempBoolIndex)';
    y = F_dff_decision(:,temptempBoolIndex,:);
    r2_6loc_delay1_timeCourse = nan(size(F_dff_decision,1),size(F_dff_decision,3));
    parfor tempi=1:size(F_dff_decision,1)
        for tempj=1:temp_timeSize
            temp_x = x;
            temp_y = y(tempi,:,tempj)'; %#ok<*PFTUSW>
            
            temp_mdl= fitglm(temp_x,temp_y,'linear','Distribution','normal','Intercept',false);
            r2_6loc_delay1_timeCourse(tempi,tempj) = temp_mdl.Rsquared.Adjusted;
        end
    end
    %r2_6loc_delay1_timeCourse_mean = mean(r2_6loc_delay1_timeCourse(r2_6loc>0.05,:),1);%0.05
    
    
    %%
    r2_rProb_timeCourse_baseline;
    r2_rProb_timeCourse_lastT;
    r2_rProb_timeCourse_delay1;
    
    r2_6loc_baseline_timeCourse;
    r2_6loc_lastT_timeCourse;
    r2_6loc_delay1_timeCourse;
    
    temp1 = 0.01;%0.01
    temp1_num = sum(p_rProb<temp1);
    r2_rProb_timeCourse_baseline_mean = mean(r2_rProb_timeCourse_baseline(p_rProb<temp1,:),1);
    r2_rProb_timeCourse_lastT_mean = mean(r2_rProb_timeCourse_lastT(p_rProb<temp1,:),1);
    r2_rProb_timeCourse_delay1_mean = mean(r2_rProb_timeCourse_delay1(p_rProb<temp1,:),1);
    
    temp_p1 = nan(size(r2_rProb_timeCourse_lastT_mean));
    for tempi=1:length(temp_p1)
        [~,temp_p1(tempi)] = ttest(r2_rProb_timeCourse_baseline_mean,r2_rProb_timeCourse_lastT_mean(tempi),'tail','left');
    end
    
    temp2 = 0.15;%0.05
    temp2_num = sum(r2_6loc>temp2);
    r2_6loc_baseline_timeCourse_mean = mean(r2_6loc_baseline_timeCourse(r2_6loc>temp2,:),1);
    r2_6loc_lastT_timeCourse_mean = mean(r2_6loc_lastT_timeCourse(r2_6loc>temp2,:),1);
    r2_6loc_delay1_timeCourse_mean = mean(r2_6loc_delay1_timeCourse(r2_6loc>temp2,:),1);
    
    temp_p2 = nan(size(r2_6loc_lastT_timeCourse_mean));
    for tempi=1:length(temp_p2)
        [~,temp_p2(tempi)] = ttest(r2_6loc_baseline_timeCourse_mean,r2_6loc_lastT_timeCourse_mean(tempi),'tail','left');
    end
    
end


%% Step 2: Tuning for trial-level memory precision (from decoder)
t_test = tic;
if if_compute_step2 == 1
    
    memoryPrecision_trialLevel_resampleIter = cell(length(trialIndex_bool_memoryCorrect),1);
    memoryPrecisionSigma_trialLevel_resampleIter = cell(length(trialIndex_bool_memoryCorrect),1);
    locCorr_trialLevel_resampleIter = cell(length(trialIndex_bool_memoryCorrect),1);
    seqDistri_trialLevel_resampleIter = cell(length(trialIndex_bool_memoryCorrect),1);
    locDistri_trialLevel_resampleIter = cell(length(trialIndex_bool_memoryCorrect),1);

    fun_sigmaBased_singleTrialPrecision_Name_v = autoGetFunName_myScripts('fun_sigmaBased_singleTrialPrecision', [targetPATH '\functions']);
    fun_sigmaBased_singleTrialPrecision = str2func(fun_sigmaBased_singleTrialPrecision_Name_v);
    
    %% Sub-step 2A: Compute memory precision on training set trials (noChoiceCorrect, choiceMemoryCorrect trials)
    %temp_trialIndex_bool = allMemoryCorrectBoolIndex;
    
    if if_trainingSet_memoryCorrect0_allMemory1 == 0
        temp_trialIndex_bool = allMemoryCorrectBoolIndex;
    elseif if_trainingSet_memoryCorrect0_allMemory1 == 1
        temp_trialIndex_bool = allMemoryBoolIndex;
        %temp_trialIndex_bool = allMemoryCorrectBoolIndex | allMemoryNonLengthErrorBoolIndex;
    end
    
    temp_trialIndex = find(temp_trialIndex_bool==1);
    temp_trialIndex_resample = temp_trialIndex(temp_trialIndex_valid_resample);
    
    for target_length=1:3
        if target_length == 1
            svm_train_lengthx_outputs = svm_train_length1_outputs;
        elseif target_length == 2
            svm_train_lengthx_outputs = svm_train_length2_outputs;
        elseif target_length == 3
            svm_train_lengthx_outputs = svm_train_length3_outputs;
        end
        temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
        
        temp_svm_resample = svm_train_lengthx_outputs.temp_svm_resample;
        for tempIter=1:resampleIterCount
            temptemp_svm = temp_svm_resample{tempIter};
            temptemp_trialIndex = temp_trialIndex_resample(tempIter,:);
            
            %temptemp_seqIndex = seqIndex(temptemp_trialIndex);
            %temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
            
            if if_trainingSet_memoryCorrect0_allMemory1 == 0
                temptemp_seqIndex = seqIndex(temptemp_trialIndex);
                temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
                temptemp_Posterior_2d_n11n = temptemp_svm.Posterior_2d_n11n;
                
            elseif if_trainingSet_memoryCorrect0_allMemory1 == 1
                
                if if_memeoryPrecision_stimuli0_response1 == 0
                    temptemp_seqIndex = seqIndex(temptemp_trialIndex);
                    temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
                    temptemp_Posterior_2d_n11n = temptemp_svm.Posterior_2d_n11n;
                    
                elseif if_memeoryPrecision_stimuli0_response1 == 1
                    % only compute allMemoryCorrectBoolIndex | allMemoryNonLengthErrorBoolIndex
                    temptemp_seqIndex = seqIndex(temptemp_trialIndex);
                    temptemp_seqIndex_response = seqIndex_response(temptemp_trialIndex);
                    
                    temptempBoolIndexA = ismember(temptemp_seqIndex,temp_range);
                    temptempBoolIndexB = ismember(temptemp_seqIndex_response,temp_range);
                    temptempBoolIndex = temptempBoolIndexA & temptempBoolIndexB;
                    
                    temptempIndexA = find(temptempBoolIndexA==true);
                    temptempIndex = find(temptempBoolIndex==true);
                    
                    temptempBoolIndexC = ismember(temptempIndexA,temptempIndex);
                    
                    temptemp_Posterior_2d_n11n = temptemp_svm.Posterior_2d_n11n;
                    temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_n11n(temptempBoolIndexC,:);
                end
            end
            
            temptemp_seqIndex_target_length = temptemp_seqIndex(temptempBoolIndex);
            temptemp_trialIndex_target_length = temptemp_trialIndex(temptempBoolIndex);
            
            temptemp_trialIndex_target_length;
            %temptemp_Posterior_2d_n11n = temptemp_svm.Posterior_2d_n11n;
            
            %temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
            %temp_p = temptemp_Posterior_2d_n11n;
            %temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
            %temp_p_production = prod(temp_p,2);
            
            
            temp_trialProb_stimuli_to_response = nan(size(temptemp_Posterior_2d_n11n,1),sum(numSeq(1:valid_length)));
            temp_trialProb_stimuli_to_response_n11n = nan(size(temptemp_Posterior_2d_n11n,1),sum(numSeq(1:valid_length)));
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                
                if isnan(temp_locDistri(1))
                    continue
                end
                for tempj=1:sum(numSeq(1:valid_length))
                    temptemp_boolIndex_location_seq = boolIndex_location_seq_T(tempj,:);
                    
                    temptemp_p = temp_locDistri;
                    temptemp_p(~temptemp_boolIndex_location_seq) = 1 - temptemp_p(~temptemp_boolIndex_location_seq);
                    temptemp_p_seq = prod(temptemp_p,2);
                    
                    temp_trialProb_stimuli_to_response(tempi,tempj) = temptemp_p_seq;
                end
                temp_trialProb_stimuli_to_response_n11n(tempi,:) = ...
                    temp_trialProb_stimuli_to_response(tempi,:)./sum(temp_trialProb_stimuli_to_response(tempi,:));
            end
            
            temptemp_p_production = nan(size(temptemp_Posterior_2d_n11n,1),1);
            temptemp_p_production_n11n = nan(size(temptemp_Posterior_2d_n11n,1),1);
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                temptemp_p_production(tempi) = ...
                    temp_trialProb_stimuli_to_response(tempi,temptemp_seqIndex_target_length(tempi));
                temptemp_p_production_n11n(tempi) = ...
                    temp_trialProb_stimuli_to_response_n11n(tempi,temptemp_seqIndex_target_length(tempi));
            end
            
            if if_pProd_trialLevel_n11n == 0
                temp_p_production = temptemp_p_production;
            elseif if_pProd_trialLevel_n11n == 1
                temp_p_production = temptemp_p_production_n11n;
            end
            
            
            
            a = 1;
            
            temp_precision = nan(length(temptemp_trialIndex_target_length),1);
            parfor tempi=1:length(temptemp_trialIndex_target_length)
                %for tempi=1:length(temptemp_trialIndex_target_length)
                
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                tempSeqIndex = temptemp_seqIndex_target_length(tempi);
                valid_length = 3;
                
                temp_score = score_stimuli_to_response(tempSeqIndex,1:sum(numSeq(1:valid_length)));
                
                temp_options = struct;
                %temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.boolIndex_location_seq_T = boolIndex_location_seq_T;
                %temp_options.Posterior_2d_n11n_mean = Posterior_2d_n11n_mean;
                temp_options.numSeq = numSeq;
                temp_options.valid_length = valid_length;
                %temp_options.score_stimuli_to_response = score_stimuli_to_response;
                temp_options.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
                
                temp_options.temp_locDistri = temp_locDistri;
                temp_options.temp_score = temp_score;
                temp_options.if_entropy = if_entropy;
                temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.tempSeqLength = target_length;
                
                
                
                temp_precision(tempi) = fun_sigmaBased_singleTrialPrecision(temp_options); %#ok<*PFBNS>
                
            end
            
            a = 1;
            temp_p_production;
            temp_precision;
            %temp1 = corr(temp_p_production,temp_precision);
            
            temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
            
            temp_r_loc = nan(size(temptemp_Posterior_2d_n11n,1),1);
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                [temp_r_loc(tempi),~] = corr(temptemp_Posterior_2d_n11n(tempi,:)',temp_boolIndex_location_seq(tempi,:)');
            end
            for tempi=1:length(temptemp_trialIndex_target_length)
                locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_r_loc(tempi)];
            end
            
            
            for tempi=1:length(temptemp_trialIndex_target_length)
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_p_production(tempi)];
                
                memoryPrecisionSigma_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecisionSigma_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_precision(tempi)];
                
                % marker here in 20250107, 14:08
                seqDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [seqDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_trialProb_stimuli_to_response_n11n(tempi,:)'];
                
                locDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [locDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temptemp_Posterior_2d_n11n(tempi,:)'];
                
                a = 1;
            end
            
        end
    end
    temptempTrialResampleNum = zeros(length(memoryPrecision_trialLevel_resampleIter),1);
    for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
        temptempBoolIndex(tempi) = ~isempty(memoryPrecision_trialLevel_resampleIter{tempi});
        temptempTrialResampleNum(tempi) = length(memoryPrecision_trialLevel_resampleIter{tempi});
    end
    aa1 = sum(seqIndex_valid<=41); %#ok<*NASGU>
    aa2 = sum(temptempTrialResampleNum>0);
    
    a1 = temptempTrialResampleNum(temptempTrialResampleNum>0);
    [M,I] = sort(a1); %#ok<*ASGLU>
    
    memoryPrecision_trialLevel_resampleMean = nan(length(memoryPrecision_trialLevel_resampleIter),1);
    memoryPrecisionSigma_trialLevel_resampleMean = nan(length(memoryPrecisionSigma_trialLevel_resampleIter),1);
    seqDistri_trialLevel_resampleMean = nan(length(seqDistri_trialLevel_resampleIter),sum(numSeq(1:valid_length)));
    locDistri_trialLevel_resampleMean = nan(length(locDistri_trialLevel_resampleIter),numFrames);    
    for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
        if temptempTrialResampleNum(tempi) == 0
            continue
        end
        memoryPrecision_trialLevel_resampleMean(tempi) = mean(memoryPrecision_trialLevel_resampleIter{tempi});
        memoryPrecisionSigma_trialLevel_resampleMean(tempi) = mean(memoryPrecisionSigma_trialLevel_resampleIter{tempi});
        seqDistri_trialLevel_resampleMean(tempi,:) = mean(seqDistri_trialLevel_resampleIter{tempi},2);
        locDistri_trialLevel_resampleMean(tempi,:) = mean(locDistri_trialLevel_resampleIter{tempi},2);
        
    end
    a1 = sum(~isnan(memoryPrecision_trialLevel_resampleMean));
    memoryPrecision_trialLevel_resampleMean;
    
    %% Sub-step 2B: Compute memory precision on noChoiceError and choiceMemoryError trials
    %temp_trialIndex_bool = allMemoryErrorBoolIndex;
    
    if if_trainingSet_memoryCorrect0_allMemory1 == 0
        temp_trialIndex_bool = allMemoryErrorBoolIndex;
    elseif if_trainingSet_memoryCorrect0_allMemory1 == 1
        temp_trialIndex_bool = allMemoryLengthErrorBoolIndex;
    end
    
    temp_trialIndex = find(temp_trialIndex_bool==1);
    for target_length=1:3
        if target_length == 1
            svm_train_lengthx_outputs = svm_train_length1_outputs;
        elseif target_length == 2
            svm_train_lengthx_outputs = svm_train_length2_outputs;
        elseif target_length == 3
            svm_train_lengthx_outputs = svm_train_length3_outputs;
        end
        temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
        
        temptemp_trialIndex = temp_trialIndex;
        
        if if_memeoryPrecision_stimuli0_response1 == 0
            temptemp_seqIndex = seqIndex(temptemp_trialIndex);
        elseif if_memeoryPrecision_stimuli0_response1 == 1
            temptemp_seqIndex = seqIndex_response(temptemp_trialIndex);
        end
        
        temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
        %temptempBoolIndex = true(1,length(temptemp_seqIndex));
        
        temptemp_seqIndex_target_length = temptemp_seqIndex(temptempBoolIndex);
        temptemp_trialIndex_target_length = temptemp_trialIndex(temptempBoolIndex);
        
        x = F_dff_decisionBin1(:,temptemp_trialIndex_target_length);
        %y = boolIndex_location_allTrial(:,temptemp_trialIndex_target_length);
        y = boolIndex_location_seq(:,temptemp_seqIndex_target_length);
        
        x_T = x';
        y_T = y';
        
        temp_svm_resample = svm_train_lengthx_outputs.temp_svm_resample;
        
        for tempIter=1:resampleIterCount
            temptemp_svm = temp_svm_resample{tempIter};
            
            temp_Mdl_CV_binary = temptemp_svm.temp_Mdl_CV_binary;
            
            temptemp_Posterior_2d_kfold = zeros(KFold_num,length(temptemp_trialIndex_target_length),numFrames);
            for temploc=1:numFrames
                for tempk=1:KFold_num
                    [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},x_T);% Very time-consuming!!!
                    if size(tempPosterior,2) == 1
                        tempPosterior_2 = tempPosterior(:,1);
                    else
                        tempPosterior_2 = tempPosterior(:,2);
                    end
                    temptemp_Posterior_2d_kfold(tempk,:,temploc) = tempPosterior_2;
                end
            end
            temptemp_Posterior_2d_kfoldMean = squeeze(mean(temptemp_Posterior_2d_kfold,1));
            temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_kfoldMean;
            
            %temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
            %temp_p = temptemp_Posterior_2d_n11n;
            %temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
            %temp_p_production = prod(temp_p,2);
            
            
            temp_trialProb_stimuli_to_response = nan(size(temptemp_Posterior_2d_n11n,1),sum(numSeq(1:valid_length)));
            temp_trialProb_stimuli_to_response_n11n = nan(size(temptemp_Posterior_2d_n11n,1),sum(numSeq(1:valid_length)));
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                
                if isnan(temp_locDistri(1))
                    continue
                end
                for tempj=1:sum(numSeq(1:valid_length))
                    temptemp_boolIndex_location_seq = boolIndex_location_seq_T(tempj,:);
                    
                    temptemp_p = temp_locDistri;
                    temptemp_p(~temptemp_boolIndex_location_seq) = 1 - temptemp_p(~temptemp_boolIndex_location_seq);
                    temptemp_p_seq = prod(temptemp_p,2);
                    
                    temp_trialProb_stimuli_to_response(tempi,tempj) = temptemp_p_seq;
                end
                temp_trialProb_stimuli_to_response_n11n(tempi,:) = ...
                    temp_trialProb_stimuli_to_response(tempi,:)./sum(temp_trialProb_stimuli_to_response(tempi,:));
            end
            
            temptemp_p_production = nan(size(temptemp_Posterior_2d_n11n,1),1);
            temptemp_p_production_n11n = nan(size(temptemp_Posterior_2d_n11n,1),1);
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                temptemp_p_production(tempi) = ...
                    temp_trialProb_stimuli_to_response(tempi,temptemp_seqIndex_target_length(tempi));
                temptemp_p_production_n11n(tempi) = ...
                    temp_trialProb_stimuli_to_response_n11n(tempi,temptemp_seqIndex_target_length(tempi));
            end
            
            if if_pProd_trialLevel_n11n == 0
                temp_p_production = temptemp_p_production;
            elseif if_pProd_trialLevel_n11n == 1
                temp_p_production = temptemp_p_production_n11n;
            end
            
            
            temp_precision = nan(length(temptemp_trialIndex_target_length),1);
            parfor tempi=1:length(temptemp_trialIndex_target_length)
                
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                tempSeqIndex = temptemp_seqIndex_target_length(tempi);
                valid_length = 3;
                
                temp_score = score_stimuli_to_response(tempSeqIndex,1:sum(numSeq(1:valid_length)));
                
                temp_options = struct;
                %temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.boolIndex_location_seq_T = boolIndex_location_seq_T;
                %temp_options.Posterior_2d_n11n_mean = Posterior_2d_n11n_mean;
                temp_options.numSeq = numSeq;
                temp_options.valid_length = valid_length;
                %temp_options.score_stimuli_to_response = score_stimuli_to_response;
                temp_options.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
                
                temp_options.temp_locDistri = temp_locDistri;
                temp_options.temp_score = temp_score;
                temp_options.if_entropy = if_entropy;
                temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.tempSeqLength = target_length;
                
                
                
                temp_precision(tempi) = fun_sigmaBased_singleTrialPrecision(temp_options); %#ok<*PFBNS>
                
            end
            
            a = 1;
            temp_p_production;
            temp_precision;
            %temp1 = corr(temp_p_production,temp_precision);
            
            
            temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
            
            temp_r_loc = nan(size(temptemp_Posterior_2d_n11n,1),1);
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                [temp_r_loc(tempi),~] = corr(temptemp_Posterior_2d_n11n(tempi,:)',temp_boolIndex_location_seq(tempi,:)');
            end
            for tempi=1:length(temptemp_trialIndex_target_length)
                locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_r_loc(tempi)];
            end
            
            
            for tempi=1:length(temptemp_trialIndex_target_length)
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_p_production(tempi)];
                
                memoryPrecisionSigma_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecisionSigma_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_precision(tempi)];
                
                seqDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [seqDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_trialProb_stimuli_to_response_n11n(tempi,:)'];
                
                locDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [locDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temptemp_Posterior_2d_n11n(tempi,:)'];
                
            end
            
        end
    end
    
    %% Sub-step 2C: Compute memory precision on choiceOffload trials
    
    temp_trialIndex_bool = choiceOffloadBoolIndex;
    temp_trialIndex = find(temp_trialIndex_bool==1);
    
    temptemp_Posterior_2d_resampleMean_length123 = zeros(length(temp_trialIndex_bool),numFrames);
    
    for target_length=1:3
        if target_length == 1
            svm_train_lengthx_outputs = svm_train_length1_outputs;
        elseif target_length == 2
            svm_train_lengthx_outputs = svm_train_length2_outputs;
        elseif target_length == 3
            svm_train_lengthx_outputs = svm_train_length3_outputs;
        end
        temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
        
        temptemp_trialIndex = temp_trialIndex;
        
        %if if_memeoryPrecision_stimuli0_response1 == 0
        %   temptemp_seqIndex = seqIndex(temptemp_trialIndex);
        %elseif if_memeoryPrecision_stimuli0_response1 == 1
        %   temptemp_seqIndex = seqIndex_response(temptemp_trialIndex);
        %end
        temptemp_seqIndex = seqIndex(temptemp_trialIndex);
        
        temptempBoolIndex = ismember(temptemp_seqIndex,temp_range);
        %temptempBoolIndex = true(1,length(temptemp_seqIndex));
        
        temptemp_seqIndex_target_length = temptemp_seqIndex(temptempBoolIndex);
        temptemp_trialIndex_target_length = temptemp_trialIndex(temptempBoolIndex);
        
        x = F_dff_decisionBin1(:,temptemp_trialIndex_target_length);
        %y = boolIndex_location_allTrial(:,temptemp_trialIndex_target_length);
        y = boolIndex_location_seq(:,temptemp_seqIndex_target_length);
        
        x_T = x';
        y_T = y';
        
        temp_svm_resample = svm_train_lengthx_outputs.temp_svm_resample;
        
        temptemp_Posterior_2d_iter = zeros(resampleIterCount,length(temptemp_trialIndex_target_length),numFrames);
        
        for tempIter=1:resampleIterCount
            temptemp_svm = temp_svm_resample{tempIter};
            
            temp_Mdl_CV_binary = temptemp_svm.temp_Mdl_CV_binary;
            
            temptemp_Posterior_2d_kfold = zeros(KFold_num,length(temptemp_trialIndex_target_length),numFrames);
            for temploc=1:numFrames
                for tempk=1:KFold_num
                    [~,~,~,tempPosterior] = predict(temp_Mdl_CV_binary{temploc}.Trained{tempk},x_T);% Very time-consuming!!!
                    if size(tempPosterior,2) == 1
                        tempPosterior_2 = tempPosterior(:,1);
                    else
                        tempPosterior_2 = tempPosterior(:,2);
                    end
                    temptemp_Posterior_2d_kfold(tempk,:,temploc) = tempPosterior_2;
                end
            end
            temptemp_Posterior_2d_kfoldMean = squeeze(mean(temptemp_Posterior_2d_kfold,1));
            temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_kfoldMean;
            temptemp_Posterior_2d_iter(tempIter,:,:) = temptemp_Posterior_2d_kfoldMean;
            
            if if_memeoryPrecision_stimuli0_response1 == 0
                temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
                
            elseif if_memeoryPrecision_stimuli0_response1 == 1
                if if_inferOffloadResponse == 1
                    % infer response seq from temptemp_Posterior_2d_n11n
                    [M,I] = sort(temptemp_Posterior_2d_n11n,2,'descend');
                    temp_boolIndex_location_seq = false(size(temptemp_Posterior_2d_n11n,1),numFrames);
                    for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                        temp_boolIndex_location_seq(tempi,I(tempi,1:target_length)) = true;
                    end
                elseif if_inferOffloadResponse == 0
                    temp_boolIndex_location_seq = boolIndex_location_seq(:,temptemp_seqIndex_target_length)';
                end
            end
            
            %temp_p = temptemp_Posterior_2d_n11n;
            %temp_p(~temp_boolIndex_location_seq) = 1 - temp_p(~temp_boolIndex_location_seq);
            %temp_p_production = prod(temp_p,2);
            
            
            temp_trialProb_stimuli_to_response = nan(size(temptemp_Posterior_2d_n11n,1),sum(numSeq(1:valid_length)));
            temp_trialProb_stimuli_to_response_n11n = nan(size(temptemp_Posterior_2d_n11n,1),sum(numSeq(1:valid_length)));
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                
                if isnan(temp_locDistri(1))
                    continue
                end
                for tempj=1:sum(numSeq(1:valid_length))
                    temptemp_boolIndex_location_seq = boolIndex_location_seq_T(tempj,:);
                    
                    temptemp_p = temp_locDistri;
                    temptemp_p(~temptemp_boolIndex_location_seq) = 1 - temptemp_p(~temptemp_boolIndex_location_seq);
                    temptemp_p_seq = prod(temptemp_p,2);
                    
                    temp_trialProb_stimuli_to_response(tempi,tempj) = temptemp_p_seq;
                end
                temp_trialProb_stimuli_to_response_n11n(tempi,:) = ...
                    temp_trialProb_stimuli_to_response(tempi,:)./sum(temp_trialProb_stimuli_to_response(tempi,:));
            end
            
            temptemp_p_production = nan(size(temptemp_Posterior_2d_n11n,1),1);
            temptemp_p_production_n11n = nan(size(temptemp_Posterior_2d_n11n,1),1);
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                temptemp_p_production(tempi) = ...
                    temp_trialProb_stimuli_to_response(tempi,temptemp_seqIndex_target_length(tempi));
                temptemp_p_production_n11n(tempi) = ...
                    temp_trialProb_stimuli_to_response_n11n(tempi,temptemp_seqIndex_target_length(tempi));
            end
            
            if if_pProd_trialLevel_n11n == 0
                temp_p_production = temptemp_p_production;
            elseif if_pProd_trialLevel_n11n == 1
                temp_p_production = temptemp_p_production_n11n;
            end
            
            
            temp_precision = nan(length(temptemp_trialIndex_target_length),1);
            parfor tempi=1:length(temptemp_trialIndex_target_length)
                
                temp_locDistri = temptemp_Posterior_2d_n11n(tempi,:);
                tempSeqIndex = temptemp_seqIndex_target_length(tempi);
                valid_length = 3;
                
                temp_score = score_stimuli_to_response(tempSeqIndex,1:sum(numSeq(1:valid_length)));
                
                temp_options = struct;
                %temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.boolIndex_location_seq_T = boolIndex_location_seq_T;
                %temp_options.Posterior_2d_n11n_mean = Posterior_2d_n11n_mean;
                temp_options.numSeq = numSeq;
                temp_options.valid_length = valid_length;
                %temp_options.score_stimuli_to_response = score_stimuli_to_response;
                temp_options.if_precision_meanProb0_sumProb1 = if_precision_meanProb0_sumProb1;
                
                temp_options.temp_locDistri = temp_locDistri;
                temp_options.temp_score = temp_score;
                temp_options.if_entropy = if_entropy;
                temp_options.tempSeqIndex = tempSeqIndex;
                temp_options.tempSeqLength = target_length;
                
                
                temp_precision(tempi) = fun_sigmaBased_singleTrialPrecision(temp_options); %#ok<*PFBNS>
                
            end
            
            a = 1;
            temp_p_production;
            temp_precision;
            %temp1 = corr(temp_p_production,temp_precision);
            
            
            temp_boolIndex_location_seq;
            
            temp_r_loc = nan(size(temptemp_Posterior_2d_n11n,1),1);
            for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                [temp_r_loc(tempi),~] = corr(temptemp_Posterior_2d_n11n(tempi,:)',temp_boolIndex_location_seq(tempi,:)');
            end
            for tempi=1:length(temptemp_trialIndex_target_length)
                locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [locCorr_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_r_loc(tempi)];
            end
            
            
            for tempi=1:length(temptemp_trialIndex_target_length)
                memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecision_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_p_production(tempi)];
                
                memoryPrecisionSigma_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [memoryPrecisionSigma_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_precision(tempi)];
                
                seqDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [seqDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temp_trialProb_stimuli_to_response_n11n(tempi,:)'];
                
                locDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)} = ...
                    [locDistri_trialLevel_resampleIter{temptemp_trialIndex_target_length(tempi)}, temptemp_Posterior_2d_n11n(tempi,:)'];
                
            end
            
        end
        temptemp_Posterior_2d_resampleMean = squeeze(mean(temptemp_Posterior_2d_iter,1));
        
        temptemp_Posterior_2d_resampleMean_length123(temptemp_trialIndex_target_length,:) = temptemp_Posterior_2d_resampleMean;
    end
    %temptemp_Posterior_2d_resampleMean_length123_B = ...
    %    temptemp_Posterior_2d_resampleMean_length123(temp_trialIndex_bool&(seqIndex<=41),:);
    
    if if_memeoryPrecision_stimuli0_response1 == 1
        for target_length=1:3
            temp_range = sum(numSeq(1:(target_length-1)))+1:sum(numSeq(1:target_length));
            
            temptempBoolIndex = ismember(seqIndex,temp_range);
            temptempBoolIndex2 = temptempBoolIndex & choiceOffloadBoolIndex;
            
            if if_inferOffloadResponse == 1
                temptemp_Posterior_2d_n11n = temptemp_Posterior_2d_resampleMean_length123(temptempBoolIndex2,:);
                
                % infer response seq from temptemp_Posterior_2d_n11n
                [M,I] = sort(temptemp_Posterior_2d_n11n,2,'descend');
                temp_boolIndex_location_seq = false(size(temptemp_Posterior_2d_n11n,1),numFrames);
                for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                    temp_boolIndex_location_seq(tempi,I(tempi,1:target_length)) = true;
                end
                
                temptempSeqIndex_response = nan(1,size(temp_boolIndex_location_seq,1));
                for tempi=1:size(temptemp_Posterior_2d_n11n,1)
                    temptempBoolIndex_seq = temp_boolIndex_location_seq(tempi,:);
                    for tempj=1:sum(numSeq(1:3))
                        temptempBoolIndex_seq2 = boolIndex_location_seq_T(tempj,:);
                        if sum(temptempBoolIndex_seq==temptempBoolIndex_seq2) == numFrames
                            temptempSeqIndex_response(tempi) = tempj;
                            break
                        end
                    end
                end
                temptempSeqIndex_response;
                seqIndex_response(temptempBoolIndex2) = temptempSeqIndex_response;
            elseif if_inferOffloadResponse == 0
                seqIndex_response(temptempBoolIndex2) = seqIndex(temptempBoolIndex2); %#ok<*SAGROW>
            end
            
        end
        
    end
    
    
    temptempTrialResampleNum = zeros(length(memoryPrecision_trialLevel_resampleIter),1);
    for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
        temptempBoolIndex(tempi) = ~isempty(memoryPrecision_trialLevel_resampleIter{tempi});
        temptempTrialResampleNum(tempi) = length(memoryPrecision_trialLevel_resampleIter{tempi});
    end
    aa1 = sum(seqIndex<=41);
    aa2 = sum(temptempTrialResampleNum>0);
    
    a1 = temptempTrialResampleNum(temptempTrialResampleNum>0);
    [M,I] = sort(a1);
    
    memoryPrecision_trialLevel_resampleMean = nan(length(memoryPrecision_trialLevel_resampleIter),1);
    memoryPrecisionSigma_trialLevel_resampleMean = nan(length(memoryPrecisionSigma_trialLevel_resampleIter),1);
    seqDistri_trialLevel_resampleMean = nan(length(seqDistri_trialLevel_resampleIter),sum(numSeq(1:valid_length)));
    locDistri_trialLevel_resampleMean = nan(length(locDistri_trialLevel_resampleIter),numFrames);    
    for tempi=1:length(memoryPrecision_trialLevel_resampleIter)
        if temptempTrialResampleNum(tempi) == 0
            continue
        end
        memoryPrecision_trialLevel_resampleMean(tempi) = mean(memoryPrecision_trialLevel_resampleIter{tempi});
        memoryPrecisionSigma_trialLevel_resampleMean(tempi) = mean(memoryPrecisionSigma_trialLevel_resampleIter{tempi});
        seqDistri_trialLevel_resampleMean(tempi,:) = mean(seqDistri_trialLevel_resampleIter{tempi},2);
        locDistri_trialLevel_resampleMean(tempi,:) = mean(locDistri_trialLevel_resampleIter{tempi},2);
        
    end
    a1 = sum(~isnan(memoryPrecision_trialLevel_resampleMean));
    memoryPrecision_trialLevel_resampleMean;
    
    a1 = sum(~isnan(memoryPrecision_trialLevel_resampleMean));
    a2 = sum(seqIndex<=41);
    %fprintf('t_step2 = %.1f secs.\n',toc(t_test));
    fprintf('valid trial num of memoryPrecision is %d, theretical trial num (stimuli) is %d.\n',a1,a2);
    
    
    locCorr_trialLevel_resampleMean = nan(length(locCorr_trialLevel_resampleIter),1);
    for tempi=1:length(locCorr_trialLevel_resampleIter)
        if temptempTrialResampleNum(tempi) == 0
            continue
        end
        locCorr_trialLevel_resampleMean(tempi) = mean(locCorr_trialLevel_resampleIter{tempi});
    end
    locCorr_trialLevel = locCorr_trialLevel_resampleMean;
    
    
end
%% Sub-step 2D: Compute tuning based on memory precision
F_dff_decisionBin1;
memoryPrecision_trialLevel_resampleMean;

if if_memoryPrecisionTuning_decoder0_behavior1 == 0
    if if_memoryPrecision_accuracy0_sigma1 == 0
        temp_y = memoryPrecision_trialLevel_resampleMean;
    elseif if_memoryPrecision_accuracy0_sigma1 == 1
        temp_y = memoryPrecisionSigma_trialLevel_resampleMean;
    end
elseif if_memoryPrecisionTuning_decoder0_behavior1 == 1
    temp_y = nan(trial_num,1);
    temp_y(allMemoryCorrectBoolIndex) = 1;
    temp_y(allMemoryErrorBoolIndex) = 0;
end

% x = F_dff_decisionBin1';
% y = temp_y;

x = temp_y;
y = F_dff_decisionBin1_zScore';


r2 = zeros(roiNum,1);
p = zeros(roiNum,1);
beta = zeros(roiNum,1);
r = zeros(roiNum,1);
parfor tempi=1:roiNum
    % for tempi=1:roiNum
    warning('off');
    %temp_mdl = fitglm(x(:,tempi),y);
    temp_mdl = fitglm(x,y(:,tempi));
    r2(tempi) = temp_mdl.Rsquared.Adjusted;
    p(tempi) = temp_mdl.Coefficients.pValue(2);
    beta(tempi) = temp_mdl.Coefficients.Estimate(2);
    %r(tempi) = corr(x(~isnan(y),tempi),y(~isnan(y)));
    r(tempi) = corr(x(~isnan(x)),y(~isnan(x),tempi));
    warning('on');
end
r2; %#ok<*VUNUS>
p;

[M,I] = sort(r2,'descend');

r2_memoryPrecision = r2;
p_memoryPrecision = p;
beta_memoryPrecision = beta;
r_memoryPrecision = r;

beta_memoryPrecision_raw = beta_memoryPrecision;
if if_beta_n11n == 1
    beta_memoryPrecision = rescale(beta_memoryPrecision_raw,-1,1);
end

fprintf('t_step2 = %.1f secs.\n',toc(t_test));

a = 1;
% end

%% Step 3: Tuning for decision (preference of choiceMemory)
a = 1;

% temp_load = load([output_path,'\','trial_para.mat']);
% trial_para = temp_load.trial_para;
% trial_para_choiceCondition_flag = trial_para.choiceCondition_flag;
%
% temp_trial_para_ifSelectOffloading = trial_para_ifSelectOffloading==1;
%
% choiceBoolIndex = trial_para_choiceCondition_flag == 2;
% allMemoryBoolIndex = ~temp_trial_para_ifSelectOffloading;
% choiceMemoryBoolIndex = choiceBoolIndex & allMemoryBoolIndex;
% choiceOffloadBoolIndex = choiceBoolIndex & temp_trial_para_ifSelectOffloading;
%
% choiceMemoryCorrectBoolIndex = choiceMemoryBoolIndex & trialIndex_bool_memoryCorrect;
% choiceMemoryErrorBoolIndex = choiceMemoryBoolIndex & ~trialIndex_bool_memoryCorrect;


choiceBoolIndex;
choiceMemoryBoolIndex;


% x_raw = F_dff_decisionBin1';
% x = x_raw(choiceBoolIndex,:);
% y_raw = choiceMemoryBoolIndex';
% y = y_raw(choiceBoolIndex);

x_raw = choiceMemoryBoolIndex';
x = x_raw(choiceBoolIndex);
y_raw = F_dff_decisionBin1_zScore';
y = y_raw(choiceBoolIndex,:);


r2 = zeros(roiNum,1);
p = zeros(roiNum,1);
beta = zeros(roiNum,1);
r = zeros(roiNum,1);
parfor tempi=1:roiNum
    % for tempi=1:roiNum
    warning('off');
    %temp_mdl = fitglm(x(:,tempi),y);
    temp_mdl = fitglm(x,y(:,tempi));
    r2(tempi) = temp_mdl.Rsquared.Adjusted;
    p(tempi) = temp_mdl.Coefficients.pValue(2);
    beta(tempi) = temp_mdl.Coefficients.Estimate(2);
    %r(tempi) = corr(x(:,tempi),y);
    r(tempi) = corr(x,y(:,tempi));
    warning('on');
end
r2; %#ok<*VUNUS>
p;

[M,I] = sort(r2,'descend');

r2_choiceMemory = r2;
p_choiceMemory = p;
beta_choiceMemory = beta;
r_choiceMemory = r;

beta_choiceMemory_raw = beta_choiceMemory;
if if_beta_n11n == 1
    beta_choiceMemory = rescale(beta_choiceMemory_raw,-1,1);
end


%% Tuning for decision (baseline period)
% temptemp_range = (baselinePeriod_interval(2)-8):baselinePeriod_interval(3);
% % temptemp_range = (baselinePeriod_interval(1)):baselinePeriod_interval(3);
% F_dff_baselineBin = mean(F_dff_baselinePeriod(:,:,temptemp_range),3);

% x_raw = F_dff_baselineBin';
% x = x_raw(choiceBoolIndex,:);
% y_raw = choiceMemoryBoolIndex';
% y = y_raw(choiceBoolIndex);

x_raw = choiceMemoryBoolIndex';
x = x_raw(choiceBoolIndex);
y_raw = F_dff_baselineBin_zScore';
y = y_raw(choiceBoolIndex,:);


r2 = zeros(roiNum,1);
p = zeros(roiNum,1);
beta = zeros(roiNum,1);
r = zeros(roiNum,1);
parfor tempi=1:roiNum
    % for tempi=1:roiNum
    warning('off');
    %temp_mdl = fitglm(x(:,tempi),y);
    temp_mdl = fitglm(x,y(:,tempi));
    r2(tempi) = temp_mdl.Rsquared.Adjusted;
    p(tempi) = temp_mdl.Coefficients.pValue(2);
    beta(tempi) = temp_mdl.Coefficients.Estimate(2);
    %r(tempi) = corr(x(:,tempi),y);
    r(tempi) = corr(x,y(:,tempi));
    warning('on');
end
r2; %#ok<*VUNUS>
p;

[M,I] = sort(r2,'descend');

r2_choiceMemory_baseline = r2;
p_choiceMemory_baseline = p;
beta_choiceMemory_baseline = beta;
r_choiceMemory_baseline = r;

beta_choiceMemory_baseline_raw = beta_choiceMemory_baseline;
if if_beta_n11n == 1
    beta_choiceMemory_baseline = rescale(beta_choiceMemory_baseline_raw,-1,1);
end

%% Others
r2_rProb;
p_rProb;

r2_memoryPrecision;
p_memoryPrecision;

r2_choiceMemory;
p_choiceMemory;

a = 1;

% threshold_p = 0.001;%0.05,0.01

a1 = sum(p_rProb<threshold_p);
a2 = sum(p_memoryPrecision<threshold_p);
a3 = sum(p_choiceMemory<threshold_p);


r2_allFeatures = [r2_rProb,r2_memoryPrecision,r2_choiceMemory];
p_allFeatures = [p_rProb,p_memoryPrecision,p_choiceMemory];

r2_allFeatures2 = r2_allFeatures;
% r2_allFeatures2(p_rProb>=threshold_p,:) = -1;

r2_allFeatures3 = r2_allFeatures;
% r2_allFeatures3 = r2_allFeatures(p_rProb<threshold_p,:);
% cellIndex_allFeatures3 = find(p_rProb<threshold_p);
% cellIndex_suite2p_allFeatures3 = cellIndex_suite2p(cellIndex_allFeatures3);
cellIndex_suite2p_allFeatures3 = cellIndex_suite2p;

[M_A,I_A] = sort(r2_allFeatures3(:,1),'descend');
[M_B,I_B] = sort(r2_allFeatures3(:,2),'descend');
[M_C,I_C] = sort(r2_allFeatures3(:,3),'descend');

r2_allFeatures_rProbSorted = r2_allFeatures3(I_A,:);
r2_allFeatures_memoryPrecisionSorted = r2_allFeatures3(I_B,:);
r2_allFeatures_choiceMemorySorted = r2_allFeatures3(I_C,:);


% p_allFeatures3 = p_allFeatures(p_rProb<threshold_p,:);
p_allFeatures3 = p_allFeatures;
r2_allFeatures4 = r2_allFeatures3;
r2_allFeatures4(p_allFeatures3(:,2)>=threshold_p,2) = -1;
r2_allFeatures4(p_allFeatures3(:,3)>=threshold_p,3) = -1;

threshold_r2_memoryPrecision = min(r2_allFeatures4(r2_allFeatures4(:,2)>-1,2));
threshold_r2_choiceMemory = min(r2_allFeatures4(r2_allFeatures4(:,3)>-1,3));

threshold_r2_choiceMemory_baseline = (threshold_r2_memoryPrecision+threshold_r2_choiceMemory)/2;


[r_AB,p_AB] = corr(r2_allFeatures3(:,1),r2_allFeatures3(:,2));
[r_AC,p_AC] = corr(r2_allFeatures3(:,1),r2_allFeatures3(:,3));
[r_BC,p_BC] = corr(r2_allFeatures3(:,2),r2_allFeatures3(:,3));

% lowThreshold_percent = 35;
% highThreshold_percent = 100-lowThreshold_percent;



x = r2_allFeatures3(:,2);
y = r2_allFeatures3(:,3);

temptempBoolIndex_A = x<=threshold_r2_memoryPrecision;
temptempBoolIndex_B = x>threshold_r2_memoryPrecision;
temptempBoolIndex_C = y<=threshold_r2_choiceMemory;
temptempBoolIndex_D = y>threshold_r2_choiceMemory;
num_xlow = sum(temptempBoolIndex_A);
num_xhigh = sum(temptempBoolIndex_B);
num_ylow = sum(temptempBoolIndex_C);
num_yhigh = sum(temptempBoolIndex_D);

temptempBoolIndex_AC = temptempBoolIndex_A & temptempBoolIndex_C;
temptempBoolIndex_AD = temptempBoolIndex_A & temptempBoolIndex_D;
temptempBoolIndex_BC = temptempBoolIndex_B & temptempBoolIndex_C;
temptempBoolIndex_BD = temptempBoolIndex_B & temptempBoolIndex_D;
num_xlow_ylow = sum(temptempBoolIndex_AC);
num_xlow_yhigh = sum(temptempBoolIndex_AD);
num_xhigh_ylow = sum(temptempBoolIndex_BC);
num_xhigh_yhigh = sum(temptempBoolIndex_BD);

a = 1;
cellIndex_suite2p_allFeatures3;

cellIndex_suite2p_xlow_ylow = cellIndex_suite2p_allFeatures3(temptempBoolIndex_AC);
cellIndex_suite2p_xlow_yhigh = cellIndex_suite2p_allFeatures3(temptempBoolIndex_AD);
cellIndex_suite2p_xhigh_ylow = cellIndex_suite2p_allFeatures3(temptempBoolIndex_BC);
cellIndex_suite2p_xhigh_yhigh = cellIndex_suite2p_allFeatures3(temptempBoolIndex_BD);

cellIndex_suite2p_xhigh_ylow;
cellIndex_suite2p_xlow_yhigh;
cellIndex_suite2p_xhigh_yhigh;

% cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow = cellIndex_suite2p_xhigh_ylow;
% cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh = cellIndex_suite2p_xlow_yhigh;
% cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh = cellIndex_suite2p_xhigh_yhigh;
% cellIndex_suite2p_memoryPrecisionLow_choiceMemoryLow = cellIndex_suite2p_xlow_ylow;


a = 1;
r2_allFeatures3_xhigh_ylow = r2_allFeatures3(temptempBoolIndex_BC,:);
r2_allFeatures3_xlow_yhigh = r2_allFeatures3(temptempBoolIndex_AD,:);
r2_allFeatures3_xhigh_yhigh = r2_allFeatures3(temptempBoolIndex_BD,:);

[M_BC,I_BC] = sort(r2_allFeatures3_xhigh_ylow(:,2),'descend');
[M_AD,I_AD] = sort(r2_allFeatures3_xlow_yhigh(:,3),'descend');
% [M_BD,I_BD] = sort(sum(r2_allFeatures3_xhigh_yhigh(:,2:3),2),'descend'); %#ok<*NBRAK>
[M_BD,I_BD] = sort(r2_allFeatures3_xhigh_yhigh(:,2),'descend'); %#ok<*NBRAK>

r2_allFeatures3_xhigh_ylow_sorted = r2_allFeatures3_xhigh_ylow(I_BC,:);
r2_allFeatures3_xlow_yhigh_sorted = r2_allFeatures3_xlow_yhigh(I_AD,:);
r2_allFeatures3_xhigh_yhigh_sorted = r2_allFeatures3_xhigh_yhigh(I_BD,:);

topX = 10;
topX = min([topX,length(I_BC),length(I_AD),length(I_BD)]);

I_BC_topX = I_BC(1:topX);
I_AD_topX = I_AD(1:topX);
I_BD_topX = I_BD(1:topX);

cellIndex_suite2p_xhigh_ylow_topX = cellIndex_suite2p_xhigh_ylow(I_BC_topX);
cellIndex_suite2p_xlow_yhigh_topX = cellIndex_suite2p_xlow_yhigh(I_AD_topX);
cellIndex_suite2p_xhigh_yhigh_topX = cellIndex_suite2p_xhigh_yhigh(I_BD_topX);

if if_exampleNeuron_trialTuning0_seqTuning1 == 0
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow = cellIndex_suite2p_xhigh_ylow;
    cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh = cellIndex_suite2p_xlow_yhigh;
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh = cellIndex_suite2p_xhigh_yhigh;
    cellIndex_suite2p_memoryPrecisionLow_choiceMemoryLow = cellIndex_suite2p_xlow_ylow;
    
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow_topX = cellIndex_suite2p_xhigh_ylow_topX;
    cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh_topX = cellIndex_suite2p_xlow_yhigh_topX;
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh_topX = cellIndex_suite2p_xhigh_yhigh_topX;
    
    cellIndex_suite2p_memoryMetaMismatch = [cellIndex_suite2p_xhigh_ylow_topX;cellIndex_suite2p_xlow_yhigh_topX];
    cellIndex_suite2p_memoryMetaMatch = cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh_topX;
    
    cellIndex_suite2p_choiceMemoryBaselineHigh = cellIndex_suite2p(r2_choiceMemory_baseline>threshold_r2_choiceMemory_baseline);
    
    
elseif if_exampleNeuron_trialTuning0_seqTuning1 == 1
    x = r2_seqPrecision;
    y = r2_rProb;
    
    temptempBoolIndex_A = x<=threshold_r2_seqPrecision;
    temptempBoolIndex_B = x>threshold_r2_seqPrecision;
    temptempBoolIndex_C = y<=threshold_r2_rProb;
    temptempBoolIndex_D = y>threshold_r2_rProb;
    num_xlow = sum(temptempBoolIndex_A);
    num_xhigh = sum(temptempBoolIndex_B);
    num_ylow = sum(temptempBoolIndex_C);
    num_yhigh = sum(temptempBoolIndex_D);
    
    temptempBoolIndex_AC = temptempBoolIndex_A & temptempBoolIndex_C;
    temptempBoolIndex_AD = temptempBoolIndex_A & temptempBoolIndex_D;
    temptempBoolIndex_BC = temptempBoolIndex_B & temptempBoolIndex_C;
    temptempBoolIndex_BD = temptempBoolIndex_B & temptempBoolIndex_D;
    num_xlow_ylow = sum(temptempBoolIndex_AC);
    num_xlow_yhigh = sum(temptempBoolIndex_AD);
    num_xhigh_ylow = sum(temptempBoolIndex_BC);
    num_xhigh_yhigh = sum(temptempBoolIndex_BD);
    
    cellIndex_suite2p_xlow_ylow = cellIndex_suite2p(temptempBoolIndex_AC);
    cellIndex_suite2p_xlow_yhigh = cellIndex_suite2p(temptempBoolIndex_AD);
    cellIndex_suite2p_xhigh_ylow = cellIndex_suite2p(temptempBoolIndex_BC);
    cellIndex_suite2p_xhigh_yhigh = cellIndex_suite2p(temptempBoolIndex_BD);
    
    temp_r2 = [x,y];
    temp_r2_xhigh_ylow = temp_r2(temptempBoolIndex_BC,:);
    temp_r2_xlow_yhigh = temp_r2(temptempBoolIndex_AD,:);
    temp_r2_xhigh_yhigh = temp_r2(temptempBoolIndex_BD,:);
    
    [M_BC,I_BC] = sort(temp_r2_xhigh_ylow(:,1),'descend');
    [M_AD,I_AD] = sort(temp_r2_xlow_yhigh(:,2),'descend');
    [M_BD,I_BD] = sort(temp_r2_xhigh_yhigh(:,1),'descend'); %#ok<*NBRAK>
    
    temp_r2_xhigh_ylow_sorted = temp_r2_xhigh_ylow(I_BC,:);
    temp_r2_xlow_yhigh_sorted = temp_r2_xlow_yhigh(I_AD,:);
    temp_r2_xhigh_yhigh_sorted = temp_r2_xhigh_yhigh(I_BD,:);
    
    topX = 100;
    topX = min([topX,length(I_BC),length(I_AD),length(I_BD)]);
    
    I_BC_topX = I_BC(1:topX);
    I_AD_topX = I_AD(1:topX);
    I_BD_topX = I_BD(1:topX);
    
    cellIndex_suite2p_xhigh_ylow_topX = cellIndex_suite2p_xhigh_ylow(I_BC_topX);
    cellIndex_suite2p_xlow_yhigh_topX = cellIndex_suite2p_xlow_yhigh(I_AD_topX);
    cellIndex_suite2p_xhigh_yhigh_topX = cellIndex_suite2p_xhigh_yhigh(I_BD_topX);
    
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow = cellIndex_suite2p_xhigh_ylow;
    cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh = cellIndex_suite2p_xlow_yhigh;
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh = cellIndex_suite2p_xhigh_yhigh;
    cellIndex_suite2p_memoryPrecisionLow_choiceMemoryLow = cellIndex_suite2p_xlow_ylow;
    
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow_topX = cellIndex_suite2p_xhigh_ylow_topX;
    cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh_topX = cellIndex_suite2p_xlow_yhigh_topX;
    cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh_topX = cellIndex_suite2p_xhigh_yhigh_topX;
    
    cellIndex_suite2p_memoryMetaMismatch = [cellIndex_suite2p_xhigh_ylow_topX;cellIndex_suite2p_xlow_yhigh_topX];
    cellIndex_suite2p_memoryMetaMatch = cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh_topX;
    
    cellIndex_suite2p_choiceMemoryBaselineHigh = cellIndex_suite2p(r2_rProb_baseline>threshold_r2_rProb_baseline);
    
end
a = 1;
%% Others
lowThreshold_percent = 50;%30
highThreshold_percent = 100-lowThreshold_percent;

% memoryPrecision_trialLevel_resampleMean;
if if_memoryPrecision_accuracy0_sigma1 == 0
    memoryPrecision_trialLevel = memoryPrecision_trialLevel_resampleMean;
elseif if_memoryPrecision_accuracy0_sigma1 == 1
    memoryPrecision_trialLevel = memoryPrecisionSigma_trialLevel_resampleMean;
    if if_entropy == 0
        memoryPrecision_trialLevel = rescale(memoryPrecision_trialLevel,0.001,0.999);
    end
end

if if_memeoryPrecision_stimuli0_response1 == 0
    memoryPrecision_trialLevel_stimuli = memoryPrecision_trialLevel;
elseif if_memeoryPrecision_stimuli0_response1 == 1
    memoryPrecision_trialLevel_response = memoryPrecision_trialLevel;
end

if if_memeoryPrecision_stimuli0_response1 == 0
    memoryPrecision_trialLevel_resampleMean_stimuli = memoryPrecision_trialLevel_resampleMean;
elseif if_memeoryPrecision_stimuli0_response1 == 1
    memoryPrecision_trialLevel_resampleMean_response = memoryPrecision_trialLevel_resampleMean;
end

seqDistri_trialLevel = seqDistri_trialLevel_resampleMean;
locDistri_trialLevel = locDistri_trialLevel_resampleMean;

if if_memeoryPrecision_stimuli0_response1 == 0
    seqDistri_trialLevel_stimuli = seqDistri_trialLevel;
    locDistri_trialLevel_stimuli = locDistri_trialLevel;    
elseif if_memeoryPrecision_stimuli0_response1 == 1
    seqDistri_trialLevel_response = seqDistri_trialLevel;
    locDistri_trialLevel_response = locDistri_trialLevel;    
end

choiceBoolIndex;
choiceMemoryBoolIndex;
choiceOffloadBoolIndex;

trial_num = length(choiceBoolIndex);

% lowThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel,lowThreshold_percent);
% highThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel,highThreshold_percent);
% lowThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),lowThreshold_percent);
% highThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),highThreshold_percent);

% temp1 = prctile(memoryPrecision_trialLevel(choiceMemoryBoolIndex),[lowThreshold_percent,highThreshold_percent]);
% lowThreshold_memoryPrecision_choiceMemory = temp1(1);
% highThreshold_memoryPrecision_choiceMemory = temp1(2);
% temp2 = prctile(memoryPrecision_trialLevel(choiceOffloadBoolIndex),[lowThreshold_percent,highThreshold_percent]);
% lowThreshold_memoryPrecision_choiceOffload = temp2(1);
% highThreshold_memoryPrecision_choiceOffload = temp2(2);

if if_memoryPrecisionThreshold_median0_optimal1 == 0
    lowThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),lowThreshold_percent);
    highThreshold_memoryPrecision = prctile(memoryPrecision_trialLevel(choiceBoolIndex),highThreshold_percent);
    
elseif if_memoryPrecisionThreshold_median0_optimal1 == 1
    
    %temptempBoolIndex_A = (~isnan(memoryPrecision_trialLevel)) & (~isnan(meta_trialLevel));
    temptempBoolIndex_A = (~isnan(memoryPrecision_trialLevel));
    
    memoryPrecision_trialLevel_choiceMemoryCorrect = memoryPrecision_trialLevel(temptempBoolIndex_A & choiceMemoryCorrectBoolIndex');
    memoryPrecision_trialLevel_choiceMemoryError = memoryPrecision_trialLevel(temptempBoolIndex_A & choiceMemoryErrorBoolIndex');
    
    
    x2 = memoryPrecision_trialLevel_choiceMemoryCorrect;
    n=100;
    n=2^ceil(log2(n)); % round up n to the next power of 2;
    [pdf2,xmesh2,bandwidth2] = ksdensity(x2','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
    temp_ratio2 = 1;
    y2 = pdf2*temp_ratio2;
    
    x3 = memoryPrecision_trialLevel_choiceMemoryError;
    n=100;
    n=2^ceil(log2(n)); % round up n to the next power of 2;
    [pdf3,xmesh3,bandwidth3] = ksdensity(x3','NumPoints',n,'Function','pdf','Support',[-0.01 1.01],'BoundaryCorrection', 'Reflection');
    temp_ratio3 = 1;
    y3 = pdf3*temp_ratio3;
    
    
    temp1 = y2;
    temp2 = y3;
    
    xmesh = xmesh3;
    
    if true
        temp_thresholdRange = 0.01:0.001:0.99;
        
        choiceMemoryCorrectError_hit_minus_falseAlarm_multi = zeros(1,length(temp_thresholdRange));
        choiceMemoryCorrectError_hit_multi = zeros(1,length(temp_thresholdRange));
        choiceMemoryCorrectError_falseAlarm_multi = zeros(1,length(temp_thresholdRange));
        
        for temptempi=1:length(temp_thresholdRange)
            
            temp_threshold = temp_thresholdRange(temptempi);
            
            temptempBoolIndex = xmesh > temp_threshold;
            temp_hit = sum(temp1(temptempBoolIndex))./sum(temp1);
            temp_falseAlarm = sum(temp2(temptempBoolIndex))./sum(temp2);
            
            hit_minus_falseAlarm = temp_hit - temp_falseAlarm;
            
            choiceMemoryCorrectError_hit_multi(temptempi) = temp_hit;
            choiceMemoryCorrectError_falseAlarm_multi(temptempi) = temp_falseAlarm;
            
            choiceMemoryCorrectError_hit_minus_falseAlarm_multi(temptempi) = hit_minus_falseAlarm;
        end
        
        [M,I] = max(choiceMemoryCorrectError_hit_minus_falseAlarm_multi);
        
        tempThreshold = temp_thresholdRange(I);
        
        choiceMemoryCorrectError_hit = choiceMemoryCorrectError_hit_multi(I);
        choiceMemoryCorrectError_falseAlarm = choiceMemoryCorrectError_falseAlarm_multi(I);
        
        lowThreshold_memoryPrecision = tempThreshold;
        highThreshold_memoryPrecision = lowThreshold_memoryPrecision;
    end
    
end






trialBoolIndex_memoryPrecisionLow_choice = false(trial_num,1);
trialBoolIndex_memoryPrecisionHigh_choice = false(trial_num,1);
for tempi=1:trial_num
    if choiceBoolIndex(tempi) == 0
        continue
    end
    if memoryPrecision_trialLevel(tempi) <= lowThreshold_memoryPrecision
        trialBoolIndex_memoryPrecisionLow_choice(tempi) = true;
    end
    if memoryPrecision_trialLevel(tempi) > highThreshold_memoryPrecision
        trialBoolIndex_memoryPrecisionHigh_choice(tempi) = true;
    end
end

% trialBoolIndex_memoryPrecisionLow_choiceMemory = trialBoolIndex_memoryPrecisionLow_choice & choiceMemoryBoolIndex';
% trialBoolIndex_memoryPrecisionLow_choiceOffload = trialBoolIndex_memoryPrecisionLow_choice & choiceOffloadBoolIndex';
% trialBoolIndex_memoryPrecisionHigh_choiceMemory = trialBoolIndex_memoryPrecisionHigh_choice & choiceMemoryBoolIndex';
% trialBoolIndex_memoryPrecisionHigh_choiceOffload = trialBoolIndex_memoryPrecisionHigh_choice & choiceOffloadBoolIndex';
% trialBoolIndex_memoryPrecisionLowError_choiceMemory = trialBoolIndex_memoryPrecisionLow_choice & choiceMemoryErrorBoolIndex';
% trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory = trialBoolIndex_memoryPrecisionHigh_choice & choiceMemoryCorrectBoolIndex';
%
% a1 = sum(trialBoolIndex_memoryPrecisionLow_choiceMemory);
% a2 = sum(trialBoolIndex_memoryPrecisionLow_choiceOffload);
% a3 = sum(trialBoolIndex_memoryPrecisionHigh_choiceMemory);
% a4 = sum(trialBoolIndex_memoryPrecisionHigh_choiceOffload);
% a5 = sum(trialBoolIndex_memoryPrecisionLowError_choiceMemory);
% a6 = sum(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory);

% match trials
trialBoolIndex_memoryPrecisionLow_choiceOffload = trialBoolIndex_memoryPrecisionLow_choice & choiceOffloadBoolIndex';
trialBoolIndex_memoryPrecisionHigh_choiceMemory = trialBoolIndex_memoryPrecisionHigh_choice & choiceMemoryBoolIndex';
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory = trialBoolIndex_memoryPrecisionHigh_choice & choiceMemoryCorrectBoolIndex';

% trialBoolIndex_memoryPrecisionLow_choiceOffloadLow
% trialBoolIndex_memoryPrecisionHigh_choiceMemoryHigh
% trialBoolIndex_memoryPrecisionHighCorrect_choiceMemoryHigh

% mismatch trials
trialBoolIndex_memoryPrecisionLow_choiceMemory = trialBoolIndex_memoryPrecisionLow_choice & choiceMemoryBoolIndex';
trialBoolIndex_memoryPrecisionLowError_choiceMemory = trialBoolIndex_memoryPrecisionLow_choice & choiceMemoryErrorBoolIndex';
trialBoolIndex_memoryPrecisionHigh_choiceOffload = trialBoolIndex_memoryPrecisionHigh_choice & choiceOffloadBoolIndex';

% trialBoolIndex_memoryPrecisionLow_choiceMemoryHigh
% trialBoolIndex_memoryPrecisionLowError_choiceMemoryHigh
% trialBoolIndex_memoryPrecisionHigh_choiceOffloadLow

a1 = sum(trialBoolIndex_memoryPrecisionLow_choiceOffload);
a2 = sum(trialBoolIndex_memoryPrecisionHigh_choiceMemory);
a3 = sum(trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory);
a4 = sum(trialBoolIndex_memoryPrecisionLow_choiceMemory);
a5 = sum(trialBoolIndex_memoryPrecisionLowError_choiceMemory);
a6 = sum(trialBoolIndex_memoryPrecisionHigh_choiceOffload);




memoryPrecision_trialLevel_choiceMemory = memoryPrecision_trialLevel(choiceMemoryBoolIndex);
memoryPrecision_trialLevel_choiceOffload = memoryPrecision_trialLevel(choiceOffloadBoolIndex);

memoryPrecision_trialLevel_choiceMemoryCorrect = memoryPrecision_trialLevel(choiceMemoryCorrectBoolIndex);
memoryPrecision_trialLevel_choiceMemoryError = memoryPrecision_trialLevel(choiceMemoryErrorBoolIndex);


% choiceMemoryBoolIndex_mildSeq = choiceMemoryBoolIndex;
% choiceOffloadBoolIndex_mildSeq = choiceOffloadBoolIndex;
memoryPrecision_trialLevel_mildSeq_choiceMemory = memoryPrecision_trialLevel(choiceMemoryBoolIndex_mildSeq);
memoryPrecision_trialLevel_mildSeq_choiceOffload = memoryPrecision_trialLevel(choiceOffloadBoolIndex_mildSeq);

aa1 = mean(memoryPrecision_trialLevel_choiceMemory,'omitnan');
aa2 = mean(memoryPrecision_trialLevel_choiceOffload,'omitnan');

aa3 = mean(memoryPrecision_trialLevel_choiceMemoryCorrect,'omitnan');
aa4 = mean(memoryPrecision_trialLevel_choiceMemoryError,'omitnan');

trialBoolIndex_memoryPrecisionLow_choiceOffload;
trialBoolIndex_memoryPrecisionHigh_choiceMemory;
trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory;
trialBoolIndex_memoryPrecisionLow_choiceMemory;
trialBoolIndex_memoryPrecisionLowError_choiceMemory;
trialBoolIndex_memoryPrecisionHigh_choiceOffload;


% memoryPrecision_trialLevel_mildSeq_allMemoryCorrect = memoryPrecision_trialLevel(allMemoryCorrectBoolIndex_mildSeq);
% memoryPrecision_trialLevel_mildSeq_allMemoryError = memoryPrecision_trialLevel(allMemoryErrorBoolIndex_mildSeq);

%% memoryMetaMismatch
% memoryMetaMismatch = struct;
if exist('memoryMetaMismatch','var') == 0
    memoryMetaMismatch = struct;
end

memoryMetaMismatch.cellIndex_suite2p_memoryMetaMismatch = cellIndex_suite2p_memoryMetaMismatch;
memoryMetaMismatch.cellIndex_suite2p_memoryMetaMatch = cellIndex_suite2p_memoryMetaMatch;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionLow_choiceOffload = trialBoolIndex_memoryPrecisionLow_choiceOffload;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionHigh_choiceMemory = trialBoolIndex_memoryPrecisionHigh_choiceMemory;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory = trialBoolIndex_memoryPrecisionHighCorrect_choiceMemory;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionLow_choiceMemory = trialBoolIndex_memoryPrecisionLow_choiceMemory;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionLowError_choiceMemory = trialBoolIndex_memoryPrecisionLowError_choiceMemory;
memoryMetaMismatch.trialBoolIndex_memoryPrecisionHigh_choiceOffload = trialBoolIndex_memoryPrecisionHigh_choiceOffload;

% if if_loadDemoData == 0
if true
    % if false
    cellIndex_suite2p_A = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1); %#ok<*UNRCH>
    % cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,2);
    cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
    
    % memoryPrecisionLow_choiceMemoryHigh
    cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh_topX = cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh_topX;
    cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryHigh_topX = nan(length(cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh_topX),1);
    for tempi=1:length(cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh_topX)
        temptempIndex = find(cellIndex_suite2p_A==cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh_topX(tempi));
        cellIndex_suite2p_A(temptempIndex);
        cellIndex_suite2p_B(temptempIndex);
        cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryHigh_topX(tempi) = cellIndex_suite2p_B(temptempIndex);
    end
    
    % memoryPrecisionHigh_choiceMemoryLow
    cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow_topX = cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow_topX;
    cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryLow_topX = nan(length(cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow_topX),1);
    for tempi=1:length(cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow_topX)
        temptempIndex = find(cellIndex_suite2p_A==cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow_topX(tempi));
        cellIndex_suite2p_A(temptempIndex);
        cellIndex_suite2p_B(temptempIndex);
        cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryLow_topX(tempi) = cellIndex_suite2p_B(temptempIndex);
    end
    
    
    % memoryPrecisionHigh_choiceMemoryHigh
    cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh_topX = cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh_topX;
    cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryHigh_topX = nan(length(cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh_topX),1);
    for tempi=1:length(cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh_topX)
        temptempIndex = find(cellIndex_suite2p_A==cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh_topX(tempi));
        cellIndex_suite2p_A(temptempIndex);
        cellIndex_suite2p_B(temptempIndex);
        cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryHigh_topX(tempi) = cellIndex_suite2p_B(temptempIndex);
    end
    
    memoryMetaMismatch.cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh_topX = cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh_topX;
    memoryMetaMismatch.cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow_topX = cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow_topX;
    memoryMetaMismatch.cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh_topX = cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh_topX;
    
    memoryMetaMismatch.cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryHigh_topX = cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryHigh_topX;
    memoryMetaMismatch.cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryLow_topX = cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryLow_topX;
    memoryMetaMismatch.cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryHigh_topX = cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryHigh_topX;
    
end

% end

%% Plot
if if_plot == 1
    close all
    
    
    %     %% Pdf of memory precision in choiceMemory trials and choiceOffload trials
    %     fig = figure('Name','memory precsion pdf','NumberTitle','off');
    %     set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %     x = memoryPrecision_trialLevel_choiceMemory;
    %     n=100;
    %     n=2^ceil(log2(n)); % round up n to the next power of 2;
    %     [pdf1,xmesh1,bandwidth1] = ksdensity(x','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
    %     temp_ratio1 = sum(choiceMemoryBoolIndex)/sum(choiceBoolIndex);
    %     %temp_ratio1 = 1;
    %     y1 = pdf1*temp_ratio1;
    %     h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',[0.4660 0.6740 0.1880]);
    %     hold on
    %
    %     x = memoryPrecision_trialLevel_choiceOffload;
    %     n=100;
    %     n=2^ceil(log2(n)); % round up n to the next power of 2;
    %     [pdf2,xmesh2,bandwidth2] = ksdensity(x','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
    %     temp_ratio2 = sum(choiceOffloadBoolIndex)/sum(choiceBoolIndex);
    %     %temp_ratio2 = 1;
    %     y2 = pdf2*temp_ratio2;
    %     h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',[0.6350 0.0780 0.1840]);
    %     hold on
    %
    %     [x1_min,x1_max] = bounds(xmesh1);
    %     [x2_min,x2_max] = bounds(xmesh2);
    %     x_min = min(x1_min,x2_min);
    %     x_max = max(x1_max,x2_max);
    %
    %
    %     [y1_min,y1_max] = bounds(y1);
    %     [y2_min,y2_max] = bounds(y2);
    %     y_min = min(y1_min,y2_min);
    %     y_max = max(y1_max,y2_max);
    %
    %     plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max],...
    %         'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
    %     hold on
    %
    %     if highThreshold_percent ~= lowThreshold_percent
    %         plot([highThreshold_memoryPrecision highThreshold_memoryPrecision],[y_min y_max],...
    %             'LineWidth',1,'color',[0.4940 0.1840 0.5560]);
    %         hold on
    %
    %         legend('choiceMemory','choiceOffload',...
    %             sprintf('%d percentile',lowThreshold_percent),...
    %             sprintf('%d percentile',highThreshold_percent),...
    %             'Location','northeast','fontsize',11)
    %     else
    %         legend('choiceMemory','choiceOffload',...
    %             sprintf('%d percentile',lowThreshold_percent),...
    %             'Location','northeast','fontsize',11)
    %     end
    %
    %
    %     xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    %     ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1
    %     set(gca, 'FontSize', 16)
    %     set(gca,'box','off');% 取消右、上边框
    %     xlabel('Memory precision', 'FontSize', 14, 'FontWeight', 'bold');
    %     ylabel('Pdf', 'FontSize', 14, 'FontWeight', 'bold');
    %     temp_title = title(sprintf('Memory precision of single trial, %s',FOVName_currentFOV2),...
    %         'FontSize',14,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    
    
    %     %% Pdf of memory precision in choiceMemoryCorrect trials and choiceMemoryError trials
    %     fig = figure('Name','memory precsion pdf 2','NumberTitle','off');
    %     set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %     x = memoryPrecision_trialLevel_choiceMemoryCorrect;
    %     n=100;
    %     n=2^ceil(log2(n)); % round up n to the next power of 2;
    %     [pdf1,xmesh1,bandwidth1] = ksdensity(x','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
    %     temp_ratio1 = sum(choiceMemoryCorrectBoolIndex)/sum(choiceMemoryBoolIndex);
    %     %temp_ratio1 = 1;
    %     y1 = pdf1*temp_ratio1;
    %     h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',[0.4660 0.6740 0.1880]);
    %     hold on
    %
    %     x = memoryPrecision_trialLevel_choiceMemoryError;
    %     n=100;
    %     n=2^ceil(log2(n)); % round up n to the next power of 2;
    %     [pdf2,xmesh2,bandwidth2] = ksdensity(x','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
    %     temp_ratio2 = sum(choiceMemoryErrorBoolIndex)/sum(choiceMemoryBoolIndex);
    %     %temp_ratio2 = 1;
    %     y2 = pdf2*temp_ratio2;
    %     h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',[0.6350 0.0780 0.1840]);
    %     hold on
    %
    %     [x1_min,x1_max] = bounds(xmesh1);
    %     [x2_min,x2_max] = bounds(xmesh2);
    %     x_min = min(x1_min,x2_min);
    %     x_max = max(x1_max,x2_max);
    %
    %
    %     [y1_min,y1_max] = bounds(y1);
    %     [y2_min,y2_max] = bounds(y2);
    %     y_min = min(y1_min,y2_min);
    %     y_max = max(y1_max,y2_max);
    %
    %     plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max],...
    %         'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
    %     hold on
    %
    %     if highThreshold_percent ~= lowThreshold_percent
    %         plot([highThreshold_memoryPrecision highThreshold_memoryPrecision],[y_min y_max],...
    %             'LineWidth',1,'color',[0.4940 0.1840 0.5560]);
    %         hold on
    %
    %         legend('choiceMemoryCorrect','choiceMemoryError',...
    %             sprintf('%d percentile',lowThreshold_percent),...
    %             sprintf('%d percentile',highThreshold_percent),...
    %             'Location','northeast','fontsize',11)
    %     else
    %         legend('choiceMemoryCorrect','choiceMemoryError',...
    %             sprintf('%d percentile',lowThreshold_percent),...
    %             'Location','northeast','fontsize',11)
    %     end
    %
    %
    %     xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    %     ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1
    %     set(gca, 'FontSize', 16)
    %     set(gca,'box','off');% 取消右、上边框
    %     xlabel('Memory precision', 'FontSize', 14, 'FontWeight', 'bold');
    %     ylabel('Pdf', 'FontSize', 14, 'FontWeight', 'bold');
    %     temp_title = title(sprintf('Memory precision of single trial, %s',FOVName_currentFOV2),...
    %         'FontSize',14,'FontWeight','bold');
    %     temp_title.Interpreter = 'none';
    
    
    %     %% Pdf of memory precision in choiceMemory_mildSeq trials and choiceOffload_mildSeq trials,
    %     fig = figure('Name','memory precsion pdf (mildSeq)','NumberTitle','off');
    %     %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %     set(gcf,'Position',[10 50 255 200*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %
    %     if_smoothHistogram = 0;
    %
    %     [~,temp_p] = ttest2(memoryPrecision_trialLevel_mildSeq_choiceMemory,...
    %         memoryPrecision_trialLevel_mildSeq_choiceOffload);
    %
    %     if if_smoothHistogram == 1
    %         x = memoryPrecision_trialLevel_mildSeq_choiceMemory;
    %         n=100;
    %         n=2^ceil(log2(n)); % round up n to the next power of 2;
    %         [pdf1,xmesh1,bandwidth1] = ksdensity(x','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
    %         %temp_ratio1 = sum(choiceMemoryBoolIndex)/sum(choiceBoolIndex);
    %         temp_ratio1 = length(memoryPrecision_trialLevel_mildSeq_choiceMemory)./...
    %             (length(memoryPrecision_trialLevel_mildSeq_choiceMemory)+length(memoryPrecision_trialLevel_mildSeq_choiceOffload));
    %         %temp_ratio1 = 1;
    %         y1 = pdf1*temp_ratio1;
    %         h1 = plot(xmesh1,y1,'LineWidth',1.5,'color',[0.4660 0.6740 0.1880]);
    %         hold on
    %
    %         x = memoryPrecision_trialLevel_mildSeq_choiceOffload;
    %         n=100;
    %         n=2^ceil(log2(n)); % round up n to the next power of 2;
    %         [pdf2,xmesh2,bandwidth2] = ksdensity(x','NumPoints',n,'Function','pdf','Support',[0 1],'BoundaryCorrection', 'Reflection');
    %         %temp_ratio2 = sum(choiceOffloadBoolIndex)/sum(choiceBoolIndex);
    %         temp_ratio2 = length(memoryPrecision_trialLevel_mildSeq_choiceOffload)./...
    %             (length(memoryPrecision_trialLevel_mildSeq_choiceMemory)+length(memoryPrecision_trialLevel_mildSeq_choiceOffload));
    %         %temp_ratio2 = 1;
    %         y2 = pdf2*temp_ratio2;
    %         h2 = plot(xmesh2,y2,'LineWidth',1.5,'color',[0.6350 0.0780 0.1840]);
    %         hold on
    %
    %         [x1_min,x1_max] = bounds(xmesh1);
    %         [x2_min,x2_max] = bounds(xmesh2);
    %
    %         x_min = min(x1_min,x2_min);
    %         x_max = max(x1_max,x2_max);
    %
    %     elseif if_smoothHistogram == 0
    %         h_NumBins = 8;%10
    %
    %         x = memoryPrecision_trialLevel_mildSeq_choiceMemory;
    %         h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
    %         hold on
    %         h1.NumBins = h_NumBins;
    %         h1.EdgeColor = color_choiceMemory;
    %
    %         x = memoryPrecision_trialLevel_mildSeq_choiceOffload;
    %         h2 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
    %         hold on
    %         h2.NumBins = h_NumBins;
    %         h2.EdgeColor = color_choiceOffload;
    %
    %
    %         y1 = h1.Values;
    %         y2 = h2.Values;
    %         x_min = 0;
    %         x_max = 1;
    %     end
    %
    %
    %
    %     [y1_min,y1_max] = bounds(y1);
    %     [y2_min,y2_max] = bounds(y2);
    %     y_min = min(y1_min,y2_min);
    %     y_max = max(y1_max,y2_max);
    %
    %     %plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max],...
    %     %    'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
    %     %hold on
    %
    %     if highThreshold_percent ~= lowThreshold_percent
    %         plot([highThreshold_memoryPrecision highThreshold_memoryPrecision],[y_min y_max],...
    %             'LineWidth',1,'color',[0.4940 0.1840 0.5560]);
    %         hold on
    %
    %         legend('choiceMemory','choiceOffload',...
    %             sprintf('%d percentile',lowThreshold_percent),...
    %             sprintf('%d percentile',highThreshold_percent),...
    %             'Location','northeast','fontsize',11)
    %     else
    %         %legend('choiceMemory','choiceOffload',...
    %         %    sprintf('%d percentile',lowThreshold_percent),...
    %         %    'Location','northeast','fontsize',11)
    %
    %         le = legend('ChoiceMemory','ChoiceOffload',...
    %             'Location','northeast','fontsize',8);
    %         le.ItemTokenSize = ones(1,2)*10;
    %
    %     end
    %
    %     %temp_p = 0;
    %
    %     tempTxt = sprintf('');
    %     if temp_p < 0.001
    %         tempTxt = sprintf('***');
    %     elseif temp_p < 0.01
    %         tempTxt = sprintf('**');
    %     elseif temp_p < 0.05
    %         tempTxt = sprintf('*');
    %     end
    %     text(0.4,y_min+(y_max-y_min)*0.7,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
    %         'HorizontalAlignment','center');
    %
    %
    %     set(gca,'linewidth',1.5)
    %     xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    %     ylim([y_min y_max+(y_max-y_min)*0.4]);%0.1
    %     %xticks([0 1]);
    %     set(gca, 'FontSize', 12)
    %     set(gca,'box','off');% 取消右、上边框
    %     xlabel('Memory precision', 'FontSize', 12, 'FontWeight', 'bold');
    %     if if_smoothHistogram == 1
    %         ylabel('Pdf', 'FontSize', 12, 'FontWeight', 'bold');
    %     elseif if_smoothHistogram == 0
    %         ylabel('Trial count', 'FontSize', 12, 'FontWeight', 'bold');
    %     end
    %
    %     temp_trialNum = length(memoryPrecision_trialLevel_mildSeq_choiceMemory)+length(memoryPrecision_trialLevel_mildSeq_choiceOffload);
    %
    %     if sum(mildSeqBoolIndex_seqLevel) > 1
    %         temp_title = title(sprintf('Memory precision of single trial (mildSeq), %s \n seqNum=%d (%.2f~%.2f), trialNum=%d',...
    %             FOVName_currentFOV2,sum(mildSeqBoolIndex_seqLevel),rProbLowThreshold_mildSeq,rProbHighThreshold_mildSeq,temp_trialNum),...
    %             'FontSize',12,'FontWeight','bold');
    %     elseif sum(mildSeqBoolIndex_seqLevel) == 1
    %         %temp_title = title(sprintf('%s \n Example seq: [%d], %d trials',FOVName_currentFOV2,...
    %         %   seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',12,'FontWeight','bold');
    %         temp_title = title(sprintf('%s \n Example seq: %d (%d trials)',FOVName_currentFOV2,...
    %             seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',12,'FontWeight','bold');
    %     end
    %     temp_title.Interpreter = 'none';
    
    
    
    %% Single neuron tuning property of memoryPrecision and choiceMemory (trial-level tuning)
    if false
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[10 50 372*0.9 290*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 50 372*0.9 290*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        x = r2_allFeatures3(:,2);
        y = r2_allFeatures3(:,3);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        
        
        %plot([threshold_r2_memoryPrecision threshold_r2_memoryPrecision],[y_min-(y_max-y_min)*0.02 y_max+(y_max-y_min)*0.02],...
        %    'LineWidth',1,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        %hold on
        %plot([x_min-(x_max-x_min)*0.02 x_max+(x_max-x_min)*0.02],[threshold_r2_choiceMemory threshold_r2_choiceMemory],...
        %    'LineWidth',1,'color',[0.25 0.25 0.25]);%[0.8500 0.3250 0.0980]
        %hold on
        plot([threshold_r2_memoryPrecision threshold_r2_memoryPrecision],[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
            'LineWidth',1,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        plot([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08],[threshold_r2_choiceMemory threshold_r2_choiceMemory],...
            'LineWidth',1,'color',[0.25 0.25 0.25]);%[0.8500 0.3250 0.0980]
        hold on
        
        % threshold_r2_memoryPrecision
        % threshold_r2_choiceMemory
        %scatter(x,y,20,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        scatter(x,y,10,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        %     scatter(x(~temptempBoolIndex_AC),y(~temptempBoolIndex_AC),...
        %         20,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        %     hold on
        
        
        
        temptempBoolIndex_A = x<=threshold_r2_memoryPrecision;
        temptempBoolIndex_B = x>threshold_r2_memoryPrecision;
        temptempBoolIndex_C = y<=threshold_r2_choiceMemory;
        temptempBoolIndex_D = y>threshold_r2_choiceMemory;
        num_xlow = sum(temptempBoolIndex_A);
        num_xhigh = sum(temptempBoolIndex_B);
        num_ylow = sum(temptempBoolIndex_C);
        num_yhigh = sum(temptempBoolIndex_D);
        
        temptempBoolIndex_AC = temptempBoolIndex_A & temptempBoolIndex_C;
        temptempBoolIndex_AD = temptempBoolIndex_A & temptempBoolIndex_D;
        temptempBoolIndex_BC = temptempBoolIndex_B & temptempBoolIndex_C;
        temptempBoolIndex_BD = temptempBoolIndex_B & temptempBoolIndex_D;
        num_xlow_ylow = sum(temptempBoolIndex_AC);
        num_xlow_yhigh = sum(temptempBoolIndex_AD);
        num_xhigh_ylow = sum(temptempBoolIndex_BC);
        num_xhigh_yhigh = sum(temptempBoolIndex_BD);
        
        x2 = x_max-(x_max-x_min)*0.1;
        y2 = y_max-(y_max-y_min)*0.1;
        text(x2-(x_max-x_min)*0.07,y2-(y_max-y_min)*0.07,sprintf('%d',num_xlow_ylow),'FontSize',12,'HorizontalAlignment','center');
        text(x2-(x_max-x_min)*0.07,y2+(y_max-y_min)*0.07,sprintf('%d',num_xlow_yhigh),'FontSize',12,'HorizontalAlignment','center');
        text(x2+(x_max-x_min)*0.07,y2-(y_max-y_min)*0.07,sprintf('%d',num_xhigh_ylow),'FontSize',12,'HorizontalAlignment','center');
        text(x2+(x_max-x_min)*0.07,y2+(y_max-y_min)*0.07,sprintf('%d',num_xhigh_yhigh),'FontSize',12,'HorizontalAlignment','center');
        
        plot([x2 x2],[y2-(y_max-y_min)*0.09 y2+(y_max-y_min)*0.09],'LineWidth',1,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        hold on
        plot([x2-(x_max-x_min)*0.09 x2+(x_max-x_min)*0.09],[y2 y2],'LineWidth',1,'color',[0.25 0.25 0.25]);%[0.8500 0.3250 0.0980]
        hold on
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('r2 of memory precision', 'FontSize', 12, 'FontWeight', 'bold');
        %ylabel('r2 of choice memory', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('r2 of meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title = title(sprintf('Single neuron tuning property, %s\n n=%d, all tuning to offloading rate (from %d)',FOVName_currentFOV2,length(x),roiNum),...
        %    'FontSize',12,'FontWeight','bold');
        %temp_title = title(sprintf('Single neuron tuning property, %s\n n=%d, all tuning to offloading rate (from %d)',FOVName_currentFOV2,length(x),roiNum),'FontSize',11);
        %temp_title = title(sprintf('%s (trial-level)\n Single neuron tuning property\n n=%d, all tuning to offloading rate (from %d)',FOVName_currentFOV2,length(x),roiNum),'FontSize',11);
        temp_title = title(sprintf('%s (trial-level)\n Single neuron tuning property, n=%d',FOVName_currentFOV2,length(x)),'FontSize',11);
        temp_title.Interpreter = 'none';
        
    end
    
    
    %% Single neuron tuning property of memoryPrecision and choiceMemory (seq-level tuning)
    %if true
    if false
        fig = figure('Name','tuning property (seq-level)','NumberTitle','off');
        %set(gcf,'Position',[400 50 372*0.9 290*0.97]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 469*0.95 379]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 445*0.8*1.02 379]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 445*0.8*1.02 379*0.7*0.9]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 445*0.8*1.02*0.965 379*0.7*0.9*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 445*0.8*1.02*0.965*0.88 379*0.7*0.9*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[50 50 445*0.8*1.02*0.965*0.88 379*0.7*0.9*0.92*0.91]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 445*0.8*1.02*0.965*0.88*1.04 379*0.7*0.9*0.92*0.91*1.47*0.95]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        x = r2_seqPrecision;
        y = r2_rProb;
        
        [temp_r, temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        
        %x_max = x_max*1.35;%1
        %y_max = y_max*1.35;
        
        %         plot([threshold_r2_seqPrecision threshold_r2_seqPrecision],[y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08],...
        %             'LineWidth',1,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        %         hold on
        %         plot([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08],[threshold_r2_rProb threshold_r2_rProb],...
        %             'LineWidth',1,'color',[0.25 0.25 0.25]);%[0.8500 0.3250 0.0980]
        %         hold on
        
        temptempBoolIndex_A = x<=threshold_r2_seqPrecision;
        temptempBoolIndex_B = x>threshold_r2_seqPrecision;
        temptempBoolIndex_C = y<=threshold_r2_rProb;
        temptempBoolIndex_D = y>threshold_r2_rProb;
        
        temptempBoolIndex_AC = temptempBoolIndex_A & temptempBoolIndex_C;% no selective
        temptempBoolIndex_AD = temptempBoolIndex_A & temptempBoolIndex_D;% pure meta
        temptempBoolIndex_BC = temptempBoolIndex_B & temptempBoolIndex_C;% pure precision
        temptempBoolIndex_BD = temptempBoolIndex_B & temptempBoolIndex_D;% mixed selective
        
        a = 1;
        
        %scatter(x,y,10,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        %hold on
        
        %scatter(x(~temptempBoolIndex_AC),y(~temptempBoolIndex_AC),10,...
        %    'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        %hold on
        
        temp_MarkerFaceAlpha = 1;%0.5
        
        scatter(x(temptempBoolIndex_AC),y(temptempBoolIndex_AC),10,...
            'filled','MarkerFaceColor',[1 1 1]*0.85,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        scatter(x(temptempBoolIndex_AD),y(temptempBoolIndex_AD),10,...
            'filled','MarkerFaceColor',[1 1 1]*0.5,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        scatter(x(temptempBoolIndex_BC),y(temptempBoolIndex_BC),10,...
            'filled','MarkerFaceColor',[1 1 1]*0.5,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        scatter(x(temptempBoolIndex_BD),y(temptempBoolIndex_BD),10,...
            'filled','MarkerFaceColor',[1 1 1]*0,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        
        
        
        num_xlow = sum(temptempBoolIndex_A);
        num_xhigh = sum(temptempBoolIndex_B);
        num_ylow = sum(temptempBoolIndex_C);
        num_yhigh = sum(temptempBoolIndex_D);
        
        temptempBoolIndex_AC = temptempBoolIndex_A & temptempBoolIndex_C;
        temptempBoolIndex_AD = temptempBoolIndex_A & temptempBoolIndex_D;
        temptempBoolIndex_BC = temptempBoolIndex_B & temptempBoolIndex_C;
        temptempBoolIndex_BD = temptempBoolIndex_B & temptempBoolIndex_D;
        num_xlow_ylow = sum(temptempBoolIndex_AC);
        num_xlow_yhigh = sum(temptempBoolIndex_AD);
        num_xhigh_ylow = sum(temptempBoolIndex_BC);
        num_xhigh_yhigh = sum(temptempBoolIndex_BD);
        
        %x2 = x_max-(x_max-x_min)*0.1;%0.1
        x2 = threshold_r2_seqPrecision+(x_max-x_min)*0.18;
        y2 = y_max-(y_max-y_min)*0.1;
        %         text(x2-(x_max-x_min)*0.07,y2-(y_max-y_min)*0.07,sprintf('%d',num_xlow_ylow),'FontSize',10,'HorizontalAlignment','center');
        %         text(x2-(x_max-x_min)*0.07,y2+(y_max-y_min)*0.07,sprintf('%d',num_xlow_yhigh),'FontSize',10,'HorizontalAlignment','center');
        %         text(x2+(x_max-x_min)*0.07,y2-(y_max-y_min)*0.07,sprintf('%d',num_xhigh_ylow),'FontSize',10,'HorizontalAlignment','center');
        %         text(x2+(x_max-x_min)*0.07,y2+(y_max-y_min)*0.07,sprintf('%d',num_xhigh_yhigh),'FontSize',10,'HorizontalAlignment','center');
        %
        %         plot([x2 x2],[y2-(y_max-y_min)*0.09 y2+(y_max-y_min)*0.09],'LineWidth',1,'color',[0.25 0.25 0.25]);%[0 0.4470 0.7410]
        %         hold on
        %         plot([x2-(x_max-x_min)*0.09 x2+(x_max-x_min)*0.09],[y2 y2],'LineWidth',1,'color',[0.25 0.25 0.25]);%[0.8500 0.3250 0.0980]
        %         hold on
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        %text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.98,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        %text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.85,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        text(x_min+(x_max-x_min)*0.03,y_min+(y_max-y_min)*0.98,sprintf('r = %.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        text(x_min+(x_max-x_min)*0.03,y_min+(y_max-y_min)*0.88,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        
        
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            %yticks([0 0.6]);
            xticks([0:0.1:0.3]);
        end
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        %xlabel('r2 of memory precision', 'FontSize', 10, 'FontWeight', 'normal');
        %ylabel('r2 of meta-memory', 'FontSize', 10, 'FontWeight', 'normal');
        xlabel('r2 (Memory precision)', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('r2 (Meta-memory)', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('%s (seq-level)\n Single neuron tuning property, n=%d',FOVName_currentFOV2,length(x)),'FontSize',9);
        %temp_title = title(sprintf('%s (seq-level)\n Regression explained variance of neurons, n=%d',FOVName_currentFOV2,length(x)),'FontSize',10);
        temp_title = title(sprintf('%s (seq-level)\n Regression explained variance, n = %d',FOVName_currentFOV2,length(x)),'FontSize',9);
        temp_title.Interpreter = 'none';
    end
    
    
    %% Single neuron tuning property of memory (location) and choiceMemory (seq-level tuning)
    if true
    %if if_loadDemoData == 0
        fig = figure('Name','tuning property (seq-level)','NumberTitle','off');
        %set(gcf,'Position',[400 50 445*0.8*1.02*0.965*0.88 379*0.7*0.9*0.92*0.91]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[400 50 308.4 200*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[400 50 308.4*0.85 200*0.92*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        x = r2_6loc;
        %x = temp_glm_r2;
        %x(x<-0.05) = -0.05;
        
        y = r2_rProb;
        
        [temp_r, temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        
        
        temptempBoolIndex_A = ~boolIndex_6loc;
        temptempBoolIndex_B = boolIndex_6loc;
        %temptempBoolIndex_A = ~temptempBoolIndex_6loc;
        %temptempBoolIndex_B = temptempBoolIndex_6loc;
        %temptempBoolIndex_A = p_anova_6loc>0.05;
        %temptempBoolIndex_B = p_anova_6loc<=0.05;
        %temptempBoolIndex_C = y<=threshold_r2_rProb;
        %temptempBoolIndex_D = y>threshold_r2_rProb;
        temptempBoolIndex_C = p_rProb>threshold_p;
        temptempBoolIndex_D =  p_rProb<=threshold_p;
        %temptempBoolIndex_C = ~temptempBoolIndex_rProb;
        %temptempBoolIndex_D = temptempBoolIndex_rProb;
        
        
        temptempBoolIndex_AC = temptempBoolIndex_A & temptempBoolIndex_C;% no selective
        temptempBoolIndex_AD = temptempBoolIndex_A & temptempBoolIndex_D;% pure meta
        temptempBoolIndex_BC = temptempBoolIndex_B & temptempBoolIndex_C;% pure location
        temptempBoolIndex_BD = temptempBoolIndex_B & temptempBoolIndex_D;% mixed selective
        
        
        temp_MarkerFaceAlpha = 1;%0.5
        temp_LineWidth = 0.75;%0.5
        
        scatter(x(temptempBoolIndex_AC),y(temptempBoolIndex_AC),10,...
            'filled','MarkerFaceColor',[1 1 1]*0.85,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        hold on
        
        scatter(x(temptempBoolIndex_AD),y(temptempBoolIndex_AD),10,...
            'filled','MarkerFaceColor',[1 1 1]*0.5,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth);
        hold on
        
        %scatter(x(temptempBoolIndex_BC),y(temptempBoolIndex_BC),10,...
        %    'filled','MarkerFaceColor',[1 1 1]*0.5,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7);
        scatter(x(temptempBoolIndex_BC),y(temptempBoolIndex_BC),10,...
            'filled','MarkerFaceColor',[1 1 1]*1,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth,'MarkerEdgeColor',[1 1 1]*0);
        hold on
        
        %scatter(x(temptempBoolIndex_BD),y(temptempBoolIndex_BD),10,...
        %    'filled','MarkerFaceColor',[1 1 1]*0,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7);
        scatter(x(temptempBoolIndex_BD),y(temptempBoolIndex_BD),10,...
            'filled','MarkerFaceColor',[1 1 1]*0,'MarkerFaceAlpha', temp_MarkerFaceAlpha, 'MarkerEdgeAlpha', 0.7,'LineWidth',temp_LineWidth,'MarkerEdgeColor',[1 1 1]*0);
        hold on
        
        
        tempTxt = 'p > 0.05';
        if temp_p < 0.05
            tempTxt = 'p < 0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p < 0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p < 0.001';
        end
        text(x_min+(x_max-x_min)*0.03,y_min+(y_max-y_min)*0.98,sprintf('r = %.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        text(x_min+(x_max-x_min)*0.03,y_min+(y_max-y_min)*0.88,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        
        
        if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
            %yticks([0 0.6]);
            %xticks([0:0.1:0.3]);
        end
        
        set(gca,'linewidth',1.5)
        xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        set(gca, 'FontSize', 8)
        set(gca,'box','off');% 取消右、上边框
        xlabel('r2 (Memory)', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('r2 (Meta-memory)', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('%s (seq-level)\n Regression explained variance, n = %d',FOVName_currentFOV2,length(x)),'FontSize',9);
        %temp_title.Interpreter = 'none';
        title(sprintf('Regression explained variance, n = %d',length(x)),'FontSize',9);
        
    end
    
    
    %if if_plot_memoryPrecision_trialLevelEvidence == 1
    %if if_plot_memoryPrecision_trialLevelEvidence == 1 && if_memeoryPrecision_stimuli0_response1 == 1
    
    % old, for sequence-based memory strength
    %if if_plot_memoryPrecision_trialLevelEvidence == 1 && if_memeoryPrecision_stimuli0_response1 == 1 ...
    %        || if_plot_memoryPrecision_trialLevelEvidence_timeCourse == 1
    
    
    % new, for entropy-based memory strength
    if if_plot_memoryPrecision_trialLevelEvidence == 1 && if_memeoryPrecision_stimuli0_response1 == 0 ...
            || if_plot_memoryPrecision_trialLevelEvidence_timeCourse == 1
        
        %% seqAccuracy_memoryPrecisionX_allMemory
        memoryPrecision_trialLevel_stimuli;
        memoryPrecision_trialLevel_response;
        
        memoryPrecision_trialLevel_mildSeq_allMemoryCorrect = memoryPrecision_trialLevel(allMemoryCorrectBoolIndex_mildSeq);
        memoryPrecision_trialLevel_mildSeq_allMemoryStimuliError = memoryPrecision_trialLevel_stimuli(allMemoryErrorBoolIndex_mildSeq);
        memoryPrecision_trialLevel_mildSeq_allMemoryResponseError = memoryPrecision_trialLevel_response(allMemoryErrorBoolIndex_mildSeq);
        
        
        memoryPrecision_trialLevel;
        allMemoryCorrectBoolIndex;
        allMemoryErrorBoolIndex;
        seqIndex;
        
        memoryPrecision_trialLevel_allSeq_allMemoryCorrectMean = nan(sum(numSeq(1:valid_length)),1);
        memoryPrecision_trialLevel_allSeq_allMemoryStimuliErrorMean = nan(sum(numSeq(1:valid_length)),1);
        memoryPrecision_trialLevel_allSeq_allMemoryResponseErrorMean = nan(sum(numSeq(1:valid_length)),1);
        
        for tempi=1:sum(numSeq(1:valid_length))
            tempSeqBoolIndex = seqIndex==tempi;
            tempSeqBoolIndex_response = seqIndex_response==tempi;
            
            memoryPrecision_trialLevel_allSeq_allMemoryCorrectMean(tempi) = ...
                mean(memoryPrecision_trialLevel(allMemoryCorrectBoolIndex&tempSeqBoolIndex));
            
            memoryPrecision_trialLevel_allSeq_allMemoryStimuliErrorMean(tempi) = ...
                mean(memoryPrecision_trialLevel_stimuli(allMemoryErrorBoolIndex&tempSeqBoolIndex));
            
            memoryPrecision_trialLevel_allSeq_allMemoryResponseErrorMean(tempi) = ...
                mean(memoryPrecision_trialLevel_response(allMemoryErrorBoolIndex&tempSeqBoolIndex_response));
        end
        d1 = memoryPrecision_trialLevel_allSeq_allMemoryCorrectMean;
        d2 = memoryPrecision_trialLevel_allSeq_allMemoryStimuliErrorMean;
        d3 = memoryPrecision_trialLevel_allSeq_allMemoryResponseErrorMean;
        
        [~,temp_p_correctErrorPrecision_stimuli] = ttest(d1,d2);
        [~,temp_p_correctErrorPrecision_response] = ttest(d1,d3);
        
        
        
        lowThreshold_memoryPrecision;
        
        
        allMemoryCorrectBoolIndex;
        allMemoryErrorBoolIndex;
        
        trialBoolIndex_memoryPrecisionLow_allMemoryCorrect = (memoryPrecision_trialLevel<=lowThreshold_memoryPrecision) & allMemoryCorrectBoolIndex';
        trialBoolIndex_memoryPrecisionHigh_allMemoryCorrect = (memoryPrecision_trialLevel>lowThreshold_memoryPrecision) & allMemoryCorrectBoolIndex';
        
        trialBoolIndex_memoryPrecisionLow_allMemoryStimuliError = (memoryPrecision_trialLevel_stimuli<=lowThreshold_memoryPrecision) & allMemoryErrorBoolIndex';
        trialBoolIndex_memoryPrecisionHigh_allMemoryStimuliError = (memoryPrecision_trialLevel_stimuli>lowThreshold_memoryPrecision) & allMemoryErrorBoolIndex';
        
        trialBoolIndex_memoryPrecisionLow_allMemoryResponseError = (memoryPrecision_trialLevel_response<=lowThreshold_memoryPrecision) & allMemoryErrorBoolIndex';
        trialBoolIndex_memoryPrecisionHigh_allMemoryResponseError = (memoryPrecision_trialLevel_response>lowThreshold_memoryPrecision) & allMemoryErrorBoolIndex';
        
        
        seqIndex;
        seqIndex_response;
        
        seqAccuracy_memoryPrecisionLow_allMemoryStimuli = nan(sum(numSeq(1:valid_length)),1);
        seqAccuracy_memoryPrecisionHigh_allMemoryStimuli = nan(sum(numSeq(1:valid_length)),1);
        
        seqAccuracy_memoryPrecisionLow_allMemoryResponse = nan(sum(numSeq(1:valid_length)),1);
        seqAccuracy_memoryPrecisionHigh_allMemoryResponse = nan(sum(numSeq(1:valid_length)),1);
        
        for tempi=1:sum(numSeq(1:valid_length))
            % Stimuli-labeled trial
            tempTrialBoolIndex_targetSeq = seqIndex'==tempi;
            
            tempTrialBoolIndex_targetSeq_low_allMemoryCorrect = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_allMemoryCorrect;
            tempTrialBoolIndex_targetSeq_low_allMemoryError = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_allMemoryStimuliError;
            
            tempAcc_low = sum(tempTrialBoolIndex_targetSeq_low_allMemoryCorrect)/...
                (sum(tempTrialBoolIndex_targetSeq_low_allMemoryCorrect)+sum(tempTrialBoolIndex_targetSeq_low_allMemoryError));
            seqAccuracy_memoryPrecisionLow_allMemoryStimuli(tempi) = tempAcc_low;
            
            tempTrialBoolIndex_targetSeq_high_allMemoryCorrect = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_allMemoryCorrect;
            tempTrialBoolIndex_targetSeq_high_allMemoryError = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_allMemoryStimuliError;
            
            tempAcc_high = sum(tempTrialBoolIndex_targetSeq_high_allMemoryCorrect)/...
                (sum(tempTrialBoolIndex_targetSeq_high_allMemoryCorrect)+sum(tempTrialBoolIndex_targetSeq_high_allMemoryError));
            seqAccuracy_memoryPrecisionHigh_allMemoryStimuli(tempi) = tempAcc_high;
            
            
            % Response-labeled trial
            tempTrialBoolIndex_targetSeq = seqIndex_response'==tempi;
            
            tempTrialBoolIndex_targetSeq_low_allMemoryCorrect = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_allMemoryCorrect;
            tempTrialBoolIndex_targetSeq_low_allMemoryError = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_allMemoryResponseError;
            
            tempAcc_low = sum(tempTrialBoolIndex_targetSeq_low_allMemoryCorrect)/...
                (sum(tempTrialBoolIndex_targetSeq_low_allMemoryCorrect)+sum(tempTrialBoolIndex_targetSeq_low_allMemoryError));
            seqAccuracy_memoryPrecisionLow_allMemoryResponse(tempi) = tempAcc_low;
            
            tempTrialBoolIndex_targetSeq_high_allMemoryCorrect = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_allMemoryCorrect;
            tempTrialBoolIndex_targetSeq_high_allMemoryError = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_allMemoryResponseError;
            
            tempAcc_high = sum(tempTrialBoolIndex_targetSeq_high_allMemoryCorrect)/...
                (sum(tempTrialBoolIndex_targetSeq_high_allMemoryCorrect)+sum(tempTrialBoolIndex_targetSeq_high_allMemoryError));
            seqAccuracy_memoryPrecisionHigh_allMemoryResponse(tempi) = tempAcc_high;
            
            
        end
        
        b1 = seqAccuracy_memoryPrecisionLow_allMemoryStimuli;
        b2 = seqAccuracy_memoryPrecisionHigh_allMemoryStimuli;
        
        b3 = seqAccuracy_memoryPrecisionLow_allMemoryResponse;
        b4 = seqAccuracy_memoryPrecisionHigh_allMemoryResponse;
        
        b1_mean = mean(b1,'omitnan');
        b2_mean = mean(b2,'omitnan');
        b3_mean = mean(b3,'omitnan');
        b4_mean = mean(b4,'omitnan');
        
        [~,temp_p_lowHigh_stimuli] = ttest(b1,b2);
        [~,temp_p_lowHigh_response] = ttest(b3,b4);
        
        
        
        
        %% seqRProb_memoryPrecisionX_choice
        memoryPrecision_trialLevel_stimuli;
        memoryPrecision_trialLevel_response;
        
        memoryPrecision_trialLevel_mildSeq_choiceOffload = memoryPrecision_trialLevel(choiceOffloadBoolIndex_mildSeq);
        memoryPrecision_trialLevel_mildSeq_choiceMemoryStimuli = memoryPrecision_trialLevel_stimuli(choiceMemoryBoolIndex_mildSeq);
        memoryPrecision_trialLevel_mildSeq_choiceMemoryResponse = memoryPrecision_trialLevel_response(choiceMemoryBoolIndex_mildSeq);
        
        
        memoryPrecision_trialLevel_allSeq_choiceOffloadMean = nan(sum(numSeq(1:valid_length)),1);
        memoryPrecision_trialLevel_allSeq_choiceMemoryStimuliMean = nan(sum(numSeq(1:valid_length)),1);
        memoryPrecision_trialLevel_allSeq_choiceMemoryResponseMean = nan(sum(numSeq(1:valid_length)),1);
        
        for tempi=1:sum(numSeq(1:valid_length))
            tempSeqBoolIndex = seqIndex==tempi;
            tempSeqBoolIndex_response = seqIndex_response==tempi;
            
            memoryPrecision_trialLevel_allSeq_choiceOffloadMean(tempi) = ...
                mean(memoryPrecision_trialLevel(choiceOffloadBoolIndex&tempSeqBoolIndex));
            
            memoryPrecision_trialLevel_allSeq_choiceMemoryStimuliMean(tempi) = ...
                mean(memoryPrecision_trialLevel_stimuli(choiceMemoryBoolIndex&tempSeqBoolIndex));
            
            memoryPrecision_trialLevel_allSeq_choiceMemoryResponseMean(tempi) = ...
                mean(memoryPrecision_trialLevel_response(choiceMemoryBoolIndex&tempSeqBoolIndex_response));
        end
        d4 = memoryPrecision_trialLevel_allSeq_choiceOffloadMean;
        d5 = memoryPrecision_trialLevel_allSeq_choiceMemoryStimuliMean;
        d6 = memoryPrecision_trialLevel_allSeq_choiceMemoryResponseMean;
        
        [~,temp_p_choicePrecision_stimuli] = ttest(d4,d5);
        [~,temp_p_choicePrecision_response] = ttest(d4,d6);
        
        
        
        lowThreshold_memoryPrecision;
        
        trialBoolIndex_memoryPrecisionLow_choiceOffload;
        trialBoolIndex_memoryPrecisionHigh_choiceOffload;
        
        
        trialBoolIndex_memoryPrecisionLow_choiceMemoryStimuli = (memoryPrecision_trialLevel_stimuli<=lowThreshold_memoryPrecision) & choiceMemoryBoolIndex';
        trialBoolIndex_memoryPrecisionHigh_choiceMemoryStimuli = (memoryPrecision_trialLevel_stimuli>lowThreshold_memoryPrecision) & choiceMemoryBoolIndex';
        
        trialBoolIndex_memoryPrecisionLow_choiceMemoryResponse = (memoryPrecision_trialLevel_response<=lowThreshold_memoryPrecision) & choiceMemoryBoolIndex';
        trialBoolIndex_memoryPrecisionHigh_choiceMemoryResponse = (memoryPrecision_trialLevel_response>lowThreshold_memoryPrecision) & choiceMemoryBoolIndex';
        
        
        seqIndex;
        seqIndex_response;
        
        seqRProb_memoryPrecisionLow_choiceStimuli = nan(sum(numSeq(1:valid_length)),1);
        seqRProb_memoryPrecisionHigh_choiceStimuli = nan(sum(numSeq(1:valid_length)),1);
        
        seqRProb_memoryPrecisionLow_choiceResponse = nan(sum(numSeq(1:valid_length)),1);
        seqRProb_memoryPrecisionHigh_choiceResponse = nan(sum(numSeq(1:valid_length)),1);
        
        for tempi=1:sum(numSeq(1:valid_length))
            % Stimuli-labeled trial
            tempTrialBoolIndex_targetSeq = seqIndex'==tempi;
            
            tempTrialBoolIndex_targetSeq_low_choiceOffload = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_choiceOffload;
            tempTrialBoolIndex_targetSeq_low_choiceMemory = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_choiceMemoryStimuli;
            
            tempRProb_low = sum(tempTrialBoolIndex_targetSeq_low_choiceOffload)/...
                (sum(tempTrialBoolIndex_targetSeq_low_choiceOffload)+sum(tempTrialBoolIndex_targetSeq_low_choiceMemory));
            seqRProb_memoryPrecisionLow_choiceStimuli(tempi) = tempRProb_low;
            
            tempTrialBoolIndex_targetSeq_high_choiceOffload = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_choiceOffload;
            tempTrialBoolIndex_targetSeq_high_choiceMemory = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_choiceMemoryStimuli;
            
            tempRProb_high = sum(tempTrialBoolIndex_targetSeq_high_choiceOffload)/...
                (sum(tempTrialBoolIndex_targetSeq_high_choiceOffload)+sum(tempTrialBoolIndex_targetSeq_high_choiceMemory));
            seqRProb_memoryPrecisionHigh_choiceStimuli(tempi) = tempRProb_high;
            
            
            % Response-labeled trial
            tempTrialBoolIndex_targetSeq = seqIndex_response'==tempi;
            
            tempTrialBoolIndex_targetSeq_low_choiceOffload = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_choiceOffload;
            tempTrialBoolIndex_targetSeq_low_choiceMemory = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionLow_choiceMemoryResponse;
            
            tempRProb_low = sum(tempTrialBoolIndex_targetSeq_low_choiceOffload)/...
                (sum(tempTrialBoolIndex_targetSeq_low_choiceOffload)+sum(tempTrialBoolIndex_targetSeq_low_choiceMemory));
            seqRProb_memoryPrecisionLow_choiceResponse(tempi) = tempRProb_low;
            
            tempTrialBoolIndex_targetSeq_high_choiceOffload = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_choiceOffload;
            tempTrialBoolIndex_targetSeq_high_choiceMemory = ...
                tempTrialBoolIndex_targetSeq & trialBoolIndex_memoryPrecisionHigh_choiceMemoryResponse;
            
            tempRProb_high = sum(tempTrialBoolIndex_targetSeq_high_choiceOffload)/...
                (sum(tempTrialBoolIndex_targetSeq_high_choiceOffload)+sum(tempTrialBoolIndex_targetSeq_high_choiceMemory));
            seqRProb_memoryPrecisionHigh_choiceResponse(tempi) = tempRProb_high;
            
            
        end
        
        c1 = seqRProb_memoryPrecisionLow_choiceStimuli;
        c2 = seqRProb_memoryPrecisionHigh_choiceStimuli;
        
        c3 = seqRProb_memoryPrecisionLow_choiceResponse;
        c4 = seqRProb_memoryPrecisionHigh_choiceResponse;
        
        c1_mean = mean(c1,'omitnan');
        c2_mean = mean(c2,'omitnan');
        c3_mean = mean(c3,'omitnan');
        c4_mean = mean(c4,'omitnan');
        
        [~,temp_p_c_lowHigh_stimuli] = ttest(c1,c2);
        [~,temp_p_c_lowHigh_response] = ttest(c3,c4);
        
        
        
        if if_plot_memoryPrecision_trialLevelEvidence_timeCourse == 1
            
            memoryPrecision_trialLevel_mildSeq_allMemoryCorrect;
            memoryPrecision_trialLevel_mildSeq_allMemoryStimuliError;
            memoryPrecision_trialLevel_crossTime_delay1;
            
            %memoryPrecision_trialLevel_mildSeq_allMemoryCorrect = memoryPrecision_trialLevel(allMemoryCorrectBoolIndex_mildSeq);
            %memoryPrecision_trialLevel_mildSeq_allMemoryStimuliError = memoryPrecision_trialLevel_stimuli(allMemoryErrorBoolIndex_mildSeq);
            
            memoryPrecision_trialLevel_time_mildSeq_allMemoryCorrect = memoryPrecision_trialLevel_crossTime_delay1(allMemoryCorrectBoolIndex_mildSeq,:);
            memoryPrecision_trialLevel_time_mildSeq_allMemoryStimuliError = memoryPrecision_trialLevel_crossTime_delay1(allMemoryErrorBoolIndex_mildSeq,:);
            
            x1_mean = mean(memoryPrecision_trialLevel_time_mildSeq_allMemoryCorrect,'all','omitnan');
            x2_mean = mean(memoryPrecision_trialLevel_time_mildSeq_allMemoryStimuliError,'all');
            
            
            %% memoryPrecision_trialLevelEvidence in stimuli-labeled trials, timeCourse
            fig = figure('Name','memory precsion mildSeq','NumberTitle','off');
            %set(gcf,'Position',[35+0 42+0 222*1.08*1.1 260*0.85/2]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[35+0 42+0 289 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[35+0 42+0 289*0.7 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[35+0 42+0 289*0.7*0.78 200*0.78]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[35+0 42+0 289*0.7*0.78*0.94 200*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[35+0 42+0 289*0.7*0.78*0.94*1.67*1.04*1.015 200*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[35+0 42+0 289*0.7*0.78*0.94*1.07*1.04 200*0.78*0.94*1.1]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[35+0 42+0 289*0.7*0.78*0.94*1.07*1.04*0.84 200*0.78*0.94*1.1*0.98]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            %t = tiledlayout(1,10,'TileSpacing','compact','Padding','loose'); %#ok<*NASGU>
            t = tiledlayout(1,9,'TileSpacing','compact','Padding','loose'); %#ok<*NASGU>
            
            %% Plot choice tuning
            
            nexttile([1 7])
            
            x1 = memoryPrecision_trialLevel_time_mildSeq_allMemoryCorrect;
            x2 = memoryPrecision_trialLevel_time_mildSeq_allMemoryStimuliError;
            
            [M,I] = sort(mean(x1,2),'descend');
            x1 = x1(I,:);
            
            [M,I] = sort(mean(x2,2),'descend');
            x2 = x2(I,:);
            
            x12 = [x1;x2];
            
            %cmin = min(x12,[],'all');
            %cmax = max(x12,[],'all');
            cmin = 0.5;%0
            cmax = 0.9;
            
            temp_ticks_decisionPeriodA = decisionPeriodA_interval;
            
            x = [temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)];
            y = [1,size(x12,1)];
            C = x12;
            imagesc(x,y,C,[cmin cmax]);
            hold on
            
            for tempi=1:length(decisionPeriodA_interval)
                plot(temp_ticks_decisionPeriodA(tempi)*[1 1],[1 size(x12,1)],...
                    '-','LineWidth',1,'Color',[1 1 1]);
                hold on
            end
            
            plot([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)],[1 1]*size(x1,1),...
                '-','LineWidth',1,'Color',[1 1 1]);
            hold on
            
            backgounrdColor = [1 1 1]*1;
            
            set(gca,'linewidth',1.5)
            set(gca,'color',backgounrdColor);
            xlim([temp_ticks_decisionPeriodA(1) temp_ticks_decisionPeriodA(end)]);
            ylim([1-1 size(x12,1)+1]);
            set(gca,'xticklabel',[])
            set(gca,'xtick',[])
            
            xticks(temp_ticks_decisionPeriodA);
            %xticklabels({'Delay1','ChoiceCue',''});
            %xticklabels({'Delay1','Decision',''});
            xticklabels({'Delay-on','Decision',''});
            
            xtickangle(0);
            set(gca, 'FontSize', 8);%10
            set(gca,'box','off');% 取消右、上边框
            
            
            tempLabelStr = string(1:2);
            tempLabelStr(1) = 'Correct';
            tempLabelStr(2) = 'Error';
            ytl=string(tempLabelStr);
            % 设置ttext的x坐标位置
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-2.4;%-2.4
            %ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-6.5;%-2.4
            ytext_xp=(xt(1))*ones(1,length(tempLabelStr))-5.5;%-3.5
            % 设置ttext的y坐标位置
            % ytext_yp=yt;
            ytext_yp = nan(1,2);
            ytext_yp(1) = size(x1,1)/2 - 0;
            ytext_yp(2) = size(x1,1) + size(x2,1)/2 - 0;
            ytext_yp = ytext_yp + 0;%0.5
            
            text(ytext_xp,ytext_yp,ytl,'HorizontalAlignment','center','rotation',90,'fontsize',8);%9
            set(gca,'yticklabel','');
            
            
            %ylabel(sprintf('Trials'), 'FontSize', 12, 'FontWeight', 'bold');
            %ylabel(sprintf('Trials'),'Position',[min(temp_ticks_decisionPeriodA)-3.5 size(x12,1)/2], 'FontSize', 10, 'FontWeight', 'normal');
            %ylabel(sprintf('Trials'),'Position',[min(temp_ticks_decisionPeriodA)-7.5 size(x12,1)/2], 'FontSize', 11, 'FontWeight', 'normal');
            
            %c = colorbar('fontsize',9);
            c = colorbar('eastoutside','FontSize',6.5);%9
            %c.Position = c.Position+[0 0 -0.015 -0.18];
            %c.Position = c.Position+[0 0 -0.0375 -0.23];
            %c.Position = c.Position+[0 0 -0.0305 -0.24];
            %c.Position = c.Position+[0.035 0 -0.0455 -0.24];
            c.Position = c.Position+[0.035 0 -0.060 -0.24];
            %c.Ticks = [0 1];
            c.Ticks = [cmin cmax];
            
            %             if if_colormap_loadEnhanced == 1
            %                 load('parula_enhanced');
            %                 colormap(parula_enhanced);
            %             elseif if_colormap_loadEnhanced == 0
            %                 colormap parula
            %             end
            
            %temp1 = gray;
            %temp1 = temp1(end:-1:1,:);
            %colormap(temp1);
            
            colormap(coolwarm());
            
            %temp1 = coolwarm(300);
            %temp2 = ((300-256)/2)+1;
            %temp3 = temp1(temp2:temp2+255,:);
            %colormap(temp3);
            
            
            
            temp_trialNum = size(x12,1);
            
            %temp_title = title(sprintf('Stimuli-labeled trial\nSeq %d (%d trials)',...
            %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
            %temp_title = title(sprintf('Precision time course\nSeq %d (%d trials)',...
            %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',9,'FontWeight','normal');
            %temp_title.Interpreter = 'none';
            
            %title(sprintf('Precision time course'),'FontSize',9,'FontWeight','normal');
            %title(sprintf('Memory strength time course'),'FontSize',9,'FontWeight','normal');
            title(sprintf('Quality time course'),'FontSize',9,'FontWeight','normal');
            %title(sprintf('Location correlation time course'),'FontSize',9,'FontWeight','normal');
            subtitle(sprintf('Seq %d (%d trials)',seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',7,'FontWeight','normal');
            
            
            a = 1;
        end
        
        %% memoryPrecision_trialLevelEvidence in stimuli-labeled trials
        fig = figure('Name','memory precsion pdf (mildSeq)','NumberTitle','off');
        
        
        if if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 0
            %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[10 450 355 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            t = tiledlayout(1,2,'TileSpacing','compact','Padding','loose');
            
            
            %nexttile([1 2])
            nexttile
            
            x1 = memoryPrecision_trialLevel_mildSeq_allMemoryCorrect;
            x2 = memoryPrecision_trialLevel_mildSeq_allMemoryStimuliError;
            %x3 = memoryPrecision_trialLevel_mildSeq_allMemoryResponseError;
            
            
            [~,temp_p] = ttest2(x1,x2);
            
            
            h_NumBins = 8;%10
            
            x = x1;
            h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h1.NumBins = h_NumBins;
            h1.EdgeColor = [0 0 0];
            
            x = x2;
            h2 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h2.NumBins = h_NumBins;
            h2.EdgeColor = [0.5 0.5 0.5];
            
            
            y1 = h1.Values;
            y2 = h2.Values;
            
            
            x_min = 0;
            x_max = 1;
            
            [y1_min,y1_max] = bounds(y1);
            [y2_min,y2_max] = bounds(y2);
            y_min = min([y1_min,y2_min]);
            y_max = max([y1_max,y2_max]);
            
            
            if if_plot_additionalSmooth == 1
                n=100;
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                
                [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
                plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',[0 0 0]);
                hold on
                
                [pdf2,xmesh2,~] = ksdensity(x2,'NumPoints',n,'Function','pdf');
                plot(xmesh2,pdf2*sum(y2)*h2.BinWidth,':','LineWidth',1.5,'color',[0.5 0.5 0.5]*0.5);
                hold on
            end
            
            
            plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max],...
                'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
            hold on
            
            le = legend('Correct','Error',...
                'Location','northeast','fontsize',7);
            le.ItemTokenSize = ones(1,2)*8;
            
            
            tempTxt = sprintf('');
            if temp_p < 0.001
                tempTxt = sprintf('***');
            elseif temp_p < 0.01
                tempTxt = sprintf('**');
            elseif temp_p < 0.05
                tempTxt = sprintf('*');
            end
            text(0.4,y_min+(y_max-y_min)*0.7,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
            ylim([y_min y_max+(y_max-y_min)*0.45]);%0.1
            %xticks([0 1]);
            set(gca, 'FontSize', 10)
            set(gca,'box','off');% 取消右、上边框
            xlabel('Memory precision', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('Trial count', 'FontSize', 11, 'FontWeight', 'bold');
            
            
            temp_trialNum = length(x1)+length(x2);
            
            %temp_title = title(sprintf('%s \n Example seq: %d (%d trials)',FOVName_currentFOV2,...
            %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',11,'FontWeight','bold');
            temp_title = title(sprintf('Stimuli-labeled trial\nSeq %d (%d trials)',...
                seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
            temp_title.Interpreter = 'none';
            
            
        elseif if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 1
            %set(gcf,'Position',[10 450 355 200*2*0.85*0.9*0.98*0.985]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
            
            %set(gcf,'Position',[10 450 355*1.5*1.05 295]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355*1.5*1.05*0.78 295*0.78]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94 295*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*1.11 295*0.78*0.94*1.09]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*1.11*0.93 295*0.78*0.94*1.09*1.16]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[10 450 355*1.5*1.05*0.78*0.94*1.11*0.93*0.79 295*0.78*0.94*1.09*1.16*0.98]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,3,'TileSpacing','compact','Padding','loose');
            
            
            
            %nexttile([1 2])
            nexttile
            
            x1 = memoryPrecision_trialLevel_mildSeq_allMemoryStimuliError;
            x2 = memoryPrecision_trialLevel_mildSeq_allMemoryCorrect;
            %x3 = memoryPrecision_trialLevel_mildSeq_allMemoryResponseError;
            
            x1_mean = mean(x1);
            x2_mean = mean(x2);
            
            [~,temp_p] = ttest2(x1,x2);
            
            temp_1 = x1;
            temp_2 = x2;
            
            if ~isempty(temp_1) && ~isempty(temp_2) == 1
                
                %temp_y_min = min([temp_1;temp_2]);
                temp_y_min = min([temp_1;temp_2;0]);
                temp_y_max = max([temp_1;temp_2]);
                
                
                temp_data = [temp_1;temp_2];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                
                temp_label = [g1;g2];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 2, 1);
                
                %violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                %    'GroupOrder',[{'A'};{'B'}]);
                
                h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'}]);
                h(1).ViolinPlot.FaceAlpha = 0.1;
                h(2).ViolinPlot.FaceAlpha = 0.1;
                
                
                tempTxt = sprintf('');
                if temp_p < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p < 0.05
                    tempTxt = sprintf('*');
                end
                text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                set(gca,'linewidth',1.5)
                %xlim([0.15 2.65])
                xlim([0.5 2.5])
                %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
                set(gca, 'FontSize', 8);%10
                
                %xtl = ["Correct"; "Error"];
                xtl = ["Error"; "Correct"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;
                %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
                xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.235;%0.56,0.4
                %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25-->0
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9
                set(gca,'xticklabel','');
                
                yticks([0 1]);
                
                set(gca,'box','off');% 取消右、上边框
                
                %ylabel('Memory precision', 'FontSize', 11, 'FontWeight', 'bold');
                ylabel('WM strength', 'FontSize', 9, 'FontWeight', 'normal');%10
                temp_trialNum = length(x1)+length(x2);
                %temp_title = title(sprintf('Stimuli-labeled trial\nSeq %d (%d trials)',...
                %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
                %temp_title = title(sprintf('Stimuli-labeled\nSeq %d (%d trials)',...
                %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',9,'FontWeight','normal');
                %temp_title.Interpreter = 'none';
                
                title(sprintf('Stimuli-labeled'),'FontSize',9,'FontWeight','normal');
                subtitle(sprintf('Seq %d (%d trials)',...
                    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',7,'FontWeight','normal');
                
                
            end
            
        end
        
        
        % Compare memory precision of error and correct trials
        nexttile
        
        temp_p = temp_p_correctErrorPrecision_stimuli;
        temptempBoolIndex = (~isnan(d2)) & (~isnan(d1));
        temp_1 = d2(temptempBoolIndex);
        temp_2 = d1(temptempBoolIndex);
        
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        %temp_y_min = min([temp_1;temp_2]);
        temp_y_min = min([temp_1;temp_2;0]);
        temp_y_max = max([temp_1;temp_2]);
        
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            
            scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
            hold on
            
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
        end
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        %xlim([0.15 2.85])
        xlim([0.5 2.5])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
        set(gca, 'FontSize', 8);%10
        
        xticks([1 2]);
        
        xtl = ["Error"; "Correct"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.325;%0.56,0.4
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.22;%0.56,0.4
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%9
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        if if_monkey_D0_Z1 == 1
            yticks([0 0.3]);
        end
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('WM strength', 'FontSize', 9, 'FontWeight', 'normal');%10
        %title(sprintf('Stimuli-labeled\nAll seqs'),'fontsize',9);
        
        title(sprintf('Stimuli-labeled'),'FontSize',9,'FontWeight','normal');
        subtitle(sprintf('All seqs'),'FontSize',7,'FontWeight','normal');
        
        
        
        % Compare accuracy of lowPrecision & highPrecision
        nexttile
        
        temp_p = temp_p_lowHigh_stimuli;
        
        temptempBoolIndex = (~isnan(b1)) & (~isnan(b2));
        temp_1 = b1(temptempBoolIndex);
        temp_2 = b2(temptempBoolIndex);
                
        if isempty(temp_1)
            temp_1 = [0.4;0.5];
            temp_2 = [0.4;0.5];
        end
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            
        end
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
        set(gca, 'FontSize', 8)
        
        xticks([1 2]);
        
        xtl = ["Low-precision"; "High-precision"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.205;%0.56,0.4
        %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',7.5);%25
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Accuracy', 'FontSize', 9, 'FontWeight', 'normal');
        %title(sprintf('Stimuli-labeled trial\nSeq accuracy'),'fontsize',9);
        %title(sprintf('Stimuli-labeled\nAll seqs'),'fontsize',9);
        
        title(sprintf('Stimuli-labeled'),'FontSize',9,'FontWeight','normal');
        subtitle(sprintf('All seqs'),'FontSize',7,'FontWeight','normal');
        
        
        
        
        %% memoryPrecision_trialLevelEvidence in response-labeled trials
        fig = figure('Name','memory precsion pdf (mildSeq)','NumberTitle','off');
        
        if if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 0
            %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[380 450 355 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            t = tiledlayout(1,2,'TileSpacing','compact','Padding','loose');
            
            
            %nexttile([1 2])
            nexttile
            
            x1 = memoryPrecision_trialLevel_mildSeq_allMemoryCorrect;
            %x2 = memoryPrecision_trialLevel_mildSeq_allMemoryStimuliError;
            x2 = memoryPrecision_trialLevel_mildSeq_allMemoryResponseError;
            
            [~,temp_p] = ttest2(x1,x2);
            
            
            h_NumBins = 8;%10
            
            x = x1;
            h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h1.NumBins = h_NumBins;
            h1.EdgeColor = [0 0 0];
            
            x = x2;
            h2 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h2.NumBins = h_NumBins;
            h2.EdgeColor = [0.5 0.5 0.5];
            
            
            y1 = h1.Values;
            y2 = h2.Values;
            
            x_min = 0;
            x_max = 1;
            
            [y1_min,y1_max] = bounds(y1);
            [y2_min,y2_max] = bounds(y2);
            y_min = min([y1_min,y2_min]);
            y_max = max([y1_max,y2_max]);
            
            if if_plot_additionalSmooth == 1
                n=100;
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                
                [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
                plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',[0 0 0]);
                hold on
                
                [pdf2,xmesh2,~] = ksdensity(x2,'NumPoints',n,'Function','pdf');
                plot(xmesh2,pdf2*sum(y2)*h2.BinWidth,':','LineWidth',1.5,'color',[0.5 0.5 0.5]*0.5);
                hold on
            end
            
            plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max],...
                'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
            hold on
            
            le = legend('Correct','Error',...
                'Location','northeast','fontsize',7);
            le.ItemTokenSize = ones(1,2)*8;
            
            
            tempTxt = sprintf('');
            if temp_p < 0.001
                tempTxt = sprintf('***');
            elseif temp_p < 0.01
                tempTxt = sprintf('**');
            elseif temp_p < 0.05
                tempTxt = sprintf('*');
            end
            text(0.4,y_min+(y_max-y_min)*0.7,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
            ylim([y_min y_max+(y_max-y_min)*0.45]);%0.1
            %xticks([0 1]);
            set(gca, 'FontSize', 10)
            set(gca,'box','off');% 取消右、上边框
            xlabel('Memory precision', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('Trial count', 'FontSize', 11, 'FontWeight', 'bold');
            
            
            temp_trialNum = length(x1)+length(x2);
            
            %temp_title = title(sprintf('%s \n Example seq: %d (%d trials)',FOVName_currentFOV2,...
            %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',11,'FontWeight','bold');
            temp_title = title(sprintf('Response-labeled trial\nSeq %d (%d trials)',...
                seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
            
            temp_title.Interpreter = 'none';
            
        elseif if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 1
            %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355 200*2*0.85*0.9*0.98*0.985]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
            
            set(gcf,'Position',[10 450 355*1.5*1.05 295]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            t = tiledlayout(2,3,'TileSpacing','compact','Padding','loose');
            
            
            
            %nexttile([1 2])
            nexttile
            
            x1 = memoryPrecision_trialLevel_mildSeq_allMemoryResponseError;
            x2 = memoryPrecision_trialLevel_mildSeq_allMemoryCorrect;
            %x2 = memoryPrecision_trialLevel_mildSeq_allMemoryStimuliError;
            
            [~,temp_p] = ttest2(x1,x2);
            
            temp_1 = x1;
            temp_2 = x2;
            
            if ~isempty(temp_1) && ~isempty(temp_2) == 1
                
                temp_y_min = min([temp_1;temp_2]);
                temp_y_max = max([temp_1;temp_2]);
                
                
                temp_data = [temp_1;temp_2];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                
                temp_label = [g1;g2];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 2, 1);
                
                %violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                %    'GroupOrder',[{'A'};{'B'}]);
                
                h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'}]);
                h(1).ViolinPlot.FaceAlpha = 0.1;
                h(2).ViolinPlot.FaceAlpha = 0.1;
                
                
                tempTxt = sprintf('');
                if temp_p < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p < 0.05
                    tempTxt = sprintf('*');
                end
                text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                set(gca,'linewidth',1.5)
                %xlim([0.15 2.65])
                xlim([0.5 2.5])
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
                set(gca, 'FontSize', 10)
                
                %xtl = ["Correct"; "Error"];
                xtl = ["Error"; "Correct"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;
                xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.335;%0.56,0.4
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
                set(gca,'xticklabel','');
                
                yticks([0 1]);
                
                set(gca,'box','off');% 取消右、上边框
                
                ylabel('Memory precision', 'FontSize', 10, 'FontWeight', 'normal');
                temp_trialNum = length(x1)+length(x2);
                temp_title = title(sprintf('Response-labeled\nSeq %d (%d trials)',...
                    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',9,'FontWeight','normal');
                temp_title.Interpreter = 'none';
                
            end
            
        end
        
        
        
        % Compare memory precision of error and correct trials
        nexttile
        
        temp_p = temp_p_correctErrorPrecision_response;
        temptempBoolIndex = (~isnan(d3)) & (~isnan(d1));
        temp_1 = d3(temptempBoolIndex);
        temp_2 = d1(temptempBoolIndex);
        
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            
        end
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
        set(gca, 'FontSize', 10)
        
        xticks([1 2]);
        
        xtl = ["Error"; "Correct"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.325;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Memory precision', 'FontSize', 10, 'FontWeight', 'normal');
        title(sprintf('Response-labeled\nAll seqs'),'fontsize',9);
        
        
        % Compare accuracy of lowPrecision & highPrecision
        nexttile
        
        temp_p = temp_p_lowHigh_response;
        temptempBoolIndex = (~isnan(b3)) & (~isnan(b4));
        temp_1 = b3(temptempBoolIndex);
        temp_2 = b4(temptempBoolIndex);
        
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            
        end
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 10)
        
        xticks([1 2]);
        
        xtl = ["Low-precision"; "High-precision"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.415;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Accuracy', 'FontSize', 10, 'FontWeight', 'normal');
        %title(sprintf('Response-labeled trial\nSeq accuracy'),'fontsize',10);
        title(sprintf('Response-labeled\nAll seqs'),'fontsize',9);
        
        
        
        
        
        
        
        
        %% memoryPrecision & meta trialLevelEvidence in stimuli-labeled trials
        fig = figure('Name','memory precsion pdf (mildSeq)','NumberTitle','off');
        
        if if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 0
            %set(gcf,'Position',[10 450 355 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[10 450 355*1.03 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            t = tiledlayout(1,2,'TileSpacing','compact','Padding','loose');
            
            
            nexttile
            
            x1 = memoryPrecision_trialLevel_mildSeq_choiceMemoryStimuli;
            %x1 = memoryPrecision_trialLevel_mildSeq_choiceMemoryResponse;
            x2 = memoryPrecision_trialLevel_mildSeq_choiceOffload;
            
            [~,temp_p] = ttest2(x1,x2);
            
            
            h_NumBins = 8;%10
            
            x = x1;
            h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h1.NumBins = h_NumBins;
            h1.EdgeColor = color_choiceMemory;
            
            x = x2;
            h2 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h2.NumBins = h_NumBins;
            h2.EdgeColor = color_choiceOffload;
            
            
            y1 = h1.Values;
            y2 = h2.Values;
            
            x_min = 0;
            x_max = 1;
            
            [y1_min,y1_max] = bounds(y1);
            [y2_min,y2_max] = bounds(y2);
            y_min = min([y1_min,y2_min]);
            y_max = max([y1_max,y2_max]);
            
            if if_plot_additionalSmooth == 1
                n=100;
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                
                [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
                plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',color_choiceMemory);
                hold on
                
                [pdf2,xmesh2,~] = ksdensity(x2,'NumPoints',n,'Function','pdf');
                plot(xmesh2,pdf2*sum(y2)*h2.BinWidth,':','LineWidth',1.5,'color',color_choiceOffload);
                hold on
            end
            
            plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max],...
                'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
            hold on
            
            le = legend('Choice-memory','Choice-offload',...
                'Location','northeast','fontsize',7);
            le.ItemTokenSize = ones(1,2)*8;
            
            
            tempTxt = sprintf('');
            if temp_p < 0.001
                tempTxt = sprintf('***');
            elseif temp_p < 0.01
                tempTxt = sprintf('**');
            elseif temp_p < 0.05
                tempTxt = sprintf('*');
            end
            text(0.4,y_min+(y_max-y_min)*0.7,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
            ylim([y_min y_max+(y_max-y_min)*0.45]);%0.1
            %xticks([0 1]);
            set(gca, 'FontSize', 10)
            set(gca,'box','off');% 取消右、上边框
            xlabel('Memory precision', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('Trial count', 'FontSize', 11, 'FontWeight', 'bold');
            
            
            temp_trialNum = length(x1)+length(x2);
            
            %temp_title = title(sprintf('%s \n Example seq: %d (%d trials)',FOVName_currentFOV2,...
            %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',11,'FontWeight','bold');
            temp_title = title(sprintf('Stimuli-labeled trial\nSeq %d (%d trials)',...
                seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
            
            temp_title.Interpreter = 'none';
            
        elseif if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 1
            
            
            %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[10 450 355 200*2*0.85*0.9*0.98*0.985]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
            
            %nexttile([1 2])
            nexttile
            
            x1 = memoryPrecision_trialLevel_mildSeq_choiceOffload;
            x2 = memoryPrecision_trialLevel_mildSeq_choiceMemoryStimuli;
            %x1 = memoryPrecision_trialLevel_mildSeq_choiceMemoryResponse;
            
            [~,temp_p] = ttest2(x1,x2);
            
            temp_1 = x1;
            temp_2 = x2;
            
            if ~isempty(temp_1) && ~isempty(temp_2) == 1
                
                temp_y_min = min([temp_1;temp_2]);
                temp_y_max = max([temp_1;temp_2]);
                
                
                temp_data = [temp_1;temp_2];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                
                temp_label = [g1;g2];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 2, 1);
                
                %violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                %    'GroupOrder',[{'A'};{'B'}]);
                
                h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'}]);
                h(1).ViolinPlot.FaceAlpha = 0.1;
                h(2).ViolinPlot.FaceAlpha = 0.1;
                
                
                tempTxt = sprintf('');
                if temp_p < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p < 0.05
                    tempTxt = sprintf('*');
                end
                text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                set(gca,'linewidth',1.5)
                xlim([0.15 2.65])
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
                set(gca, 'FontSize', 10)
                
                %xtl = ["Choice-memory"; "Choice-offload"];
                xtl = ["Choice-offload"; "Choice-memory"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;
                xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25,10
                set(gca,'xticklabel','');
                
                
                set(gca,'box','off');% 取消右、上边框
                
                ylabel('Memory precision', 'FontSize', 11, 'FontWeight', 'bold');
                temp_trialNum = length(x1)+length(x2);
                temp_title = title(sprintf('Stimuli-labeled trial\nSeq %d (%d trials)',...
                    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
                temp_title.Interpreter = 'none';
                
            end
            
        end
        
        
        nexttile
        
        temp_p = temp_p_c_lowHigh_stimuli;
        temptempBoolIndex = (~isnan(c1)) & (~isnan(c2));
        temp_1 = c1(temptempBoolIndex);
        temp_2 = c2(temptempBoolIndex);
        
        
        if isempty(temp_1)
            temp_1 = [0.4;0.5];
            temp_2 = [0.4;0.5];
        end
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        %violinplot([temp_1 temp_2],[],'ViolinColor',[1 1 1]*0.5,'BoxColor',[1 1 1]*0.2);
        
        temp_data = [temp_1;temp_2];
        
        g1 = repmat({'A'},length(temp_1),1);
        g2 = repmat({'B'},length(temp_2),1);
        
        temp_label = [g1;g2];
        
        temptemp_color1 = [1 1 1]*0.5;
        temptemp_color2 = repmat(temptemp_color1, 2, 1);
        
        h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
            'GroupOrder',[{'A'};{'B'}]);
        h(1).ViolinPlot.FaceAlpha = 0.1;
        h(2).ViolinPlot.FaceAlpha = 0.1;
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.65])
        %ylim([temp_y_min-(temp_y_max-temp_y_min)*0.25 temp_y_max+(temp_y_max-temp_y_min)*0.25]);
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 10)
        
        %set(gca,'XTickLabel', ["Low"; "High"],'FontSize', 12);%给坐标加标签
        %set(gca,'XTickLabel', ["LowPrecision"; "HighPrecision"],'FontSize', 10);%给坐标加标签
        %xtickangle(20);
        
        %xtl = ["LowPrecision"; "HighPrecision"];
        xtl = ["Low-precision"; "High-precision"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.435;%0.56,0.4
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',10);%25
        set(gca,'xticklabel','');
        
        
        %yticks([0.7 0.8]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Offloading rate', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('Stimuli-labeled trial\nSeq offloading rate'),'fontsize',10);
        
        
        
        
        
        
        %% memoryPrecision & meta trialLevelEvidence in response-labeled trials
        fig = figure('Name','memory precsion pdf (mildSeq)','NumberTitle','off');
        
        if if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 0
            set(gcf,'Position',[380 450 355*1.03 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            t = tiledlayout(1,2,'TileSpacing','compact','Padding','loose');
            
            
            nexttile
            
            %x1 = memoryPrecision_trialLevel_mildSeq_choiceMemoryStimuli;
            x1 = memoryPrecision_trialLevel_mildSeq_choiceMemoryResponse;
            x2 = memoryPrecision_trialLevel_mildSeq_choiceOffload;
            
            [~,temp_p] = ttest2(x1,x2);
            
            
            h_NumBins = 8;%10
            
            x = x1;
            h1 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h1.NumBins = h_NumBins;
            h1.EdgeColor = color_choiceMemory;
            
            x = x2;
            h2 = histogram(x,'FaceAlpha',1,'Normalization','count','DisplayStyle','stairs','LineWidth',1.5);
            hold on
            h2.NumBins = h_NumBins;
            h2.EdgeColor = color_choiceOffload;
            
            
            y1 = h1.Values;
            y2 = h2.Values;
            
            x_min = 0;
            x_max = 1;
            
            [y1_min,y1_max] = bounds(y1);
            [y2_min,y2_max] = bounds(y2);
            y_min = min([y1_min,y2_min]);
            y_max = max([y1_max,y2_max]);
            
            if if_plot_additionalSmooth == 1
                n=100;
                n=2^ceil(log2(n)); % round up n to the next power of 2;
                
                [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
                plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',color_choiceMemory);
                hold on
                
                [pdf2,xmesh2,~] = ksdensity(x2,'NumPoints',n,'Function','pdf');
                plot(xmesh2,pdf2*sum(y2)*h2.BinWidth,':','LineWidth',1.5,'color',color_choiceOffload);
                hold on
            end
            
            
            plot([lowThreshold_memoryPrecision lowThreshold_memoryPrecision],[y_min y_max],...
                'LineWidth',1,'color',[0.9290 0.6940 0.1250]);
            hold on
            
            le = legend('Choice-memory','Choice-offload',...
                'Location','northeast','fontsize',7);
            le.ItemTokenSize = ones(1,2)*8;
            
            
            tempTxt = sprintf('');
            if temp_p < 0.001
                tempTxt = sprintf('***');
            elseif temp_p < 0.01
                tempTxt = sprintf('**');
            elseif temp_p < 0.05
                tempTxt = sprintf('*');
            end
            text(0.4,y_min+(y_max-y_min)*0.7,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            
            set(gca,'linewidth',1.5)
            xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
            ylim([y_min y_max+(y_max-y_min)*0.45]);%0.1
            %xticks([0 1]);
            set(gca, 'FontSize', 10)
            set(gca,'box','off');% 取消右、上边框
            xlabel('Memory precision', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('Trial count', 'FontSize', 11, 'FontWeight', 'bold');
            
            
            temp_trialNum = length(x1)+length(x2);
            
            %temp_title = title(sprintf('%s \n Example seq: %d (%d trials)',FOVName_currentFOV2,...
            %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',11,'FontWeight','bold');
            temp_title = title(sprintf('Response-labeled trial\nSeq %d (%d trials)',...
                seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
            
            temp_title.Interpreter = 'none';
            
            
        elseif if_trialEvidenceExampleSeq_histogram0_violinPlot1 == 1
            %set(gcf,'Position',[10 50 550 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %set(gcf,'Position',[10 450 355 200*2*0.85*0.9*0.98*0.985]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            %t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
            
            set(gcf,'Position',[10 450 355*1.5*1.05 295]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            t = tiledlayout(2,3,'TileSpacing','compact','Padding','loose');
            
            
            %nexttile([1 2])
            nexttile
            
            x1 = memoryPrecision_trialLevel_mildSeq_choiceOffload;
            %x1 = memoryPrecision_trialLevel_mildSeq_choiceMemoryStimuli;
            x2 = memoryPrecision_trialLevel_mildSeq_choiceMemoryResponse;
            
            [~,temp_p] = ttest2(x1,x2);
            
            temp_1 = x1;
            temp_2 = x2;
            
            if ~isempty(temp_1) && ~isempty(temp_2) == 1
                
                temp_y_min = min([temp_1;temp_2]);
                temp_y_max = max([temp_1;temp_2]);
                
                
                temp_data = [temp_1;temp_2];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                
                temp_label = [g1;g2];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 2, 1);
                
                %violinplot(temp_data,temp_label,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                %    'GroupOrder',[{'A'};{'B'}]);
                
                h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'}]);
                h(1).ViolinPlot.FaceAlpha = 0.1;
                h(2).ViolinPlot.FaceAlpha = 0.1;
                
                
                tempTxt = sprintf('');
                if temp_p < 0.001
                    tempTxt = sprintf('***');
                elseif temp_p < 0.01
                    tempTxt = sprintf('**');
                elseif temp_p < 0.05
                    tempTxt = sprintf('*');
                end
                text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                
                set(gca,'linewidth',1.5)
                xlim([0.15 2.85])
                ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
                set(gca, 'FontSize', 10)
                
                %xtl = ["Choice-memory"; "Choice-offload"];
                xtl = ["Choice-offload"; "Choice-memory"];
                xt=get(gca,'XTick');
                yt=get(gca,'YTick');
                xtext_xp=xt;
                xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.425;%0.56,0.4
                text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25,10
                set(gca,'xticklabel','');
                
                yticks([0 1]);
                
                set(gca,'box','off');% 取消右、上边框
                
                ylabel('Memory precision', 'FontSize', 10, 'FontWeight', 'normal');
                temp_trialNum = length(x1)+length(x2);
                %temp_title = title(sprintf('Response-labeled trial\nSeq %d (%d trials)',...
                %    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',10,'FontWeight','normal');
                temp_title = title(sprintf('Response-labeled\nSeq %d (%d trials)',...
                    seqSet_inOne_inOne(exampleSeq),temp_trialNum),'FontSize',9,'FontWeight','normal');
                temp_title.Interpreter = 'none';
                
            end
            
        end
        
        
        % Compare memory precision of choice-offload and choice-memory trials
        nexttile
        
        temp_p = temp_p_choicePrecision_response;
        temptempBoolIndex = (~isnan(d4)) & (~isnan(d6));
        temp_1 = d4(temptempBoolIndex);
        temp_2 = d6(temptempBoolIndex);
        
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            
        end
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.85])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
        set(gca, 'FontSize', 10)
        
        xticks([1 2]);
        
        xtl = ["Choice-offload"; "Choice-memory"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.355;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Memory precision', 'FontSize', 10, 'FontWeight', 'normal');
        title(sprintf('Response-labeled\nAll seqs'),'fontsize',9);
        
        
        
        % Compare offloading rate of lowPrecision & highPrecision
        nexttile
        
        temp_p = temp_p_c_lowHigh_response;
        temptempBoolIndex = (~isnan(c3)) & (~isnan(c4));
        temp_1 = c3(temptempBoolIndex);
        temp_2 = c4(temptempBoolIndex);
        
        
        temp1_SEM = std(temp_1)/sqrt(length(temp_1));
        temp2_SEM = std(temp_2)/sqrt(length(temp_2));
        
        
        temp_y_min = min([temp_1;temp_2]);
        temp_y_max = max([temp_1;temp_2]);
        
        
        if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
            temp_data = [temp_1;temp_2];
            
            g1 = repmat({'A'},length(temp_1),1);
            g2 = repmat({'B'},length(temp_2),1);
            
            temp_label = [g1;g2];
            
            temptemp_color1 = [1 1 1]*0.5;
            temptemp_color2 = repmat(temptemp_color1, 2, 1);
            
            h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                'GroupOrder',[{'A'};{'B'}]);
            h(1).ViolinPlot.FaceAlpha = 0.1;
            h(2).ViolinPlot.FaceAlpha = 0.1;
            
        elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
            
            plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
            hold on
            
            for tempi=1:length(temp_1)
                plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                hold on
            end
            
            scatter(1*ones(1,length(temp_1)),temp_1,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            scatter(2*ones(1,length(temp_2)),temp_2,30,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);
            hold on
            
        end
        
        
        
        tempTxt = sprintf('');
        if temp_p < 0.001
            tempTxt = sprintf('***');
        elseif temp_p < 0.01
            tempTxt = sprintf('**');
        elseif temp_p < 0.05
            tempTxt = sprintf('*');
        end
        text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.08,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','center');
        
        set(gca,'linewidth',1.5)
        xlim([0.15 2.65])
        ylim([temp_y_min-(temp_y_max-temp_y_min)*0.15 temp_y_max+(temp_y_max-temp_y_min)*0.15]);
        set(gca, 'FontSize', 10)
        
        xticks([1 2]);
        xtl = ["Low-precision"; "High-precision"];
        xt=get(gca,'XTick');
        yt=get(gca,'YTick');
        xtext_xp=xt;
        xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.405;%0.56,0.4
        text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
        set(gca,'xticklabel','');
        
        
        yticks([0 1]);
        
        
        set(gca,'box','off');% 取消右、上边框
        ylabel('Offloading rate', 'FontSize', 10, 'FontWeight', 'normal');
        %title(sprintf('Response-labeled trial\nSeq offloading rate'),'fontsize',10);
        title(sprintf('Response-labeled\nAll seqs'),'fontsize',9);
        
        
        
        
    end
    
    
    %% Single neuron tuning property of memoryPrecision and choiceMemory (trial-level tuning)
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        tempStr = 'trial-level';
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        tempStr = 'seq-level';
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        tempStr = 'mix-level';
    end
    if if_plot_3d_beta == true
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        %set(gcf,'Position',[10 450 400 400]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[10 150 600 600]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        t = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');
        
        %     x = r2_memoryPrecision;
        %     y = r2_choiceMemory;
        %     z = r2_choiceMemory_baseline;
        
        %     x = r_memoryPrecision;
        %     y = r_choiceMemory;
        %     z = r_choiceMemory_baseline;
        
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_memoryPrecision;
            y = beta_choiceMemory;
            z = beta_choiceMemory_baseline;
            tempStr = 'trial-level';
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_seqPrecision;
            y = beta_gProb;
            z = beta_gProb_baseline;
            tempStr = 'seq-level';
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_seqPrecision;
            y = beta_gProb;
            z = beta_choiceMemory_baseline;
            tempStr = 'mix-level';
        end
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        [z_min,z_max] = bounds(z);
        
        tempSize = 10;%10
        
        %scatter(x,y,10,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        scatter3(x,y,z,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        axis equal
        
        set(gca,'linewidth',1.5)
        %     xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %     ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        %     zlim([z_min-(z_max-z_min)*0.08 z_max+(z_max-z_min)*0.08]);
        set(gca, 'FontSize', 12)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Memory precision', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
        zlabel('Baseline', 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title = title(sprintf('%s (trial-level)\n Single neuron tuning property, n=%d',FOVName_currentFOV2,length(x)),'FontSize',11);
        temp_title = title(sprintf('%s (%s)\n Single neuron tuning property, n=%d',FOVName_currentFOV2,tempStr,length(x)),'FontSize',11);
        temp_title.Interpreter = 'none';
        
    end
    
    
    
    %% Single neuron tuning property of memoryPrecision and choiceMemory (trial-level tuning)
    if if_plot_2d_3sub_beta == true
    %if if_plot_2d_3sub_beta == true && if_loadDemoData == 0
        fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
        %set(gcf,'Position',[510 450 335 281]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[510 450 240*3 240]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        %set(gcf,'Position',[510 450 240*3 240*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        set(gcf,'Position',[510 450 240*3*0.67 240*1.05*0.67*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
        
        %t = tiledlayout(1,3,'TileSpacing','compact','Padding','loose');
        t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        
        %t.Title.String = sprintf('%s (trial-level)\n Single neuron tuning property, n=%d',FOVName_currentFOV2,length(x));
        %t.Title.String = sprintf('%s (%s)\n Single neuron tuning property, n=%d',FOVName_currentFOV2,tempStr,length(x));
        t.Title.String = sprintf('%s (%s)\n Regression coefficients of neurons, n=%d',FOVName_currentFOV2,tempStr,length(x));
        t.Title.FontSize = 9;%11
        t.Title.Interpreter = 'none';
        
        tempSize = 4;%10
        
        nexttile
        
        %x = beta_memoryPrecision;
        %y = beta_choiceMemory;
        
        if if_tuning_location0_precision1 == 1
            if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
                x = beta_memoryPrecision;
                y = beta_choiceMemory;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
                x = beta_seqPrecision;
                y = beta_gProb;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
                x = beta_seqPrecision;
                y = beta_gProb;
            end
            
        elseif if_tuning_location0_precision1 == 0
            if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
                x = beta_6loc_peak;
                y = beta_choiceMemory;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
                x = beta_6loc_peak;
                y = beta_gProb;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
                x = beta_6loc_peak;
                y = beta_gProb;
            end
        end
        
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        
        tempTxt = 'p>0.05';
        if temp_p < 0.05
            tempTxt = 'p<0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p<0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p<0.001';
        end
        text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.18,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
        text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.05,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        
        
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        xticks([-1 0 1]);
        yticks([-1 0 1]);
        
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        if if_tuning_location0_precision1 == 1
            xlabel('Memory precision', 'FontSize', 9, 'FontWeight', 'normal');
        elseif if_tuning_location0_precision1 == 0
            xlabel('Memory', 'FontSize', 9, 'FontWeight', 'normal');
        end
        ylabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title.Interpreter = 'none';
        
        
        nexttile
        %x = beta_memoryPrecision;
        %y = beta_choiceMemory_baseline;
        
        if if_tuning_location0_precision1 == 1
            if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
                x = beta_memoryPrecision;
                y = beta_choiceMemory_baseline;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
                x = beta_seqPrecision;
                y = beta_gProb_baseline;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
                x = beta_seqPrecision;
                y = beta_choiceMemory_baseline;
            end
        elseif if_tuning_location0_precision1 == 0
            if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
                x = beta_6loc_peak;
                y = beta_choiceMemory_baseline;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
                x = beta_6loc_peak;
                y = beta_gProb_baseline;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
                x = beta_6loc_peak;
                y = beta_choiceMemory_baseline;
            end
        end
        
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        
        tempTxt = 'p>0.05';
        if temp_p < 0.05
            tempTxt = 'p<0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p<0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p<0.001';
        end
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.98,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
            text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.85,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.18,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
            text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.05,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        end
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        xticks([-1 0 1]);
        yticks([-1 0 1]);
        
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        if if_tuning_location0_precision1 == 1
            xlabel('Memory precision', 'FontSize', 9, 'FontWeight', 'normal');
        elseif if_tuning_location0_precision1 == 0
            xlabel('Memory', 'FontSize', 9, 'FontWeight', 'normal');
        end
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title.Interpreter = 'none';
        
        nexttile
        
        %x = beta_choiceMemory;
        %y = beta_choiceMemory_baseline;
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            x = beta_choiceMemory;
            y = beta_choiceMemory_baseline;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
            x = beta_gProb;
            y = beta_gProb_baseline;
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            x = beta_gProb;
            y = beta_choiceMemory_baseline;
        end
        
        
        [temp_r,temp_p] = corr(x,y);
        
        [x_min,x_max] = bounds(x);
        [y_min,y_max] = bounds(y);
        xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
        xy_min = -1*xy_max;
        
        
        scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        
        plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on
        
        
        tempTxt = 'p>0.05';
        if temp_p < 0.05
            tempTxt = 'p<0.05';
        end
        if temp_p < 0.01
            tempTxt = 'p<0.01';
        end
        if temp_p < 0.001
            tempTxt = 'p<0.001';
        end
        if if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
            text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.98,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
            text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.85,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
            text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.18,sprintf('r=%.3f',temp_r), 'fontsize',8,'FontWeight','normal');
            text(x_min+(x_max-x_min)*0.68,y_min+(y_max-y_min)*0.05,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        end
        
        
        %axis equal
        
        set(gca,'linewidth',1.5)
        %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
        %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
        
        xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
        
        xticks([-1 0 1]);
        yticks([-1 0 1]);
        
        set(gca, 'FontSize', 9)
        set(gca,'box','off');% 取消右、上边框
        xlabel('Meta-memory', 'FontSize', 9, 'FontWeight', 'normal');
        ylabel('Baseline', 'FontSize', 9, 'FontWeight', 'normal');
        %temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
        %temp_title.Interpreter = 'none';
        
    end
    
end


if if_plot == 1
    %% Compare memory precision of choice-offload and choice-memory trials
    if exist('temp_p_choicePrecision_stimuli','var') == 1
        if if_plot_memoryPrecision_trialLevelEvidence == 1
            
            fig = figure('Name','WM strength in two decisions','NumberTitle','off');
            %set(gcf,'Position',[10 450 355*1.5*1.05*0.97*0.84*0.74 295*1.06*0.88]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            set(gcf,'Position',[10 450 337.1*0.4*0.9*0.9 295*1.06*0.88]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
            
            t = tiledlayout(2,1,'TileSpacing','loose','Padding','loose');
            
            
            nexttile
            
            temp_p = temp_p_choicePrecision_stimuli;
            temptempBoolIndex = (~isnan(d4)) & (~isnan(d5));
            temp_1 = d4(temptempBoolIndex);
            temp_2 = d5(temptempBoolIndex);
            
            
            temp1_SEM = std(temp_1)/sqrt(length(temp_1));
            temp2_SEM = std(temp_2)/sqrt(length(temp_2));
            
            
            temp_y_min = min([temp_1;temp_2]);
            
            %temp_y_min = min([temp_y_min,0.06]);
            
            temp_y_max = max([temp_1;temp_2]);
            
            temp_y_min = min([temp_y_min,0]);
            %temp_y_max = max([temp_y_max,1]);
            temp_y_max = ceil(temp_y_max*10)/10;
            
            if if_trialEvidenceAllSeq_violinplot0_pairline1 == 0
                temp_data = [temp_1;temp_2];
                
                g1 = repmat({'A'},length(temp_1),1);
                g2 = repmat({'B'},length(temp_2),1);
                
                temp_label = [g1;g2];
                
                temptemp_color1 = [1 1 1]*0.5;
                temptemp_color2 = repmat(temptemp_color1, 2, 1);
                
                h = violinplot(temp_data,temp_label,'ViolinAlpha',0.3,'ViolinColor',temptemp_color2,'BoxColor',[1 1 1]*0.2,...
                    'GroupOrder',[{'A'};{'B'}]);
                h(1).ViolinPlot.FaceAlpha = 0.1;
                h(2).ViolinPlot.FaceAlpha = 0.1;
                
            elseif if_trialEvidenceAllSeq_violinplot0_pairline1 == 1
                
                plot([1-0.35 1+0.35],[1 1]*mean(temp_1),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                plot([2-0.35 2+0.35],[1 1]*mean(temp_2),'Color',[0.6350 0.0780 0.1840],'LineWidth',4);
                hold on
                
                for tempi=1:length(temp_1)
                    plot([1 2],[temp_1(tempi) temp_2(tempi)],'Color',[1 1 1]*0.4);
                    hold on
                end
                
                scatter(1*ones(1,length(temp_1)),temp_1,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
                hold on
                scatter(2*ones(1,length(temp_2)),temp_2,15,'filled','MarkerFaceColor',[1 1 1]*0.05,...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.7);%30
                hold on
                
            end
            
            
            tempTxt = sprintf('');
            if temp_p < 0.001
                tempTxt = sprintf('***');
            elseif temp_p < 0.01
                tempTxt = sprintf('**');
            elseif temp_p < 0.05
                tempTxt = sprintf('*');
            end
            text(1.5,temp_y_max+(temp_y_max-temp_y_min)*0.12,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
                'HorizontalAlignment','center');
            
            set(gca,'linewidth',1.5)
            xlim([0.15 2.85])
            ylim([temp_y_min-(temp_y_max-temp_y_min)*0.125 temp_y_max+(temp_y_max-temp_y_min)*0.20]);
            set(gca, 'FontSize', 8)
            
            xticks([1 2]);
            
            %xtl = ["Choice-offload"; "Choice-memory"];
            xtl = ["Offload"; "Memory"];
            xt=get(gca,'XTick');
            yt=get(gca,'YTick');
            xtext_xp=xt;
            %xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.345;%0.56,0.4
            xtext_yp=(yt(1))*ones(1,length(xt))-(temp_y_max-temp_y_min)*0.19;%0.56,0.4
            %text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','center','rotation',25,'fontsize',9);%25
            text(xtext_xp,xtext_yp,xtl,'HorizontalAlignment','right','rotation',25,'fontsize',8);%25
            set(gca,'xticklabel','');
            
            
            %yticks([0 1]);
            
            yticks([0 temp_y_max]);
            
            set(gca,'box','off');% 取消右、上边框
            ylabel('WM strength', 'FontSize', 9, 'FontWeight', 'normal');
            title(sprintf('Stimuli-labeled\nAll seqs'),'fontsize',9);
            
        end
    end
end


% if false
% if if_loadDemoData == 0
if true
    
    cellIndex_suite2p_A = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1); %#ok<*UNRCH>
    % cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,2);
    cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
    cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow = cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryLow;
    cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh = cellIndex_suite2p_memoryPrecisionLow_choiceMemoryHigh;
    cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh = cellIndex_suite2p_memoryPrecisionHigh_choiceMemoryHigh;
    cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryLow = cellIndex_suite2p_memoryPrecisionLow_choiceMemoryLow;
    
    cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryLow = nan(length(cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow),1);
    for tempi=1:length(cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow)
        temptempIndex = find(cellIndex_suite2p_A==cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryLow(tempi));
        cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryLow(tempi) = cellIndex_suite2p_B(temptempIndex);
    end
    
    cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryHigh = nan(length(cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh),1);
    for tempi=1:length(cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh)
        temptempIndex = find(cellIndex_suite2p_A==cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryHigh(tempi));
        cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryHigh(tempi) = cellIndex_suite2p_B(temptempIndex);
    end
    
    cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryHigh = nan(length(cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh),1);
    for tempi=1:length(cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh)
        temptempIndex = find(cellIndex_suite2p_A==cellIndex_suite2p_A_memoryPrecisionHigh_choiceMemoryHigh(tempi));
        cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryHigh(tempi) = cellIndex_suite2p_B(temptempIndex);
    end
    
    cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryLow = nan(length(cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryLow),1);
    for tempi=1:length(cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryLow)
        temptempIndex = find(cellIndex_suite2p_A==cellIndex_suite2p_A_memoryPrecisionLow_choiceMemoryLow(tempi));
        cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryLow(tempi) = cellIndex_suite2p_B(temptempIndex);
    end
    
    
    cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryLow;
    cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryHigh;
    cellIndex_suite2p_B_memoryPrecisionHigh_choiceMemoryHigh;
    cellIndex_suite2p_B_memoryPrecisionLow_choiceMemoryLow;
    
    cellIndex_suite2p_A_choiceMemoryBaselineHigh = cellIndex_suite2p_choiceMemoryBaselineHigh;
    cellIndex_suite2p_B_choiceMemoryBaselineHigh = nan(length(cellIndex_suite2p_A_choiceMemoryBaselineHigh),1);
    for tempi=1:length(cellIndex_suite2p_A_choiceMemoryBaselineHigh)
        temptempIndex = find(cellIndex_suite2p_A==cellIndex_suite2p_A_choiceMemoryBaselineHigh(tempi));
        cellIndex_suite2p_B_choiceMemoryBaselineHigh(tempi) = cellIndex_suite2p_B(temptempIndex);
    end
    % end
    
end


%% ePRIRS
if if_ePAIRS == 1
    %x = beta_memoryPrecision;
    %y = beta_choiceMemory;
    %z = beta_choiceMemory_baseline;
    
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_memoryPrecision;
        y = beta_choiceMemory;
        z = beta_choiceMemory_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_seqPrecision;
        y = beta_gProb;
        z = beta_gProb_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_seqPrecision;
        y = beta_gProb;
        z = beta_choiceMemory_baseline;
    end
    
    temp_vectors = [x,y,z];
    [clusteriness, p, dists] = pairsClusterTest_elliptical(temp_vectors);
    
    temp_angles = dists.data;
    temp_angles_bootstrap = dists.bootstrap;
    % temp_angles_bootstrap_median = median(temp_angles_bootstrap,1)';
    %temp_angles_bootstrap_median = median(temp_angles_bootstrap,2);
    % temp_angles_bootstrap_median = temp_angles_bootstrap(1,:)';
    temp_angles_bootstrap_median = reshape(temp_angles_bootstrap,[],1);
    
    fig = figure('Name','asd','NumberTitle','off');
    %set(gcf,'Position',[10 450 355/2*2 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    % set(gcf,'Position',[10 450 177*1.3 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    % set(gcf,'Position',[10 450 177*2 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[10 450 412-10 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[10 450 402*0.43*1.15*1.05 200*0.53*1.08]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    
    nexttile
    
    x1 = temp_angles;
    %x1 = temp_angles_bootstrap_median;
    x2 = temp_angles_bootstrap_median;
    
    [~,temp_p] = ttest2(x1,x2);
    
    
    h_NumBins = 50;%10-->100
    
    x = x1;
    h1 = histogram(x,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
    hold on
    %h1.NumBins = h_NumBins;
    h1.EdgeColor = [0.4660 0.6740 0.1880];%[0.25 0.75 0.25]
    
    x = x2;
    h2 = histogram(x,'FaceAlpha',1,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1.5);
    hold on
    %h2.NumBins = h_NumBins;
    h2.EdgeColor = [0 0 0];
    
    
    x_min = min([x1;x2]);
    %x_max = max([x1;x2]);
    x_max = prctile([x1;x2],99);%97.5
    
    h1.BinLimits = [x_min,x_max];
    h1.NumBins = h_NumBins;
    
    h2.BinLimits = [x_min,x_max];
    h2.NumBins = h_NumBins;
    
    y1 = h1.Values;
    y2 = h2.Values;
    
    
    [y1_min,y1_max] = bounds(y1);
    [y2_min,y2_max] = bounds(y2);
    y_min = min([y1_min,y2_min]);
    y_max = max([y1_max,y2_max]);
    
    
    if if_plot_additionalSmooth == 1
        n=1000;
        n=2^ceil(log2(n)); % round up n to the next power of 2;
        
        [pdf1,xmesh1,~] = ksdensity(x1,'NumPoints',n,'Function','pdf');
        
        %plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',1.5,'color',[0 0 0]);
        plot(xmesh1,pdf1*sum(y1)*h1.BinWidth,':','LineWidth',2,'color',[0.4660 0.6740 0.1880]*0.7);
        hold on
        
        %[pdf2,xmesh2,~] = ksdensity(x2,'NumPoints',n,'Function','pdf');
        %plot(xmesh2,pdf2*sum(y2)*h2.BinWidth,':','LineWidth',1.5,'color',[0.5 0.5 0.5]*0.5);
        %hold on
    end
    
    le = legend('Observed','Bootstrap',...
        'Location','northeast','fontsize',7);
    le.ItemTokenSize = ones(1,2)*8;
    
    
    tempTxt = sprintf('');
    if temp_p < 0.001
        tempTxt = sprintf('***');
    elseif temp_p < 0.01
        tempTxt = sprintf('**');
    elseif temp_p < 0.05
        tempTxt = sprintf('*');
    end
    text((x_min+x_max)*0.3,y_min+(y_max-y_min)*0.65,tempTxt,'Color','black','FontSize',12,'FontWeight','bold',...
        'HorizontalAlignment','center');
    
    
    set(gca,'linewidth',1.5)
    %xlim([x_min-(x_max-x_min)*0.1 x_max+(x_max-x_min)*0.1]);
    xlim([x_min-(x_max-x_min)*0.02 x_max+(x_max-x_min)*0.02]);
    ylim([y_min y_max+(y_max-y_min)*0.1]);%0.1
    %xticks([0 1]);
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Nearest-neighbour angle', 'FontSize', 9, 'FontWeight', 'normal');
    %ylabel('Frequency', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Pdf', 'FontSize', 9, 'FontWeight', 'normal');
    
    temp_title = title(sprintf('Distribution tests: ePAIRS'),'FontSize',9,'FontWeight','normal');
    temp_title.Interpreter = 'none';
    
end


%% Tuning clustering
if if_ePAIRS == 1
    %x = beta_memoryPrecision;
    %y = beta_choiceMemory;
    %z = beta_choiceMemory_baseline;
    
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_memoryPrecision;
        y = beta_choiceMemory;
        z = beta_choiceMemory_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_seqPrecision;
        y = beta_gProb;
        z = beta_gProb_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_seqPrecision;
        y = beta_gProb;
        z = beta_choiceMemory_baseline;
    end
    
    temp_vectors = [x,y,z];
    
    
    KList = 1:6;
    % criterion = 'CalinskiHarabasz';
    % criterion = 'DaviesBouldin';
    criterion = 'gap';
    % criterion = 'silhouette';
    
    % temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList);
    temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','sqEuclidean');
    %temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','cosine');
    % temp_eva = evalclusters(temp_vectors,'gmdistribution',criterion,'KList',KList,'Distance','sqEuclidean');
    
    CriterionValues_tuning = temp_eva.CriterionValues;
    clusterID_tuning = temp_eva.OptimalY;
    clusterNum_tuning = temp_eva.OptimalK;
    
    
    multi_rgbColor = ...
        [228,26,28;
        55,126,184;
        77,175,74;
        152,78,163;
        255,127,0;
        [255,255,51].*0.8]/255;
    
    
    
    %% Plot clusterNum vs CriterionValues
    fig = figure('Name',' ','NumberTitle','off');
    %set(gcf,'Position',[600 520 264 260*0.92]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[600 520 308 200]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[600 520 308*0.62*1.005 200*0.44*1.05*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile
    
    h_line = [];
    
    x = KList;
    y = CriterionValues_tuning;
    h_line = [h_line plot(x,y,'LineWidth',1,'color',[0.25 0.25 0.25])];
    hold on
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    
    plot([1 1].*clusterNum_tuning,[y_min-(y_max-y_min)*0.8 CriterionValues_tuning(clusterNum_tuning)],':','color',[0.25 0.25 0.25],'linewidth',1);
    hold on
    
    plot([x_min-(x_max-x_min)*0.08 clusterNum_tuning],[1 1].*CriterionValues_tuning(clusterNum_tuning),':','color',[0.25 0.25 0.25],'linewidth',1);
    hold on
    
    xticks(KList);
    
    set(gca,'linewidth',1.5)
    xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
    ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Cluster number', 'FontSize', 9, 'FontWeight', 'normal');
    %ylabel('Clustering scores', 'FontSize', 9, 'FontWeight', 'normal');
    ylabel('Scores', 'FontSize', 9, 'FontWeight', 'normal');
    temp_title = title(sprintf('Optimal cluster number'),'FontSize',9);
    temp_title.Interpreter = 'none';
    
    
    %% Single neuron tuning property of memoryPrecision and choiceMemory (trial-level tuning)
    fig = figure('Name','tuning property (trial-level)','NumberTitle','off');
    set(gcf,'Position',[510 450 240*3 240*1.05]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    
    %t.Title.String = sprintf('%s (trial-level)\n Single neuron tuning property, n=%d',FOVName_currentFOV2,length(x));
    t.Title.String = sprintf('%s (%s)\n Single neuron tuning property, n=%d',FOVName_currentFOV2,tempStr,length(x));
    t.Title.FontSize = 11;
    t.Title.Interpreter = 'none';
    
    tempSize = 4;%10
    
    nexttile
    
    %x = beta_memoryPrecision;
    %y = beta_choiceMemory;
    
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_memoryPrecision;
        y = beta_choiceMemory;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_seqPrecision;
        y = beta_gProb;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_seqPrecision;
        y = beta_gProb;
    end
    
    
    [temp_r,temp_p] = corr(x,y);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
    xy_min = -1*xy_max;
    
    %scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    %hold on
    
    for tempi=1:clusterNum_tuning
        temptempBoolIndex = clusterID_tuning==tempi;
        temptempColor = multi_rgbColor(tempi,:);
        scatter(x(temptempBoolIndex),y(temptempBoolIndex),tempSize,...
            'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
    end
    
    
    plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    
    
    %axis equal
    
    set(gca,'linewidth',1.5)
    %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
    %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
    
    xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Memory precision', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title.Interpreter = 'none';
    
    
    nexttile
    %x = beta_memoryPrecision;
    %y = beta_choiceMemory_baseline;
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_memoryPrecision;
        y = beta_choiceMemory_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_seqPrecision;
        y = beta_gProb_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_seqPrecision;
        y = beta_choiceMemory_baseline;
    end
    
    
    [temp_r,temp_p] = corr(x,y);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
    xy_min = -1*xy_max;
    
    
    %scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    %hold on
    
    for tempi=1:clusterNum_tuning
        temptempBoolIndex = clusterID_tuning==tempi;
        temptempColor = multi_rgbColor(tempi,:);
        scatter(x(temptempBoolIndex),y(temptempBoolIndex),tempSize,...
            'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
    end
    
    plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    
    %axis equal
    
    set(gca,'linewidth',1.5)
    
    %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
    %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
    
    xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Memory precision', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Baseline', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title.Interpreter = 'none';
    
    nexttile
    
    %x = beta_choiceMemory;
    %y = beta_choiceMemory_baseline;
    
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_choiceMemory;
        y = beta_choiceMemory_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_gProb;
        y = beta_gProb_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_gProb;
        y = beta_choiceMemory_baseline;
    end
    
    
    [temp_r,temp_p] = corr(x,y);
    
    [x_min,x_max] = bounds(x);
    [y_min,y_max] = bounds(y);
    xy_max = max([abs(x_min) abs(x_max) abs(y_min) abs(y_max)]);
    xy_min = -1*xy_max;
    
    
    %scatter(x,y,tempSize,'filled','MarkerFaceColor',[0.25 0.25 0.25],'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
    %hold on
    
    for tempi=1:clusterNum_tuning
        temptempBoolIndex = clusterID_tuning==tempi;
        temptempColor = multi_rgbColor(tempi,:);
        scatter(x(temptempBoolIndex),y(temptempBoolIndex),tempSize,...
            'filled','MarkerFaceColor',temptempColor,'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
    end
    
    plot([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],[0 0],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    plot([0 0],[xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04],'--','LineWidth',1,'color',[0.5 0.5 0.5]);
    hold on
    
    %axis equal
    
    set(gca,'linewidth',1.5)
    %xlim([x_min-(x_max-x_min)*0.08 x_max+(x_max-x_min)*0.08]);
    %ylim([y_min-(y_max-y_min)*0.08 y_max+(y_max-y_min)*0.08]);
    
    xlim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    ylim([xy_min-(xy_max-xy_min)*0.04 xy_max+(xy_max-xy_min)*0.04]);
    
    set(gca, 'FontSize', 12)
    set(gca,'box','off');% 取消右、上边框
    xlabel('Meta-memory', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Baseline', 'FontSize', 12, 'FontWeight', 'bold');
    temp_title = title(sprintf('r=%.3f, p=%.3f',temp_r,temp_p), 'FontSize', 12, 'FontWeight', 'bold');
    temp_title.Interpreter = 'none';
    
end

if false
    if if_monkey_D0_Z1 == 0 && currentSessionIndex_AB == 8
        beta_6loc;
        
        
        %if if_ePAIRS == 1
        if true
            if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
                x = beta_memoryPrecision;
                y = beta_choiceMemory;
                z = beta_choiceMemory_baseline;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
                x = beta_seqPrecision;
                y = beta_gProb;
                z = beta_gProb_baseline;
            elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
                x = beta_seqPrecision;
                y = beta_gProb;
                z = beta_choiceMemory_baseline;
            end
            
            %l = beta_6loc;
            l = std_beta_6loc;
            
            temp_vectors = [x,y,z,l];
            [clusteriness, p, dists] = pairsClusterTest_elliptical(temp_vectors);
            
            temp_angles = dists.data;
            temp_angles_bootstrap = dists.bootstrap;
            % temp_angles_bootstrap_median = median(temp_angles_bootstrap,1)';
            %temp_angles_bootstrap_median = median(temp_angles_bootstrap,2);
            % temp_angles_bootstrap_median = temp_angles_bootstrap(1,:)';
            temp_angles_bootstrap_median = reshape(temp_angles_bootstrap,[],1);
            
            
            %temp_vectors;
            %temp_vectors = [x,y,z,l];
            temp_vectors = [x,y,z];
            
            KList = 1:6;%1:6
            % criterion = 'CalinskiHarabasz';
            % criterion = 'DaviesBouldin';
            criterion = 'gap';
            % criterion = 'silhouette';
            
            % temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList);
            temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','sqEuclidean');
            %temp_eva = evalclusters(temp_vectors,'kmeans',criterion,'KList',KList,'Distance','cosine');
            % temp_eva = evalclusters(temp_vectors,'gmdistribution',criterion,'KList',KList,'Distance','sqEuclidean');
            
            CriterionValues_tuning = temp_eva.CriterionValues;
            clusterID_tuning = temp_eva.OptimalY;
            clusterNum_tuning = temp_eva.OptimalK;
            
        end
    end
end


%% PCA
if false
    if if_plotBeta_trialTuning0_seqTuning1_mix2 == 0
        x = beta_memoryPrecision;
        y = beta_choiceMemory;
        z = beta_choiceMemory_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 1
        x = beta_seqPrecision;
        y = beta_gProb;
        z = beta_gProb_baseline;
    elseif if_plotBeta_trialTuning0_seqTuning1_mix2 == 2
        x = beta_seqPrecision;
        y = beta_gProb;
        z = beta_choiceMemory_baseline;
    end
    
    %     temp1 = max(beta_6loc,[],2);
    %     temp2 = min(beta_6loc,[],2);
    %     temp3 = nan(length(beta_6loc),1);
    %     for tempi=1:length(beta_6loc)
    %         if abs(temp1(tempi)) > abs(temp2(tempi))
    %             temp3(tempi) = temp1(tempi);
    %         else
    %             temp3(tempi) = temp2(tempi);
    %         end
    %     end
    %     beta_6loc_peak = temp3;
    %
    %     l = beta_6loc_peak;
    
    %     temp_vectors = [x,y,z,l];
    %     [coeff4,score4,latent4,tsquared4,explained4,mu4] = pca(temp_vectors);
    
    temp_vectors = [x,y,z];
    [coeff3,score3,latent3,tsquared3,explained3,mu3] = pca(temp_vectors);
    
    
    temp_vectors = [x,y];
    [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(temp_vectors);
    
end




%% Correlation between memory strength and behavior across seqs, using all trials (correct + error)
memoryPrecision_trialLevel;
seqIndex;
if_memoryPrecision_accuracy0_sigma1;

temp_seqPrecision_neuron = nan(sum(numSeq(1:valid_length)),1);
for tempi=1:sum(numSeq(1:valid_length))
    temptempBoolIndex = (seqIndex == tempi) & (~isnan(memoryPrecision_trialLevel))';
    
    temp1 = memoryPrecision_trialLevel(temptempBoolIndex);
    temp_seqPrecision_neuron(tempi) = mean(temp1,'omitnan');
end


if true
    fig = figure('Name','All trials','NumberTitle','off'); %#ok<*NASGU>
    set(gcf,'Position',[50+0 50+0 226*0.75*0.78*0.94*1.25 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    temp_x = temp_seqPrecision_neuron;
    %temp_y = gAcc_target_collapsed_inOne(1:41)';
    temp_y = seqPrecision_behavior;
    tempBoolIndex = ~isnan(temp_x);
    
    x = temp_x(~isnan(temp_x));
    y = temp_y(1:sum(numSeq(1:valid_length)));
    y = y(~isnan(temp_x));
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    h = [];
    for tempi=1:valid_length
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        temp_size = 10;
        temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h];
    end
    
    %[temp_xmin,temp_xmax] = bounds(x);
    %[temp_ymin,temp_ymax] = bounds(y);
    
    temp_xmin = 0;
    temp_xmax = 1;
    temp_ymin = 0;
    temp_ymax = 1;
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    h = [h temp_h];
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    text(temp_xmin+(temp_xmax-temp_xmin)*0.55,0.22,sprintf('r = %.3f',r), 'fontsize',9,'FontWeight','normal');
    text(temp_xmin+(temp_xmax-temp_xmin)*0.55,0.07,sprintf('%s',tempTxt), 'fontsize',9,'FontWeight','normal');
    
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    
    xticks([0:0.2:1]);
    yticks([0:0.2:1]);
    
    xtickangle(0);
    
    if if_memoryPrecision_accuracy0_sigma1 == 0
        xlabel(sprintf('Seq probability'), 'FontSize', 9, 'FontWeight', 'normal');
    elseif if_memoryPrecision_accuracy0_sigma1 == 1
        xlabel(sprintf('Memory strength'), 'FontSize', 9, 'FontWeight', 'normal');
    end
    
    temp_title = title(sprintf('Across seqs'), 'FontSize', 9, 'FontWeight', 'normal');
    temp_title.Interpreter = 'none';
    
    if if_memoryPrecision_accuracy0_sigma1_forBehavior == 0
        ylabel(sprintf('Recall accuracy'), 'FontSize', 9, 'FontWeight', 'normal');
    elseif if_memoryPrecision_accuracy0_sigma1_forBehavior == 1
        ylabel(sprintf('Recall variability'), 'FontSize', 9, 'FontWeight', 'normal');
    end
    
end



%% Correlation between memory strength and behavior across seqs, using all noChoice trials (correct + error)
memoryPrecision_trialLevel;
seqIndex;
if_memoryPrecision_accuracy0_sigma1;

temp2_seqPrecision_neuron = nan(sum(numSeq(1:valid_length)),1);
for tempi=1:sum(numSeq(1:valid_length))
    temptempBoolIndex = (seqIndex == tempi) & (~isnan(memoryPrecision_trialLevel))' & (~choiceBoolIndex);
    
    temp1 = memoryPrecision_trialLevel(temptempBoolIndex);
    temp2_seqPrecision_neuron(tempi) = mean(temp1,'omitnan');
end


%% Plot correlation between WM strength and Recall variability, in noChoice trials
if true
    fig = figure('Name','All noChoice trials','NumberTitle','off'); %#ok<*NASGU>
    %set(gcf,'Position',[250+0 50+0 226*0.75*0.78*0.94*1.25 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[250+0 50+0 226*0.75*0.78*0.94*1.25*0.87 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[250+0 50+0 226*0.75*0.78*0.94*1.25*0.87*0.9*1.03 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    temp_x = temp2_seqPrecision_neuron;
    %temp_y = gAcc_target_collapsed_inOne(1:41)';
    temp_y = seqPrecision_behavior;
    tempBoolIndex = ~isnan(temp_x);
    
    x = temp_x(~isnan(temp_x));
    y = temp_y(1:sum(numSeq(1:valid_length)));
    y = y(~isnan(temp_x));
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    h = [];
    for tempi=1:valid_length
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        %temp_size = 10;      
        temp_size = (((tempi+1).^3)*2 + 3) ./ 5;
                        
        temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h];
    end
    
    
    %[temp_xmin,temp_xmax] = bounds(x);
    %[temp_ymin,temp_ymax] = bounds(y);
    
    temp_xmin = 0;
    temp_xmax = 1;
    temp_ymin = 0;
    temp_ymax = 1;
    
    if if_monkey_D0_Z1 == 1 && currentSessionIndex_AB == 9
        temp_xmin = 0;
        temp_xmax = 0.4;
    end
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    h = [h temp_h];
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    if if_entropy_behaviorInverted == 1
        text(temp_xmin+(temp_xmax-temp_xmin)*0.46,0.22,sprintf('r = %.3f',r), 'fontsize',8,'FontWeight','normal');
        text(temp_xmin+(temp_xmax-temp_xmin)*0.46,0.07,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
    elseif if_entropy_behaviorInverted == 0
        text(temp_xmin+(temp_xmax-temp_xmin)*0.46,0.95,sprintf('r = %.3f',r), 'fontsize',8,'FontWeight','normal');
        text(temp_xmin+(temp_xmax-temp_xmin)*0.46,0.80,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
    end
    
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    
    %xticks([0:0.5:1]);
    xticks([0:1]);
    %yticks([0:0.2:1]);
    yticks([0:1]);
    
    if if_monkey_D0_Z1 == 1 && currentSessionIndex_AB == 9
        xticks([temp_xmin,temp_xmax]);
    end
    
    xtickangle(0);
    
    if if_memoryPrecision_accuracy0_sigma1 == 0
        xlabel(sprintf('Seq probability'), 'FontSize', 9, 'FontWeight', 'normal');
    elseif if_memoryPrecision_accuracy0_sigma1 == 1
        xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    end
    
    if if_memoryPrecision_accuracy0_sigma1_forBehavior == 0
        ylabel(sprintf('Recall accuracy'), 'FontSize', 9, 'FontWeight', 'normal');
    elseif if_memoryPrecision_accuracy0_sigma1_forBehavior == 1
        ylabel(sprintf('Recall variability'), 'FontSize', 9, 'FontWeight', 'normal');
    end
    
    temp_title = title(sprintf('Across seqs'), 'FontSize', 9, 'FontWeight', 'normal');
    temp_title.Interpreter = 'none';
    
end


%% Plot correlation between WM strength and Recall variability, in noChoice trials, without outliers

temptempBoolIndex = isoutlier(temp2_seqPrecision_neuron,'median');
% temptempBoolIndex = isoutlier(temp2_seqPrecision_neuron,'mean');
% temptempBoolIndex = isoutlier(temp2_seqPrecision_neuron,'quartiles');
% temptempBoolIndex = isoutlier(temp2_seqPrecision_neuron,'grubbs');
% temptempBoolIndex = isoutlier(temp2_seqPrecision_neuron,'gesd');

% find(temptempBoolIndex');

if true
    fig = figure('Name','All noChoice trials','NumberTitle','off'); %#ok<*NASGU>
    set(gcf,'Position',[450+0 50+0 226*0.75*0.78*0.94*1.25*0.87*0.9*1.03 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    temp_x = temp2_seqPrecision_neuron;
    
    temp_x(temptempBoolIndex) = nan;
    
    temp_y = seqPrecision_behavior;
    tempBoolIndex = ~isnan(temp_x);
    
    x = temp_x(~isnan(temp_x));
    y = temp_y(1:sum(numSeq(1:valid_length)));
    y = temp_y(~isnan(temp_x));
    
    [r,p] = corr(x,y);
    %p
    
    mdl = fitglm(x,y);
    
    
    h = [];
    for tempi=1:valid_length
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        %temp_size = 10;
        temp_size = (((tempi+1).^3)*2 + 3) ./ 5;
        temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h];
    end
    
    %[temp_xmin,temp_xmax] = bounds(x);
    %[temp_ymin,temp_ymax] = bounds(y);
    
    temp_xmin = 0;
    temp_xmax = 1;
    temp_ymin = 0;
    temp_ymax = 1;
    
    if if_monkey_D0_Z1 == 1 && currentSessionIndex_AB == 9
        %temp_xmin = 0;
        %temp_xmax = 0.4;%0.4
        
        [temp_xmin,temp_xmax] = bounds(x);
        %[temp_ymin,temp_ymax] = bounds(y);
    end
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    h = [h temp_h];
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    if if_entropy_behaviorInverted == 1
        text(temp_xmin+(temp_xmax-temp_xmin)*0.46,0.22,sprintf('r = %.3f',r), 'fontsize',8,'FontWeight','normal');
        text(temp_xmin+(temp_xmax-temp_xmin)*0.46,0.07,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
    elseif if_entropy_behaviorInverted == 0
        if if_monkey_D0_Z1 == 1 && currentSessionIndex_AB == 9
            text(temp_xmin+(temp_xmax-temp_xmin)*0.39,temp_ymax+(temp_ymax-temp_ymin)*0.145,sprintf('r = %.3f',r), 'fontsize',8,'FontWeight','normal');
            text(temp_xmin+(temp_xmax-temp_xmin)*0.515,temp_ymax-(temp_ymax-temp_ymin)*0.03,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');            
        else
            text(temp_xmin+(temp_xmax-temp_xmin)*0.46,0.95,sprintf('r = %.3f',r), 'fontsize',8,'FontWeight','normal');
            text(temp_xmin+(temp_xmax-temp_xmin)*0.46,0.80,sprintf('%s',tempTxt), 'fontsize',8,'FontWeight','normal');
        end
    end
    
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
    
    if if_monkey_D0_Z1 == 1 && currentSessionIndex_AB == 9
        xlim([temp_xmin-(temp_xmax-temp_xmin)*0.25 temp_xmax+(temp_xmax-temp_xmin)*0.25]);
        ylim([temp_ymin-(temp_ymax-temp_ymin)*0.25 temp_ymax+(temp_ymax-temp_ymin)*0.25]);
    end
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    
    %xticks([0:0.5:1]);
    %xticks([0:1]);
    %yticks([0:0.2:1]);
    %yticks([0:1]);
    
    if if_monkey_D0_Z1 == 1 && currentSessionIndex_AB == 9
        %xticks([temp_xmin,temp_xmax]);
        yticks([0:1]);
    else
        xticks([0:1]);
        yticks([0:1]);
    end
    
    xtickangle(0);
    
    if if_memoryPrecision_accuracy0_sigma1 == 0
        xlabel(sprintf('Seq probability'), 'FontSize', 9, 'FontWeight', 'normal');
    elseif if_memoryPrecision_accuracy0_sigma1 == 1
        xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    end
    
    if if_memoryPrecision_accuracy0_sigma1_forBehavior == 0
        ylabel(sprintf('Recall accuracy'), 'FontSize', 9, 'FontWeight', 'normal');
    elseif if_memoryPrecision_accuracy0_sigma1_forBehavior == 1
        ylabel(sprintf('Recall variability'), 'FontSize', 9, 'FontWeight', 'normal');
    end
    
    temp_title = title(sprintf('Across seqs'), 'FontSize', 9, 'FontWeight', 'normal');
    temp_title.Interpreter = 'none';
    
end

%% Plot correlation between WM strength and offloading rate, in noChoice trials
if true
    fig = figure('Name','All noChoice trials','NumberTitle','off'); %#ok<*NASGU>
    %set(gcf,'Position',[250+0 50+0 226*0.75*0.78*0.94*1.25 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    %set(gcf,'Position',[250+0 50+0 226*0.75*0.78*0.94*1.25*0.87 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    set(gcf,'Position',[450+0 50+0 226*0.75*0.78*0.94*1.25*0.87*0.9*1.03 186*0.78*0.94]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    
    temp_x = temp2_seqPrecision_neuron;
    %temp_y = gAcc_target_collapsed_inOne(1:41)';
    %temp_y = seqPrecision_behavior;
    temp_y = offloadingProb_inOne(1:41)';
    tempBoolIndex = ~isnan(temp_x);
    
    x = temp_x(~isnan(temp_x));
    y = temp_y(1:sum(numSeq(1:valid_length)));
    y = y(~isnan(temp_x));
    
    [r,p] = corr(x,y);
    
    mdl = fitglm(x,y);
    
    
    h = [];
    for tempi=1:valid_length
        temp_range = (sum(numSeq(1:tempi-1))+1):sum(numSeq(1:tempi));
        tempIndex = find(tempBoolIndex(temp_range)==true);
        temp_range2 = temp_range(tempIndex);
        
        temp_size = 10;
        temp_h = scatter(temp_x(temp_range2), temp_y(temp_range2), ...
            temp_size, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);
        hold on
        h = [h temp_h];
    end
    
    %[temp_xmin,temp_xmax] = bounds(x);
    %[temp_ymin,temp_ymax] = bounds(y);
    
    temp_xmin = 0;
    temp_xmax = 1;
    temp_ymin = 0;
    temp_ymax = 1;
    
    if if_monkey_D0_Z1 == 1 && currentSessionIndex_AB == 9
        temp_xmin = 0;
        temp_xmax = 0.4;
    end
    
    x_fit = temp_xmin:0.001:temp_xmax;
    y_fit = predict(mdl,x_fit')';
    
    temp_h = plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.35 0.35 0.35 0.7]);
    hold on
    h = [h temp_h];
    
    
    tempTxt = 'p > 0.05';
    if p < 0.05
        tempTxt = 'p < 0.05';
    end
    if p < 0.01
        tempTxt = 'p < 0.01';
    end
    if p < 0.001
        tempTxt = 'p < 0.001';
    end
    text(temp_xmin+(temp_xmax-temp_xmin)*0.6,1.05,sprintf('r = %.3f',r), 'fontsize',6.5,'FontWeight','normal');
    text(temp_xmin+(temp_xmax-temp_xmin)*0.6,0.92,sprintf('%s',tempTxt), 'fontsize',6.5,'FontWeight','normal');
    
    
    xlim([temp_xmin-(temp_xmax-temp_xmin)*0.1 temp_xmax+(temp_xmax-temp_xmin)*0.1]);
    ylim([temp_ymin-(temp_ymax-temp_ymin)*0.1 temp_ymax+(temp_ymax-temp_ymin)*0.1]);
    
    set(gca,'linewidth',1.5)
    set(gca, 'FontSize', 8)
    set(gca,'box','off');% 取消右、上边框
    
    %xticks([0:0.5:1]);
    xticks([0:0.5:1]);
    %yticks([0:0.2:1]);
    yticks([0:1]);
    
    if if_monkey_D0_Z1 == 1 && currentSessionIndex_AB == 9
        xticks([temp_xmin,temp_xmax]);
    end
    
    xtickangle(0);
    
    if if_memoryPrecision_accuracy0_sigma1 == 0
        xlabel(sprintf('Seq probability'), 'FontSize', 9, 'FontWeight', 'normal');
    elseif if_memoryPrecision_accuracy0_sigma1 == 1
        xlabel(sprintf('WM strength'), 'FontSize', 9, 'FontWeight', 'normal');
    end
    
    ylabel(sprintf('Offloading rate'), 'FontSize', 9, 'FontWeight', 'normal');
    
    
    temp_title = title(sprintf('Across seqs'), 'FontSize', 9, 'FontWeight', 'normal');
    temp_title.Interpreter = 'none';
    
end

a = 1;


%% Angle between tuning strength vecters
if true
    % if false
    
    if_fitglm0_lasso1 = 0;%0
    temp_resampleNum = 16;%16
    
    fprintf('(Unit view of angle) if_fitglm0_lasso1=%d,temp_resampleNum=%d.\n',...
        if_fitglm0_lasso1,temp_resampleNum);
    
    t_angle = tic;
    beta_6loc;
    beta_gProb;
    beta_choiceMemory_baseline;
    
    lasso_linearAxis_jjb_Name_v = autoGetFunName_myScripts('lasso_linearAxis_jjb', [targetPATH '\functions']);
    fun_lasso_linearAxis_jjb = str2func(lasso_linearAxis_jjb_Name_v);
    
    
    beta_6loc_norm = nan(size(beta_6loc,1),1);
    for tempi=1:size(beta_6loc,1)
        temp1 = beta_6loc(tempi,:);
        beta_6loc_norm(tempi) = norm(temp1);
    end
    
    beta_gProb_norm = abs(beta_gProb);
    beta_choiceMemory_baseline_norm = abs(beta_choiceMemory_baseline);
    
    
    temp_vectors = [beta_6loc_norm,beta_gProb_norm,beta_choiceMemory_baseline_norm];
    
    a = temp_vectors(:,1);
    b = temp_vectors(:,2);
    c = temp_vectors(:,3);
    
    temp_degree_ab = subspace(a,b)*180/pi;
    temp_degree_ac = subspace(a,c)*180/pi;
    temp_degree_bc = subspace(b,c)*180/pi;
    
    
    %%% Compute shuffle level of degree within features
    
    %% Resample beta_6loc in no choice correct trials (from half trials)
    x_raw = boolIndex_location_seq_T(seqIndex,:);
    x_raw2 = x_raw(allMemoryCorrectBoolIndex&(~choiceBoolIndex),:);
    x = (x_raw2-mean(x_raw2,1))./std(x_raw2,0,1);
    
    y_raw = F_dff_decisionBin1_zScore';
    y_raw2 = y_raw(allMemoryCorrectBoolIndex&(~choiceBoolIndex),:);
    y = (y_raw2-mean(y_raw2,1))./std(y_raw2,0,1);
    
    
    beta_6loc_halfA_norm_resampled = cell(temp_resampleNum,1);
    beta_6loc_halfB_norm_resampled = cell(temp_resampleNum,1);
    
    %for tempIter=1:temp_resampleNum
    parfor tempIter=1:temp_resampleNum
        
        temptemp_trialNum = size(x,1);
        temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
        
        temptempBoolIndex = false(temptemp_trialNum,1);
        temptempBoolIndex(temptempIndex) = true;
        temptempBoolIndex_A = temptempBoolIndex;
        temptempBoolIndex_B = ~temptempBoolIndex;
        
        beta_A = zeros(roiNum,numFrames);
        beta_B = zeros(roiNum,numFrames);
        
        %parfor tempi=1:roiNum
        for tempi=1:roiNum
            warning('off');
                        
            if if_fitglm0_lasso1 == 0
                temp_mdl_A = fitglm(x(temptempBoolIndex_A,:),y(temptempBoolIndex_A,tempi));
                beta_A(tempi,:) = temp_mdl_A.Coefficients.Estimate(2:end);
                temp_mdl_B = fitglm(x(temptempBoolIndex_B,:),y(temptempBoolIndex_B,tempi));
                beta_B(tempi,:) = temp_mdl_B.Coefficients.Estimate(2:end);
            elseif if_fitglm0_lasso1 == 1
                [beta_A(tempi,:),~] = fun_lasso_linearAxis_jjb(x(temptempBoolIndex_A,:),y(temptempBoolIndex_A,tempi));
                [beta_B(tempi,:),~] = fun_lasso_linearAxis_jjb(x(temptempBoolIndex_B,:),y(temptempBoolIndex_B,tempi));
            end            
            warning('on');
        end
        temp_beta_6loc_halfA = beta_A;
        temp_beta_6loc_halfB = beta_B;
        
        
        temp_beta_6loc_halfA_norm = nan(size(beta_6loc,1),1);
        temp_beta_6loc_halfB_norm = nan(size(beta_6loc,1),1);
        for tempi=1:size(beta_6loc,1)
            temp1 = temp_beta_6loc_halfA(tempi,:);
            temp2 = temp_beta_6loc_halfB(tempi,:);
            
            temp_beta_6loc_halfA_norm(tempi) = norm(temp1);
            temp_beta_6loc_halfB_norm(tempi) = norm(temp2);
        end
        
        
        beta_6loc_halfA_norm_resampled{tempIter} = temp_beta_6loc_halfA_norm;
        beta_6loc_halfB_norm_resampled{tempIter} = temp_beta_6loc_halfB_norm;
        
    end
    
    beta_6loc_halfA_norm_resampled;
    beta_6loc_halfB_norm_resampled;
    
    betaAngle_6loc_halfAB_norm_resampled = nan(temp_resampleNum,1);
    for tempi=1:temp_resampleNum
        a = beta_6loc_halfA_norm_resampled{tempi};
        b = beta_6loc_halfB_norm_resampled{tempi};
        betaAngle_6loc_halfAB_norm_resampled(tempi) = subspace(a,b)*180/pi;
    end
    
    
    
    %% Resample beta_gProb in all trials (from half trials)
    
    beta_gProb_halfA_norm_resampled = cell(temp_resampleNum,1);
    beta_gProb_halfB_norm_resampled = cell(temp_resampleNum,1);
    
    %for tempIter=1:temp_resampleNum
    parfor tempIter=1:temp_resampleNum
        
        temptemp_trialNum = size(F_dff_decisionBin1_zScore,2);
        temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
        
        temptempBoolIndex = false(temptemp_trialNum,1);
        temptempBoolIndex(temptempIndex) = true;
        temptempBoolIndex_A = temptempBoolIndex;
        temptempBoolIndex_B = ~temptempBoolIndex;
        
        temp_dff_seqMerged_A = nan(roiNum,sum(numSeq));
        temp_dff_seqMerged_B = nan(roiNum,sum(numSeq));
        for tempi=1:sum(numSeq)
            
            temp_dff = F_dff_decisionBin1_zScore(:,(seqIndex==tempi) & temptempBoolIndex_A');
            temp_dff_seqMerged_A(:,tempi) = mean(temp_dff,2);
            
            temp_dff = F_dff_decisionBin1_zScore(:,(seqIndex==tempi) & temptempBoolIndex_B');
            temp_dff_seqMerged_B(:,tempi) = mean(temp_dff,2);
        end
        
        
        beta_A = nan(roiNum,1);
        beta_B = nan(roiNum,1);
        for tempi=1:roiNum
            x = (1-offloadingProb_inOne)';
            y_A = temp_dff_seqMerged_A(tempi,:)';
            y_B = temp_dff_seqMerged_B(tempi,:)';
            
            if if_fitglm0_lasso1 == 0
                temp_mdl_A = fitglm(x,y_A);
                beta_A(tempi) = temp_mdl_A.Coefficients.Estimate(2);
                temp_mdl_B = fitglm(x,y_B);
                beta_B(tempi) = temp_mdl_B.Coefficients.Estimate(2);
            elseif if_fitglm0_lasso1 == 1
                [beta_A(tempi),~] = fun_lasso_linearAxis_jjb(x,y_A);
                [beta_B(tempi),~] = fun_lasso_linearAxis_jjb(x,y_B);
            end
        end
        
        temp_beta_gProb_halfA_norm = abs(beta_A);
        temp_beta_gProb_halfB_norm = abs(beta_B);
        
        beta_gProb_halfA_norm_resampled{tempIter} = temp_beta_gProb_halfA_norm;
        beta_gProb_halfB_norm_resampled{tempIter} = temp_beta_gProb_halfB_norm;
        
    end
    
    beta_gProb_halfA_norm_resampled;
    beta_gProb_halfB_norm_resampled;
    
    betaAngle_gProb_halfAB_norm_resampled = nan(temp_resampleNum,1);
    for tempi=1:temp_resampleNum
        a = beta_gProb_halfA_norm_resampled{tempi};
        b = beta_gProb_halfB_norm_resampled{tempi};
        betaAngle_gProb_halfAB_norm_resampled(tempi) = subspace(a,b)*180/pi;
    end
    
    
    %% Resample beta_choiceMemory_baseline in all trials (from half trials)
    
    x_raw = choiceMemoryBoolIndex';
    x = x_raw(choiceBoolIndex);
    y_raw = F_dff_baselineBin_zScore';
    y = y_raw(choiceBoolIndex,:);
    
    beta_choiceMemory_baseline_halfA_norm_resampled = cell(temp_resampleNum,1);
    beta_choiceMemory_baseline_halfB_norm_resampled = cell(temp_resampleNum,1);
    
    %for tempIter=1:temp_resampleNum
    parfor tempIter=1:temp_resampleNum
        
        temptemp_trialNum = size(x,1);
        temptempIndex = sort(randperm(temptemp_trialNum,round(temptemp_trialNum/2)));
        
        temptempBoolIndex = false(temptemp_trialNum,1);
        temptempBoolIndex(temptempIndex) = true;
        temptempBoolIndex_A = temptempBoolIndex;
        temptempBoolIndex_B = ~temptempBoolIndex;
        
        
        beta_A = zeros(roiNum,1);
        beta_B = zeros(roiNum,1);
        for tempi=1:roiNum     
            
            if if_fitglm0_lasso1 == 0
                temp_mdl_A = fitglm(x(temptempBoolIndex_A),y(temptempBoolIndex_A,tempi));
                beta_A(tempi) = temp_mdl_A.Coefficients.Estimate(2);
                temp_mdl_B = fitglm(x(temptempBoolIndex_B),y(temptempBoolIndex_B,tempi));
                beta_B(tempi) = temp_mdl_B.Coefficients.Estimate(2);
            elseif if_fitglm0_lasso1 == 1
                [beta_A(tempi),~] = fun_lasso_linearAxis_jjb(x(temptempBoolIndex_A),y(temptempBoolIndex_A,tempi));
                [beta_B(tempi),~] = fun_lasso_linearAxis_jjb(x(temptempBoolIndex_B),y(temptempBoolIndex_B,tempi));
            end
            
        end
        
        temp_beta_choiceMemory_baseline_halfA_norm = abs(beta_A);
        temp_beta_choiceMemory_baseline_halfB_norm = abs(beta_B);
        
        beta_choiceMemory_baseline_halfA_norm_resampled{tempIter} = temp_beta_choiceMemory_baseline_halfA_norm;
        beta_choiceMemory_baseline_halfB_norm_resampled{tempIter} = temp_beta_choiceMemory_baseline_halfB_norm;
        
    end
    
    beta_choiceMemory_baseline_halfA_norm_resampled;
    beta_choiceMemory_baseline_halfB_norm_resampled;
    
    betaAngle_choiceMemory_baseline_halfAB_norm_resampled = nan(temp_resampleNum,1);
    for tempi=1:temp_resampleNum
        a = beta_choiceMemory_baseline_halfA_norm_resampled{tempi};
        b = beta_choiceMemory_baseline_halfB_norm_resampled{tempi};
        betaAngle_choiceMemory_baseline_halfAB_norm_resampled(tempi) = subspace(a,b)*180/pi;
    end
    
    
    %% Summary of halfAB
    temp_degree_ab;
    temp_degree_ac;
    temp_degree_bc;
    
    betaAngle_6loc_halfAB_norm_resampled;
    betaAngle_gProb_halfAB_norm_resampled;
    betaAngle_choiceMemory_baseline_halfAB_norm_resampled;
    
    temp_degree_aa = mean(betaAngle_6loc_halfAB_norm_resampled);
    temp_degree_aa_std = std(betaAngle_6loc_halfAB_norm_resampled);
    
    temp_degree_bb = mean(betaAngle_gProb_halfAB_norm_resampled);
    temp_degree_bb_std = std(betaAngle_gProb_halfAB_norm_resampled);
    
    temp_degree_cc = mean(betaAngle_choiceMemory_baseline_halfAB_norm_resampled);
    temp_degree_cc_std = std(betaAngle_choiceMemory_baseline_halfAB_norm_resampled);
    
    
    %fprintf('Degree between beta vectors (unit view): \nbetween [%.1f, %.1f, %.1f], within [%.1f, %.1f, %.1f].\n',...
    %    temp_degree_ab,temp_degree_ac,temp_degree_bc,temp_degree_aa,temp_degree_bb,temp_degree_cc);
    fprintf('Degree between beta vectors (unit view): within [%.1f, %.1f, %.1f], between [%.1f, %.1f, %.1f].\n',...
       temp_degree_cc,temp_degree_aa,temp_degree_bb,temp_degree_ac,temp_degree_bc,temp_degree_ab);
    
    fprintf('t_angle = %.1f secs.\n',toc(t_angle));
    
    
    betaAngle_6loc_halfAB_norm_resampled_unitView = betaAngle_6loc_halfAB_norm_resampled;
    betaAngle_gProb_halfAB_norm_resampled_unitView = betaAngle_gProb_halfAB_norm_resampled;
    betaAngle_choiceMemory_baseline_halfAB_norm_resampled_unitView = betaAngle_choiceMemory_baseline_halfAB_norm_resampled;

end


%% End