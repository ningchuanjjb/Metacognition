% Chuan's 1st script (20251214)
% This script: To compute location & offloading rate selectivity of each neuron from multi-FOVs.
%% Initialization
close all


if_compute = 0;

if if_compute == 1
    clear
    if_compute = 1;
    if_load_beta = 0;
else
    if_load_beta = 1;
end

% if_load_beta = 1;

if_load_selectiveIndex = 1;

if_save = 0;


% if_monkey_D0_Z1 = 0;%1


if_glm_12B0_24B1_6B2 = 2;
if_singleFOV = 0;



order_glm = 0;
plot_lengthFlag = 0;%0 mean all length

if_load_allTrial0_memoryCorrect1 = 0;%0

if_strictHit = 1;
if_allSDT_withOnlyHit = 1;

if if_strictHit == 1
    if_allSDT_withOnlyHit = 1;
end

if if_allSDT_withOnlyHit == 1
    if_glm_12B0_24B1_6B2 = 2;
end

if_seqSubspace_onlyhit0_allSDT1 = 1;
if if_seqSubspace_onlyhit0_allSDT1 == 1
    order_glm = 0;
    if_load_allTrial0_memoryCorrect1 = 0;
end

if_seqTrialCountWeighted = 1;

%%
targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)


if if_load_allTrial0_memoryCorrect1 == 0
    searchName_glmData = 'allTrial.mat';
elseif if_load_allTrial0_memoryCorrect1 == 1
    searchName_glmData = 'memoryCorrect.mat';
end

if if_monkey_D0_Z1 == 0
    path_glm = 'C:\ASDROOT\STUDY\temp\data\Ding';
elseif if_monkey_D0_Z1 == 1
    path_glm = 'C:\ASDROOT\STUDY\temp\data\Zelku';
end
%path_glm = 'C:\ASDROOT\STUDY\temp\data';

numFrames = 6;

t0 = tic;

%% if_compute
if if_compute == 1
    if if_seqSubspace_onlyhit0_allSDT1 == 0
        
    elseif if_seqSubspace_onlyhit0_allSDT1 == 1
        if if_glm_12B0_24B1_6B2 == 1
            fprintf('Now is all SDT mode with 24B (P:Hit, 1-P:Miss, P0:False alarm, 1-P0:Correct rejection).\n');
        elseif if_glm_12B0_24B1_6B2 == 0
            fprintf('Now is all SDT mode with 12B (P:Hit, 1-P:Miss, P0:False alarm, 1-P0:Correct rejection).\n');
        elseif if_glm_12B0_24B1_6B2 == 2
            fprintf('Now is all SDT mode with 6B (P:Hit, 1-P:Miss, P0:False alarm, 1-P0:Correct rejection).\n');
        end
    end
    
    if if_strictHit == 0
        fprintf('Loose hit.\n');
    elseif if_strictHit == 1
        fprintf('Strict hit (only work in all SDT with only Hit).\n');
    end
    if if_allSDT_withOnlyHit == 1
        fprintf('All SDT with only Hit (only work in 6B mode).\n');
    end
    
    % Load all sessions' data files
    [glm_file_load, glm_file_name] = load_glmData_multiFOV_v1(searchName_glmData,path_glm);
    fileSize = length(glm_file_name);
    
    
    glm_beta_lengthx_delay1Bin_multiFOV = [];
    glm_beta_lengthx_delay2Bin_multiFOV = [];
    glm_beta_lengthx_baselineBin_multiFOV = [];
    
    glm_r2_lengthx_delay1Bin_multiFOV = [];
    glm_r2_lengthx_delay2Bin_multiFOV = [];
    glm_r2_lengthx_baselineBin_multiFOV = [];
    
    lasso_jjb_Name_v = autoGetFunName_myScripts('lasso_jjb', [targetPATH '\functions']);
    fun_lasso_jjb = str2func(lasso_jjb_Name_v);
    lasso_repeat_jjb_Name_v = autoGetFunName_myScripts('lasso_repeat_jjb', [targetPATH '\functions']);
    fun_lasso_repeat_jjb = str2func(lasso_repeat_jjb_Name_v);
    neuronGLM_eachLength_par_Name_v = autoGetFunName_myScripts('fun_glm_neuron_eachLength_par', [targetPATH '\functions']);
    fun_neuronGLM_eachLength_par = str2func(neuronGLM_eachLength_par_Name_v);
    glm_preparation_Name_v = autoGetFunName_myScripts('fun_glm_preparation', [targetPATH '\functions']);
    fun_glm_preparation = str2func(glm_preparation_Name_v);
    glm_dataLoad_Name_v = autoGetFunName_myScripts('fun_glm_dataLoad', [targetPATH '\functions']);
    fun_glm_dataLoad = str2func(glm_dataLoad_Name_v);
    
    
    %fileSize = 3;
    for tempIndex=1:fileSize
        %if true
        %tempIndex = 20;
        % tempIndex = 1;
        
        if if_singleFOV == 1
            tempIndex = 20; %#ok<*FXSET>
        end
        
        tempstr = sprintf('file%d',tempIndex);
        tempData = glm_file_load.(tempstr).glmData;
        
        a = 1;
        
        %% fun_glm_dataLoad
        glm_dataLoad_options = struct;
        glm_dataLoad_options.tempData = tempData;
        glm_dataLoad_options.if_load_allTrial0_memoryCorrect1 = if_load_allTrial0_memoryCorrect1;
        glm_dataLoad_options.plot_lengthFlag = plot_lengthFlag;
        
        glm_dataLoad_output = fun_glm_dataLoad(glm_dataLoad_options);
        
        F_dff_lengthx_delay1Bin = glm_dataLoad_output.F_dff_lengthx_delay1Bin;
        F_dff_lengthx_delay2Bin = glm_dataLoad_output.F_dff_lengthx_delay2Bin;
        F_dff_lengthx_baselineBin = glm_dataLoad_output.F_dff_lengthx_baselineBin;
        sequence_lengthx_onehot = glm_dataLoad_output.sequence_lengthx_onehot;
        responseSeq_lengthx_onehot = glm_dataLoad_output.responseSeq_lengthx_onehot;
        
        a = 1;
        
        %% fun_glm_preparation
        glm_prep_options = struct;
        glm_prep_options.order_glm = order_glm;
        glm_prep_options.plot_lengthFlag = plot_lengthFlag;
        glm_prep_options.sequence_lengthx_onehot = sequence_lengthx_onehot;
        glm_prep_options.numFrames = numFrames;
        
        glm_prep = fun_glm_preparation(glm_prep_options);
        
        order_glm_valid = glm_prep.order_glm_valid;
        sequence_lengthx_onehot_oneOrder = glm_prep.sequence_lengthx_onehot_oneOrder;
        sequence_lengthx_onehot_order = glm_prep.sequence_lengthx_onehot_order;
        temp_locValidBool = glm_prep.temp_locValidBool;
        temp_locValidBool_real = glm_prep.temp_locValidBool_real;
        
        a = 1;
        %% get designMatrix_allSDT
        designMatrix_allSDT_P = sequence_lengthx_onehot & responseSeq_lengthx_onehot; % Hit
        designMatrix_allSDT_oneMinusP = sequence_lengthx_onehot & ~responseSeq_lengthx_onehot; % Miss
        designMatrix_allSDT_P0 = ~sequence_lengthx_onehot & responseSeq_lengthx_onehot; % False alarm
        designMatrix_allSDT_oneMinusP0 = ~sequence_lengthx_onehot & ~responseSeq_lengthx_onehot; % Correct rejection
        
        %tempBool_or = designMatrix_allSDT_P | designMatrix_allSDT_oneMinusP | designMatrix_allSDT_P0 | designMatrix_allSDT_oneMinusP0;
        %tempBool_and = designMatrix_allSDT_P & designMatrix_allSDT_oneMinusP & designMatrix_allSDT_P0 & designMatrix_allSDT_oneMinusP0;
        %tempBool_or_sum = sum(tempBool_or,'all');
        %tempBool_and_sum = sum(tempBool_and,'all');
        
        a = 1;
        
        designMatrix_allSDT = [designMatrix_allSDT_P,designMatrix_allSDT_oneMinusP,designMatrix_allSDT_P0,designMatrix_allSDT_oneMinusP0];
        
        %% fun_neuronGLM_eachLength_par, compure beta from GLM
        
        glm_options = struct;
        glm_options.order_glm_valid = order_glm_valid;
        glm_options.numFrames = numFrames;
        glm_options.plot_lengthFlag = plot_lengthFlag;
        glm_options.sequence_lengthx_onehot = sequence_lengthx_onehot;
        glm_options.sequence_lengthx_onehot_order = sequence_lengthx_onehot_order;
        glm_options.fun_lasso_jjb = fun_lasso_jjb;
        glm_options.fun_lasso_repeat_jjb = fun_lasso_repeat_jjb;
        glm_options.designMatrix_allSDT = designMatrix_allSDT;
        glm_options.if_seqSubspace_onlyhit0_allSDT1 = if_seqSubspace_onlyhit0_allSDT1;
        glm_options.if_strictHit = if_strictHit;
        glm_options.if_allSDT_withOnlyHit = if_allSDT_withOnlyHit;
        glm_options.if_seqTrialCountWeighted = if_seqTrialCountWeighted;
        glm_options.if_glm_12B0_24B1_6B2 = if_glm_12B0_24B1_6B2;
        
        a = 1;
        [glm_beta_lengthx_delay1Bin,glm_r2_lengthx_delay1Bin] = fun_neuronGLM_eachLength_par(F_dff_lengthx_delay1Bin,glm_options);
        [glm_beta_lengthx_delay2Bin,glm_r2_lengthx_delay2Bin] = fun_neuronGLM_eachLength_par(F_dff_lengthx_delay2Bin,glm_options);
        [glm_beta_lengthx_baselineBin,glm_r2_lengthx_baselineBin] = fun_neuronGLM_eachLength_par(F_dff_lengthx_baselineBin,glm_options);
        
        glm_beta_lengthx_delay1Bin_multiFOV = [glm_beta_lengthx_delay1Bin_multiFOV; glm_beta_lengthx_delay1Bin];  %#ok<*AGROW>
        glm_beta_lengthx_delay2Bin_multiFOV = [glm_beta_lengthx_delay2Bin_multiFOV; glm_beta_lengthx_delay2Bin];  %#ok<*AGROW>
        glm_beta_lengthx_baselineBin_multiFOV = [glm_beta_lengthx_baselineBin_multiFOV; glm_beta_lengthx_baselineBin];  %#ok<*AGROW>
        
        glm_r2_lengthx_delay1Bin_multiFOV = [glm_r2_lengthx_delay1Bin_multiFOV; glm_r2_lengthx_delay1Bin];  %#ok<*AGROW>
        glm_r2_lengthx_delay2Bin_multiFOV = [glm_r2_lengthx_delay2Bin_multiFOV; glm_r2_lengthx_delay2Bin];  %#ok<*AGROW>
        glm_r2_lengthx_baselineBin_multiFOV = [glm_r2_lengthx_baselineBin_multiFOV; glm_r2_lengthx_baselineBin];  %#ok<*AGROW>
        
        a = 1;
        
        fprintf('Progress %d/%d, t=%.1f secs.\n',tempIndex,fileSize,toc(t0));
        
        if if_singleFOV == 1
            break
        end
    end
    
end
a = 1; %#ok<*NASGU>



% numFrames = 6;
pointKindsNum = 4;
target_seqSet = get_target_seqSet_v1(numFrames,pointKindsNum);
numSeq = zeros(1,pointKindsNum);
for tempi=1:pointKindsNum
    numSeq(tempi) = length(target_seqSet{tempi});
end

target_seqSet_inOne = [];
for tempi=1:size(target_seqSet,2)
    target_seqSet_inOne = [target_seqSet_inOne target_seqSet{tempi}']; %#ok<*AGROW>
end

boolIndex_location_seq = false(numFrames,length(target_seqSet_inOne));
for tempi=1:length(target_seqSet_inOne)
    currentSequence = target_seqSet_inOne{tempi};
    boolIndex_location_seq(currentSequence,tempi) = true;
end

glm_dataLoad_Name_v = autoGetFunName_myScripts('fun_glm_dataLoad', [targetPATH '\functions']);
fun_glm_dataLoad = str2func(glm_dataLoad_Name_v);

a = 1;

if if_load_selectiveIndex == 1
    % Load all sessions' data files
    [glm_file_load, glm_file_name] = load_glmData_multiFOV_v1(searchName_glmData,path_glm);
    fileSize = length(glm_file_name);
    
    cellIndex_multiFOV = [];
    cellIndex_suite2p_multiFOV = [];
    cell_FOVIndex_multiFOV = [];
    FOVcellNum_multiFOV = [];
    FOVName_multiFOV = [];
    
    selectiveCellBoolIndex_length1_multiFOV = [];
    selectiveCellBoolIndex_length2_multiFOV = [];
    selectiveCellBoolIndex_length3_multiFOV = [];
    selectiveCellBoolIndex_length4_multiFOV = [];
    
    selectiveCellBoolIndex_length1_delay1Bin_multiFOV = [];
    selectiveCellBoolIndex_length2_delay1Bin_multiFOV = [];
    selectiveCellBoolIndex_length3_delay1Bin_multiFOV = [];
    selectiveCellBoolIndex_length4_delay1Bin_multiFOV = [];
    
    selectiveCellBoolIndex_length1_delay2Bin_multiFOV = [];
    selectiveCellBoolIndex_length2_delay2Bin_multiFOV = [];
    selectiveCellBoolIndex_length3_delay2Bin_multiFOV = [];
    selectiveCellBoolIndex_length4_delay2Bin_multiFOV = [];
    
    selectiveCellBoolIndex_length1_T_multiFOV = [];
    selectiveCellBoolIndex_length2_T_multiFOV = [];
    selectiveCellBoolIndex_length3_T_multiFOV = [];
    selectiveCellBoolIndex_length4_T_multiFOV = [];
    
    F_dff_lengthx_delay1Bin_seq_multiFOV = [];
    F_dff_lengthx_delay1Bin_seq_trial_multiFOV = [];
    
    selectiveCellBoolIndex_rProb_glm_baseline_multiFOV = [];
    selectiveCellBoolIndex_rProb_glm_delay1_multiFOV = [];
    selectiveCellBoolIndex_rProb_glm_delay2_multiFOV = [];
    selectiveCellBoolIndex_precision_glm_baseline_multiFOV = [];
    selectiveCellBoolIndex_precision_glm_delay1_multiFOV = [];
    selectiveCellBoolIndex_precision_glm_delay2_multiFOV = [];
    
    dff_seqMean_baseline_multiFOV = [];
    dff_seqMean_delay1_multiFOV = [];
    dff_seqMean_delay2_multiFOV = [];
    
    dff_seqFrame_baseline_multiFOV = [];
    dff_seqFrame_delay1_multiFOV = [];
    dff_seqFrame_delay2_multiFOV = [];
    
    r2_rProb_glm_baseline_multiFOV = [];
    r2_rProb_glm_delay1_multiFOV = [];
    r2_rProb_glm_delay2_multiFOV = [];
    
    p_rProb_glm_baseline_multiFOV = [];
    p_rProb_glm_delay1_multiFOV = [];
    p_rProb_glm_delay2_multiFOV = [];
    
    r2_precision_glm_baseline_multiFOV = [];
    r2_precision_glm_delay1_multiFOV = [];
    r2_precision_glm_delay2_multiFOV = [];
    
    p_precision_glm_baseline_multiFOV = [];
    p_precision_glm_delay1_multiFOV = [];
    p_precision_glm_delay2_multiFOV = [];
    
    for tempIndex=1:fileSize
        %if true
        %tempIndex = 20;
        % tempIndex = 1;
        
        if if_singleFOV == 1
            tempIndex = 20; %#ok<*FXSET>
        end
        
        tempstr = sprintf('file%d',tempIndex);
        tempData = glm_file_load.(tempstr).glmData;
        
        
        glm_dataLoad_options = struct;
        glm_dataLoad_options.tempData = tempData;
        glm_dataLoad_options.if_load_allTrial0_memoryCorrect1 = if_load_allTrial0_memoryCorrect1;
        glm_dataLoad_options.plot_lengthFlag = plot_lengthFlag;
        
        glm_dataLoad_output = fun_glm_dataLoad(glm_dataLoad_options);
        
        F_dff_lengthx_delay1Bin = glm_dataLoad_output.F_dff_lengthx_delay1Bin;
        sequence_lengthx_onehot = glm_dataLoad_output.sequence_lengthx_onehot;
        responseSeq_lengthx_onehot = glm_dataLoad_output.responseSeq_lengthx_onehot;
        
        [temp_roiNum,temp_trialNum] = size(F_dff_lengthx_delay1Bin);
        
        boolIndex_location_seq;
        boolIndex_location_seq_T = boolIndex_location_seq';
        
        seqNum = size(boolIndex_location_seq_T,1);
        seqIndex = zeros(temp_trialNum,1);
        for tempi=1:temp_trialNum
            temp_seq_onehot = sequence_lengthx_onehot(tempi,:);
            for tempj=1:seqNum
                if sum(temp_seq_onehot == boolIndex_location_seq_T(tempj,:)) == numFrames
                    break
                end
            end
            seqIndex(tempi) = tempj;
        end
        isCorrect = sum(sequence_lengthx_onehot == responseSeq_lengthx_onehot,2) == numFrames;
        
        F_dff_lengthx_delay1Bin_seq = zeros(temp_roiNum,seqNum,'single');
        F_dff_lengthx_delay1Bin_seq_trial = cell(temp_roiNum,seqNum);
        for tempi=1:seqNum
            %temptempBoolIndex = seqIndex==tempi;
            temptempBoolIndex = (seqIndex==tempi) & isCorrect;
            if sum(temptempBoolIndex) > 0
                F_dff_lengthx_delay1Bin_seq(:,tempi) = mean(F_dff_lengthx_delay1Bin(:,temptempBoolIndex),2);
                for tempj=1:temp_roiNum
                    F_dff_lengthx_delay1Bin_seq_trial{tempj,tempi} = F_dff_lengthx_delay1Bin(tempj,temptempBoolIndex);
                end
            end
        end
        a = 1;
        
        
        
        
        cellIndex = tempData.cellIndex;
        cellIndex_suite2p = tempData.cellIndex_suite2p;
        cell_FOVIndex = tempIndex*ones(length(cellIndex),1);
        FOVcellNum = length(cellIndex);
        
        selectiveCellBoolIndex_length1 = tempData.selectiveCellBoolIndex_length1;
        selectiveCellBoolIndex_length2 = tempData.selectiveCellBoolIndex_length2;
        selectiveCellBoolIndex_length3 = tempData.selectiveCellBoolIndex_length3;
        selectiveCellBoolIndex_length4 = tempData.selectiveCellBoolIndex_length4;
        
        selectiveCellBoolIndex_length1_delay1Bin = tempData.selectiveCellBoolIndex_length1_delay1Bin;
        selectiveCellBoolIndex_length2_delay1Bin = tempData.selectiveCellBoolIndex_length2_delay1Bin;
        selectiveCellBoolIndex_length3_delay1Bin = tempData.selectiveCellBoolIndex_length3_delay1Bin;
        selectiveCellBoolIndex_length4_delay1Bin = tempData.selectiveCellBoolIndex_length4_delay1Bin;
        
        selectiveCellBoolIndex_length1_delay2Bin = tempData.selectiveCellBoolIndex_length1_delay2Bin;
        selectiveCellBoolIndex_length2_delay2Bin = tempData.selectiveCellBoolIndex_length2_delay2Bin;
        selectiveCellBoolIndex_length3_delay2Bin = tempData.selectiveCellBoolIndex_length3_delay2Bin;
        selectiveCellBoolIndex_length4_delay2Bin = tempData.selectiveCellBoolIndex_length4_delay2Bin;
        
        selectiveCellBoolIndex_length1_T = tempData.selectiveCellBoolIndex_length1_T;
        selectiveCellBoolIndex_length2_T = tempData.selectiveCellBoolIndex_length2_T;
        selectiveCellBoolIndex_length3_T = tempData.selectiveCellBoolIndex_length3_T;
        selectiveCellBoolIndex_length4_T = tempData.selectiveCellBoolIndex_length4_T;
        
        [~,FOVName,~] = fileparts(glm_file_name{tempIndex});
        
        cellIndex_multiFOV = [cellIndex_multiFOV; cellIndex];
        cellIndex_suite2p_multiFOV = [cellIndex_suite2p_multiFOV; cellIndex_suite2p];
        cell_FOVIndex_multiFOV = [cell_FOVIndex_multiFOV; cell_FOVIndex];
        FOVcellNum_multiFOV = [FOVcellNum_multiFOV; FOVcellNum];
        FOVName_multiFOV = [FOVName_multiFOV; FOVName];
        
        selectiveCellBoolIndex_length1_multiFOV = [selectiveCellBoolIndex_length1_multiFOV; selectiveCellBoolIndex_length1];
        selectiveCellBoolIndex_length2_multiFOV = [selectiveCellBoolIndex_length2_multiFOV; selectiveCellBoolIndex_length2];
        selectiveCellBoolIndex_length3_multiFOV = [selectiveCellBoolIndex_length3_multiFOV; selectiveCellBoolIndex_length3];
        selectiveCellBoolIndex_length4_multiFOV = [selectiveCellBoolIndex_length4_multiFOV; selectiveCellBoolIndex_length4];
        
        selectiveCellBoolIndex_length1_delay1Bin_multiFOV = [selectiveCellBoolIndex_length1_delay1Bin_multiFOV; selectiveCellBoolIndex_length1_delay1Bin];
        selectiveCellBoolIndex_length2_delay1Bin_multiFOV = [selectiveCellBoolIndex_length2_delay1Bin_multiFOV; selectiveCellBoolIndex_length2_delay1Bin];
        selectiveCellBoolIndex_length3_delay1Bin_multiFOV = [selectiveCellBoolIndex_length3_delay1Bin_multiFOV; selectiveCellBoolIndex_length3_delay1Bin];
        selectiveCellBoolIndex_length4_delay1Bin_multiFOV = [selectiveCellBoolIndex_length4_delay1Bin_multiFOV; selectiveCellBoolIndex_length4_delay1Bin];
        
        selectiveCellBoolIndex_length1_delay2Bin_multiFOV = [selectiveCellBoolIndex_length1_delay2Bin_multiFOV; selectiveCellBoolIndex_length1_delay2Bin];
        selectiveCellBoolIndex_length2_delay2Bin_multiFOV = [selectiveCellBoolIndex_length2_delay2Bin_multiFOV; selectiveCellBoolIndex_length2_delay2Bin];
        selectiveCellBoolIndex_length3_delay2Bin_multiFOV = [selectiveCellBoolIndex_length3_delay2Bin_multiFOV; selectiveCellBoolIndex_length3_delay2Bin];
        selectiveCellBoolIndex_length4_delay2Bin_multiFOV = [selectiveCellBoolIndex_length4_delay2Bin_multiFOV; selectiveCellBoolIndex_length4_delay2Bin];
        
        selectiveCellBoolIndex_length1_T_multiFOV = [selectiveCellBoolIndex_length1_T_multiFOV; selectiveCellBoolIndex_length1_T];
        selectiveCellBoolIndex_length2_T_multiFOV = [selectiveCellBoolIndex_length2_T_multiFOV; selectiveCellBoolIndex_length2_T];
        selectiveCellBoolIndex_length3_T_multiFOV = [selectiveCellBoolIndex_length3_T_multiFOV; selectiveCellBoolIndex_length3_T];
        selectiveCellBoolIndex_length4_T_multiFOV = [selectiveCellBoolIndex_length4_T_multiFOV; selectiveCellBoolIndex_length4_T];
        
        F_dff_lengthx_delay1Bin_seq_multiFOV = [F_dff_lengthx_delay1Bin_seq_multiFOV; F_dff_lengthx_delay1Bin_seq];
        F_dff_lengthx_delay1Bin_seq_trial_multiFOV = [F_dff_lengthx_delay1Bin_seq_trial_multiFOV; F_dff_lengthx_delay1Bin_seq_trial];
        
        
        
        
        rProb_glm_output_all = tempData.rProb_glm_output_all;
        selectiveCellBoolIndex_rProb_glm_baseline = rProb_glm_output_all.baseline.selectiveCellBoolIndex_rProb_glm;
        selectiveCellBoolIndex_rProb_glm_delay1 = rProb_glm_output_all.delay1.selectiveCellBoolIndex_rProb_glm;
        selectiveCellBoolIndex_rProb_glm_delay2 = rProb_glm_output_all.delay2.selectiveCellBoolIndex_rProb_glm;
        
        selectiveCellBoolIndex_rProb_glm_baseline_multiFOV = [selectiveCellBoolIndex_rProb_glm_baseline_multiFOV; selectiveCellBoolIndex_rProb_glm_baseline];
        selectiveCellBoolIndex_rProb_glm_delay1_multiFOV = [selectiveCellBoolIndex_rProb_glm_delay1_multiFOV; selectiveCellBoolIndex_rProb_glm_delay1];
        selectiveCellBoolIndex_rProb_glm_delay2_multiFOV = [selectiveCellBoolIndex_rProb_glm_delay2_multiFOV; selectiveCellBoolIndex_rProb_glm_delay2];
        
        selectiveCellBoolIndex_precision_glm_baseline = rProb_glm_output_all.baseline.selectiveCellBoolIndex_precision_glm;
        selectiveCellBoolIndex_precision_glm_delay1 = rProb_glm_output_all.delay1.selectiveCellBoolIndex_precision_glm;
        selectiveCellBoolIndex_precision_glm_delay2 = rProb_glm_output_all.delay2.selectiveCellBoolIndex_precision_glm;
        
        selectiveCellBoolIndex_precision_glm_baseline_multiFOV = [selectiveCellBoolIndex_precision_glm_baseline_multiFOV; selectiveCellBoolIndex_precision_glm_baseline];
        selectiveCellBoolIndex_precision_glm_delay1_multiFOV = [selectiveCellBoolIndex_precision_glm_delay1_multiFOV; selectiveCellBoolIndex_precision_glm_delay1];
        selectiveCellBoolIndex_precision_glm_delay2_multiFOV = [selectiveCellBoolIndex_precision_glm_delay2_multiFOV; selectiveCellBoolIndex_precision_glm_delay2];
        
        
        dff_seqMean_baseline = rProb_glm_output_all.baseline.temp_F_dff_mean;
        dff_seqMean_delay1 = rProb_glm_output_all.delay1.temp_F_dff_mean;
        dff_seqMean_delay2 = rProb_glm_output_all.delay2.temp_F_dff_mean;
        
        dff_seqMean_baseline_multiFOV = [dff_seqMean_baseline_multiFOV; dff_seqMean_baseline];
        dff_seqMean_delay1_multiFOV = [dff_seqMean_delay1_multiFOV; dff_seqMean_delay1];
        dff_seqMean_delay2_multiFOV = [dff_seqMean_delay2_multiFOV; dff_seqMean_delay2];
        
        dff_seqFrame_baseline = rProb_glm_output_all.baseline.temp_F_dff_frame;
        dff_seqFrame_delay1 = rProb_glm_output_all.delay1.temp_F_dff_frame;
        dff_seqFrame_delay2 = rProb_glm_output_all.delay2.temp_F_dff_frame;
        
        dff_seqFrame_baseline_multiFOV = [dff_seqFrame_baseline_multiFOV; dff_seqFrame_baseline];
        dff_seqFrame_delay1_multiFOV = [dff_seqFrame_delay1_multiFOV; dff_seqFrame_delay1];
        dff_seqFrame_delay2_multiFOV = [dff_seqFrame_delay2_multiFOV; dff_seqFrame_delay2];
        
        
        r2_rProb_glm_baseline = rProb_glm_output_all.baseline.r2_rProb_glm;
        r2_rProb_glm_delay1 = rProb_glm_output_all.delay1.r2_rProb_glm;
        r2_rProb_glm_delay2 = rProb_glm_output_all.delay2.r2_rProb_glm;
        
        r2_rProb_glm_baseline_multiFOV = [r2_rProb_glm_baseline_multiFOV; r2_rProb_glm_baseline];
        r2_rProb_glm_delay1_multiFOV = [r2_rProb_glm_delay1_multiFOV; r2_rProb_glm_delay1];
        r2_rProb_glm_delay2_multiFOV = [r2_rProb_glm_delay2_multiFOV; r2_rProb_glm_delay2];
        
        p_rProb_glm_baseline = rProb_glm_output_all.baseline.p_rProb_glm;
        p_rProb_glm_delay1 = rProb_glm_output_all.delay1.p_rProb_glm;
        p_rProb_glm_delay2 = rProb_glm_output_all.delay2.p_rProb_glm;
        
        p_rProb_glm_baseline_multiFOV = [p_rProb_glm_baseline_multiFOV; p_rProb_glm_baseline];
        p_rProb_glm_delay1_multiFOV = [p_rProb_glm_delay1_multiFOV; p_rProb_glm_delay1];
        p_rProb_glm_delay2_multiFOV = [p_rProb_glm_delay2_multiFOV; p_rProb_glm_delay2];
        
        
        r2_precision_glm_baseline = rProb_glm_output_all.baseline.r2_precision_glm;
        r2_precision_glm_delay1 = rProb_glm_output_all.delay1.r2_precision_glm;
        r2_precision_glm_delay2 = rProb_glm_output_all.delay2.r2_precision_glm;
        
        r2_precision_glm_baseline_multiFOV = [r2_precision_glm_baseline_multiFOV; r2_precision_glm_baseline];
        r2_precision_glm_delay1_multiFOV = [r2_precision_glm_delay1_multiFOV; r2_precision_glm_delay1];
        r2_precision_glm_delay2_multiFOV = [r2_precision_glm_delay2_multiFOV; r2_precision_glm_delay2];
        
        p_precision_glm_baseline = rProb_glm_output_all.baseline.p_precision_glm;
        p_precision_glm_delay1 = rProb_glm_output_all.delay1.p_precision_glm;
        p_precision_glm_delay2 = rProb_glm_output_all.delay2.p_precision_glm;
        
        p_precision_glm_baseline_multiFOV = [p_precision_glm_baseline_multiFOV; p_precision_glm_baseline];
        p_precision_glm_delay1_multiFOV = [p_precision_glm_delay1_multiFOV; p_precision_glm_delay1];
        p_precision_glm_delay2_multiFOV = [p_precision_glm_delay2_multiFOV; p_precision_glm_delay2];
        
    end
end

selectiveCellBoolIndex_length1_multiFOV = selectiveCellBoolIndex_length1_multiFOV==1;
selectiveCellBoolIndex_length2_multiFOV = selectiveCellBoolIndex_length2_multiFOV==1;
selectiveCellBoolIndex_length3_multiFOV = selectiveCellBoolIndex_length3_multiFOV==1;
selectiveCellBoolIndex_length4_multiFOV = selectiveCellBoolIndex_length4_multiFOV==1;

selectiveCellBoolIndex_length1_delay1Bin_multiFOV = selectiveCellBoolIndex_length1_delay1Bin_multiFOV==1;
selectiveCellBoolIndex_length2_delay1Bin_multiFOV = selectiveCellBoolIndex_length2_delay1Bin_multiFOV==1;
selectiveCellBoolIndex_length3_delay1Bin_multiFOV = selectiveCellBoolIndex_length3_delay1Bin_multiFOV==1;
selectiveCellBoolIndex_length4_delay1Bin_multiFOV = selectiveCellBoolIndex_length4_delay1Bin_multiFOV==1;

selectiveCellBoolIndex_length1_delay2Bin_multiFOV = selectiveCellBoolIndex_length1_delay2Bin_multiFOV==1;
selectiveCellBoolIndex_length2_delay2Bin_multiFOV = selectiveCellBoolIndex_length2_delay2Bin_multiFOV==1;
selectiveCellBoolIndex_length3_delay2Bin_multiFOV = selectiveCellBoolIndex_length3_delay2Bin_multiFOV==1;
selectiveCellBoolIndex_length4_delay2Bin_multiFOV = selectiveCellBoolIndex_length4_delay2Bin_multiFOV==1;

selectiveCellBoolIndex_length1_T_multiFOV = selectiveCellBoolIndex_length1_T_multiFOV==1;
selectiveCellBoolIndex_length2_T_multiFOV = selectiveCellBoolIndex_length2_T_multiFOV==1;
selectiveCellBoolIndex_length3_T_multiFOV = selectiveCellBoolIndex_length3_T_multiFOV==1;
selectiveCellBoolIndex_length4_T_multiFOV = selectiveCellBoolIndex_length4_T_multiFOV==1;

selectiveCellBoolIndex_all_multiFOV = ...
    selectiveCellBoolIndex_length1_multiFOV |...
    selectiveCellBoolIndex_length2_multiFOV |...
    selectiveCellBoolIndex_length3_multiFOV |...
    selectiveCellBoolIndex_length4_multiFOV |...
    selectiveCellBoolIndex_length1_T_multiFOV |...
    selectiveCellBoolIndex_length2_T_multiFOV |...
    selectiveCellBoolIndex_length3_T_multiFOV |...
    selectiveCellBoolIndex_length4_T_multiFOV;

selectiveCellBoolIndex_length1_all_multiFOV = selectiveCellBoolIndex_length1_multiFOV | selectiveCellBoolIndex_length1_T_multiFOV;
selectiveCellBoolIndex_length2_all_multiFOV = selectiveCellBoolIndex_length2_multiFOV | selectiveCellBoolIndex_length2_T_multiFOV;
selectiveCellBoolIndex_length3_all_multiFOV = selectiveCellBoolIndex_length3_multiFOV | selectiveCellBoolIndex_length3_T_multiFOV;
selectiveCellBoolIndex_length4_all_multiFOV = selectiveCellBoolIndex_length4_multiFOV | selectiveCellBoolIndex_length4_T_multiFOV;

selectiveCellBoolIndex_rProb_glm_baseline_multiFOV = selectiveCellBoolIndex_rProb_glm_baseline_multiFOV==1;
selectiveCellBoolIndex_rProb_glm_delay1_multiFOV = selectiveCellBoolIndex_rProb_glm_delay1_multiFOV==1;
selectiveCellBoolIndex_rProb_glm_delay2_multiFOV = selectiveCellBoolIndex_rProb_glm_delay2_multiFOV==1;

selectiveCellBoolIndex_precision_glm_baseline_multiFOV = selectiveCellBoolIndex_precision_glm_baseline_multiFOV==1;
selectiveCellBoolIndex_precision_glm_delay1_multiFOV = selectiveCellBoolIndex_precision_glm_delay1_multiFOV==1;
selectiveCellBoolIndex_precision_glm_delay2_multiFOV = selectiveCellBoolIndex_precision_glm_delay2_multiFOV==1;
        
FOVcellNum_multiFOV;
FOVAllCellRange_multiFOV = zeros(size(FOVcellNum_multiFOV,1),2);
for tempi=1:size(FOVcellNum_multiFOV,1)
    FOVAllCellRange_multiFOV(tempi,1) = sum(FOVcellNum_multiFOV(1:tempi-1))+1;
    FOVAllCellRange_multiFOV(tempi,2) = sum(FOVcellNum_multiFOV(1:tempi));
end
a = 1;

selevtivity_multiFOV = struct;
selevtivity_multiFOV.cellIndex_multiFOV = cellIndex_multiFOV;
selevtivity_multiFOV.cellIndex_suite2p_multiFOV = cellIndex_suite2p_multiFOV;
selevtivity_multiFOV.cell_FOVIndex_multiFOV = cell_FOVIndex_multiFOV;
selevtivity_multiFOV.FOVcellNum_multiFOV = FOVcellNum_multiFOV;
selevtivity_multiFOV.FOVAllCellRange_multiFOV = FOVAllCellRange_multiFOV;
selevtivity_multiFOV.FOVName_multiFOV = FOVName_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_all_multiFOV = selectiveCellBoolIndex_all_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length1_all_multiFOV = selectiveCellBoolIndex_length1_all_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length2_all_multiFOV = selectiveCellBoolIndex_length2_all_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length3_all_multiFOV = selectiveCellBoolIndex_length3_all_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length4_all_multiFOV = selectiveCellBoolIndex_length4_all_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length1_multiFOV = selectiveCellBoolIndex_length1_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length2_multiFOV = selectiveCellBoolIndex_length2_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length3_multiFOV = selectiveCellBoolIndex_length3_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length4_multiFOV = selectiveCellBoolIndex_length4_multiFOV;

selevtivity_multiFOV.selectiveCellBoolIndex_length1_delay1Bin_multiFOV = selectiveCellBoolIndex_length1_delay1Bin_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length2_delay1Bin_multiFOV = selectiveCellBoolIndex_length2_delay1Bin_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length3_delay1Bin_multiFOV = selectiveCellBoolIndex_length3_delay1Bin_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length4_delay1Bin_multiFOV = selectiveCellBoolIndex_length4_delay1Bin_multiFOV;

selevtivity_multiFOV.selectiveCellBoolIndex_length1_delay2Bin_multiFOV = selectiveCellBoolIndex_length1_delay2Bin_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length2_delay2Bin_multiFOV = selectiveCellBoolIndex_length2_delay2Bin_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length3_delay2Bin_multiFOV = selectiveCellBoolIndex_length3_delay2Bin_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length4_delay2Bin_multiFOV = selectiveCellBoolIndex_length4_delay2Bin_multiFOV;

selevtivity_multiFOV.selectiveCellBoolIndex_length1_T_multiFOV = selectiveCellBoolIndex_length1_T_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length2_T_multiFOV = selectiveCellBoolIndex_length2_T_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length3_T_multiFOV = selectiveCellBoolIndex_length3_T_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_length4_T_multiFOV = selectiveCellBoolIndex_length4_T_multiFOV;

selevtivity_multiFOV.selectiveCellBoolIndex_rProb_glm_baseline_multiFOV = selectiveCellBoolIndex_rProb_glm_baseline_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_rProb_glm_delay1_multiFOV = selectiveCellBoolIndex_rProb_glm_delay1_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_rProb_glm_delay2_multiFOV = selectiveCellBoolIndex_rProb_glm_delay2_multiFOV;

selevtivity_multiFOV.selectiveCellBoolIndex_precision_glm_baseline_multiFOV = selectiveCellBoolIndex_precision_glm_baseline_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_precision_glm_delay1_multiFOV = selectiveCellBoolIndex_precision_glm_delay1_multiFOV;
selevtivity_multiFOV.selectiveCellBoolIndex_precision_glm_delay2_multiFOV = selectiveCellBoolIndex_precision_glm_delay2_multiFOV;

selevtivity_multiFOV.dff_seqMean_baseline_multiFOV = dff_seqMean_baseline_multiFOV;
selevtivity_multiFOV.dff_seqMean_delay1_multiFOV = dff_seqMean_delay1_multiFOV;
selevtivity_multiFOV.dff_seqMean_delay2_multiFOV = dff_seqMean_delay2_multiFOV;

selevtivity_multiFOV.dff_seqFrame_baseline_multiFOV = dff_seqFrame_baseline_multiFOV;
selevtivity_multiFOV.dff_seqFrame_delay1_multiFOV = dff_seqFrame_delay1_multiFOV;
selevtivity_multiFOV.dff_seqFrame_delay2_multiFOV = dff_seqFrame_delay2_multiFOV;

selevtivity_multiFOV.r2_rProb_glm_baseline_multiFOV = r2_rProb_glm_baseline_multiFOV;
selevtivity_multiFOV.r2_rProb_glm_delay1_multiFOV = r2_rProb_glm_delay1_multiFOV;
selevtivity_multiFOV.r2_rProb_glm_delay2_multiFOV = r2_rProb_glm_delay2_multiFOV;

selevtivity_multiFOV.p_rProb_glm_baseline_multiFOV = p_rProb_glm_baseline_multiFOV;
selevtivity_multiFOV.p_rProb_glm_delay1_multiFOV = p_rProb_glm_delay1_multiFOV;
selevtivity_multiFOV.p_rProb_glm_delay2_multiFOV = p_rProb_glm_delay2_multiFOV;

selevtivity_multiFOV.r2_precision_glm_baseline_multiFOV = r2_precision_glm_baseline_multiFOV;
selevtivity_multiFOV.r2_precision_glm_delay1_multiFOV = r2_precision_glm_delay1_multiFOV;
selevtivity_multiFOV.r2_precision_glm_delay2_multiFOV = r2_precision_glm_delay2_multiFOV;

selevtivity_multiFOV.p_precision_glm_baseline_multiFOV = p_precision_glm_baseline_multiFOV;
selevtivity_multiFOV.p_precision_glm_delay1_multiFOV = p_precision_glm_delay1_multiFOV;
selevtivity_multiFOV.p_precision_glm_delay2_multiFOV = p_precision_glm_delay2_multiFOV;

a = 1;

%% get B_fileName
if if_glm_12B0_24B1_6B2 == 1
    if if_singleFOV == 1
        B_fileName = 'glm_beta_singleFOV_allTrial_eachLength_full_weighted_24B.mat';
    elseif if_singleFOV == 0
        B_fileName = 'glm_beta_multiFOV_allTrial_eachLength_full_weighted_24B.mat';
    end
elseif if_glm_12B0_24B1_6B2 == 0
    if if_singleFOV == 1
        B_fileName = 'glm_beta_singleFOV_allTrial_eachLength_full_weighted_12plus12B.mat';
    elseif if_singleFOV == 0
        B_fileName = 'glm_beta_multiFOV_allTrial_eachLength_full_weighted_12plus12B.mat';
    end
elseif if_glm_12B0_24B1_6B2 == 2
    if if_singleFOV == 1
        B_fileName = 'glm_beta_singleFOV_allTrial_eachLength_full_weighted_6B.mat';
    elseif if_singleFOV == 0
        B_fileName = 'glm_beta_multiFOV_allTrial_eachLength_full_weighted_6B.mat';
    end
end
fprintf('B_fileName=%s.\n',B_fileName);


%% load
if if_load_beta == 1
    %load(['C:\ASDROOT\STUDY\temp','\data\',B_fileName]);
    load([path_glm,'\',B_fileName]);    
end

%% save
if if_save == 1
    
    if exist('glm_r2_lengthx_delay1Bin_multiFOV','var') == 0
        %save(['C:\ASDROOT\STUDY\temp','\data\',B_fileName],...
        save([path_glm,'\',B_fileName],...        
            'glm_beta_lengthx_delay1Bin_multiFOV',...
            'glm_beta_lengthx_delay2Bin_multiFOV',...
            'glm_beta_lengthx_baselineBin_multiFOV');
    else
        %save(['C:\ASDROOT\STUDY\temp','\data\',B_fileName],...
        save([path_glm,'\',B_fileName],...        
            'glm_beta_lengthx_delay1Bin_multiFOV',...
            'glm_beta_lengthx_delay2Bin_multiFOV',...
            'glm_beta_lengthx_baselineBin_multiFOV',...
            'glm_r2_lengthx_delay1Bin_multiFOV',...
            'glm_r2_lengthx_delay2Bin_multiFOV',...
            'glm_r2_lengthx_baselineBin_multiFOV');
    end
    
end


%%

r2Valid_Name_v = autoGetFunName_myScripts('fun_r2Valid', [targetPATH '\functions']);
fun_r2Valid = str2func(r2Valid_Name_v);

r2_threshold = 0.00;%0.05-->0.1
% r2Valid_boolIndex_inOne = fun_r2Valid(if_glm_12B0_24B1_6B2,r2_threshold,if_singleFOV);

a = 1;
glm_beta_lengthx_delay2Bin_multiFOV;
glm_r2_lengthx_delay2Bin_multiFOV;

FOV_num = fileSize;


%% End