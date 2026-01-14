function [F_designMatrixMean]=...
    fun_FmeanProj2PCA_v7(if_glm_12B0_24B1_6B2,r2Valid_boolIndex_inOne,if_singleFOV,if_strictHit)
%% Initialization

% if_singleFOV = 0;

order_glm = 0;
plot_lengthFlag = 0;%0 mean all length

if_load_allTrial0_memoryCorrect1 = 0;

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

numFrames = 6;


%% get B_fileName
if if_glm_12B0_24B1_6B2 == 1
    if if_singleFOV == 1
        B_fileName = 'glm_beta_singleFOV_allTrial_allLength_full_weighted_24B.mat';
    elseif if_singleFOV == 0
        B_fileName = 'glm_beta_multiFOV_allTrial_allLength_full_weighted_24B.mat';
    end
elseif if_glm_12B0_24B1_6B2 == 0
    if if_singleFOV == 1
        B_fileName = 'glm_beta_singleFOV_allTrial_allLength_full_weighted_12plus12B.mat';
    elseif if_singleFOV == 0
        B_fileName = 'glm_beta_multiFOV_allTrial_allLength_full_weighted_12plus12B.mat';
    end
elseif if_glm_12B0_24B1_6B2 == 2
    if if_singleFOV == 1
        B_fileName = 'glm_beta_singleFOV_allTrial_allLength_full_weighted_6B.mat';
    elseif if_singleFOV == 0
        B_fileName = 'glm_beta_multiFOV_allTrial_allLength_full_weighted_6B.mat';
    end
end
% fprintf('B_fileName=%s.\n',B_fileName);


%% load
load(['C:\ASDROOT\STUDY\temp','\data\',B_fileName]); %#ok<*LOAD>

glm_beta_lengthx_baselineBin_multiFOV;
glm_beta_lengthx_delay1Bin_multiFOV;
glm_beta_lengthx_delay2Bin_multiFOV;

a = 1; %#ok<*NASGU>
[num_roi,num_B] = size(glm_beta_lengthx_baselineBin_multiFOV);


%% Compute PCA based on GLM, in hit subspace
pca_score = cell(1,3);
pca_coeff = cell(1,3);
for tempIndex=1:3
    if tempIndex == 1
        glm_beta = glm_beta_lengthx_baselineBin_multiFOV;
    elseif tempIndex == 2
        glm_beta = glm_beta_lengthx_delay1Bin_multiFOV;
    elseif tempIndex == 3
        glm_beta = glm_beta_lengthx_delay2Bin_multiFOV;
    end
    
    glm_beta_raw = glm_beta;
    glm_beta = glm_beta_raw(r2Valid_boolIndex_inOne,:);
    
    a = 1;
    
    B_P = glm_beta(:,(1:numFrames)+numFrames*0);
    
    x = B_P';
    
    [temp_coeff,temp_score,temp_latent,temp_tsquared,temp_explained,temp_mu] = pca(x); %#ok<*ASGLU>
    pca_score{tempIndex} = temp_score;
    pca_coeff{tempIndex} = temp_coeff;
end
a = 1;

t0 = tic;
%% Compute
% if if_compute == 1
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

% Load all sessions' data files
[glm_file_load, glm_file_name] = load_glmData_multiFOV_v1(searchName_glmData,'C:\ASDROOT\STUDY\temp\data');
fileSize = length(glm_file_name);


F_designMatrixMean_baselineBin_multiFOV = [];
F_designMatrixMean_delay1Bin_multiFOV = [];
F_designMatrixMean_delay2Bin_multiFOV = [];



lasso_jjb_Name_v = autoGetFunName_myScripts('lasso_jjb', [targetPATH '\functions']);
fun_lasso_jjb = str2func(lasso_jjb_Name_v);
lasso_repeat_jjb_Name_v = autoGetFunName_myScripts('lasso_repeat_jjb', [targetPATH '\functions']);
fun_lasso_repeat_jjb = str2func(lasso_repeat_jjb_Name_v);
neuronGLM_par_Name_v = autoGetFunName_myScripts('fun_glm_neuron_par', [targetPATH '\functions']);
fun_neuronGLM_par = str2func(neuronGLM_par_Name_v);
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
    
    %% fun_neuronGLM_par, compure beta from GLM
    
    glm_options = struct;
    glm_options.order_glm_valid = order_glm_valid;
    %glm_options.lasso_repeatNum = lasso_repeatNum;
    glm_options.numFrames = numFrames;
    glm_options.plot_lengthFlag = plot_lengthFlag;
    glm_options.sequence_lengthx_onehot = sequence_lengthx_onehot;
    glm_options.sequence_lengthx_onehot_order = sequence_lengthx_onehot_order;
    glm_options.fun_lasso_jjb = fun_lasso_jjb;
    glm_options.fun_lasso_repeat_jjb = fun_lasso_repeat_jjb;
    glm_options.designMatrix_allSDT = designMatrix_allSDT;
    glm_options.if_seqSubspace_onlyhit0_allSDT1 = if_seqSubspace_onlyhit0_allSDT1;
    glm_options.if_seqTrialCountWeighted = if_seqTrialCountWeighted;
    glm_options.if_glm_12B0_24B1_6B2 = if_glm_12B0_24B1_6B2; %#ok<*STRNU>
    
    a = 1;
    
    
    F_dff_lengthx_baselineBin;
    F_dff_lengthx_delay1Bin;
    F_dff_lengthx_delay2Bin;
    
    designMatrix_allSDT; %#ok<*VUNUS>
    
    F_dff_lengthx_baselineBin_zscore = ...
        (F_dff_lengthx_baselineBin-mean(F_dff_lengthx_baselineBin,2))./std(F_dff_lengthx_baselineBin,0,2);
    F_dff_lengthx_delay1Bin_zscore = ...
        (F_dff_lengthx_delay1Bin-mean(F_dff_lengthx_delay1Bin,2))./std(F_dff_lengthx_delay1Bin,0,2);
    F_dff_lengthx_delay2Bin_zscore = ...
        (F_dff_lengthx_delay2Bin-mean(F_dff_lengthx_delay2Bin,2))./std(F_dff_lengthx_delay2Bin,0,2);

    
    num_roi;
    [temp_num_roi,temp_num_trial] = size(F_dff_lengthx_baselineBin);

    num_B = size(designMatrix_allSDT,2);
    
    F_designMatrixMean_baselineBin = zeros(temp_num_roi,num_B);
    F_designMatrixMean_delay1Bin = zeros(temp_num_roi,num_B);
    F_designMatrixMean_delay2Bin = zeros(temp_num_roi,num_B);
    
    a = 1;
    
    for temptempi=1:numFrames
        if if_strictHit == 0
            tempBoolIndex = designMatrix_allSDT(:,temptempi) == 1;
        elseif if_strictHit == 1
            tempBoolIndex1 = designMatrix_allSDT(:,temptempi) == 1;
            tempBoolIndex2 = sum(designMatrix_allSDT(:,(numFrames+1):(numFrames*3)),2) == 0;
            tempBoolIndex = tempBoolIndex1 & tempBoolIndex2;
        end
        F_designMatrixMean_baselineBin(:,temptempi) = mean(F_dff_lengthx_baselineBin_zscore(:,tempBoolIndex),2);
        F_designMatrixMean_delay1Bin(:,temptempi) = mean(F_dff_lengthx_delay1Bin_zscore(:,tempBoolIndex),2);
        F_designMatrixMean_delay2Bin(:,temptempi) = mean(F_dff_lengthx_delay2Bin_zscore(:,tempBoolIndex),2);
    end    
    for temptempi=(numFrames+1):num_B
        tempBoolIndex = designMatrix_allSDT(:,temptempi) == 1;
        F_designMatrixMean_baselineBin(:,temptempi) = mean(F_dff_lengthx_baselineBin_zscore(:,tempBoolIndex),2);
        F_designMatrixMean_delay1Bin(:,temptempi) = mean(F_dff_lengthx_delay1Bin_zscore(:,tempBoolIndex),2);
        F_designMatrixMean_delay2Bin(:,temptempi) = mean(F_dff_lengthx_delay2Bin_zscore(:,tempBoolIndex),2);
    end
    
    F_designMatrixMean_baselineBin_multiFOV = [F_designMatrixMean_baselineBin_multiFOV; F_designMatrixMean_baselineBin]; %#ok<*AGROW>
    F_designMatrixMean_delay1Bin_multiFOV = [F_designMatrixMean_delay1Bin_multiFOV; F_designMatrixMean_delay1Bin];
    F_designMatrixMean_delay2Bin_multiFOV = [F_designMatrixMean_delay2Bin_multiFOV; F_designMatrixMean_delay2Bin];
    
    a = 1;
    
    if if_singleFOV == 1
        break
    end
end

% end
fprintf('Progress %d/%d, t=%.1f secs.\n',tempIndex,fileSize,toc(t0));
a = 1;

F_designMatrixMean_baselineBin_multiFOV_raw = F_designMatrixMean_baselineBin_multiFOV;
F_designMatrixMean_delay1Bin_multiFOV_raw = F_designMatrixMean_delay1Bin_multiFOV;
F_designMatrixMean_delay2Bin_multiFOV_raw = F_designMatrixMean_delay2Bin_multiFOV;

F_designMatrixMean_baselineBin_multiFOV = F_designMatrixMean_baselineBin_multiFOV(r2Valid_boolIndex_inOne,:);
F_designMatrixMean_delay1Bin_multiFOV = F_designMatrixMean_delay1Bin_multiFOV(r2Valid_boolIndex_inOne,:);
F_designMatrixMean_delay2Bin_multiFOV = F_designMatrixMean_delay2Bin_multiFOV(r2Valid_boolIndex_inOne,:);

F_designMatrixMean = struct;
F_designMatrixMean.F_designMatrixMean_baselineBin_multiFOV = F_designMatrixMean_baselineBin_multiFOV;
F_designMatrixMean.F_designMatrixMean_delay1Bin_multiFOV = F_designMatrixMean_delay1Bin_multiFOV;
F_designMatrixMean.F_designMatrixMean_delay2Bin_multiFOV = F_designMatrixMean_delay2Bin_multiFOV;


%% End