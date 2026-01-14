%% Initialization
close all


if_compute = 1;

if if_compute == 1
    clear
    if_compute = 1;
    if_load_beta = 0;
else
    if_load_beta = 1;
end

% if_load_beta = 1;

if_save = 0;

if_glm_12B0_24B1_6B2 = 2;
if_singleFOV = 0;


% 1 B_P;
% 2 B_oneMinusP;
% 3 B_P0;
% 4 B_oneMinusP0;

x1_flag = 1;
x1_proj_flag = 2;

x2_flag = 1;
x2_proj_flag = 3;


if_plot_x_B0_Fmean1 = 1;

if_plot_proj_B0_Fmean1_trial2 = 1;% -1 means don't proj

if_plot_projTrial_hit1_miss2_FA3_CR4 = 1;



order_glm = 0;
plot_lengthFlag = 0;%0 mean all length

if_load_allTrial0_memoryCorrect1 = 0;%0

% lasso_repeatNum = 1;%10-->5-->3-->1

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

if_fit_PC = 0;
if_plotPCA = 1;
if_plotPCA_2d0_3d1 = 1;

if_profile = 0;


if if_profile == 1
    profile on
end

%%
targetPATH = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts';
cd(targetPATH)




path_behavior = 'D:\twoPhotonData_motionCorrected\behavior';

searchName_gAcc = 'from23-04-26to23-06-12_D_gAcc_1';
searchName_rProb = 'from23-04-26to23-06-12_D_offloadingProb_1';

% Load other processed results
load_gAcc = loadMat_single(searchName_gAcc, path_behavior);
gAcc_noChoice_collapsed_inOne = load_gAcc.gAcc_noChoice_collapsed_inOne;

% Load other processed results
load_rProb = loadMat_single(searchName_rProb, path_behavior);
offloadingProb_collapsed = load_rProb.offloadingProb_all;

offloadingProb_inOne = [];
pointKindsNum = 4;
for tempi=1:pointKindsNum
    offloadingProb_inOne = [offloadingProb_inOne offloadingProb_collapsed{tempi}']; %#ok<*AGROW>
end

cd(targetPATH)


if if_load_allTrial0_memoryCorrect1 == 0
    searchName_glmData = 'allTrial.mat';
elseif if_load_allTrial0_memoryCorrect1 == 1
    searchName_glmData = 'memoryCorrect.mat';
end

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
    [glm_file_load, glm_file_name] = load_glmData_multiFOV_v1(searchName_glmData,'C:\ASDROOT\STUDY\temp\data');
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
        glm_options.if_strictHit = if_strictHit;
        glm_options.if_allSDT_withOnlyHit = if_allSDT_withOnlyHit;                
        glm_options.if_seqTrialCountWeighted = if_seqTrialCountWeighted;
        glm_options.if_glm_12B0_24B1_6B2 = if_glm_12B0_24B1_6B2;
        
        a = 1;
        [glm_beta_lengthx_delay1Bin,glm_r2_lengthx_delay1Bin] = fun_neuronGLM_par(F_dff_lengthx_delay1Bin,glm_options);
        [glm_beta_lengthx_delay2Bin,glm_r2_lengthx_delay2Bin] = fun_neuronGLM_par(F_dff_lengthx_delay2Bin,glm_options);
        [glm_beta_lengthx_baselineBin,glm_r2_lengthx_baselineBin] = fun_neuronGLM_par(F_dff_lengthx_baselineBin,glm_options);
        
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
a = 1;

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
fprintf('B_fileName=%s.\n',B_fileName);


%% load
if if_load_beta == 1
    load(['C:\ASDROOT\STUDY\temp','\data\',B_fileName]);
    
    temp_locValidBool = true(1,numFrames);
    temp_locValidBool_real = true(1,numFrames);    
end

%% save
if if_save == 1
        
    if exist('glm_r2_lengthx_delay1Bin_multiFOV','var') == 0
        save(['C:\ASDROOT\STUDY\temp','\data\',B_fileName],...
            'glm_beta_lengthx_delay1Bin_multiFOV',...
            'glm_beta_lengthx_delay2Bin_multiFOV',...
            'glm_beta_lengthx_baselineBin_multiFOV');
    else
        save(['C:\ASDROOT\STUDY\temp','\data\',B_fileName],...
            'glm_beta_lengthx_delay1Bin_multiFOV',...
            'glm_beta_lengthx_delay2Bin_multiFOV',...
            'glm_beta_lengthx_baselineBin_multiFOV',...
            'glm_r2_lengthx_delay1Bin_multiFOV',...
            'glm_r2_lengthx_delay2Bin_multiFOV',...
            'glm_r2_lengthx_baselineBin_multiFOV');
    end
    
end


%% Get shuffled beta
num_roi = size(glm_beta_lengthx_baselineBin_multiFOV,1);
num_betaGroup = size(glm_beta_lengthx_baselineBin_multiFOV,2)/6;

shuffleBeta_Name_v = autoGetFunName_myScripts('shuffleBeta', [targetPATH '\functions']);
fun_shuffleBeta = str2func(shuffleBeta_Name_v);

glm_beta_lengthx_baselineBin_multiFOV_shuffled = ...
    fun_shuffleBeta(glm_beta_lengthx_baselineBin_multiFOV,num_roi,numFrames,num_betaGroup);
glm_beta_lengthx_delay1Bin_multiFOV_shuffled = ...
    fun_shuffleBeta(glm_beta_lengthx_delay1Bin_multiFOV,num_roi,numFrames,num_betaGroup);
glm_beta_lengthx_delay2Bin_multiFOV_shuffled = ...
    fun_shuffleBeta(glm_beta_lengthx_delay2Bin_multiFOV,num_roi,numFrames,num_betaGroup);

a = 1;


%% Compute PCA based on GLM
a = 1;

r2Valid_Name_v = autoGetFunName_myScripts('fun_r2Valid', [targetPATH '\functions']);
fun_r2Valid = str2func(r2Valid_Name_v);
trialProj2PCA_Name_v = autoGetFunName_myScripts('fun_trialProj2PCA', [targetPATH '\functions']);
fun_trialProj2PCA = str2func(trialProj2PCA_Name_v);
FmeanProj2PCA_Name_v = autoGetFunName_myScripts('fun_FmeanProj2PCA', [targetPATH '\functions']);
fun_FmeanProj2PCA = str2func(FmeanProj2PCA_Name_v);

r2_threshold = 0.00;%0.05-->0.1
% r2Valid_boolIndex_inOne = fun_r2Valid(if_glm_12B0_24B1_6B2,r2_threshold,if_singleFOV);
r2Valid_boolIndex_inOne = fun_r2Valid(if_glm_12B0_24B1_6B2,r2_threshold,if_singleFOV,B_fileName);

a = 1;

fprintf('if_plot_proj_B0_Fmean1_trial2 = %d.\n',if_plot_proj_B0_Fmean1_trial2);
if if_plot_proj_B0_Fmean1_trial2 == 2
    [F_dff_PCA,designMatrix_allSDT_multiFOV] = fun_trialProj2PCA(if_glm_12B0_24B1_6B2,r2Valid_boolIndex_inOne);
    
    F_dff_PCA_lengthx_baselineBin_multiFOV = F_dff_PCA.F_dff_PCA_lengthx_baselineBin_multiFOV;
    F_dff_PCA_lengthx_delay1Bin_multiFOV = F_dff_PCA.F_dff_PCA_lengthx_delay1Bin_multiFOV;
    F_dff_PCA_lengthx_delay2Bin_multiFOV = F_dff_PCA.F_dff_PCA_lengthx_delay2Bin_multiFOV;
    
    temp_range = (1:numFrames) + numFrames*(if_plot_projTrial_hit1_miss2_FA3_CR4-1);
    trialBoolIndex_location_proj = designMatrix_allSDT_multiFOV(:,temp_range) == 1;   
    
elseif if_plot_proj_B0_Fmean1_trial2 == 1
    F_designMatrixMean = fun_FmeanProj2PCA(if_glm_12B0_24B1_6B2,r2Valid_boolIndex_inOne,if_singleFOV,if_strictHit);
    %F_designMatrixMean = fun_FmeanProj2PCA(if_glm_12B0_24B1_6B2,r2Valid_boolIndex_inOne,if_singleFOV,0);
    
    F_designMatrixMean_baselineBin_multiFOV = F_designMatrixMean.F_designMatrixMean_baselineBin_multiFOV;
    F_designMatrixMean_delay1Bin_multiFOV = F_designMatrixMean.F_designMatrixMean_delay1Bin_multiFOV;
    F_designMatrixMean_delay2Bin_multiFOV = F_designMatrixMean.F_designMatrixMean_delay2Bin_multiFOV;    
end


if if_plotPCA == 1 && if_seqSubspace_onlyhit0_allSDT1 == 1
    close all
    
    fig1 = figure('Name','Fig1','NumberTitle','off');
    set(gcf,'Position',[35+000 35+000 1230 800]);%set figure, 原点的位置x，原点的位置y，宽，高，其坐标为point
    t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    
    score = cell(2,3);
    score_shuffled = cell(2,3);
    score_proj = cell(2,3);
    score_new = cell(2,3);
    explained = cell(2,3);
    pca_mu = cell(2,3);
    for trueTempIndex=1:2
        temp_x_min = zeros(1,3);
        temp_x_max = zeros(1,3);
        temp_y_min = zeros(1,3);
        temp_y_max = zeros(1,3);
        for tempIndex=1:3
            if tempIndex == 1
                glm_beta = glm_beta_lengthx_baselineBin_multiFOV;
                glm_beta_shuffled = glm_beta_lengthx_baselineBin_multiFOV_shuffled;
            elseif tempIndex == 2
                glm_beta = glm_beta_lengthx_delay1Bin_multiFOV;
                glm_beta_shuffled = glm_beta_lengthx_delay1Bin_multiFOV_shuffled;
            elseif tempIndex == 3
                glm_beta = glm_beta_lengthx_delay2Bin_multiFOV;
                glm_beta_shuffled = glm_beta_lengthx_delay2Bin_multiFOV_shuffled; 
            end       
            
            if if_plot_proj_B0_Fmean1_trial2 == 2
                if tempIndex == 1
                    F_dff_PCA_trial = F_dff_PCA_lengthx_baselineBin_multiFOV;
                elseif tempIndex == 2
                    F_dff_PCA_trial = F_dff_PCA_lengthx_delay1Bin_multiFOV;
                elseif tempIndex == 3
                    F_dff_PCA_trial = F_dff_PCA_lengthx_delay2Bin_multiFOV;
                end
                                
                F_dff_PCA_trial_location = cell(1,numFrames);
                F_dff_PCA_trial_location_mean = zeros(size(F_dff_PCA_trial,1),numFrames);
                for temptempi=1:numFrames
                    tempBoolIndex = trialBoolIndex_location_proj(:,temptempi);
                    temp_F_dff_PCA_trial_location = F_dff_PCA_trial(:,tempBoolIndex);
                    F_dff_PCA_trial_location_mean(:,temptempi) = mean(temp_F_dff_PCA_trial_location,2);
                    F_dff_PCA_trial_location{temptempi} = temp_F_dff_PCA_trial_location;
                end
                F_dff_PCA_trial_location;
                F_dff_PCA_trial_location_mean;
                
            elseif if_plot_proj_B0_Fmean1_trial2 == 1
                F_designMatrixMean_baselineBin_multiFOV;
                F_designMatrixMean_delay1Bin_multiFOV;
                F_designMatrixMean_delay2Bin_multiFOV;
                
                a = 1;
                
                if tempIndex == 1
                    F_mean = F_designMatrixMean_baselineBin_multiFOV;
                elseif tempIndex == 2
                    F_mean = F_designMatrixMean_delay1Bin_multiFOV;
                elseif tempIndex == 3
                    F_mean = F_designMatrixMean_delay2Bin_multiFOV;
                end      
                a = 1;
                
                F_P = F_mean(:,(1:numFrames)+numFrames*0);
                F_oneMinusP = F_mean(:,(1:numFrames)+numFrames*1);
                F_P0 = F_mean(:,(1:numFrames)+numFrames*2);
                F_oneMinusP0 = F_mean(:,(1:numFrames)+numFrames*3); 
            end
            

            glm_beta_raw = glm_beta;
            glm_beta = glm_beta_raw(r2Valid_boolIndex_inOne,:);
            glm_beta_shuffled_raw = glm_beta_shuffled;
            glm_beta_shuffled = glm_beta_shuffled_raw(r2Valid_boolIndex_inOne,:);

            a = 1;
            
            B_P = glm_beta(:,(1:numFrames)+numFrames*0);
            B_oneMinusP = glm_beta(:,(1:numFrames)+numFrames*1);
            B_P0 = glm_beta(:,(1:numFrames)+numFrames*2);
            B_oneMinusP0 = glm_beta(:,(1:numFrames)+numFrames*3);
            
            B_P_shuffled = glm_beta_shuffled(:,(1:numFrames)+numFrames*0);
            B_oneMinusP_shuffled = glm_beta_shuffled(:,(1:numFrames)+numFrames*1);
            B_P0_shuffled = glm_beta_shuffled(:,(1:numFrames)+numFrames*2);
            B_oneMinusP0_shuffled = glm_beta_shuffled(:,(1:numFrames)+numFrames*3);
            
            a = 1;
            
            if trueTempIndex == 1
                x_flag = x1_flag;
                x_proj_flag = x1_proj_flag;
            elseif trueTempIndex == 2
                x_flag = x2_flag;
                x_proj_flag = x2_proj_flag;
            end
            
            a = 1;
            if x_flag == 1
                x = B_P';
                x_shuffled = B_P_shuffled';
            elseif x_flag == 2
                x = B_oneMinusP';
                %x_shuffled = B_P_shuffled';
                x_shuffled = B_oneMinusP_shuffled';
            elseif x_flag == 3
                x = B_P0';
                %x_shuffled = B_P_shuffled';
                x_shuffled = B_P0_shuffled';
            elseif x_flag == 4
                x = B_oneMinusP0';
                %x_shuffled = B_P_shuffled';
                x_shuffled = B_oneMinusP0_shuffled';
            end
            
            a = 1;
            
            % mark here
            if if_plot_proj_B0_Fmean1_trial2 ~= 1
                if x_proj_flag == 1
                    x_proj = B_P';
                elseif x_proj_flag == 2
                    x_proj = B_oneMinusP';
                elseif x_proj_flag == 3
                    x_proj = B_P0';
                elseif x_proj_flag == 4
                    x_proj = B_oneMinusP0';
                end
            elseif if_plot_proj_B0_Fmean1_trial2 == 1
                if x_proj_flag == 1
                    x_proj = F_P';
                elseif x_proj_flag == 2
                    x_proj = F_oneMinusP';
                elseif x_proj_flag == 3
                    x_proj = F_P0';
                elseif x_proj_flag == 4
                    x_proj = F_oneMinusP0';
                end                
            end
            
            a = 1;
            if_centered = 1;
            
            if if_centered == 0
                [temp_coeff,temp_score,temp_latent,temp_tsquared,temp_explained,temp_mu] = pca(x,'Centered',false);                
            elseif if_centered == 1
                [temp_coeff,temp_score,temp_latent,temp_tsquared,temp_explained,temp_mu] = pca(x);
            end
            score{trueTempIndex,tempIndex} = temp_score;
            explained{trueTempIndex,tempIndex} = temp_explained;
            pca_mu{trueTempIndex,tempIndex} = temp_mu;
            
            if if_centered == 0
                temp_score_shuffled = x_shuffled(1,:) * temp_coeff;
                temp_score_proj = x_proj * temp_coeff;
            elseif if_centered == 1
                temp_score_shuffled = (x_shuffled(1,:)-temp_mu) * temp_coeff;
                temp_score_proj = (x_proj-mean(x_proj,1)) * temp_coeff;
            end
            
            score_shuffled{trueTempIndex,tempIndex} = temp_score_shuffled;
            score_proj{trueTempIndex,tempIndex} = temp_score_proj;

            
            if if_plot_x_B0_Fmean1 == 1 && if_plot_proj_B0_Fmean1_trial2 == 1
                if x_flag == 1
                    x_new = F_P';
                elseif x_flag == 2
                    x_new = F_oneMinusP';
                elseif x_flag == 3
                    x_new = F_P0';
                elseif x_flag == 4
                    x_new = F_oneMinusP0';
                end
                
                if if_centered == 0
                    temp_score_new = x_new * temp_coeff;
                elseif if_centered == 1
                    temp_score_new = (x_new-mean(x_new,1)) * temp_coeff;
                end
                score_new{trueTempIndex,tempIndex} = temp_score_new;
                
                score_raw = score;
                score = score_new;
                temp_score = temp_score_new;
            end
            
            
            [temp_x_min(tempIndex),temp_x_max(tempIndex)] = bounds(temp_score(:,1));
            [temp_y_min(tempIndex),temp_y_max(tempIndex)] = bounds(temp_score(:,2));
        end
        
        x_min = min(temp_x_min);
        x_max = max(temp_x_max);
        y_min = min(temp_y_min);
        y_max = max(temp_y_max);
        
        delta_x = x_max-x_min;
        delta_y = y_max-y_min;
        if delta_x > delta_y
            y_min = y_min - (delta_x-delta_y)/2;
            y_max = y_max + (delta_x-delta_y)/2;
        elseif delta_x <= delta_y
            x_min = x_min - (delta_y-delta_x)/2;
            x_max = x_max + (delta_y-delta_x)/2;
        end
        
        for tempIndex=1:3
            
            nexttile
            
            multi_rgbColor = ...
                [228,26,28;
                55,126,184;
                77,175,74;
                152,78,163;
                255,127,0;
                255,255,51]/255;            
            
            backgounrdColor = [1 1 1]*0.825;%0.875
                        
            if if_plotPCA_2d0_3d1 == 0
                plot(score{trueTempIndex,tempIndex}(temp_locValidBool,1),score{trueTempIndex,tempIndex}(temp_locValidBool,2),'k','lineWidth',1.5);
            elseif if_plotPCA_2d0_3d1 == 1
                plot3(score{trueTempIndex,tempIndex}(temp_locValidBool,1),score{trueTempIndex,tempIndex}(temp_locValidBool,2),score{trueTempIndex,tempIndex}(temp_locValidBool,3),'k','lineWidth',1.5);
            end
            hold on
            
            
            % plot shuffle
            if if_plotPCA_2d0_3d1 == 0
                scatter(score_shuffled{trueTempIndex,tempIndex}(1),score_shuffled{trueTempIndex,tempIndex}(2),60,[0.25 0.25 0.25],'filled');%36
            elseif if_plotPCA_2d0_3d1 == 1
                scatter3(score_shuffled{trueTempIndex,tempIndex}(1),score_shuffled{trueTempIndex,tempIndex}(2),score_shuffled{trueTempIndex,tempIndex}(3),60,[0.25 0.25 0.25],'filled');
            end
            hold on
                    
            h_plot = [];
            tempj = 0;
            for tempi=1:numFrames
                tempj = tempj + 1;                              
                
                temph = [];
                if temp_locValidBool_real(tempi) == false
                elseif temp_locValidBool_real(tempi) == true
                    if if_plotPCA_2d0_3d1 == 0
                        temph = scatter(score{trueTempIndex,tempIndex}(tempj,1),score{trueTempIndex,tempIndex}(tempj,2),120,multi_rgbColor(tempi,:),'filled');%36
                    elseif if_plotPCA_2d0_3d1 == 1
                        temph = scatter3(score{trueTempIndex,tempIndex}(tempj,1),score{trueTempIndex,tempIndex}(tempj,2),score{trueTempIndex,tempIndex}(tempj,3),120,multi_rgbColor(tempi,:),'filled');
                    end
                end
                hold on
                h_plot = [h_plot temph];
            end
            
            if if_plot_proj_B0_Fmean1_trial2 == 0 || if_plot_proj_B0_Fmean1_trial2 == 1
                tempj = 0;
                for tempi=1:numFrames
                    tempj = tempj + 1;
                    if temp_locValidBool_real(tempi) == false
                    elseif temp_locValidBool_real(tempi) == true
                        if if_plotPCA_2d0_3d1 == 0
                            scatter(score_proj{trueTempIndex,tempIndex}(tempj,1),score_proj{trueTempIndex,tempIndex}(tempj,2),30,multi_rgbColor(tempi,:),'filled');%36
                        elseif if_plotPCA_2d0_3d1 == 1
                            scatter3(score_proj{trueTempIndex,tempIndex}(tempj,1),score_proj{trueTempIndex,tempIndex}(tempj,2),score_proj{trueTempIndex,tempIndex}(tempj,3),30,multi_rgbColor(tempi,:),'filled');
                        end
                    end
                    hold on
                end
            elseif if_plot_proj_B0_Fmean1_trial2 == 2
                F_dff_PCA_trial_location_mean;
                for tempi=1:numFrames
                        if if_plotPCA_2d0_3d1 == 0
                            scatter(F_dff_PCA_trial_location_mean(1,tempi),F_dff_PCA_trial_location_mean(2,tempi),30,multi_rgbColor(tempi,:),'filled');%36
                        elseif if_plotPCA_2d0_3d1 == 1
                            scatter3(F_dff_PCA_trial_location_mean(1,tempi),F_dff_PCA_trial_location_mean(2,tempi),...
                                F_dff_PCA_trial_location_mean(3,tempi),30,multi_rgbColor(tempi,:),'filled');
                        end
                        
                    hold on   
                end
            end
            
            grid on
            
            temp_labels = {'1','2','3','4','5','6'};
            
            xlim([x_min,x_max]);
            ylim([y_min,y_max]);
            if if_plotPCA_2d0_3d1 == 1
                zlim([mean([x_min y_min]),mean([x_max y_max])]);
            end
            set(gca,'color',backgounrdColor);
            %axis equal
            xlabel(sprintf('1st Principal Component (%.1f%%)',explained{trueTempIndex,tempIndex}(1)),'fontsize',16)
            ylabel(sprintf('2nd Principal Component (%.1f%%)',explained{trueTempIndex,tempIndex}(2)),'fontsize',16)
            if if_plotPCA_2d0_3d1 == 1
                zlabel(sprintf('3rd Principal Component (%.1f%%)',explained{trueTempIndex,tempIndex}(3)))
            end
            

            if tempIndex == 1
                tempstr3 = sprintf('baseline');
            elseif tempIndex == 2
                tempstr3 = sprintf('delay1');
            elseif tempIndex == 3
                tempstr3 = sprintf('delay2');
            end
            
            if order_glm == 0
                title(tempstr3,'FontSize',14);
            else
                tempstr2 = sprintf('length%d,order%d',plot_lengthFlag,order_glm);
                title(sprintf('%s,%s',tempstr2,tempstr3),'FontSize',14);
            end
            
            
        end
        
    end
    
end

a = 1;



if if_profile == 1
    profile viewer
end

%% End